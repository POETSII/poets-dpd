#ifndef affinity_queue_hpp
#define affinity_queue_hpp

#include <pthread.h>
#include <cassert>

#include <thread>
#include <mutex>
#include <vector>
#include <stdexcept>
#include <functional>
#include <atomic>
#include <condition_variable>
#include <iostream>

#include <sys/types.h>

#include "cache_info.hpp"

class AffinityQueue
{
public:
    AffinityQueue(const AffinityQueue &) = delete;
    AffinityQueue &operator=(const AffinityQueue &) = delete;

    AffinityQueue()
        : AffinityQueue(cache_info::get_system_cache_info().get_llc_nodes())
    {}

    AffinityQueue(unsigned threads_per_affinity_set)
        : AffinityQueue(make_cpu_sets(threads_per_affinity_set))
    {}

    AffinityQueue(const std::vector<cpu_set_obj> &affinity_sets)
    {
        //std::cerr<<" affinity_sets.size()=="<<affinity_sets.size()<<"\n";
        for(unsigned i=0; i<affinity_sets.size(); i++){
            m_queues.push_back( std::make_unique<TaskQueue>(*this, i, affinity_sets[i]) );
        }
    }

    ~AffinityQueue()
    {
        // Request quit ahead.
        // Destructor will then join them all
        for(unsigned i=0; i<m_queues.size(); i++){
            m_queues[i]->request_quit();
        }
    }

    unsigned NumQueues() const
    { return m_queues.size(); }

    // Not thread safe
    intptr_t AddTask(unsigned queue_id, std::function<void()> execute, const std::vector<intptr_t> &parents)
    {
        if(m_sentinel){
            throw std::runtime_error("AddTask after queue is sealed.");
        }

        intptr_t id=m_tasks.size();

        auto &queue=*m_queues.at(queue_id);

        std::unique_ptr<Task> t=std::make_unique<Task>();
        t->task_id=id;
        t->queue_id=queue_id;
        t->blockers=0;
        t->blockers_reset_value=0;
        t->execute=execute;

        int num_parents=0;
        for(auto pid : parents){
            if(pid==-1){
                continue; // We allow -1 to represent invalid handle
            }
            Task *pt=m_tasks.at(pid);
            pt->following.push_back(t.get());
            t->blockers_reset_value += 1;
            t->blockers += 1;
            num_parents += 1;
        }

        m_tasks.push_back(t.get());
        queue.tasks.push_back(std::move(t));

        return id;
    }

    //! Adds a dependency saying parent must finish before child
    /*! Warning: this can create cyclic dependencies, and they
        will not automatically be detected.
    */
    void AddDependency(intptr_t parent, intptr_t child)
    {
        if(m_sentinel){
            throw std::runtime_error("Tasks are sealed.");
        }

        if(parent==-1 || child==-1){
            return;
        }

        assert(child);

        auto *p=m_tasks.at(parent);
        auto *c=m_tasks.at(child);

        p->following.push_back(c);
        c->blockers += 1;
        c->blockers_reset_value += 1;
    }

    // Must be called before RunTasks.
    // No more tasks can be added after
    // Performs a linear scan of all tasks...
    void SealTasks()
    {
        // Create a sentinel task
        m_sentinel=std::make_unique<Task>();
        m_sentinel->task_id=-1;
        m_sentinel->queue_id=-1;
        m_sentinel->blockers=0;
        m_sentinel->blockers_reset_value=0;
        m_sentinel->execute=[&](){
            std::unique_lock<std::mutex> lk(m_global_mutex);
            m_running=false;
            m_global_cond.notify_one();
        };

        //std::cerr<<"Sealing\n";
        for(auto tt : m_tasks){
            //std::cerr<<"  Seal "<<tt->task_id<<", followers.size()=="<<tt->following.size()<<"\n";
            if(tt->following.empty()){
                //std::cerr<<"    leaf task\n";
                tt->following.push_back(m_sentinel.get());
                m_sentinel->blockers += 1;
                m_sentinel->blockers_reset_value += 1;
            }
            if(tt->blockers_reset_value==0){
                m_queues[tt->queue_id]->ready_reset_value.push_back(tt);
            }
        }
    }

    void RunTasks()
    {
        if(!m_sentinel){
            SealTasks();
        }

        assert(!m_running);
        for(auto t : m_tasks){
            assert(t->blockers == t->blockers_reset_value);
        }
        for(auto &pqueue : m_queues){
            assert(pqueue->ready_complete==0);
            assert(pqueue->ready.empty());
        }
        assert(m_sentinel->blockers==m_sentinel->blockers_reset_value);

        // We do this outside of locks, as worker threads will not look
        // at queues until m_running is true
        for(auto &pqueue : m_queues){
            auto &queue=*pqueue;
            queue.ready_complete=0;
            queue.ready.insert(queue.ready.end(), queue.ready_reset_value.begin(), queue.ready_reset_value.end());
        }

        m_running=true;

        for(auto &pqueue : m_queues){
            auto &queue=*pqueue;
            std::unique_lock<std::mutex> lk(queue.mutex);
            queue.cond.notify_all();
        }

        std::unique_lock<std::mutex> lk(m_global_mutex);
        m_global_cond.wait(lk, [&](){
            return !m_running;
        });
        m_sentinel->blockers=m_sentinel->blockers_reset_value;


        for(auto &pqueue : m_queues){
            pqueue->ready_complete=0;
            pqueue->ready.clear();
        }

        assert(!m_running);
        for(auto t : m_tasks){
            assert(t->blockers == t->blockers_reset_value);
        }
        for(auto &pqueue : m_queues){
            assert(pqueue->ready_complete==0);
            assert(pqueue->ready.empty());
        }
        assert(m_sentinel->blockers==m_sentinel->blockers_reset_value);
    }

private:
    struct Task
    {
        intptr_t task_id;
        int queue_id;
        std::atomic<int> blockers;
        int blockers_reset_value;
        std::vector<Task*> following;
        std::function<void()> execute;
    };

    struct TaskQueue
    {
        TaskQueue(AffinityQueue &parent, int _id, cpu_set_obj _affinity)
            : queue_id(_id)
            , affinity(_affinity)
            , waiters(0)
        {
            thread_count=_affinity.size();
            for(unsigned i=0; i<thread_count; i++){
                workers.emplace_back([&](){
                    parent.task_worker(*this);
                });
            }
        }

        void request_quit()
        {
            if(!quit){
                std::unique_lock<std::mutex> lk(mutex);
                quit=true;
                cond.notify_all();
            }
        }

        ~TaskQueue()
        {
            request_quit();

            for(unsigned i=0; i<thread_count; i++){
                workers[i].join();
            }
        }

        int queue_id;
        cpu_set_obj affinity;
        unsigned thread_count;

        std::mutex mutex;
        std::condition_variable cond;

        unsigned ready_complete = 0;
        std::vector<Task*> ready;

        bool quit = false;
        int waiters = 0;

        std::vector<Task*> ready_reset_value; // Tasks with no predecessors
        std::vector<std::unique_ptr<Task>> tasks;

        std::vector<std::thread> workers;
    };

    std::vector<std::unique_ptr<TaskQueue>> m_queues;
    std::vector<Task*> m_tasks;

    std::vector<std::thread> m_threads;

    // Used to protext the queue and detect when task graph has finished
    std::mutex m_global_mutex;
    std::condition_variable m_global_cond;
    bool m_running=false;                    // Set before running; cleared by the sentinel task
    std::unique_ptr<Task> m_sentinel;  // Sentinel task that follows all leaf (childless) tasks

    static std::vector<cpu_set_obj> make_cpu_sets(unsigned threads_per_affinity_set)
    {
        // TODO : This is not stable, and not the right approach. Should be using
        // libnuma or equivalent.
        unsigned nthreads=std::thread::hardware_concurrency();
        if(nthreads%threads_per_affinity_set){
            throw std::runtime_error("threads_per_affinity_set does not divide pthread_getconcurrency");
        }
        //std::cerr<<"nthreads="<<nthreads<<"\n";

        std::vector<cpu_set_obj> sets;
        for(unsigned i=0; i<nthreads; i += threads_per_affinity_set){
            cpu_set_t s;
            CPU_ZERO(&s);
            for(unsigned j=i; j<i+threads_per_affinity_set; j++){
                CPU_SET(j, &s);
            }
            sets.push_back(s);
        }
        return sets;
    }

    void task_worker(TaskQueue &queue)
    {
        std::vector<Task*> local_ready;

        auto notify_child=[&](Task *child)
        {
            unsigned count=child->blockers.fetch_sub(1);
            //std::cerr<<" notify_child "<<child->task_id<<", new count="<<count<<"\n";
            if(count==1){
                // This thread was the last to decrement, so it
                // must "own" the thread. No other thread can touch anything in the
                // task until blockers drops to 0.
                assert(child->blockers.load()==0);

                if(child->queue_id==queue.queue_id){
                    local_ready.push_back(child);
                }else if(child->queue_id>=0){
                    auto &other_queue=*m_queues[child->queue_id];

                    std::unique_lock<std::mutex> lk(other_queue.mutex);
                    other_queue.ready.push_back(child);
                    if(other_queue.waiters){
                        other_queue.cond.notify_one();
                    }
                }else{
                    //std::cerr<<"Executing sentinel task\n";
                    assert(child->queue_id==-1);
                    child->execute();
                }
            }
        };

        pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &queue.affinity.set);

        ///////////////////////////////////////////
        // Enter critical section
        std::unique_lock<std::mutex> lk(queue.mutex);

        while(!queue.quit){        
            if(!m_running || queue.ready_complete==queue.ready.size()){
                queue.waiters++;
                queue.cond.wait(lk);
                queue.waiters--;
                continue;
            }

            assert(queue.ready_complete < queue.ready.size());

            //std::cerr<<"  Thread "<<&lk<<" : ready size="<<(queue.ready.size()-queue.ready_complete)<<"\n";

            Task *t=queue.ready[queue.ready_complete];
            queue.ready_complete++;
            assert(t);

            //std::cerr<<"  Thread "<<&lk<<" : got task "<<t->task_id<<"\n";


            lk.unlock();
            // Exit critical section
            ///////////////////////////////////

            assert(t->blockers==0);
            t->execute();
            t->blockers = t->blockers_reset_value;

            // Do this outside the lock, as there are atomics and sometimes we have to notify other queues
            //std::cerr<<"  t->following.size()=="<<t->following.size()<<"\n";
            for(auto child : t->following){
                // This populates the local_ready queue
                notify_child(child);
            }

            ///////////////////////////////////
            // Enter critical section
            lk.lock();

            //std::cerr<<"Entered lock, local_ready.size()="<<local_ready.size()<<"\n";

            for(unsigned i=0; i<local_ready.size(); i++){
                assert(local_ready[i]->blockers==0);
            }

            if(!local_ready.empty()){
                queue.ready.insert(queue.ready.end(), local_ready.begin(), local_ready.end());
                local_ready.clear();
            }
        }
    }
};

#endif
