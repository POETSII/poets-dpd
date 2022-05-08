#ifndef async_host_link_hpp
#define async_host_link_hpp

#ifndef PDPD_TINSEL

#include <thread>
#include <mutex>
#include <condition_variable>
#include <cassert>
#include <vector>
#include <queue>
#include <atomic>
#include <cstring>

/* All public methods should be considered non thread safe.
    They cannot be called in parallel with themselves.
*/
template<class TImpl>
class AsyncHostLink
{
public:
    using HostLinkParams = typename TImpl::HostLinkParams;
    using HostLink = typename TImpl::HostLink;

private:
    enum State
    {
        None,
        Aquiring,
        Aquired
    };

    HostLinkParams m_params;

    std::mutex m_lock;
    std::condition_variable m_cond;

    State m_state = None;
    std::thread m_aquire_thread;

    std::shared_ptr<HostLink> m_hostlink;

    void aquire_thread()
    {
        std::unique_lock lk(m_lock);
        assert(m_state==Aquiring);
        if(m_hostlink){
            fprintf(stderr, "AsyncHostLink - hostlink already aquired.");
            exit(1);
        }
        lk.unlock();

        //fprintf(stderr, "AsyncHostLink::beginAquire worker started.\n");
        auto hostlink=std::make_shared<HostLink>(m_params);
        //fprintf(stderr, "AsyncHostLink::beginAquire worker finished.\n");

        lk.lock();
        m_state=Aquired;
        m_hostlink=hostlink;
        m_cond.notify_all();
    }

    void beginAquisition(std::unique_lock<std::mutex> &lk)
    {
        assert(lk.owns_lock());
        if(m_state!=Aquired){
            if(m_state!=Aquiring){
                m_state=Aquiring;
                m_aquire_thread=std::thread([&](){
                    aquire_thread();
                });
            }
        }
    }

    void completeAquisition(std::unique_lock<std::mutex> &lk)
    {
        assert(m_state != None);
        if(m_state==Aquiring){
            fprintf(stderr, "AsyncHostLink::completeAquisition - waiting for worker thread.\n");
            m_cond.wait(lk, [&]() -> bool {return m_state==Aquired;} );

            assert(m_state==Aquired);
        }

        if(m_aquire_thread.joinable()){
            m_aquire_thread.join();
        }
    }

public:
    AsyncHostLink()
    {}

    AsyncHostLink(const HostLinkParams &params)
        : m_params(params)
    {}

    // I think they dont exist anyway due to members not having them.
    AsyncHostLink(const AsyncHostLink &) = delete;
    AsyncHostLink &operator=(const AsyncHostLink &) = delete;

    ~AsyncHostLink()
    {
        reset();
    }

    void SetParams(const HostLinkParams &params)
    {
        std::unique_lock lk(m_lock);
        if(m_state!=None){
            throw std::runtime_error("Attempt to change params after aquire started.");
        }
        m_params=params;
    }

    void beginAquisition()
    {
        std::unique_lock lk(m_lock);
        beginAquisition(lk);
    }

    //! Return true if the hostlink has been aquired and initialised
    bool isAquired() const
    {
        std::unique_lock lk(m_lock);
        return m_state==Aquired;
    }

    //! Blocks until the hostlink has been aquired and initialised
    void completeAquisition()
    {
        std::unique_lock lk(m_lock);
        if(m_state==Aquired){
            // no-op
        }else if(m_state==Aquiring){
            completeAquisition(lk);
        }else{
            m_hostlink=std::make_shared<HostLink>(m_params);
            m_state=Aquired;
        }
    }

    HostLink *get()
    {
        HostLink *p=m_hostlink.get();
        if(p){
            return p;
        }

        completeAquisition();
        
        return m_hostlink.get();
    }

    HostLink *operator->()
    {
        return m_hostlink.get();
    }

    void reset()
    {
        std::unique_lock lk(m_lock);
        if(m_aquire_thread.joinable()){
            lk.unlock();
            m_aquire_thread.join();
            lk.lock();
        }
        m_state=None;
        m_hostlink.reset();
    }

};

template<class TImpl>
class AsyncMessageMover
{
public:
    using HostLink = typename TImpl::HostLink;
    const int BytesPerMsg = 1<<TImpl::LogBytesPerMsg;
private:
    std::mutex m_mutex;
    std::condition_variable m_cond;
    std::thread m_worker_thread;

    HostLink *m_hostlink = nullptr;

    std::atomic<bool> m_reader_waiting;
    std::atomic<bool> m_reader_quitting;
    std::queue<std::vector<char>> m_queue;
    std::vector<std::vector<char>> m_free;
    int m_page_low_water_mark  = 1<<12;
    int m_page_high_water_mark = 3<<12;
    int m_page_size            = 4<<12;
    int m_queue_max_bytes      = 1<<24;
    int m_queue_total_bytes=0;

    void worker()
    {
        std::vector<char> page(m_page_size);
        int page_valid=0;

        uint64_t blocked_read_count=0;
        while(!m_reader_quitting.load()){
            assert(page.size()==(size_t)m_page_size);

            if(m_queue_total_bytes >= m_queue_max_bytes){
                std::unique_lock lk(m_mutex);
                m_cond.wait(lk, [&](){
                    return m_reader_quitting.load() || m_queue_total_bytes < m_queue_max_bytes;
                });
                if(m_reader_quitting.load()){
                    break;
                }
            }

            int todo=m_page_size - page_valid;

            int new_page_valid=page_valid;
            m_hostlink->recvBulkNonBlock(
                page.size(),
                &new_page_valid,
                &page[0]
            );

            if(new_page_valid==page_valid){
                blocked_read_count++;
            }else{
                blocked_read_count=0;
                fprintf(stderr, "  Got %u bytes\n", new_page_valid-page_valid);
            }
            page_valid=new_page_valid;

            bool reader_waiting=m_reader_waiting.load();

            if(
                (page_valid>=m_page_high_water_mark)
                || (reader_waiting && (page_valid>=m_page_low_water_mark) )
                || (reader_waiting && (page_valid>=BytesPerMsg) && blocked_read_count)
            ){
                // Flush the page
                std::unique_lock lk(m_mutex);
                // TODO : this is too much work to do under the lock

                std::vector<char> next_page;
                if(!m_free.empty()){
                    next_page=std::move(m_free.back());
                    m_free.pop_back();
                }
                next_page.resize(m_page_size);
                int next_page_valid=0;

                int full=page_valid/BytesPerMsg;
                int partial=page_valid%BytesPerMsg;
                assert(full>0);

                fprintf(stderr, "fLUSH : full=%d, partial=%d\n", full, partial);

                if(partial){
                    memmove(&next_page[0], &page[page_valid-partial], partial);
                    next_page_valid=partial;
                    page.resize(full*BytesPerMsg);
                }
                    
                m_queue_total_bytes+=page.size();
                assert(page.size()>0);
                m_queue.push(std::move(page));
                m_cond.notify_one();

                std::swap(next_page, page);
                std::swap(next_page_valid, page_valid);
            }else if(blocked_read_count<4){
                std::this_thread::yield();
            }else{
                // TODO : Arghh! This is horrible. Should switch to blocking read on host-link socket
                usleep(100);
            }
        }
    }
public:
    AsyncMessageMover(HostLink *hostlink)
        : m_hostlink(hostlink)
    {
        m_reader_waiting=false;
        m_reader_quitting=false;
        m_worker_thread=std::thread([&](){ worker(); });
    }

    ~AsyncMessageMover()
    {
        assert(m_worker_thread.joinable());
        m_reader_quitting=true;
        m_worker_thread.join();
    }

    void pull_page(std::vector<char> &page)
    {
        std::unique_lock lk(m_mutex);
        fprintf(stderr, "pull_page\n");

        if(page.capacity()>0){
            m_free.push_back(std::move(page)); // TODO: can this grow without bounds?
        }

        if(m_queue.empty()){
            m_reader_waiting=true;
            m_cond.wait(lk, [&](){ return !m_queue.empty(); });
            m_reader_waiting=false;
        }

        page=std::move(m_queue.front());
        m_queue.pop();
        m_queue_total_bytes -= page.size();

        assert(page.size()>0);
        int full=page.size()/BytesPerMsg;
        fprintf(stderr, "  Pulled %d\n", full);
    }
};

#endif

#endif
