#include "dpd/core/affinity_queue.hpp"

#include <iostream>
#include <random>

#define REQUIRE(cond) \
    { assert(cond); if(!(cond)){ fprintf(stderr, "Fail : " # cond); exit(1);  }}


std::mt19937_64 rng(1);

void t1()
{
    std::cerr<<"t1\n";

    AffinityQueue aq(2);

    aq.AddTask(0, [=](){
    }, {} );

    aq.SealTasks();

    aq.RunTasks();
}

void t2()
{
    std::cerr<<"t2\n";

    AffinityQueue aq(2);

    volatile int a=0, b=0;
    int n=1000000;

    aq.AddTask(0, [&](){
        for(int i=0; i<n; i++){
            a += 1;
        }
    }, {} );
    aq.AddTask(1, [&](){
        for(int i=0; i<n; i++){
            b += 1;
        }
    }, {} );

    aq.SealTasks();

    std::cerr<<"Run1\n";
    aq.RunTasks();
    REQUIRE(a==n);
    REQUIRE(b==n);
    std::cerr<<"Run2\n";
    aq.RunTasks();
    REQUIRE(a==2*n);
    REQUIRE(b==2*n);
        std::cerr<<"Fin\n";
}

void t3()
{
    std::cerr<<"t3\n";

    AffinityQueue aq(2);

    int n=1000000;
    volatile int a=0;

    int t1 = aq.AddTask(0, [&](){
        for(int i=0; i<n; i++){
            a += 1;
        }
    }, {} );
    
    aq.AddTask(0, [&](){
        for(int i=0; i<n; i++){
            a += 2;
        }
    }, {t1} );

    aq.SealTasks();

    aq.RunTasks();
    REQUIRE(a == 3*n);
    aq.RunTasks();
    REQUIRE(a == 6*n);
}

void t4()
{
    std::cerr<<"t4\n";

    AffinityQueue aq(2);

    int n=1000000;
    volatile int a=0;

    int t1 = aq.AddTask(0, [&](){
        for(int i=0; i<n; i++){
            a += 1;
        }
    }, {} );
    
    aq.AddTask(1, [&](){
        for(int i=0; i<n; i++){
            a += 2;
        }
    }, {t1} );

    aq.SealTasks();

    aq.RunTasks();
    REQUIRE(a == 3*n);
    aq.RunTasks();
    REQUIRE(a == 6*n);
}

void t5()
{
    std::cerr<<"t5\n";

    AffinityQueue aq(2);

    int n=1000;
    volatile int a=0;

    int tprev=aq.AddTask(0, [&](){
        a += 1;
    }, {} );

    for(int i=1; i<n; i++){
        tprev=aq.AddTask(0, [&](){
            a += 1;
        }, {tprev} );
    }

    aq.SealTasks();

    aq.RunTasks();
    REQUIRE(a == n);
    aq.RunTasks();
    REQUIRE(a == 2*n);
}

void t6()
{
    std::cerr<<"t6\n";

    AffinityQueue aq(2);

    int n=1000;
    volatile int a=0;

    int tprev=aq.AddTask(0, [&](){
        a += 1;
    }, {} );

    for(int i=1; i<n; i++){
        tprev=aq.AddTask((i%2), [&](){
            a += 1;
        }, {tprev} );
    }

    aq.SealTasks();

    aq.RunTasks();
    REQUIRE(a == n);
    aq.RunTasks();
    REQUIRE(a == 2*n);
}

void t7()
{
    std::cerr<<"t7\n";

    AffinityQueue aq(2);

    int n=1000;
    
    std::atomic<int> a;
    a=0;

    for(int i=0; i<n; i++){
        aq.AddTask( (i%aq.NumQueues()), [&](){
            a.fetch_add(1);
        }, {});
    }

    aq.SealTasks();

    aq.RunTasks();
    REQUIRE(a.load() == n);
    aq.RunTasks();
    REQUIRE(a.load() == 2*n);
}

void t8()
{
    std::cerr<<"t8\n";

    AffinityQueue aq(2);

    int n=1000;
    
    std::atomic<int> a;
    a=0;

    std::vector<intptr_t> prev;

    for(int i=0; i<n; i++){
        int tid=aq.AddTask( (i%aq.NumQueues()), [&](){
            a.fetch_add(1);
        }, prev);
        prev.push_back(tid);
    }

    aq.SealTasks();

    aq.RunTasks();
    REQUIRE(a.load() == n);
    aq.RunTasks();
    REQUIRE(a.load() == 2*n);
}


void t9()
{
    std::cerr<<"t9\n";

    AffinityQueue aq(2);

    std::atomic<int> count;
    count=0;

    std::function<void(int d, intptr_t parent)> recurse;
    recurse=[&](int d, intptr_t parent)
    {
        intptr_t me=aq.AddTask( rng()%aq.NumQueues(), [=,&count](){
            //std::cerr<<" count="<<count.load()<<"\n";
            count.fetch_add(1);
        }, {parent});
        if(d < 10){
            recurse(d+1, me);
            recurse(d+1, me);
        }
    };

    recurse(0, -1);

    //std::cerr<<"Sealing\n";
    aq.SealTasks();

    //std::cerr<<"Running\n";
    aq.RunTasks();
    REQUIRE(count.load() == 2047);
    aq.RunTasks();
    REQUIRE(count.load() == 2*2047);
}


int main()
{
    for(int i=0; i<100; i++){
        t1();
        t2();
        t3();
        t4();
        t5();
        t6();
        t7();
        t8();
        t9();
    }

    std::cerr<<"Success\n";
}
