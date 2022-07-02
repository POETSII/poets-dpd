#ifndef non_conflict_grid_hpp
#define non_conflict_grid_hpp

#include <cassert>
#include <vector>
#include <functional>
#include <tbb/parallel_for.h>
#include <tbb/flow_graph.h>
#include <tbb/flow_graph_abstractions.h>
#include <tbb/blocked_range.h>

struct graph_grid_by_parfor
{
    std::array<std::vector<std::function<void()>>,8> tasks;

    template<class F>
    graph_grid_by_parfor(unsigned w, unsigned h, unsigned d, F make_task)
    {
        construct(w,h,d,make_task);
    }

    template<class F>
    void construct(unsigned w, unsigned h, unsigned d, F make_task)
    {
        assert(tasks[0].empty());

        auto to_linear=[&](unsigned x, unsigned y, unsigned z)
        { return x + y*w + z*w*h; };

        for(int x=0; x<w; x++){
            for(int y=0; y<h; y++){
                for(int z=0; z<d; z++){
                    int ox=x&1, oy=y&1, oz=z&1;

                    tasks[ox+2*oy+4*oz].push_back(make_task(x,y,z));
                }
            }
        }
    }

    void execute()
    {
       for(int i=0; i<8; i++){
        tbb::blocked_range<unsigned> r(0, tasks[i].size());
        tbb::parallel_for(r, [&](const tbb::blocked_range<unsigned> &rr){
            for(unsigned j=rr.begin(); j<rr.end(); j++){
                tasks[i][j]();
            }
        });
       }
    }
};

/*
    Schedules tasks in a 3d grid such that no touching cubes
    execute at the same time. This includes touching via
    faces, edges, and corners.
*/
struct graph_grid_by_flow
{
    typedef tbb::flow::continue_node< tbb::flow::continue_msg > node_t;
    typedef const tbb::flow::continue_msg &msg_t;

    tbb::flow::graph g;
    tbb::flow::broadcast_node< tbb::flow::continue_msg > input;
    std::vector<node_t> nodes;

    template<class F>
    graph_grid_by_flow(unsigned w, unsigned h, unsigned d, F make_task)
        : input(g)
    {
        auto to_linear=[&](unsigned x, unsigned y, unsigned z)
        { return x + y*w + z*w*h; };

        auto get_node=[&](int x, int y, int z) -> node_t &
        { return nodes.at(to_linear( (x+w)%w, (y+h)%h, (z+d)%d )); };

        nodes.reserve(w*h*d);

        for(int z=0; z<d; z++){        
            for(int y=0; y<h; y++){
                for(int x=0; x<w; x++){
                    auto t=make_task(x,y,z);
                    nodes.emplace_back(g, [=](msg_t){
                        t();
                    });
                }
            }
        }

        auto for_each_with_off=[&](int ox, int oy, int oz, auto f){
            for(int x=ox; x<w; x+=2){
                for(int y=oy; y<h; y+=2){
                    for(int z=oz; z<d; z+=2){
                        f(x,y,z);
                    }
                }
            }
        };

        for_each_with_off(0, 0, 0, [&](int x, int y, int z){
            auto &n=get_node(x,y,z);
            make_edge(input, n );
            make_edge(n, get_node(x-1, y, z) );
            make_edge(n, get_node(x+1, y, z) );
        });

        for_each_with_off(1, 0, 0, [&](int x, int y, int z){
            auto &n=get_node(x,y,z);
            make_edge(n, get_node(x, y-1, z) );
            make_edge(n, get_node(x, y+1, z) );
        });

        for_each_with_off(1, 1, 0, [&](int x, int y, int z){
            auto &n=get_node(x,y,z);
            make_edge(n, get_node(x-1, y, z) );
            make_edge(n, get_node(x+1, y, z) );
        });

        for_each_with_off(0, 1, 0, [&](int x, int y, int z){
            auto &n=get_node(x,y,z);
            make_edge(input, get_node(x,y,z-1) );
            make_edge(input, get_node(x,y,z+1) );
        });

        for_each_with_off(0, 1, 1, [&](int x, int y, int z){
            auto &n=get_node(x,y,z);
            make_edge(input, get_node(x,y-1,z) );
            make_edge(input, get_node(x,y+1,z) );
        });

        for_each_with_off(0, 0, 1, [&](int x, int y, int z){
            auto &n=get_node(x,y,z);
            make_edge(input, get_node(x-1,y,z) );
            make_edge(input, get_node(x+1,y,z) );
        });

        for_each_with_off(1, 0, 1, [&](int x, int y, int z){
            auto &n=get_node(x,y,z);
            make_edge(input, get_node(x,y-1,z) );
            make_edge(input, get_node(x,y+1,z) );
        });

        // All the (1,1,1) nodes are leaves.
    }

    void execute()
    {
        input.try_put(tbb::flow::continue_msg());
        g.wait_for_all();
    }
};

#endif

