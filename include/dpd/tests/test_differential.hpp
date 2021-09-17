#ifndef test_differential_hpp
#define test_differential_hpp

#include "test_base.hpp"
#include "dpd/core/dpd_engine.hpp"
#include "dpd/core/dpd_state_io.hpp"

std::pair<bool,std::string> test_differential(
    TestBase &test,
    DPDEngine &engine1,
    DPDEngine &engine2
){
    WorldState state1, state2;

    try{
        state1 = test.create_world();
        validate(state1);

        state2 = state1;

        std::string no_support=engine2.CanSupport(&state2);
        if(!no_support.empty()){
            return {true, "SKIP : "+no_support};
        }

        engine1.Attach(&state1);
        engine2.Attach(&state2);

        unsigned steps_done=0;

        while(1){
            if(steps_done > 100000){
                return {false, "More than 100000 steps done without test finishing."};
            }
            unsigned todo=test.get_advance_count(state1);
            if(todo > 100000){
                return {false, "Test asked to do more than 100000 steps."};
            }
            if(todo==0){
                return {true, {}};
            }

            for(unsigned i=0; i<state1.beads.size(); i++){
                state2.beads[i] = state1.beads[i];
            }

            engine1.Run(todo);
            engine2.Run(todo);

            validate(state1);
            validate(state2);

            if(state1.t != state2.t){
                std::cerr<<"Ref is at time  "<<state1.t<<", dut is at time "<<state2.t<<".\n";
            }

            //std::cerr<<"t="<<state1.t<<"\n";
            

            double maxDiff= 0.001;

            double xdiff=0;
            int xdiff_index=-1;
            double vdiff=0;
            double fdiff=0;
            int fdiff_index=-1;
            for(unsigned i=0; i<state1.beads.size(); i++){
                //std::cerr<<" ref="<<state1.beads[i].x<<", "<<state2.beads[i].x<<"\n";
                vec3r_t dx=state1.beads[i].x-state2.beads[i].x;
                for(int j=0; j<3; j++){
                    if(dx[j] > state1.box[j]/2){
                        dx[j] -= state1.box[j];
                    }
                    if(dx[j] < -state1.box[j]/2){
                        dx[j] += state1.box[j];
                    }
                }
                double xdiffnow=(dx).l2_norm();
                if(xdiffnow > xdiff){
                    xdiff=xdiffnow;
                    xdiff_index=i;
                }
                vdiff=std::max(vdiff, (state1.beads[i].v-state2.beads[i].v).l2_norm() );
                double fdiffnow=(state1.beads[i].f-state2.beads[i].f).l2_norm();
                if(fdiffnow > fdiff){
                    fdiff=fdiffnow;
                    fdiff_index=i;
                }
            }
            if(fdiff>10*maxDiff){ // Forces can be large, so allow more slack
                std::cerr<<"At time "<<state1.t<<" the f values have diverged.\n";
                

                write_world_state(std::cerr, state1);
                write_world_state(std::cerr, state2);

                std::cerr<<"  fdiff_max at bead "<<fdiff_index<<", diff="<<fdiff<<"\n";
                std::cerr<<"      ref="<<state1.beads[fdiff_index].f<<"\n";
                std::cerr<<"      got="<<state2.beads[fdiff_index].f<<"\n";

                return {false, "Force divergence."};
            }

            if(xdiff>maxDiff){
                std::cerr<<"At time "<<state1.t<<" the x values have diverged.\n";

                write_world_state(std::cerr, state1);
                write_world_state(std::cerr, state2);

                std::cerr<<"  xdiff_max at bead "<<xdiff_index<<", diff="<<xdiff<<"\n";
                std::cerr<<"      ref="<<state1.beads[xdiff_index].x<<"\n";
                std::cerr<<"      got="<<state2.beads[xdiff_index].x<<"\n";

                return {false, "Position divergence."};
            }
            if(vdiff>maxDiff){
                std::cerr<<"At time "<<state1.t<<" the v values have diverged.\n";

                write_world_state(std::cerr, state1);
                write_world_state(std::cerr, state2);

                return {false, "Velocity divergence."};
            }

            validate(state2);

            steps_done += todo;
        }
    }catch(TestFailedException &e){
        write_world_state(std::cerr, state1);
        write_world_state(std::cerr, state2);
        return {false, e.what()};
    }catch(ValidationFailedError &e){
        write_world_state(std::cerr, state1);
        write_world_state(std::cerr, state2);
        return {false, e.what()};
    }catch(std::exception &e){
        return {false, std::string("Exception : ")+e.what()};
    }catch(...){
        return {false, "Unknown exception occured."};
    }
}

#endif
