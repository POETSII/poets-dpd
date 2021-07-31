#ifndef test_differential_hpp
#define test_differential_hpp

#include "test_base.hpp"
#include "../dpd_engine.hpp"
#include "../dpd_state_io.hpp"

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
            validate(state1);

            engine2.Run(todo);
            

            double maxDiff= 0.0001;

            double xdiff=0;
            double vdiff=0;
            double fdiff=0;
            for(unsigned i=0; i<state1.beads.size(); i++){
                xdiff=std::max(xdiff, (state1.beads[i].x-state2.beads[i].x).l2_norm() );
                vdiff=std::max(xdiff, (state1.beads[i].v-state2.beads[i].v).l2_norm() );
                fdiff=std::max(fdiff, (state1.beads[i].f-state2.beads[i].f).l2_norm() );
            }
            if(fdiff>10*maxDiff){ // Forces can be large, so allow more slack
                std::cerr<<"At time "<<state1.t<<" the f values have diverged.\n";

                write_world_state(std::cerr, state1);
                write_world_state(std::cerr, state2);

                return {false, "Force divergence."};
            }

            if(xdiff>maxDiff){
                std::cerr<<"At time "<<state1.t<<" the x values have diverged.\n";

                write_world_state(std::cerr, state1);
                write_world_state(std::cerr, state2);

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
