#ifndef test_runner_hpp
#define test_runner_hpp

#include "test_base.hpp"
#include "dpd/core/dpd_engine.hpp"
#include "dpd/core/dpd_state_io.hpp"

std::pair<bool,std::string> test_runner(
    TestBase &test,
    DPDEngine &engine
){
    WorldState state;

    try{
        std::string name=test.name();

        state = test.create_world();
        validate(state);

        std::string no_support=engine.CanSupport(&state);
        if(!no_support.empty()){
            return {true, "SKIP : "+no_support};
        }

        engine.Attach(&state);

        unsigned steps_done=0;

        while(1){
            //write_world_state(std::cout, state);
            //std::cout<<"\n";

            if(steps_done > 100000){
                return {false, "More than 100000 steps done without test finishing."};
            }
            unsigned todo=test.get_advance_count(state);
            if(todo > 100000){
                return {false, "Test asked to do more than 100000 steps."};
            }
            if(todo==0){
                return {true, {}};
            }

            engine.Run(todo);
            validate(state);

            steps_done += todo;
        }
    }catch(TestFailedException &e){
        write_world_state(std::cerr, state);
        return {false, e.what()};
    }catch(ValidationFailedError &e){
        write_world_state(std::cerr, state);
        return {false, e.what()};
    }catch(std::exception &e){
        return {false, std::string("Exception : ")+e.what()};
    }catch(...){
        return {false, "Unknown exception occured."};
    }
}

#endif
