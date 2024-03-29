
#include "dpd/core/dpd_engine.hpp"

#include "dpd/core/dpd_state_io.hpp"
#include "dpd/core/dpd_state_validator.hpp"

#include <fstream>

void usage()
{
    fprintf(stderr, "change_world : [--dt dt] [--t t] [--rng-seed seed] src-file output-file \n");
    exit(1);
}

void print_exception(const std::exception& e, int level =  0)
{
    std::cerr << std::string(level, ' ') << "exception: " << e.what() << '\n';
    try {
        std::rethrow_if_nested(e);
    } catch(const std::exception& e) {
        print_exception(e, level+1);
    } catch(...) {}
}
 

double now()
{
    timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9*ts.tv_nsec;
}


int main(int argc, const char *argv[])
{
    try{
        double tStart=now();

        std::vector<std::pair<std::string,std::string>> changes;
        std::vector<std::string> left_args;

        int i=1;
        while(i<argc){
            std::string key=argv[i];
            ++i;
            if(key.size()>2 && key.substr(0,2)=="--"){
                if(i==argc){
                    std::cerr<<"Missing parameter to "<<key<<"\n";
                    exit(1);
                }
                std::string value=argv[i];
                ++i;
                changes.push_back({key,value});
            }else{
                left_args.push_back(key);
            }
        }

        if(left_args.size()!=2){
            std::cerr<<"Expecting two positional arguments.\n";
            std::cerr<<"Leftovers = \n";
            for(int i=0; i<left_args.size(); i++){
                std::cerr<<"  "<<left_args[i]<<"\n";
            }
            std::cerr<<"Key Vals\n";
            for(auto kv : changes){
                std::cerr<<"  "<<kv.first<<", "<<kv.second<<"\n";
            }
            usage();
        }

        std::string srcFile=left_args[0];
        std::string dstFile=left_args[1];

        WorldState state=read_world_state(srcFile);

        validate(state, 100);

        for(auto kv : changes){
            if(kv.first=="--dt"){
                double dt=std::stod(kv.second);
                state.dt=dt;
                fprintf(stderr, "Setting dt to %g\n", dt);
            }else if(kv.first=="--t"){
                int t=std::stoi(kv.second);
                state.t=t;
                fprintf(stderr, "Setting t to %d\n", t);
            }else if(kv.first=="--rng-seed"){
                unsigned seed=std::abs(std::stoi(kv.second));
                state.seed=seed;
                fprintf(stderr, "Setting RNG seed to %u\n", seed);
            }else{
                std::cerr<<"Didn't understand option "<<kv.first<<"\n";
                exit(1);
            }
        }

        write_world_state(dstFile, state);

        fprintf(stderr, "Successfully changed properties.\n");

    }catch(const std::exception &e){
        print_exception(e);
        exit(1);
    }

    return 0;
}
