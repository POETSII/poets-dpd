#include "dpd/tests/test_ordered_mesh.hpp"
#include "dpd/tests/test_bonded_bead_angle.hpp"
#include "dpd/tests/test_bonded_bead_triple.hpp"
#include "dpd/tests/test_bonded_bead_pair.hpp"
#include "dpd/tests/test_colliding_bead_pair_diss.hpp"
#include "dpd/tests/test_colliding_bead_pair.hpp"
#include "dpd/tests/test_isolated_moving_bead.hpp"
#include "dpd/tests/test_multiple_moving_beads.hpp"
#include "dpd/engines/naive/naive_dpd_engine.hpp"

#include "dpd/tests/test_runner.hpp"

#include "dpd/core/logging.hpp"
#include "dpd/core/logging_impl.hpp"

#define TBB_PREVIEW_GLOBAL_CONTROL 1
#include <tbb/global_control.h>

#include <random>

void usage()
{
    fprintf(stderr, "test_engine : engine-name [filter]\n");
    fprintf(stderr, "  engine names:\n");
    for(auto s : DPDEngineFactory::ListFactories()){
        fprintf(stderr, "    %s\n", s.c_str());
    }
    fprintf(stderr,"     PDPD_LOG=log-path : do full force logging to given file.\n");
    fprintf(stderr,"     PDPD_NUM_THREADS=n : Limit TBB to using n threads.\n");
    exit(1);
}

int main(int argc, const char *argv[])
{
    int max_parallelism=tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
    if(getenv("PDPD_NUM_THREADS")){
	    max_parallelism=std::atoi(getenv("PDPD_NUM_THREADS"));
    }
    std::cerr<<"TBB is using "<<max_parallelism<<" threads.\n";
    tbb::global_control tbb_control_threads(tbb::global_control::max_allowed_parallelism, max_parallelism);

    
    std::string engine_name;
    if(argc>1){
        engine_name=argv[1];
    }else{
        usage();
    }

    std::string prefix;
    if(argc>2){
        prefix=argv[2];
    }

    std::string log_dst;
    if(getenv("PDPD_LOG")){
        log_dst=getenv("PDPD_LOG");
        ForceLogging::set_logger(new FileLogger(log_dst));
        fprintf(stderr, "Logging to %s\n", log_dst.c_str());
    }

    TestOrderedMesh::register_tests();
    TestBondedBeadAngle::register_tests();
    TestBondedBeadTriple::register_tests();
    TestBondedBeadPair::register_tests();
    TestCollidingBeadPair::register_tests();
    TestCollidingBeadPairDiss::register_tests();
    TestIsolatedMovingBead::register_tests();
    TestMultipleMovingBeads::register_tests();

    std::shared_ptr<DPDEngine> pengine = DPDEngineFactory::create(engine_name);

    fprintf(stdout, "TAP version 13\n");

    auto tests=TestBase::get_test_factories_ordered_by_complexity();

    int test_num=1;
    int failed=0;
    for(const auto & ttt : tests){
        if(!prefix.empty() && prefix != std::get<0>(ttt).substr(0,prefix.size())){
            continue;
        }

        auto test=std::get<2>(ttt)();

        auto r=test_runner(*test, *pengine);
        if(r.first){
            fprintf(stdout, "ok %d - %s :  %s\n", test_num, test->name().c_str(), r.second.c_str());
        }else{
            fprintf(stdout, "not ok %d - %s :%s\n", test_num, test->name().c_str(), r.second.c_str());
            ++failed;
        }
        test_num++;
    }
    return failed>0;
}
