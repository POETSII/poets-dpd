#include "tests/test_ordered_mesh.hpp"
#include "tests/test_bonded_bead_angle.hpp"
#include "tests/test_bonded_bead_triple.hpp"
#include "tests/test_bonded_bead_pair.hpp"
#include "tests/test_colliding_bead_pair_diss.hpp"
#include "tests/test_colliding_bead_pair.hpp"
#include "tests/test_isolated_moving_bead.hpp"
#include "tests/test_multiple_moving_beads.hpp"
#include "naive_dpd_engine.hpp"

#include "tests/test_differential.hpp"

#include <random>

void usage()
{
    fprintf(stderr, "test_dpd_engine : engine-name [filter]\n");
    fprintf(stderr, "  engine names:\n");
    for(auto s : DPDEngineFactory::ListFactories()){
        fprintf(stderr, "    %s\n", s.c_str());
    }
    exit(1);
}

int main(int argc, const char *argv[])
{
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

    TestOrderedMesh::register_tests();
    TestBondedBeadAngle::register_tests();
    TestBondedBeadTriple::register_tests();
    TestBondedBeadPair::register_tests();
    TestCollidingBeadPair::register_tests();
    TestCollidingBeadPairDiss::register_tests();
    TestIsolatedMovingBead::register_tests();
    TestMultipleMovingBeads::register_tests();

    NaiveDPDEngine<false> engine1;
    std::shared_ptr<DPDEngine> pengine2 = DPDEngineFactory::create(engine_name);

    fprintf(stdout, "TAP version 13\n");

    auto tests=TestBase::get_test_factories_ordered_by_complexity();

    int test_num=1;
    int failed=0;
    for(const auto & ttt : tests){
        if(!prefix.empty() && prefix != std::get<0>(ttt).substr(0,prefix.size())){
            continue;
        }

        auto test=std::get<2>(ttt)();

        auto r=test_differential(*test, engine1, *pengine2);
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
