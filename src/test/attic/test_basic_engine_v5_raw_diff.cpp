#include "dpd/tests/test_ordered_mesh.hpp"
#include "dpd/tests/test_bonded_bead_angle.hpp"
#include "dpd/tests/test_bonded_bead_triple.hpp"
#include "dpd/tests/test_bonded_bead_pair.hpp"
#include "dpd/tests/test_colliding_bead_pair_diss.hpp"
#include "dpd/tests/test_colliding_bead_pair.hpp"
#include "dpd/tests/test_isolated_moving_bead.hpp"
#include "dpd/tests/test_multiple_moving_beads.hpp"
#include "dpd/engines/naive/naive_dpd_engine.hpp"
#include "dpd/engines/basic/basic_dpd_engine_v5_raw.hpp"

#include "dpd/tests/test_differential.hpp"

#include <random>

int main(int argc, const char *argv[])
{
    std::string prefix;
    if(argc>1){
        prefix=argv[1];
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
    BasicDPDEngineV5Raw engine2;

    fprintf(stdout, "TAP version 13\n");

    auto tests=TestBase::get_test_factories_ordered_by_complexity();

    int test_num=1;
    int failed=0;
    for(const auto & ttt : tests){
        if(!prefix.empty() && prefix != std::get<0>(ttt).substr(0,prefix.size())){
            continue;
        }

        auto test=std::get<2>(ttt)();

        auto r=test_differential(*test, engine1, engine2);
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
