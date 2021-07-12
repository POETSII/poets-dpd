#include "tests/test_ordered_mesh.hpp"
#include "tests/test_bonded_bead_angle.hpp"
#include "tests/test_bonded_bead_triple.hpp"
#include "tests/test_bonded_bead_pair.hpp"
#include "tests/test_colliding_bead_pair_diss.hpp"
#include "tests/test_colliding_bead_pair.hpp"
#include "tests/test_isolated_moving_bead.hpp"
#include "tests/test_multiple_moving_beads.hpp"
#include "naive_dpd_engine.hpp"

#include "tests/test_runner.hpp"

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

    NaiveDPDEngine<true> engine;

    fprintf(stdout, "TAP version 13\n");

    auto tests=TestBase::get_test_factories_ordered_by_complexity();

    int test_num=1;
    for(const auto & ttt : tests){
        if(!prefix.empty() && prefix != std::get<0>(ttt).substr(0,prefix.size())){
            continue;
        }

        auto test=std::get<2>(ttt)();

        auto r=test_runner(*test, engine);
        if(r.first){
            fprintf(stdout, "ok %d - %s\n", test_num, test->name().c_str());
        }else{
            fprintf(stdout, "not ok %d - %s :%s\n", test_num, test->name().c_str(), r.second.c_str());
        }
        test_num++;
    }
}
