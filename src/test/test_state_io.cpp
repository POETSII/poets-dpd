#include "dpd/core/dpd_state.hpp"
#include "dpd/core/dpd_state_io.hpp"

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


#include <sstream>

class RoundTripFailedException
    : public std::runtime_error
{
public:
    RoundTripFailedException(const std::string &msg)
        : std::runtime_error(msg)
    {}
};

template<class T>
void require_same(const T &a, const T &b, const std::string &part)
{
    if(a!=b){
        throw RoundTripFailedException("Mismatch on "+part);
    }
}

void require_close(const double &a, const double &b, double tol, const std::string &part)
{
    if( std::abs(a-b) > tol){
        throw RoundTripFailedException("Mismatch on "+part);
    }
}

void check_diff(const WorldState &ref, const WorldState &got)
{
#define SAME(member) \
    require_same(ref. member, got. member, #member)
#define CLOSE(member) \
    require_close(ref. member, got. member, 1e-6, #member)

    SAME(origin);
    SAME(box);
    CLOSE(lambda);
    SAME(t);
    CLOSE(dt);
    SAME(seed);

    SAME(interactions.size());
    for(unsigned i=0; i<ref.interactions.size(); i++){
        CLOSE(interactions[i].conservative);
        CLOSE(interactions[i].dissipative);
    }

    SAME(bead_types.size());
    for(unsigned i=0; i<ref.bead_types.size(); i++){
        SAME(bead_types[i].id);
        SAME(bead_types[i].name);
        CLOSE(bead_types[i].r);
    }

    SAME(polymer_types.size());
    for(unsigned p=0; p<ref.polymer_types.size(); p++){
        SAME(polymer_types[p].bead_types);
        SAME(polymer_types[p].bonds.size());
        for(unsigned i=0; i<ref.polymer_types[p].bonds.size(); i++){
            SAME(polymer_types[p].bonds[i].bead_offset_head);
            SAME(polymer_types[p].bonds[i].bead_offset_tail);
            CLOSE(polymer_types[p].bonds[i].kappa);
            CLOSE(polymer_types[p].bonds[i].r0);
        }

        SAME(polymer_types[p].bond_pairs.size());
        for(unsigned i=0; i<ref.polymer_types[p].bond_pairs.size(); i++){
            SAME(polymer_types[p].bond_pairs[i].bond_offset_head);
            SAME(polymer_types[p].bond_pairs[i].bond_offset_tail);
            CLOSE(polymer_types[p].bond_pairs[i].kappa);
            CLOSE(polymer_types[p].bond_pairs[i].theta0);
        }
    }

    SAME(beads.size());
    for(unsigned i=0; i<ref.beads.size(); i++){
        SAME(beads[i].bead_id);
        SAME(beads[i].polymer_id);
        SAME(beads[i].polymer_offset);
        SAME(beads[i].bead_type);
        SAME(beads[i].polymer_type);
        SAME(beads[i].is_monomer);
        
        auto dx=ref.beads[i].x - got.beads[i].x;
        double dxm=vec3_wrapped_distance(ref.beads[i].x, got.beads[i].x, ref.box);
        if(dxm > 1e-4){
            std::stringstream tmp;
            tmp<<"x component: xref="<<ref.beads[i].x<<", xgot="<<got.beads[i].x<<", dist="<<dxm;
            throw RoundTripFailedException(tmp.str());
        }

        // Check the velocity is sufficiently close
        vec3r_t rv=ref.beads[i].v, gv=got.beads[i].v;
        auto rvm=rv.l2_norm(), gvm=gv.l2_norm();
        if(rvm==0){
            if(gvm!=0){
                throw RoundTripFailedException("v of 0 does notmatvch.");
            }
        }else{
            double max_mag=-DBL_MAX;
            for(int d=0; d<3; d++){
                max_mag=std::max(max_mag, std::abs(rv[d]));
            }

            for(int d=0; d<3; d++){
                bool ok=false;

                double delta=std::abs(rv[d]-gv[d]);
                if( (delta / max_mag) >= (1.0/64) && delta >= ldexp(1.0,-14)  ){
                    std::stringstream tmp;
                    tmp<<"c component: vref="<<rv<<", vgot="<<gv<<", max_mag="<<max_mag<<", delta["<<d<<"]="<<delta<<", delta/max_mag="<<(delta / max_mag);
                    throw RoundTripFailedException(tmp.str());
                }
            }
        }

        // Check the force is sufficiently close
        vec3r_t rf=ref.beads[i].f, gf=got.beads[i].f;
        auto rfm=rf.l2_norm(), gfm=gf.l2_norm();
        if(rfm==0){
            if(gfm!=0){
                throw RoundTripFailedException("f of 0 does notmatvch.");
            }
        }else{
            double max_mag=-DBL_MAX;
            for(int d=0; d<3; d++){
                max_mag=std::max(max_mag, std::abs(rf[d]));
            }

            for(int d=0; d<3; d++){
                bool ok=false;

                double delta=std::abs(rf[d]-gf[d]);
                if( (delta / max_mag) >= (1.0/64) && delta >= ldexp(1.0,-14)  ){
                    std::stringstream tmp;
                    tmp<<"f component: fref="<<rf<<", fgot="<<gf<<", max_mag="<<max_mag<<", delta["<<d<<"]="<<delta<<", delta/max_mag="<<(delta / max_mag);
                    throw RoundTripFailedException(tmp.str());
                }
            }
        }
    }
}

void check_round_trip(const WorldState &orig)
{
    {
        std::stringstream tmp;
        write_world_state(tmp, orig, false);
        int line_no=0;
        WorldState got=read_world_state(tmp, line_no);

        check_diff(orig, got);
    }

    {
        std::stringstream tmp;
        write_world_state(tmp, orig, true);
        int line_no=0;
        WorldState got=read_world_state(tmp, line_no);

        check_diff(orig, got);
    }
}

std::pair<bool,std::string> test_runner_intermediate_states(
    TestBase &test,
    DPDEngine &engine
){
    WorldState state;

    try{
        std::string name=test.name();

        state = test.create_world();
        validate(state);

        engine.Attach(&state);

        unsigned steps_done=0;

        while(1){
            check_round_trip(state);

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
    }catch(const RoundTripFailedException &e){
        return {false, e.what()};
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

int main()
{
    TestOrderedMesh::register_tests();
    TestBondedBeadAngle::register_tests();
    TestBondedBeadTriple::register_tests();
    TestBondedBeadPair::register_tests();
    TestCollidingBeadPair::register_tests();
    TestCollidingBeadPairDiss::register_tests();
    TestIsolatedMovingBead::register_tests();
    TestMultipleMovingBeads::register_tests();

    std::shared_ptr<DPDEngine> pengine = std::make_shared<NaiveDPDEngine<>>();
    std::string prefix;

    fprintf(stdout, "TAP version 13\n");

    auto tests=TestBase::get_test_factories_ordered_by_complexity();
    int test_num=1;
    int failed=0;
    for(const auto & ttt : tests){
        if(!prefix.empty() && prefix != std::get<0>(ttt).substr(0,prefix.size())){
            continue;
        }

        auto test=std::get<2>(ttt)();

        auto r=test_runner_intermediate_states(*test, *pengine);
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