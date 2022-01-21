#include "dpd/engines/basic/basic_dpd_engine_v5_raw_orch.hpp"
#include "dpd/core/dpd_state_builder.hpp"
#include "dpd/core/dpd_state_io.hpp"

#include "dpd/engines/naive/naive_dpd_engine.hpp"

#include "dpd/core/struct_to_c.hpp"

#include "dpd/core/with_optional_gzip_stream.hpp"


#include <random>
#include <regex>
#include <cmath>
#include <vector>

struct BeadBondingInfo{
public:
    enum BeadBondingFlags : uint32_t{
        IsAngleBondCentre = 1,
        IsAngleBondHead = 2,
        IsAngleBondTail = 4,
        AngleBondMask = IsAngleBondCentre | IsAngleBondHead | IsAngleBondTail,
        IsHookeanHead = 8,
        IsHookeanTail = 16,
        HookeanMask = IsHookeanHead | IsHookeanTail
    };
private:
    unsigned nPolymerTypes;
    std::vector<BeadBondingFlags> m_interestingness;

    uint32_t get_polymer_offset_index(unsigned polymer_type, unsigned polymer_offset) const
    {
        return polymer_offset*nPolymerTypes+polymer_type;
    }

    uint32_t get_polymer_offset_index(const Bead &b) const
    {
        return get_polymer_offset_index(b.polymer_type, b.polymer_offset);
    }

    void add_polymer_offset_index_flags(unsigned polymer_type, unsigned polymer_offset, BeadBondingFlags flags)
    {
        unsigned index=get_polymer_offset_index(polymer_type, polymer_offset);
        if(index >= m_interestingness.size()){
            m_interestingness.resize(index+1, BeadBondingFlags(0));
        }
        m_interestingness[index] = BeadBondingFlags(m_interestingness[index] | flags);
    }

    void build_polymer_index_to_interest(const WorldState &state)
    {
        nPolymerTypes=state.polymer_types.size();

        for(const PolymerType &pt : state.polymer_types){
            for(unsigned i=0; i<pt.bead_types.size(); i++){
                add_polymer_offset_index_flags(pt.polymer_id, i, BeadBondingFlags(0));
            }

            for(const Bond &b : pt.bonds){
                add_polymer_offset_index_flags(pt.polymer_id,  b.bead_offset_head, BeadBondingFlags::IsHookeanHead);
                add_polymer_offset_index_flags(pt.polymer_id,  b.bead_offset_tail, BeadBondingFlags::IsHookeanTail);
            }
            for(const BondPair &bp : pt.bond_pairs){
                const auto head=pt.bonds.at(bp.bond_offset_head).bead_offset_head;
                const auto centre=pt.bonds.at(bp.bond_offset_head).bead_offset_tail;
                const auto tail=pt.bonds.at(bp.bond_offset_tail).bead_offset_tail;

                add_polymer_offset_index_flags(pt.polymer_id, head, BeadBondingFlags::IsAngleBondHead);
                add_polymer_offset_index_flags(pt.polymer_id, centre, BeadBondingFlags::IsAngleBondCentre);
                add_polymer_offset_index_flags(pt.polymer_id, tail, BeadBondingFlags::IsAngleBondTail);
            }
        }
    }
public:
    BeadBondingInfo(const WorldState &state)
    {
        build_polymer_index_to_interest(state);
    }


    BeadBondingFlags get_bonding_info(const Bead &b) const
    {
        unsigned index=get_polymer_offset_index(b);
        return m_interestingness.at(index);
    }
};

int main(int argc, const char *argv[])
{
    unsigned nSteps=1;
    if(argc>1){
        nSteps=std::atoi(argv[1]);
    }

    int line_no=0;
    WorldState state=read_world_state(std::cin, line_no);

    double bead_steps=nSteps*state.beads.size();

    bool do_check_beads = nSteps <= 1000;

    std::vector<Bead> check_beads;
    if(do_check_beads){
        WorldState state2=state;

        BasicDPDEngineV5Raw native;
        native.Attach(&state2);

        native.Run(nSteps);

        unsigned nCheckSumBeads=std::min<unsigned>(16u, state.beads.size());

        std::vector<unsigned> selected;

        if(nCheckSumBeads<16){
            for(const auto &b : state.beads){
                selected.push_back(b.bead_id);
            }
        }else{
            BeadBondingInfo bi(state);

            std::vector<std::pair<double,unsigned>> distances;
            std::vector<std::pair<double,unsigned>> speeds;
            std::vector<unsigned> angle_centre, hookean_only;
            for(const auto &b : state2.beads){
                distances.push_back({vec3_wrapped_distance(b.x, state.beads.at(b.bead_id).x, state.box), b.bead_id});
                speeds.push_back({b.v.l2_norm(), b.bead_id});

                auto bonding=bi.get_bonding_info(b);
                if( bonding & (BeadBondingInfo::IsAngleBondCentre) ){
                    angle_centre.push_back(b.bead_id);
                }
                if( (bonding & BeadBondingInfo::HookeanMask) & !(bonding & BeadBondingInfo::AngleBondMask)){
                    hookean_only.push_back(b.bead_id);
                }
            }

            std::sort(distances.begin(), distances.end());
            std::sort(speeds.begin(), speeds.end());

            for(unsigned i=0; i<4; i++){
                selected.push_back( (distances.end()-i-1)->second );
                selected.push_back( (speeds.end()-i-1)->second );
            }

            std::mt19937_64 urng;

            std::sample(angle_centre.begin(), angle_centre.end(), std::back_inserter(selected), 4, urng);
            std::sample(hookean_only.begin(), hookean_only.end(), std::back_inserter(selected), 4, urng);
            
            while(selected.size() < nCheckSumBeads){
                selected.push_back(urng() % state.beads.size());
            }

            for(unsigned i=0; i<selected.size(); i++){
                check_beads.push_back( state2.beads.at(selected[i]));
            }
        }

        double max_dist=0, sum_dist=0;
        for(unsigned i=0; i<check_beads.size(); i++){
            const auto &begin_pos=state.beads.at(check_beads[i].bead_id).x;
            const auto &end_pos=check_beads[i].x;
            double dist=vec3_wrapped_distance(begin_pos, end_pos, state.box);

            max_dist=std::max(max_dist, dist);
            sum_dist += dist;
        }
        fprintf(stderr, "Max check bead movement = %g\n", max_dist);
        fprintf(stderr, "Mean check bead movement = %g\n", sum_dist/check_beads.size());
    }

    BasicDPDEngineV5RawOrch engine;
    engine.Attach(&state);

    engine.write_xml(std::cout, nSteps, check_beads);
}
