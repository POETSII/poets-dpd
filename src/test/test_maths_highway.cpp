#include "dpd/maths/dpd_maths_highway.hpp"

#include "dpd/maths/dpd_maths_core_half_step_raw.hpp"

#include <random>

#include "dpd/core/vec3.hpp"

namespace it = dpd_maths_highway::HWY_NAMESPACE;
namespace hn = it::hn;

std::mt19937_64 rng;
std::uniform_real_distribution<> udist;
double urng(){ return udist(rng); };

using dpd_maths_highway::packed_bead;
using it::soa_packed_beads;
using it::vector_packed_beads;

packed_bead random_bead (vec3i_t origin)
{
    packed_bead res;
    res.hash=rng() & 0xFFFFFFFFul;
    for(int d=0; d<3; d++){
        res.x[d] = floor(origin[d]+urng());
        res.v[d] = urng() * 8 - 4;
        res.f[d] = urng() * 8 - 4;
    }
    return res;
};

void check_soa_vector()
{
    std::vector<packed_bead> ref;
    soa_packed_beads packed;

    vec3i_t origin{1, 2, 3};

    for(int i=0; i<1000; i++){
        auto bit=rng()%2;
        if( (bit || ref.empty()) && ref.size() < soa_packed_beads::MAX_BEADS ){
            // add
            auto b= random_bead(origin);
            ref.push_back(b);
            packed.insert(b);
        }
        if( (!bit || ref.size()==soa_packed_beads::MAX_BEADS) && !ref.empty() ){
            unsigned index=rng() % ref.size();
            auto rb=ref[index];
            auto gb=packed.erase(index);

            fprintf(stderr, "Erase %u of %u\n", (unsigned)index, (unsigned)ref.size());

            if(memcmp(&rb,&gb,sizeof(rb))){
                throw std::runtime_error("Bead mismatch.");
            }

            if(index+1 < ref.size()){
                std::swap(ref[index], ref.back());
            }
            ref.resize(ref.size()-1);
        }

        for(unsigned j=0; j<ref.size(); j++){
            auto rb=ref[j];
            auto gb=packed.get_bead(j);
            if(memcmp(&rb,&gb,sizeof(rb))){
                throw std::runtime_error("Bead mismatch.");
            }
        }
    }
}

void check_packed_vector()
{
    std::vector<vec3i_t> origins;
    std::vector<soa_packed_beads*> cells;
    std::vector<packed_bead> ref;

    origins.resize(14);
    cells.resize(14);
    for(unsigned i=0; i<cells.size(); i++){
        cells[i]=new soa_packed_beads();
        origins[i]=vec3i_t(i/9,(i/3)%3, i%3);
    }

    vector_packed_beads vpb;

    for(unsigned i=0; i<cells.size(); i++){
        cells[i]->clear();

        unsigned n=rng()%(soa_packed_beads::MAX_BEADS+1);
        for(unsigned j=0; j<n; j++){
            auto b=random_bead(origins[i]);
            ref.push_back(b);
            cells[i]->insert(b);

            auto rb=ref.back();
            auto gb=cells[i]->get_bead(j);
            if(memcmp(&rb,&gb,sizeof(rb))){
                throw std::runtime_error("Bead mismatch.");
            }
        }
    }

    unsigned offset=0;
    for(auto c : cells){
        for(unsigned j=0; j<c->n; j++){
            auto rb=ref[offset];
            auto gb=c->get_bead(j);
            if(memcmp(&rb,&gb,sizeof(rb))){
                throw std::runtime_error("Bead mismatch.");
            }
            ++offset;
        }
    }

    vpb.gather_beads(cells.size(), &cells[0]);

    for(unsigned i=0; i<ref.size(); i++){
        auto rb=ref[i];
        auto gb=vpb.get_bead(i);
        if(memcmp(&rb,&gb,sizeof(rb))){
            throw std::runtime_error("Bead mismatch.");
        }
    }

    for(unsigned i=0; i<ref.size(); i++){
        for(int d=0; d<3; d++){
            ref[i].f[d] += 1+d;
            vpb.f[d][i] += 1+d;
        }
    }

    vpb.scatter_forces(cells.size(), &cells[0]);

    offset=0;
    for(auto c : cells){
        for(unsigned j=0; j<c->n; j++){
            auto rb=ref[offset];
            auto gb=c->get_bead(j);
            if(memcmp(&rb,&gb,sizeof(rb))){
                throw std::runtime_error("Bead mismatch.");
            }
            ++offset;
        }
    }

    offset=0;
    for(auto &c : cells)
    {
        unsigned ci=&c - &cells[0];
        float dt=0.1;
        std::vector<uint8_t> outgoing;
        float origin[3];
        origins[ci].extract(origin);
        c->newton(origin, dt, outgoing);

        for(unsigned j=0; j<c->n; j++){
            assert( c->get_f_vec(0)[j] == 0 );
            auto bnow=c->get_bead(j);

            auto bref=ref[offset];
            fprintf(stderr, "+hash=%u,x[0]=%f,v[0]=%f,f[0]=%f\n", bref.hash, bref.x[0], bref.v[0], bref.f[0]);

            dpd_maths_core_half_step_raw::update_mom(
                dt, bref
            );
            dpd_maths_core_half_step_raw::update_pos_no_wrap(
                dt, bref
            );
            fprintf(stderr, " .   x[0]=%f,v[0]=%f,f[0]=%f\n", bref.x[0], bref.v[0], bref.f[0]);

            for(int d=0; d<3; d++){
                assert(bnow.f[d]==bref.f[d]);
                assert(fabs( bnow.x[d] - bref.x[d] ) < 0.001 );
                assert(fabs( bnow.v[d] - bref.v[d] ) < 0.001 );
            }

            ++offset;
        }
    }

    for(unsigned i=0; i<cells.size(); i++){
        delete cells[i];
    }
}

int main()
{
    fprintf(stderr, "Here\n");
    
    check_soa_vector();

    for(int i=0;i<100; i++){
        check_packed_vector();
    }

    fprintf(stderr, "Done\n");
}