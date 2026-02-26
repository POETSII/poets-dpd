
#include "dpd/core/dpd_engine.hpp"
#include "dpd/core/dpd_state_io.hpp"
#include "dpd/core/dpd_state_validator.hpp"
#include "dpd/core/dpd_state_to_vtk.hpp"
#include <vector>
#include <utility>
#include <set>

#include <fstream>

void usage()
{
    fprintf(stderr, "world_state_to_snapshot_v2 : output-dir frac_bits smoothing state-input-files...\n");
    fprintf(stderr, "world_state_to_snapshot_v2 : output-dir frac_bits smoothing --vtk state-template vtk-input-files...\n");
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

int16_t make_diff(uint16_t curr, uint16_t target, uint16_t bound)
{
    int diff1=(int)target-(int)curr;
    if(diff1 > bound/2){
        diff1 -=bound;
    }
    if(diff1 < -bound/2){
        diff1 +=bound;
    }
    assert(INT16_MIN <= diff1 && diff1 <= INT16_MAX);
    return (int16_t)diff1;
}

uint16_t to_zig_zag(int16_t x)
{
    return (x << 1) ^ (x >> (16 - 1));
}


void make_delta(
    const std::array<uint16_t,3> &bounds,
    const std::vector<uint16_t> &curr,   // Where the beads are now
    const std::vector<uint16_t> &target, // Where we want them to be
    std::vector<int16_t> &delta,         // Delta to get there
    uint16_t *max_abs
){
    uint16_t m=0;
    for(unsigned i=0; i<curr.size(); i+=3){
        for(int d=0; d<3; d++){
            delta[i+d]=make_diff(curr[i+d], target[i+d], bounds[d]);
            m=std::max<uint16_t>(m, std::abs(delta[i+d]));
        }
    }
    if(max_abs){
        *max_abs=m;
    }
}

void make_delta_zigzag(
    const std::vector<int16_t> &delta, // Where we want them to be
    std::vector<uint16_t> &delta_zz
){
    for(unsigned i=0; i<delta.size(); i+=3){
        for(int d=0; d<3; d++){
            delta_zz[i+d]=to_zig_zag(delta[i+d]);
        }
    }
}



struct coded_form
{
    std::string type; // Goes in the "Type" field.
    std::vector<std::pair<std::string,int>> properties; // Any extra key-frame properties to go in
    std::vector<uint8_t> bytes; // Raw bytes to go in file
};

size_t calc_delta_zz_8_coded_size(
    const std::vector<uint16_t> &delta_zz
){
    size_t res=delta_zz.size();
    for(auto v : delta_zz){
        if(v >= 255){
            res += 2;
        }
    }
    return res;
}

coded_form make_delta_zz_8_coded_form(
    const std::vector<uint16_t> &delta_zz // The delta_zz compressed data
){
    coded_form res;
    res.type="delta_zz_8";
    res.bytes.reserve(delta_zz.size()*1.2);

    for(unsigned i=0; i<delta_zz.size(); i++){
        uint16_t val=delta_zz[i];
        if(val < 255){
            res.bytes.push_back(val);
        }else{
            res.bytes.push_back(255);
            res.bytes.push_back( val&0xFF );
            res.bytes.push_back( val>>8 );
        }
    }
    return res;
}

uint16_t smooth(
    uint16_t bound,
    uint16_t prev,
    uint16_t curr,
    uint16_t next
){
    double c=curr, p=prev, n=next;
    if(curr < bound/4){
        if( prev > bound-bound/4 ){
            c -= bound;
        }
        if( next > bound-bound/4 ){
            n -= bound;
        }
    }else if(curr > bound-bound/4){
        if( prev < bound/4 ){
            c += bound;
        }
        if( next < bound/4 ){
            n += bound;
        }
    }

    int x=round( 0.25*p + 0.5*c + 0.25*n );
    if(x < 0){
        x += bound;
    }else if(x >= bound){
        x -= bound;
    }
    return x;
}

std::vector<uint16_t> smooth(
    const std::array<uint16_t,3> &bounds,
    const std::vector<uint16_t> &prev,
    const std::vector<uint16_t> &curr,
    const std::vector<uint16_t> &next
){
    std::vector<uint16_t> res(curr.size());

    for(unsigned i=0; i<curr.size(); i++){
        res[i] = smooth(bounds[i%3], prev[i], curr[i], next[i]);
    }
    return res;
}

std::vector<std::pair<int,std::vector<uint16_t>>> smooth(
    const std::array<uint16_t,3> &bounds,
    std::vector<std::pair<int,std::vector<uint16_t>>> & slices
){
    std::vector<std::pair<int,std::vector<uint16_t>>> res;

    res.push_back(slices.front());
    for(unsigned i=1; i+1<slices.size(); i++){
        const auto &prev=slices[i-1].second;
        const auto &curr=slices[i].second;
        const auto &next=slices[i+1].second;
        res.push_back({ slices[i].first, smooth(bounds,prev,curr,next) });
    }
    return res;
}

int main(int argc, const char *argv[])
{
    try{
        if(argc<4){
            usage();
        }

        std::string dst_dir=argv[1];
        std::ofstream dst(dst_dir+"/dpd_snapshot.json");
        if(!dst.is_open()){
            std::cerr<<"Couldn't open "<<dst_dir<<"/dpd_snapshot.json\n";
            exit(1);
        }

        if(argc < 4){
            std::cerr<<"No files\n";
            exit(1);
        }

        int frac_bits=atoi(argv[2]);
        int smoothing=atoi(argv[3]);

        bool vtk=false;
        std::string prototype_name=argv[4];
        if(prototype_name=="--vtk"){
            std::cerr<<"vtk mode\n";
            vtk=true;
            prototype_name=argv[5];
        }

        WorldState prototype=read_world_state(prototype_name);

        int filter_bead_index=-1;
        for(const auto & b : prototype.bead_types){
            if(b.name=="W" || b.name=="w" || b.name=="Water"){
                filter_bead_index=b.id;
            }
        }

        std::map<unsigned,unsigned> source_beads;
        std::set<unsigned> source_polymers;
        for(unsigned i=0; i<prototype.beads.size(); i++){
            const auto &b=prototype.beads[i];
            if(b.bead_type!=filter_bead_index){
                source_beads.insert({ i, source_beads.size() });
                source_polymers.insert(prototype.beads[i].polymer_id);
            }
        }

        std::vector<std::pair<int,std::vector<uint16_t>>> slices;

        int max_dim=std::max(prototype.box[0], std::max(prototype.box[1], prototype.box[2]));
        double log2_max_dim=std::ceil(std::log2(max_dim));
        
        if(frac_bits + log2_max_dim > 16 ){
            fprintf(stderr, "This system is too large for this codec's assumptions (16-bit is enough for positions). frac_bits=%d, log2_max_dim=%g\n", frac_bits, log2_max_dim);
            exit(1);
        }
        
        double frac_bits_scale = ldexp(1, frac_bits);
        fprintf(stderr, "Frac_bits_scale=%g\n", frac_bits_scale);

        if(!vtk){
            for(int i=4; i<argc; i++){
                WorldState state1=read_world_state(argv[i]);
                if(state1.beads.size()!=prototype.beads.size()){
                    throw std::runtime_error("Bead count mis-match.");
                }
                validate(state1, 2.0);

                std::vector<uint16_t> positions(source_beads.size()*3);
                unsigned off=0;
                for(auto [index,dst] : source_beads){
                    const auto &bb=state1.beads[index];
                    positions[dst*3+0]=floor( bb.x[0] * frac_bits_scale );
                    positions[dst*3+1]=floor( bb.x[1] * frac_bits_scale );
                    positions[dst*3+2]=floor( bb.x[2] * frac_bits_scale );
                    //fprintf(stderr, "%d,%d,%d\n", positions[i*3])
                    ++off;
                }
                slices.push_back({state1.t, positions});
            }
        }else{
            // HACK  : to deal with old state files
            validate_no_check_hashes() = true;

            for(int i=6; i<argc; i++){
                WorldState state1=read_from_vtk(prototype, argv[i]);
                std::cerr<<"  loading "<<argv[i]<<"\n";
                if(state1.beads.size()!=prototype.beads.size()){
                    throw std::runtime_error("Bead count mis-match.");
                }
                state1.t=slices.size();

                std::vector<uint16_t> positions(source_beads.size()*3);
                unsigned off=0;
                for(auto [index,dst] : source_beads){
                    const auto &bb=state1.beads[index];
                    positions[dst*3+0]=floor( bb.x[0] * frac_bits_scale );
                    positions[dst*3+1]=floor( bb.x[1] * frac_bits_scale );
                    positions[dst*3+2]=floor( bb.x[2] * frac_bits_scale );
                    //fprintf(stderr, "%d,%d,%d\n", positions[i*3])
                    ++off;
                }
                slices.push_back({state1.t, positions});
            }
        }

        std::array<uint16_t,3> bounds{
            (uint16_t)(prototype.box[0] * frac_bits_scale), (uint16_t)(prototype.box[1] * frac_bits_scale), (uint16_t)(prototype.box[2] * frac_bits_scale)
        };


        std::sort(slices.begin(), slices.end());

        for(int i=0; i<smoothing; i++){
            slices=smooth(bounds, slices);
        }    

        auto q=[](const std::string &s) -> std::string
        { return "\""+s+"\""; };

        dst<<"{\"type\":\"dpd-snapshot-v1-indirect\", \"numSpecies\":"<<prototype.bead_types.size()<<", \"numBeads\":"<<source_beads.size()<<",\n";
        dst<<" \"numPolymers\":"<<source_polymers.size()<<",\n";
        dst<<" \"lower_bounds\":["<<prototype.origin[0]<<","<<prototype.origin[1]<<","<<prototype.origin[2]<<"],\n";
        dst<<" \"upper_bounds\":["<<prototype.origin[0]+prototype.box[0]<<","<<prototype.origin[1]+prototype.box[1]<<","<<prototype.origin[2]+prototype.box[2]<<"],\n";
        dst<<" \"bead_species\":[";
        bool first=true;
        for(auto [b,i] : source_beads){
            if(!first){
                dst<<",";
            }else{
                first=false;
            }
            dst<<(int)prototype.beads[b].bead_type;
        }
        dst<<"],\n";
        
        first=true;
        dst<<"  \"polymers\" : [\n";
        for(auto pid : source_polymers){
            if(!first){
                dst<<",";
            }else{
                first=false;
            }
            dst<<"[";
            const auto &polymer=prototype.polymers[pid];
            bool ffirst=true;
            for(const auto bid : polymer.bead_ids){
                if(source_beads.find(bid)==source_beads.end()){
                    continue;
                }
                if(ffirst){
                    ffirst=false;
                }else{
                    dst<<",";
                }
                dst<<source_beads[bid];
            }
            dst<<"]\n";
        }
        dst<<"],\n";

        dst<<"  \"frac_bits\":"<<frac_bits<<",\n";
        dst<<"  \"slices\":[\n";

        fprintf(stderr, "Box=%f, %f, %f\n", prototype.box[0], prototype.box[1], prototype.box[2]);
        fprintf(stderr, "Bounds=%u, %u, %u\n", bounds[0], bounds[1], bounds[2]);


        std::vector<int16_t> delta;
        std::vector<uint16_t> delta_zz;
        std::vector<uint16_t> prev;
        for(unsigned i=0; i<slices.size(); i++){
            if(i!=0){
                dst<<",";
            }

            dst<<"{ \"seq\":"<<i<<", \"t\":"<<i<<", \"frac_bits\":"<<frac_bits<<",";

            std::string path=dst_dir+"/data."+std::to_string(i)+".bin";
            std::ofstream bindst(path, std::ios_base::binary);

            const auto &sl=slices[i];
            const auto &target=sl.second;

            coded_form form;

            form.type="key16";
            form.bytes.assign( (const uint8_t*)&target[0], (const uint8_t*)&target[target.size()] );

            if(i!=0 && (i%60)){
                delta.resize(target.size());
                delta_zz.resize(target.size());

                uint16_t max_abs_delta=0;
                make_delta(bounds, prev, target, delta, &max_abs_delta);
                if(max_abs_delta < 128){
                    form.type="delta8";
                    form.bytes.assign( delta.begin(), delta.end() );
                }else{
                    make_delta_zigzag(delta, delta_zz);

                    auto pform=make_delta_zz_8_coded_form(delta_zz);
                    if(pform.bytes.size() < form.bytes.size()){
                        form=pform;
                    }
                }
            }
            dst<<q("type")<<":"<<q(form.type)<<",\n";
            for(const auto &kv : form.properties){
                dst<<q(kv.first)<<":"<<kv.second<<",\n";
            }

            bindst.write((char*)&form.bytes[0], form.bytes.size());

            std::cerr<<"  "<<i<<" : "<<form.type<<", "<<form.bytes.size()<<"\n";

            prev=target;
            
            dst<<"\"url\":\""<<"data."<<std::to_string(i)<<".bin\"";
            dst<<"}\n";
        }
        dst<<"  ]\n";
        dst<<"}\n";

    }catch(const std::exception &e){
        print_exception(e);
        exit(1);
    }

    return 0;
}
