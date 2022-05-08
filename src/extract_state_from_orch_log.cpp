#include "dpd/core/dpd_state.hpp"
#include "dpd/core/dpd_state_io.hpp"
#include <regex>
#include <fstream>
#include <iostream>

#include "dpd/maths/dpd_maths_core_half_step.hpp"

std::regex reBegin(R"(^.*a8f23f7c94100804_BEGIN\[t=([0-9a-fA-F]+),numBeads=([0-9a-fA-F]+)\].*$)");
std::regex reBead(R"(^.*a8f23f7c94100804_BEAD\[id=([0-9a-fA-F]+),t=([0-9a-fA-F]+),x=(\[[^\]]+\]),v=(\[[^\]]+\]),f=(\[[^\]]+\])\].*$)");
std::regex reEnd(R"(^.*a8f23f7c94100804_END\[checksum=([0-9a-fA-F]+)\].*$)");

void require(bool cond, const char *msg)
{
    if(!cond){
        throw std::runtime_error(msg);
    }
}

vec3r_t raw_to_vec(const std::string &x)
{
    char *p=const_cast<char*>(x.c_str());

    require(*p=='[', "Raw vector did not start with [");
    ++p;

    uint32_t rx[3];

    rx[0]=std::strtoul(p, &p, 16);
    require(*p==',', "Raw vector missing first comma.");
    ++p;

    rx[1]=std::strtoul(p, &p, 16);
    require(*p==',', "Raw vector missing second comma.");
    ++p;

    rx[2]=std::strtoul(p, &p, 16);
    require(*p==']', "Raw vector does not end with ].");
    
    vec3r_t res;
    for(int i=0; i<3; i++){
        union{ uint32_t i; float f; } u;
        u.i=rx[i];
        res[i]=u.f;
    }
    return res;
}

int main(int argc, char *argv[])
{
    std::string world_state_path=argv[1];
    
    std::cerr<<"Reading world state from "<<world_state_path<<"\n";
    int line_no=0;
    std::ifstream world_state_in(world_state_path);
    if(!world_state_in.is_open()){
        throw std::runtime_error("Could not open "+world_state_path);
    }
    WorldState state=read_world_state(world_state_in, line_no);

    std::unordered_map<BeadHash,Bead> received;

    std::istream &in=std::cin;

    unsigned numBeads, time;

    std::cerr<<"Reading log from stdin\n";
    std::string line;
    while(1){
        if(!std::getline(in, line)){
            throw std::runtime_error("Didn't find header.");
        }

        std::smatch m;
        if(std::regex_match(line, m, reBegin)){
            time=std::stoi(m[1],0, 16);
            numBeads=std::stoi(m[2], 0, 16);
            break;
        }
    }

    if(numBeads!=state.beads.size()){
        throw std::runtime_error("Num beads in state does not match num beads in log.");
    }
    if(time<state.t){
        throw std::runtime_error("Bead batch is from a time "+std::to_string(time)+" earlier than state "+std::to_string(state.t)+".");
    }

    std::cerr<<"FOund header, numBeads="<<numBeads<<"=0x"<<std::hex<<numBeads<<std::dec<<"\n";

    uint32_t checksum;
    while(1){
        if(!std::getline(in, line)){
            throw std::runtime_error("Didn't find footer.");
        }

        std::smatch m;
        if(std::regex_match(line, m, reEnd)){
            std::cerr<<"Got end\n";
            checksum=std::stoul(m[1], 0, 16);
            break;
        }else if(std::regex_match(line, m, reBead)){
            BeadHash id{std::stoul(m[1], 0, 16)};
            uint32_t t=std::stoul(m[2], 0, 16);
            vec3r_t x=raw_to_vec(m[3]);
            vec3r_t v=raw_to_vec(m[4]);
            vec3r_t f=raw_to_vec(m[5]);
            BeadHash bh{id};

            if(t!=time){
                throw std::runtime_error("Bead time does not match reported batch time from header.");
            }

            Bead b;
            b.set_hash_code(id);
            b.x=x;
            b.v=v;
            b.f=f;
            
            auto it=received.insert({id,b});
            if(!it.second){
                throw std::runtime_error("Duplicated bead.");
            }
        }else{
            std::cerr<<"Line = "<<line<<"\n";
        }
    }

    std::cerr<<"received.size()="<<received.size()<<", numBeads="<<numBeads<<"\n";
    if(received.size()!=numBeads){
        throw std::runtime_error("Not all beads received before end tag (out of order printing?)");
    }

    std::cerr<<"Applying half-step mom correction.\n";
    for(Bead &b : state.beads){

        auto it=received.find(b.get_hash_code());
        if(it==received.end()){
            throw std::runtime_error("Missing bead output.");
        }

        b.x=it->second.x;
        b.v=it->second.v;
        b.f=it->second.f;

        dpd_maths_core_half_step::update_mom<float>((float)state.dt, b);
    }

    write_world_state(std::cout, state);
}