#ifndef cache_info_hpp
#define cache_info_hpp

#include <vector>
#include <filesystem>
#include <cstdint>
#include <sys/types.h>
#include <pthread.h>
#include <fstream>
#include <iostream>
#include <set>
#include <cassert>

struct cpu_set_obj
{
    cpu_set_t set;

    cpu_set_obj()
    {
        CPU_ZERO(&set);
    }

    cpu_set_obj(const cpu_set_t &o)
        :set(o)
    {}

    bool empty() const
    { return CPU_COUNT(&set)==0; }

    unsigned size() const
    { return CPU_COUNT(&set); }

    bool operator[](unsigned index) const
    { return CPU_ISSET(index, &set); }

    struct proxy
    {
        cpu_set_t &dst;
        unsigned index;

        void operator=(bool x)
        {
            if(x){
                CPU_SET(index, &dst);
            }else{
                CPU_CLR(index, &dst);
            }
        }
    };

    proxy operator[](unsigned index)
    { return proxy{set, index}; }

    bool operator==(const cpu_set_obj &o) const
    { return CPU_EQUAL(&set, &o.set); }

    bool operator<(const cpu_set_obj &o) const
    {
        unsigned i=0;
        unsigned todo_left=size(), todo_right=o.size();
        while(1){
            if(todo_left==0 && todo_right!=0 ) return true;
            if(todo_left!=0 && todo_right==0 ) return false;
            if(todo_left==0 && todo_right==0 ) return false;
            bool left=(*this)[i], right=o[i];
            if(!left && right) return true;
            if(left && !right) return false;
            todo_left -= left;
            todo_right -= right;
            ++i;
        }
        assert(false);
        throw std::logic_error("Invariants failed");
    }

    bool is_disjoint(const cpu_set_obj &o) const
    {
        cpu_set_t i;
        CPU_AND(&i, &set, &o.set);
        return CPU_COUNT(&i)==0;
    }
};

    
inline std::ostream &operator<<(std::ostream &dst, const cpu_set_obj &s)
{
    unsigned bits=s.size();
    unsigned i=0;
    bool first=true;
    dst<<"[";
    int begin_range = -1;
    int end_range=-1;

    auto flush=[&]()
    {
        if(!first){
            dst<<",";
        }
        if(begin_range!=end_range){
            dst<<begin_range<<"-"<<end_range;
        }else{
            dst<<begin_range;
        }
        begin_range=-1;
        end_range=-1;
        first=false;
    };

    while(bits){
        if(s[i]){
            if(begin_range==-1){
                begin_range=i;
            }
            end_range=i;
            bits -= 1;
        }else{
            if(begin_range!=-1){
                flush();
            }
        }
        ++i;
    }
    if(begin_range!=-1){
        flush();
    }
    dst<<"]";
    return dst;
}

namespace cache_info
{


    enum CacheType : uint32_t
    {
        Invalid = 0,
        Data = 1 ,
        Instruction = 2,
        Unified = Data | Instruction,
        NUMA = 4
    };

    struct cache_cluster
    {
        int level; // level of heirarchy. 1=L1, 2=L2, ... Not valid for NUMA
        CacheType type;
        cpu_set_obj cpu_set;
        size_t size;

        bool same_level(const cache_cluster &o) const
        { return level==o.level && type==o.type; }

        bool operator<(const cache_cluster &o) const
        {
            if(level < o.level) return true;
            if(level > o.level) return false;
            if(type < o.type) return true;
            if(type > o.type) return false;
            return cpu_set < o.cpu_set;
        }

        bool operator==(const cache_cluster &o) const
        {
            if(level!=o.level) return false;
            if(type!=o.type) return false;
            if(cpu_set.is_disjoint(o.cpu_set)) return false;
            if( !(cpu_set == o.cpu_set)){
                throw std::runtime_error("Inconistent CPU sets.");
            }
            if(size!=o.size){
                throw std::runtime_error("Inconsistent CPU memory sizes.");
            }
            return true;
        }
    };

    struct cache_level
    {
        int level; // level of heirarchy. 1=L1, 2=L2, ... Not valid for NUMA
        CacheType type;
        size_t node_size;
        size_t total_size;

        std::vector<cpu_set_obj> nodes;
    };

    struct system_cache_info
    {
        std::vector<cache_level> levels;
        int llc_index = -1;

        size_t get_llc_node_size() const
        { return levels.at(llc_index).node_size; }

        size_t get_llc_total_size() const
        { return levels.at(llc_index).total_size; }

        const std::vector<cpu_set_obj> &get_llc_nodes() const
        { return levels.at(llc_index).nodes; }
    };

    auto read_line( const std::filesystem::path &p ) -> std::string
    {
        std::ifstream src(p);
        if(!src.is_open()){
            throw std::runtime_error("Couldn't open path "+p.string()+" for reading.");
        }
        std::string res;
        std::getline(src, res);
        return res;
    };

    auto read_lines( const std::filesystem::path &p ) -> std::vector<std::string>
    {
        std::ifstream src(p);
        if(!src.is_open()){
            throw std::runtime_error("Couldn't open path "+p.string()+" for reading.");
        }
        std::vector<std::string> res;
        std::string s;
        while(std::getline(src, s)){
            res.push_back(s);
        }
        return res;
    };

    auto read_node_meminfo_int( const std::filesystem::path &p, const std::string &key ) -> int64_t
    {
        auto lines=read_lines(p);
        for(auto s : lines){
            if(s.find(key)){
                std::stringstream ss(s);
                std::string lit_node, index, size, postfix;
                ss >> lit_node >> index >> size >> postfix;
                if(lit_node!="Node"){ throw std::runtime_error("Expecting literal Node"); }
                if(postfix!="kB"){ throw std::runtime_error("Expecting postfix kB"); }
                return std::stoi(size) * 1024;
            }
        }
        throw std::runtime_error("Couldn't find key "+key);
    };

    auto read_cache_int( const std::filesystem::path &p ) -> int64_t
    {
        std::string tmp=read_line(p);
        size_t last=0;
        int64_t res=std::stoll(tmp, &last);
        if(last < tmp.size()){
            if(tmp[last]=='K'){
                res*=1024;
            }else if(tmp[last]=='M'){
                res*=1024*1024;
            }else{
                throw std::runtime_error("Don't know how to parse suffix on "+tmp);
            }
        }
        return res;
    };

    cache_cluster read_cpu_cache_index(const std::filesystem::path &cache_base)
    {
        cache_cluster res;

        std::string type=read_line(cache_base/"type");
        if(type=="Instruction"){
            res.type=CacheType::Instruction;
        }else if(type=="Data"){
            res.type=CacheType::Data;
        }else if(type=="Unified"){
            res.type=CacheType::Unified;
        }else{
            throw std::runtime_error("Didn't understand cache index type string "+type);
        }

        res.level=read_cache_int(cache_base/"level");
        res.size=read_cache_int(cache_base/"size");

        std::string map=read_line(cache_base/"shared_cpu_map");
        CPU_ZERO(&res.cpu_set);
        for(unsigned i=0; i<map.size(); i++){
            char ch=map[i];
            if(!isxdigit(ch)){
                throw std::runtime_error("Unknown digit in cpu map "+map);
            }
            int digit=std::stoi(map.substr(i,1), 0, 16);
            for(int j=0; j<4; j++){
                if( (digit>>j)&1 ){
                    res.cpu_set[4*i+j] = true;
                }
            }
        }
        std::cerr<<"  "<<res.level<<"\n";

        return res;
    }

    cpu_set_t read_cpu_list(const std::string &s)
    {
        cpu_set_t res;
        CPU_ZERO(&res);
        int pos=0;
        int range_begin=-1; // Beginning of current range (if we are in a range)
        int number=0;  // Number we are currently building
        for(unsigned pos=0; pos<s.size(); pos++){
            char ch=s[pos];
            if(ch==','){
                if(range_begin==-1){
                    CPU_SET(number, &res);
                }else{
                    if(range_begin > number){
                        throw std::runtime_error("Invalid cpu list.");
                    }
                    for(int i=range_begin; i<number; i++){
                        CPU_SET(i, &res);
                    }
                }
                range_begin=-1;
                number=0;
            }else if(ch=='-'){
                if(range_begin!=-1){
                    throw std::runtime_error("Invalid cpu list");
                }
                range_begin=number;
                number=0;
            }else if(isdigit(ch)){
                number=number*10 + (ch-'0');
            }else{
                throw std::runtime_error("Invalid char in cpu list");
            }
        }
        return res;
    }

    cache_cluster read_numa_node(const std::filesystem::path &numa_base)
    {
        cache_cluster res;
        res.type=CacheType::NUMA;
        res.size=read_node_meminfo_int(numa_base/"meminfo", "MemTotal");
        res.level=-1;
        res.cpu_set=read_cpu_list(numa_base/"cpulist");
        return res;
    }

    // Parse information about one CPU from /sys/devices/system/cpu/cpuN/
    std::vector<cache_cluster> build_cpu_cache_info(std::filesystem::path base)
    {
        base /= "cache";

        std::vector<cache_cluster> res;

        int cache_index=-1;
        while(1){
            cache_index += 1; // First is 0

            auto cache_base= base / ("index"+std::to_string(cache_index) );
            std::cerr<<"cache_base="<<cache_base<<"\n";
            if(!exists(cache_base)){
                break;
            }

            auto c=read_cpu_cache_index(cache_base);
            res.push_back(c);
            
            std::cerr<<"  c.level="<<c.level<<", back.level="<<res.back().level<<"\n";
        }

        int num_nodes=0;
        for(auto d : std::filesystem::directory_iterator(base)){
            auto p=d.path().filename().string();
            if(p.substr(0, 4)=="node"){
                if(num_nodes>0){
                    throw std::runtime_error("CPU is in two NUMA nodes...?");
                }
                res.push_back(read_numa_node(p));
            }
        }

        for(auto c : res){
            std::cerr<<"  c.level="<<c.level<<", back.level="<<res.back().level<<"\n";
        }

        return res;
    }

    std::vector<cache_cluster> build_cpu_cache_info(unsigned cpu_index)
    {
        std::filesystem::path base("/sys/devices/system/cpu");
        base /= "cpu"+std::to_string(cpu_index);
        return build_cpu_cache_info(base);
    }

    static bool cpu_set_t_equal(const cpu_set_t &a, const cpu_set_t &b)
    {
        return CPU_EQUAL(&a, &b);
    }

    system_cache_info build_cpu_sets()
    {
        std::set<cache_cluster> all;
        for(auto d : std::filesystem::directory_iterator("/sys/devices/system/cpu")){
            auto sd=d.path().filename().string();
            if(sd.substr(0, 3)=="cpu" && isdigit(sd[3])){
                auto clusters=build_cpu_cache_info(d.path());
                for(auto cu : clusters){
                    std::cerr<<" cu :  type="<<cu.type<<", level="<<cu.level<<", cpus="<<cu.cpu_set<<"\n";
                    auto itb=all.insert(cu);
                    std::cerr<<"  did_insert="<<itb.second<<"\n";
                }
            }
        }

        if(all.empty()){
            throw std::runtime_error("This system apparently has no CPUs or caches!");
        }

        for(auto s : all){
            std::cerr<<"   type="<<s.type<<", level="<<s.level<<", cpus="<<s.cpu_set<<"\n";
        }

        system_cache_info res;
        auto it=all.begin();
        res.levels.push_back({it->level, it->type, it->size, it->size, {it->cpu_set}});
        ++it;

        while(it!=all.end()){
            if(res.levels.back().level!=it->level || res.levels.back().type!=it->type){
                res.levels.push_back({it->level, it->type, it->size, it->size, {it->cpu_set}});
            }else{
                res.levels.back().total_size += it->size;
                res.levels.back().nodes.push_back(it->cpu_set);
            }
            ++it;
        }

        if(res.levels.empty()){
            throw std::runtime_error("Couldn't find any cache levels.");
        }

        res.llc_index=-1;
        for(unsigned i=0; i<res.levels.size(); i++){
            if(res.levels[i].type & Data){
                res.llc_index=i;
            }
        }

        return res;
    }

    const system_cache_info &get_system_cache_info()
    {
        static system_cache_info info=build_cpu_sets();
        return info;
    }

}


#endif
