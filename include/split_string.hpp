#ifndef split_string_hpp
#define split_string_hpp

#include <cstdint>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <cstring>

// Ideally we would use charconv, but it's taking a while to get through
//#include <charconv>



struct split_string
{
    std::string data;            // String with spaces replaced by nulls
    std::vector<std::pair<uint16_t,uint16_t>> parts; // Contains the beginning index of each null-terminated sub-part

    std::string string_at(unsigned p) const
    {
        auto e=parts.at(p);
        return {data.begin()+e.first, data.begin()+e.first+e.second};
    }

    /*
    template<class T>
    T T_at(unsigned index) const
    {
        auto e=x.parts.at(index);
        T res;
        const char *begin=x.data.begin()+e.first;
        const char *end=x.data.begin()+e.second;
        auto c=std::from_chars(begin, end, res);
        if(c.ptr!=end){
            throw std::runtime_error("Not all chars consumed while parsing '"+string_at(index)+"'");
        }
        if(c.ec!=std::errc()){
            throw std::runtime_error("Couldn't parse '"+string_at(index)+"'");
        }
        return res;
    }

    double double_at(unsigned index) const
    { return this->template T_at<double>(index); }

    unsigned unsigned_at(unsigned index) const
    { return this->template T_at<unsigned>(index); }
    */

    double double_at(unsigned index) const
    {
        auto e=parts.at(index);
        char *ep;
        const char *begin=&data[0]+e.first;
        const char *end=&data[0]+e.first+e.second;
        double r=strtod(begin, &ep);
        if(ep!=end){
            throw std::runtime_error("Not all chars consumed while parsing '"+string_at(index)+"'");
        }
        return r;
    }

    uint64_t unsigned_at(unsigned index) const
    {
        static_assert(sizeof(unsigned long long)==sizeof(uint64_t));
        auto e=parts.at(index);
        char *ep;
        const char *begin=&data[0]+e.first;
        const char *end=&data[0]+e.first+e.second;
        uint64_t r=strtoull(begin, &ep, 10);
        if(ep!=end){
            throw std::runtime_error("Not all chars consumed while parsing '"+string_at(index)+"'");
        }
        return r;
    }

    // Rebuilds the line, mainly for error purposes
    std::string full() const
    {
        std::string tmp(data);
        for(unsigned i=0; i<tmp.size(); i++){
            if(tmp[i]==0){
                tmp[i]=' ';
            }
        }
        return tmp;
    }
};


split_string read_prefixed_line_and_split_on_space(std::istream &src, const std::string &prefix, int &line_no)
{
    split_string res;
    std::string &s=res.data;
    while(1){
        if(!std::getline(src, s)){
            throw std::runtime_error("Expecting prefix '"+prefix+"' but hit end of file.\n");
        }
        line_no++;

        if(s.size() > 10000){ // We use 16-bit indices in return
            throw std::runtime_error("This method is not expecting to deal with very long lines.");
        }
        
        int start=-1;
        unsigned curr=0;

        while(curr<=s.size()){
            if(isspace(s[curr]) || s[curr]==0){ // C++11 - s[s.size()] is null, and we are allowed to write it (but only with null)
                if(start!=-1){
                    assert( curr-start > 0 );
                    res.parts.push_back({start, curr-start});
                }
                start=-1;
                s[curr]=0;
            }else{
                if(start==-1){
                    start=curr;
                }
            }
            ++curr;
        }

        if(res.parts.empty()){
            continue;
        }

        
        if(res.data[res.parts[0].first]=='#'){
            continue;  // This is slow, but assume comments are quite infrequent
        }

        if(strcmp(&res.data[res.parts[0].first],prefix.c_str())){
            throw std::runtime_error("Expecting line prefix '"+prefix+"' but got '"+res.string_at(0)+"' at line "+std::to_string(line_no)+", value ='"+res.full()+"'");
        }

        return res;
    }
}

split_string read_prefixed_line_and_split_on_space(std::istream &src, const std::string &prefix, unsigned num_elements, int &line_no)
{
    split_string res=read_prefixed_line_and_split_on_space(src, prefix, line_no);
    if(res.parts.size()!=num_elements){
        throw std::runtime_error("Expecting "+std::to_string(num_elements)+" for prefix "+prefix+" on line "+std::to_string(line_no)+", value ='"+res.full()+"'");
    }
    return res;
}

#endif