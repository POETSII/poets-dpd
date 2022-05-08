#ifndef logging_impl_hpp
#define logging_impl_hpp

#include "logging.hpp"

#include <fstream>

class FileLogger
    : public ForceLogging
{
public:
    FileLogger(const std::string &name)
        : m_dst(name)
        , m_t(0)
    {
        fprintf(stderr, "dst=%s\n", name.c_str());
        if(!m_dst.is_open()){
            fprintf(stderr, "Couldnt open force logging file %s\n", name.c_str());
        }
    }

    ~FileLogger()
    {
        m_dst.flush();
    }

    void print(double val)
    {
        if(std::round(val)==val && 0<=val && val < ldexp(1, 32)){
            m_dst<<(uint32_t)val;
        }else{
            m_dst<<val;
        }
    }

    virtual void SetTime(long t)
    { m_t=t; }

    virtual void SetPrefix(const std::string &prefix)
    {
        m_prefix=prefix;
    }

    virtual void LogProperty(const char *name, int dims, const double *x)
    {
        m_dst<<m_prefix<<"Prop,"<<m_t<<",,,,"<<name;
        for(int i=0; i<dims; i++){
            m_dst<<",";
            print(x[i]);
        }
        for(int i=dims;i<3;i++){
            m_dst<<",";
        }
        m_dst<<"\n";
    }

    virtual void LogBeadProperty(BeadHash bead_id, const char *name, int dims, const double *x)
    {
        m_dst<<m_prefix<<"Prop,"<<m_t<<","<<bead_id.reduced_hash().hash<<",,,"<<name;
        for(int i=0; i<dims; i++){
            m_dst<<",";
            print(x[i]);
        }
        for(int i=dims;i<3;i++){
            m_dst<<",";
        }
        m_dst<<"\n";
    }

    virtual void LogBeadPairProperty(BeadHash bead_id0,BeadHash bead_id1, const char *name, int dims, const double *x)
    {
        m_dst<<m_prefix<<"Prop,"<<m_t<<","<<bead_id0.reduced_hash().hash<<","<<bead_id1.reduced_hash().hash<<",,"<<name;
        for(int i=0; i<dims; i++){
            m_dst<<",";
            print(x[i]);
        }
        for(int i=dims;i<3;i++){
            m_dst<<",";
        }
        m_dst<<"\n";
    }

    virtual void LogBeadTripleProperty(BeadHash bead_id0, BeadHash bead_id1, BeadHash bead_id2, const char *name, int dims, const double *x)
    {
        m_dst<<m_prefix<<"Prop,"<<m_t<<","<<bead_id0.reduced_hash().hash<<","<<bead_id1.reduced_hash().hash<<","<<bead_id2.reduced_hash().hash<<","<<name;
        for(int i=0; i<dims; i++){
            m_dst<<",";
            print(x[i]);
        }
        for(int i=dims;i<3;i++){
            m_dst<<",";
        }
        m_dst<<"\n";
    }
private:
    std::ofstream m_dst;
    long m_t;
    std::string m_prefix;
};

#endif
