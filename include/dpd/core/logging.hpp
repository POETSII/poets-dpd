#ifndef logging_hpp
#define logging_hpp

#include <functional>
#include <cstdint>

#include "dpd/core/dpd_state.hpp"

class ForceLogging
{
public:
    virtual ~ForceLogging()
    {}

    virtual void SetTime(long t)=0;
    virtual void LogProperty(const char *name, int dims, const double *x)=0;
    virtual void LogBeadProperty(BeadHash b1, const char *name, int dims, const double *x)=0;
    virtual void LogBeadPairProperty(BeadHash bead_id0,BeadHash bead_id1, const char *name, int dims, const double *x)=0;
    virtual void LogBeadTripleProperty(BeadHash bead_id0, BeadHash bead_id1, BeadHash bead_id2, const char *name, int dims, const double *x)=0;

    static ForceLogging *&logger()
    {
        static ForceLogging *_logger=0;
        return _logger;
    }

    static void set_logger(ForceLogging *s)
    {
        assert(!logger());
        logger()=s;
    }

    template<class T>
    void LogProperty(const char *name, int dims, const T *x)
    {
        double tmp[3];
        assert(dims<=3);
        std::copy(x, x+dims, tmp);
        LogProperty(name, dims, tmp);
    }

    template<class T>
    void LogBeadProperty(BeadHash bead_id, const char *name, int dims, const T *x)
    {
        double tmp[3];
        assert(dims<=3);
        std::copy(x, x+dims, tmp);
        LogBeadProperty(bead_id, name, dims, tmp);
    }

    template<class T>
    void LogBeadPairProperty(BeadHash bead_id0,BeadHash bead_id1, const char *name, int dims, const T *x)
    {
        double tmp[3];
        assert(dims<=3);
        std::copy(x, x+dims, tmp);
        LogBeadPairProperty(bead_id0, bead_id1, name, dims, tmp);
    }

    template<class T>
    void LogBeadTripleProperty(BeadHash bead_id0, BeadHash bead_id1, BeadHash bead_id2, const char *name, int dims, const T *x)
    {
        double tmp[3];
        assert(dims<=3);
        std::copy(x, x+dims, tmp);
        LogBeadTripleProperty(bead_id0, bead_id1, bead_id2, name, dims, tmp);
    }


    template<class T>
    void LogProperty(const char *name, const vec3g_t<T> &x)
    {
        double tmp[3]={x[0],x[1],x[2]};
        LogProperty(name, 3, tmp);
    }

    template<class T>
    void LogBeadProperty(BeadHash bead_id, const char *name, const vec3g_t<T> &x)
    {
        double tmp[3]={x[0],x[1],x[2]};
        LogBeadProperty(bead_id, name, 3, tmp);
    }

    template<class T>
    void LogBeadPairProperty(BeadHash bead_id0,BeadHash bead_id1, const char *name, const vec3g_t<T> &x)
    {
        double tmp[3]={x[0],x[1],x[2]};
        LogBeadPairProperty(bead_id0, bead_id1, name, 3, tmp);
    }

    template<class T>
    void LogBeadTripleProperty(BeadHash bead_id0, BeadHash bead_id1, BeadHash bead_id2, const char *name, const vec3g_t<T> &x)
    {
        double tmp[3]={x[0],x[1],x[2]};
        LogBeadTripleProperty(bead_id0, bead_id1, bead_id2, name, 3, tmp);
    }

};

#endif