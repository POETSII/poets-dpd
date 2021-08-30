#ifndef logging_hpp
#define logging_hpp

#include <functional>
#include <cstdint>

class ForceLogging
{
public:
    virtual ~ForceLogging()
    {}

    virtual void SetTime(long t)=0;
    virtual void LogProperty(const char *name, int dims, const double *x)=0;
    virtual void LogBeadProperty(long bead_id, const char *name, int dims, const double *x)=0;
    virtual void LogBeadPairProperty(long bead_id0,long bead_id1, const char *name, int dims, const double *x)=0;
    virtual void LogBeadTripleProperty(long bead_id0, long bead_id1, long bead_id2, const char *name, int dims, const double *x)=0;

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
    void LogBeadProperty(long bead_id, const char *name, int dims, const T *x)
    {
        double tmp[3];
        assert(dims<=3);
        std::copy(x, x+dims, tmp);
        LogBeadProperty(bead_id, name, dims, tmp);
    }

    template<class T>
    void LogBeadPairProperty(long bead_id0,long bead_id1, const char *name, int dims, const T *x)
    {
        double tmp[3];
        assert(dims<=3);
        std::copy(x, x+dims, tmp);
        LogBeadPairProperty(bead_id0, bead_id1, name, dims, tmp);
    }

    template<class T>
    void LogBeadTripleProperty(long bead_id0, long bead_id1, long bead_id2, const char *name, int dims, const T *x)
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
    void LogBeadProperty(long bead_id, const char *name, const vec3g_t<T> &x)
    {
        double tmp[3]={x[0],x[1],x[2]};
        LogBeadProperty(bead_id, name, 3, tmp);
    }

    template<class T>
    void LogBeadPairProperty(long bead_id0,long bead_id1, const char *name, const vec3g_t<T> &x)
    {
        double tmp[3]={x[0],x[1],x[2]};
        LogBeadPairProperty(bead_id0, bead_id1, name, 3, tmp);
    }

    template<class T>
    void LogBeadTripleProperty(long bead_id0, long bead_id1, long bead_id2, const char *name, const vec3g_t<T> &x)
    {
        double tmp[3]={x[0],x[1],x[2]};
        LogBeadTripleProperty(bead_id0, bead_id1, bead_id2, name, 3, tmp);
    }



    using bead_hash_to_id_t = std::function<uint32_t(long)>;

    static bead_hash_to_id_t &bead_hash_to_id()
    {
        static bead_hash_to_id_t _map;
        return _map;
    } 

    template<class ...TArgs>
    void LogBeadHashProperty(uint32_t bead_hash0, const char *name,  TArgs ...args)
    {
        bead_hash_to_id_t &map=bead_hash_to_id();
        assert(map);
        LogBeadProperty(map(bead_hash0), name, args...);
    }

    template<class ...TArgs>
    void LogBeadHashPairProperty(uint32_t bead_hash0, uint32_t bead_hash1, const char *name, TArgs ...args)
    {
        bead_hash_to_id_t &map=bead_hash_to_id();
        assert(map);
        LogBeadPairProperty(map(bead_hash0), map(bead_hash1), name, args...);
    }

    template<class ...TArgs>
    void LogBeadHashTripleProperty(uint32_t bead_hash0, uint32_t bead_hash1, uint32_t bead_hash2, const char *name, TArgs ...args)
    {
        bead_hash_to_id_t &map=bead_hash_to_id();
        assert(map);
        LogBeadTripleProperty(map(bead_hash0), map(bead_hash1), map(bead_hash2), name, args...);
    }
};

#endif