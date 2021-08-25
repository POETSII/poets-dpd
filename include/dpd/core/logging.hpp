#ifndef logging_hpp
#define logging_hpp

enum ForceLoggingFlags
{
    Asymmetric,       // need to record property seperately for (b0,b1) and (b1,b0)
    Symmetric,        // Same value for (b0,b1) and (b1,b0)
    SymmetricFlipped  // Value for (b0,b1) is negative of value for (b1,b0)
};

class ForceLogging
{
public:
    virtual ~ForceLogging()
    {}

    virtual void SetTime(long t)=0;
    virtual void LogProperty(const char *name, int dims, const double *x)=0;
    virtual void LogBeadProperty(long bead_id, const char *name, int dims, const double *x)=0;
    virtual void LogBeadPairProperty(long bead_id0,long bead_id1, const char *name, ForceLoggingFlags flags, int dims, const double *x)=0;
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
};

#endif