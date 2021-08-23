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
    void SetTime(long t);
    void LogBeadProperty(long bead_id, const char *name, int dims, const double *x);
    void LogBeadPairProperty(long bead_id0,long bead_id1, const char *name, ForceLoggingFlags flags, int dims, const double *x);
    void LogBeadTripleProperty(long bead_id0, long bead_id1, long bead_id2, const char *name, int dims, const double *x);

    static ForceLogging *logger = 0;
private:
    long m_t;
};

#endif