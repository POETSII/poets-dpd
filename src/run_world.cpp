
#include "dpd/core/dpd_engine.hpp"
#include "dpd/core/dpd_state_builder.hpp"

#include "dpd/core/dpd_state_io.hpp"
#include "dpd/core/dpd_state_to_vtk.hpp"
#include "dpd/core/dpd_state_validator.hpp"

#include "dpd/core/logging.hpp"

#include <random>
#include <fstream>

void usage()
{
    fprintf(stderr, "run_world : engine-name src-file output-base-name interval_count interval_size \n");
    fprintf(stderr, "  engine names:\n");
    for(auto s : DPDEngineFactory::ListFactories()){
        fprintf(stderr, "    %s\n", s.c_str());
    }
    fprintf(stderr,"  env:\n");
    fprintf(stderr,"     PDPD_LOG=log-path : do full force logging to given file.\n");
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

class FileLogger
    : public ForceLogging
{
public:
    FileLogger(const std::string &name)
        : m_dst(name)
        , m_t(0)
    {
        if(!m_dst.is_open()){
            fprintf(stderr, "Couldnt open force logging file %s\n", name.c_str());
        }
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

    virtual void LogProperty(const char *name, int dims, const double *x)
    {
        m_dst<<"Prop,"<<m_t<<",,,,"<<name;
        for(int i=0; i<dims; i++){
            m_dst<<",";
            print(x[i]);
        }
        for(int i=dims;i<3;i++){
            m_dst<<",";
        }
        m_dst<<"\n";
    }

    virtual void LogBeadProperty(long bead_id, const char *name, int dims, const double *x)
    {
        m_dst<<"Prop,"<<m_t<<","<<bead_id+1<<",,,"<<name;
        for(int i=0; i<dims; i++){
            m_dst<<",";
            print(x[i]);
        }
        for(int i=dims;i<3;i++){
            m_dst<<",";
        }
        m_dst<<"\n";
    }

    virtual void LogBeadPairProperty(long bead_id0,long bead_id1, const char *name, int dims, const double *x)
    {
        m_dst<<"Prop,"<<m_t<<","<<bead_id0+1<<","<<bead_id1+1<<",,"<<name;
        for(int i=0; i<dims; i++){
            m_dst<<",";
            print(x[i]);
        }
        for(int i=dims;i<3;i++){
            m_dst<<",";
        }
        m_dst<<"\n";
    }

    virtual void LogBeadTripleProperty(long bead_id0, long bead_id1, long bead_id2, const char *name, int dims, const double *x)
    {
        m_dst<<"Prop,"<<m_t<<","<<bead_id0+1<<","<<bead_id1+1<<","<<bead_id2+1<<","<<name;
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
};

int main(int argc, const char *argv[])
{
    try{
        std::string engine_name;
        if(argc>1){
            engine_name=argv[1];
        }else{
            usage();
        }

        std::string srcFile;
        if(argc>2){
            srcFile=argv[2];
        }else{
            usage();
        }

        std::string baseName;
        if(argc>3){
            baseName=argv[3];
        }else{
            usage();
        }

        std::string log_dst;
        if(getenv("PDPD_LOG")){
            log_dst=getenv("PDPD_LOG");
            ForceLogging::set_logger(new FileLogger(log_dst));
            fprintf(stderr, "Logging to %s\n", log_dst.c_str());
        }

        int interval_count=1000;
        if(argc>4){
            interval_count=std::stoi(argv[4]);
        }

        int interval_size=1;
        if(argc>5){
            interval_size=std::stoi(argv[5]);
        }

        std::ifstream input(srcFile.c_str());
        if(!input.is_open()){
            fprintf(stderr, "Couldnt open %s\n", srcFile.c_str());
            exit(1);
        }

        int line_no=0;
        WorldState state=read_world_state(input, line_no);

        validate(state);

        std::shared_ptr<DPDEngine> engine = DPDEngineFactory::create(engine_name);

        int volume=state.box[0]*state.box[1]*state.box[2];

        double t0=now();

        engine->Attach(&state);

        double t1=now();

        unsigned done=0;
        unsigned slice_i=0;
        engine->Run(interval_count, interval_size, [&](){
            double t2=now();
            done+=interval_size;

            fprintf(stdout, "%s,%s,%d,%d,%d,%g,%g,%g,%g\n",
                engine_name.c_str(), srcFile.c_str(), volume, (int)state.beads.size(), done,
                t1-t0, t2-t1, state.beads.size()/(t2-t0)*done, state.beads.size()/(t2-t1)*interval_size
            );
            t1=t2;

            std::vector<char> tmp(baseName.size()+100);
            snprintf(&tmp[0], tmp.size()-1, "%s.%06d.state", baseName.c_str(), slice_i);
            std::ofstream output(&tmp[0]);
            if(!output.is_open()){
                fprintf(stderr, "Couldnt create file %s\n", &tmp[0]);
                exit(1);
            }
            write_world_state(output, state);
            output.close();

            snprintf(&tmp[0], tmp.size()-1, "%s.%06d.vtk", baseName.c_str(), slice_i);
            output.open(&tmp[0]);
            if(!output.is_open()){
                fprintf(stderr, "Couldnt create file %s\n", &tmp[0]);
                exit(1);
            }
            write_to_vtk(output, state);
            output.close();


            ++slice_i;
            return true;
        });

        if(ForceLogging::logger()){
            delete ForceLogging::logger();
            ForceLogging::logger()=0;
        }

    }catch(const std::exception &e){
        print_exception(e);
        exit(1);
    }

    return 0;
}
