#include "dpd/tests/test_ordered_mesh.hpp"
#include "dpd/tests/test_bonded_bead_angle.hpp"
#include "dpd/tests/test_bonded_bead_triple.hpp"
#include "dpd/tests/test_bonded_bead_pair.hpp"
#include "dpd/tests/test_colliding_bead_pair_diss.hpp"
#include "dpd/tests/test_colliding_bead_pair.hpp"
#include "dpd/tests/test_isolated_moving_bead.hpp"
#include "dpd/tests/test_multiple_moving_beads.hpp"
#include "dpd/engines/naive/naive_dpd_engine.hpp"

#include "dpd/tests/test_differential.hpp"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/skewness.hpp>

#include <random>

using acc_t = boost::accumulators::accumulator_set<
    double,
    boost::accumulators::stats<
        boost::accumulators::tag::min,
        boost::accumulators::tag::mean,
        boost::accumulators::tag::max,
        boost::accumulators::tag::variance,
        boost::accumulators::tag::skewness
    >
>;

std::ostream &operator<<(std::ostream &dst, const acc_t &x)
{
    using namespace boost::accumulators;
    return dst<<min(x)<<","<<mean(x)<<","<<max(x)<<","<<sqrt(variance(x))<<","<<skewness(x);
}

struct header
{
    std::string name;
};

std::ostream &operator<<(std::ostream &dst, const header &h)
{
    std::string x=h.name;
    return dst<<"min("<<x<<"),mean("<<x<<"),max("<<x<<"),stddev("<<x<<"),skewness("<<x<<")";
}

void usage()
{
    fprintf(stderr, "engine_diff : engine1 engine2 state-file interval_count interval_size\n");
    fprintf(stderr, "  engine names:\n");
    for(auto s : DPDEngineFactory::ListFactories()){
        fprintf(stderr, "    %s\n", s.c_str());
    }
    exit(1);
}

int main(int argc, const char *argv[])
{
    std::string engine1_name;
    if(argc>1){
        engine1_name=argv[1];
    }else{
        usage();
    }

    std::string engine2_name;
    if(argc>2){
        engine2_name=argv[2];
    }else{
        usage();
    }

    std::string state_file;
    if(argc>3){
        state_file=argv[3];
    }else{
        usage();
    }

    unsigned interval_count=10;
    if(argc>3){
        interval_count=std::stoi(argv[4]);
    }else{
        usage();
    }

    unsigned interval_size=1;
    if(argc>4){
        interval_size=std::stoi(argv[5]);
    }else{
        usage();
    }


    std::shared_ptr<DPDEngine> pengine1 = DPDEngineFactory::create(engine1_name);
    std::shared_ptr<DPDEngine> pengine2 = DPDEngineFactory::create(engine2_name);

    WorldState s=read_world_state(state_file);

    size_t state_size = 2*sizeof(Bead)*s.beads.size()*interval_count;
    std::cerr<<"Total size of state = "<<(state_size*pow(2.0,-20))<<"MB\n";
    uint64_t bead_steps = 2*s.beads.size()*interval_count;
    std::cerr<<"Total mega-bead-steps = "<<(bead_steps*1e-6)<<" MBeadSteps\n";

    auto run_rep=[&](std::shared_ptr<DPDEngine> engine, const WorldState &s) -> std::vector<std::vector<Bead>>
    {
        WorldState state=s;
        std::vector<std::vector<Bead>> res;

        std::string ok=engine->CanSupport(&state);
        if(!ok.empty()){
            fprintf(stderr, "Engine can't support world : %s\n", ok.c_str());
            exit(1);
        }

        engine->Attach(&state);

        unsigned t=state.t;
        unsigned done=engine->Run(interval_count, interval_size, [&]() -> bool {
            if(t+interval_size!=state.t){
                throw std::runtime_error("Mis-matched time value.");
            }

            res.push_back(state.beads);

            t=state.t;

            return true;
        });

        return res;
    };

    std::vector<std::vector<Bead>> beads1=run_rep(pengine1,s);
    std::vector<std::vector<Bead>> beads2=run_rep(pengine2,s);

    std::cout<<"t,"<<header{"x"}<<","<<header{"v"}<<","<<header{"f"}<<"\n";
    for(unsigned i=0; i<interval_count; i++){
        const auto &s1=beads1[i], &s2=beads2[i];
        acc_t x_diff,v_diff,f_diff;
        acc_t v_mag_mean1, v_mag_mean2;
        acc_t f_mag_mean1, f_mag_mean2;
        double worst_f_err=0;
        for(unsigned j=0; j<s.beads.size(); j++){
            vec3r_t dx=(s1[j].x-s2[j].x);
            for(int d=0;d<3;d++){
                if(dx[d] > s.box[d]/2){
                    dx[d]-=s.box[d];
                }
                if(dx[d] < -s.box[d]/2){
                    dx[d]+=s.box[d];
                }
            }

            x_diff(dx.l2_norm());
            v_diff((s1[j].v-s2[j].v).l2_norm());
            double f_err=(s1[j].f-s2[j].f).l2_norm();
            f_diff(f_err);
            if(f_err > worst_f_err){
                std::cout.flush();
                std::cout<<"# Worst f_err "<<f_err<<", b1="<<s1[j].get_hash_code().reduced_hash().hash<<", "<<s1[j].x<<","<<s1[j].v<<","<<s1[j].f<<", b1="<<s1[j].get_hash_code().reduced_hash().hash<<", "<<s2[j].x<<","<<s2[j].v<<","<<s2[j].f<<"\n";
                worst_f_err=f_err;
            }

            v_mag_mean1( s1[j].v.l2_norm() );
            v_mag_mean2( s2[j].v.l2_norm() );
            f_mag_mean1( s1[j].f.l2_norm() );
            f_mag_mean2( s2[j].f.l2_norm() );
        }

        using boost::accumulators::mean;
        std::cout<<s.t+i*interval_size<<",  "<<x_diff<<",  "<<v_diff<<",  "<<f_diff<<",  "<<mean(v_mag_mean1)<<","<<mean(v_mag_mean2)<<", "<<mean(f_mag_mean1)<<","<<mean(f_mag_mean2)<<"\n";
    }
}
