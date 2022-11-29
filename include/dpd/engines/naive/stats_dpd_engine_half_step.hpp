#ifndef stats_dpd_engine_half_step_hpp
#define stats_dpd_engine_half_step_hpp

#include "dpd/core/dpd_engine.hpp"
#include "dpd/engines/naive/naive_dpd_engine_half_step.hpp"

#include "dpd/maths/dpd_maths_core_half_step.hpp"

#include "dpd/core/vec3.hpp"
#include "dpd/core/hash.hpp"

#include <cassert>
#include <cmath>
#include <array>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <cstdint>
#include <fstream>

#include <dpd/core/json_helper.hpp>


#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/skewness.hpp>
#include <boost/accumulators/statistics/kurtosis.hpp>
#include <boost/accumulators/statistics/min.hpp>

#include <cmath>

namespace dpd_stats
{
    using namespace boost::accumulators;

    struct count_stats
    {
        std::vector<uint64_t> counts;

        void clear()
        {
            counts.clear();
            counts.push_back(0);
        }

        void add(unsigned i)
        {
            if(i>=counts.size()){
                counts.resize(i+1, 0);
            }
            counts[i] += 1;
        }

        void add(const count_stats &s)
        {
            if(s.counts.size() > counts.size()){
                counts.resize(s.counts.size());
            }
            for(unsigned i=0; i<s.counts.size(); i++){
                counts[i] += s.counts[i];
            }
        }

        JSON to_json(rapidjson::Document &doc)
        {
            JSON res(doc);

            res.AddMember("type","histogram");
            uint64_t count=0;
            double sum=0;
            for(unsigned i=0; i<counts.size(); i++){
                count += counts[i];
                sum += counts[i] * i;
            }
            res.AddMember("count", count);
            if(count!=0){
                double mean=sum/count;
                double sum_sqr=0; 
                for(unsigned i=0; i<counts.size(); i++){
                    double xx=counts[i] - mean;
                    sum_sqr += counts[i]*xx*xx;
                }
                double stddev=sqrt( sum_sqr / count );
                res.AddMember("mean", mean);
                res.AddMember("stddev", stddev);
                std::vector<double> probs(counts.size());
                for(unsigned i=0; i<counts.size(); i++){
                    probs[i] = counts[i] / (double)count;
                }
                res.AddMember("probs", probs);
            }
            res.AddMember("counts", counts);
            return res;
        }
    };

    struct float_stats{
        uint64_t count;
        double minval, maxval;
        double sum;
        double sum_sqr; // Unstable sum. Life is too short.
        double abs_smallest[2], abs_biggest[2];
        
        uint64_t zero_count, nan_count, inf_count[2];
        std::array<uint64_t, 4*(DBL_MAX_EXP-DBL_MIN_EXP+1)> exp_histogram[2];

        float_stats()
        {
            clear();
        }

        void clear()
        {
            memset(this, 0, sizeof(*this));
            abs_smallest[0]=DBL_MAX;
            abs_smallest[1]=DBL_MAX;
            minval=DBL_MAX;
            maxval=-DBL_MAX;
        }

        void add(const float_stats &stats)
        {
            count += stats.count;
            sum += stats.sum;
            sum_sqr += stats.sum_sqr;
            zero_count += stats.zero_count;
            nan_count += stats.nan_count;
            for(int sign=0; sign<2; sign++){
                abs_smallest[sign] = std::min(abs_smallest[sign], stats.abs_smallest[sign]);
                abs_biggest[sign] = std::max(abs_biggest[sign], stats.abs_biggest[sign]);
                inf_count[sign] += stats.inf_count[sign];
                for(int i=0; i<exp_histogram[0].size(); i++){
                    exp_histogram[sign][i] += stats.exp_histogram[sign][i];
                }
            }
            minval=std::min(minval, stats.minval);
            maxval=std::max(maxval, stats.maxval);
        }
        
        void add(double x)
        {
            count += 1;
            if(isnan(x)){
                nan_count++;
                return;
            }
            minval=std::min(minval, x);
            maxval=std::max(maxval, x);
            
            if(x==0){
                zero_count++;
                return;
            }
        
            double ax=x;
            int sign=0;
            if(x<0){
                sign=1;
                ax=-x;
            }
            if(isinf(x)){
                inf_count[sign]++;
            }else{
                sum += x;
                sum_sqr += x*x;

                abs_smallest[sign]=std::min(abs_smallest[sign], ax );
                abs_biggest[sign]=std::max(abs_biggest[sign], ax);
                
                static const double split01 = std::pow(2, -0.75);
                static const double split10 = std::pow(2, -0.50);
                static const double split11 = std::pow(2, -0.25);
                int e;
                double f=frexp(x, &e);
                assert(DBL_MIN_EXP <= e && e <= DBL_MAX_EXP);
                assert(0.5<=f && f<1.0);
                unsigned index=4*e-4*DBL_MIN_EXP;
                // This is so dumb. Just take the MSBs already...
                if(f>=split10){
                    index += f>=split11 ? 3 : 2;
                }else{
                    index += f>=split01;
                }
                exp_histogram[sign].at(index) ++;

                //assert(std::accumulate(exp_histogram[sign].begin(), exp_histogram[sign].end(), 0) > 0);
            }
            
        }

        JSON to_json(rapidjson::Document &doc)
        {
            JSON res(doc);
            res.AddMember("type","float");
            res.AddMember("count", count);

            if(count>0){
                res.AddMember("pos_biggest", abs_biggest[0]);
                res.AddMember("pos_smallest", abs_smallest[0]);
                res.AddMember("neg_biggest", -abs_biggest[1]);
                res.AddMember("neg_smallest", -abs_smallest[1]);
                res.AddMember("min", minval);
                res.AddMember("max", maxval);
                double mean=sum / count;
                double stddev=sqrt(sum_sqr / count - mean*mean);
                res.AddMember("mean", mean);
                res.AddMember("stddev", stddev);

                std::vector<uint64_t> c;
                std::vector<double> v, probs;

                auto get_index_val=[]( int index ) -> double
                {
                    static const double fs[4]={std::pow(2.0, -1),std::pow(2.0, -0.75),std::pow(2.0, -0.5),std::pow(2.0, -0.25)};
                    return ldexp( fs[index%4], index/4+DBL_MIN_EXP);
                };

                auto nneg=std::accumulate(exp_histogram[1].begin(), exp_histogram[1].end(), 0ul);

                int first_non_zero=INT_MAX, last_non_zero=-1;
                for(int i=0; i<exp_histogram[1].size(); i++){
                    if(exp_histogram[1][i]>0){
                        first_non_zero=std::min(first_non_zero, i);
                        last_non_zero=std::max(last_non_zero, i);
                    }
                }    
                
                std::cerr<<"d=1, count="<<count<<", fnz="<<first_non_zero<<", lnz="<<last_non_zero<<"\n";
                if(first_non_zero <= last_non_zero){
                    for(int i=last_non_zero; i>=first_non_zero; i--){
                        c.push_back(exp_histogram[1][i]);
                        v.push_back(-get_index_val(i));
                        probs.push_back(exp_histogram[1][i] / (double)count);
                    }
                }
                if(nneg>0){
                    assert(first_non_zero <= last_non_zero);
                }

                c.push_back(zero_count);
                v.push_back(0.0);
                probs.push_back(zero_count / (double)count);

                auto npos=std::accumulate(exp_histogram[0].begin(), exp_histogram[0].end(), 0ul);
                first_non_zero=INT_MAX, last_non_zero=-1;
                for(int i=0; i<exp_histogram[0].size(); i++){
                    if(exp_histogram[0][i]>0){
                        first_non_zero=std::min(first_non_zero, i);
                        last_non_zero=std::max(last_non_zero, i);
                    }
                }    
                if(npos>0){
                    assert(first_non_zero <= last_non_zero);
                }

                std::cerr<<"d=0, count="<<count<<", fnz="<<first_non_zero<<", lnz="<<last_non_zero<<"\n";
                if(first_non_zero <= last_non_zero){
                    for(int i=first_non_zero; i<=last_non_zero; i++){
                        c.push_back(exp_histogram[0][i]);
                        v.push_back(get_index_val(i));  
                        probs.push_back(exp_histogram[0][i] / (double)count);
                    }
                }

                assert(npos+nneg+zero_count+nan_count==count);

                res.AddMember("counts", c);
                res.AddMember("count_lower_bounds", v);
            }

            return res;
        }
    };
};

/*
    This method must have lambda=0.5

    The process to get from t to t+dt is:
    - UpdatePos : calculate position at time t+dt using values only from time t
        // Pre: x(t),v(t),f(t) all live
        x(t+dt)   = x(t) + v(t)*dt + f(t)*(dt^2/2)
        v(t+dt/2) = v(t) + f(t)*(dt/2),  This assumes lambda=0.5
        // Post: x(t+dt),v(t+dt) live,   x(t) is dead, f(t) is dead

    - UpdateForce : calculate forces using positions at time t+dt and velocities at t+dt/2
        // Pre: x(t+dt),v(t+dt/2)
        f(t+dt) = (function of x(t+dt),v(t+dt/2))
        // Post: f(t+dt),x(t+dt),f(t),v(t) all live

    - UpdateMom : calculate velocity at time t+dt
        v(t+dt) = v(t+dt/2) + dt*f(t+dt)/2
        // Post: f(t+dt),v(t+dt),x(t+dt) all live, v(t+dt/2) dead

    To keep state clean, the approach taken here is not to pre-calculate
    v'(t), and just to recalculate on the fly. This engine is already
    pretty inefficient anyway.
*/
class StatsDPDEngineHalfStep
    : public NaiveDPDEngineHalfStep
{
    
public:
    enum InteractionDir
    {
        All,
        Centre,
        Face,
        Edge,
        Corner,
        DirCount
    };

    static const char *get_interaction_dir_name(InteractionDir dir)
    {
        switch(dir){
            case All: return "All";
            case Centre: return "Centre";
            case Face: return "Face";
            case Edge: return "Edge";
            case Corner: return "Corner";
            default: assert(0); throw std::invalid_argument("Invalid direction.");
        }
    }

    static InteractionDir calc_interaction_dir(bool dx, bool dy, bool dz)
    {
        switch(dx+dy+dz){
        case 0: return Centre;
        case 1: return Face;
        case 2: return Edge;
        case 3: return Corner;
        default: assert(0); return Centre;
        }
    }

    static InteractionDir calc_interaction_dir(vec3i_t a, vec3i_t b)
    {
        return calc_interaction_dir(
            a.x[0]!=b.x[0], a.x[1]!=b.x[1], a.x[2]!=b.x[2]
        );
    }

private:
    struct interaction_stats
    {
        uint64_t n;
        uint64_t hit;
        dpd_stats::float_stats r;

        void clear()
        {
            n=0;
            hit=0;
            r.clear();
        }

        void add(const interaction_stats &o)
        {
            n += o.n;
            hit += o.hit;
            r.add(o.r);
        }

        JSON to_json(rapidjson::Document &doc)
        {
            JSON res(r.to_json(doc));
            res.AddMember("hit_count", hit);
            if(n>0){
                res.AddMember("hit_prob", hit/(double)n);
            }
            return res;
        }
    };

    struct round_stats
    {
        uint64_t t;
        uint64_t t_count;
        dpd_stats::count_stats migrations_total;

        dpd_stats::count_stats beads_per_cell;
        std::array<dpd_stats::count_stats,5> migrations_in_per_cell_by_dir;
        std::array<dpd_stats::count_stats,5> migrations_out_per_cell_by_dir;

        dpd_stats::float_stats bead_velocity;
        dpd_stats::float_stats bead_force;
        dpd_stats::float_stats global_velocity;
        dpd_stats::float_stats global_temperature;        

        interaction_stats interactions_all;
        std::array<interaction_stats,5> interactions_by_dir;

        void clear()
        {
            t=0;
            t_count=0;
            migrations_total.clear();
            beads_per_cell.clear();
            for(int i=0; i<DirCount; i++){
                migrations_in_per_cell_by_dir[i].clear();
                migrations_out_per_cell_by_dir[i].clear();
                interactions_by_dir[i].clear();
            }
            interactions_all.clear();
            bead_velocity.clear();
            bead_force.clear();
            global_velocity.clear();
            global_temperature.clear();
        }

        void add(const round_stats &o)
        {
            if(t_count==0){
                t=o.t;
                t_count=o.t_count;
            }else{
                if(t+t_count!=o.t){
                    throw std::runtime_error("Non-contiguous t.");
                }
                t_count += o.t_count;
            }
            migrations_total.add(o.migrations_total);
            beads_per_cell.add(o.beads_per_cell);
            for(int i=0; i<DirCount; i++){
                migrations_in_per_cell_by_dir[i].add(o.migrations_in_per_cell_by_dir[i]);
                migrations_out_per_cell_by_dir[i].add(o.migrations_out_per_cell_by_dir[i]);
                interactions_by_dir[i].add(o.interactions_by_dir[i]);
            }
            interactions_all.add(o.interactions_all);
            bead_velocity.add(o.bead_velocity);
            bead_force.add(o.bead_force);
            global_velocity.add(o.global_velocity);
            global_temperature.add(o.global_temperature);
        }

        JSON to_json(rapidjson::Document &doc)
        {
            JSON res(doc);
            res.AddMember("t", t);
            res.AddMember("t_count", t_count);
            res.AddMember("migrations_total", migrations_total.to_json(doc));
            res.AddMember("beads_per_cell", beads_per_cell.to_json(doc));
            res.AddMember("bead_velocity", bead_velocity.to_json(doc));
            res.AddMember("bead_force", bead_force.to_json(doc));
            res.AddMember("global_velocity", global_velocity.to_json(doc));
            res.AddMember("global_temperature", global_temperature.to_json(doc));
            JSON interactions(doc), migrations_in(doc), migrations_out(doc);
            for(int dir =0; dir < DirCount; dir++ ){
                std::string dn=get_interaction_dir_name((InteractionDir)dir);
                interactions.AddMember(dn, interactions_by_dir[dir].to_json(doc));
                migrations_in.AddMember(dn, migrations_in_per_cell_by_dir[dir].to_json(doc));
                migrations_out.AddMember(dn, migrations_in_per_cell_by_dir[dir].to_json(doc));
            }
            res.AddMember("interactions", interactions.get());
            res.AddMember("migrations_in", migrations_in.get());
            res.AddMember("migrations_out", migrations_out.get());
            return res;
        }

    };

    round_stats m_curr_stats;
    round_stats m_total_stats;

    std::vector<std::array<uint64_t,DirCount>> m_migrate_in_by_dir;
    std::vector<std::array<uint64_t,DirCount>> m_migrate_out_by_dir;

    bool HasStepHooks() const override
    { return true; }

    void on_begin_step() override
    {
        m_curr_stats.clear();
        m_curr_stats.t=m_state->t;
        m_curr_stats.t_count=1;

        m_migrate_in_by_dir.assign(m_cells.size(), {0,0,0,0,0});
        m_migrate_out_by_dir.assign(m_cells.size(), {0,0,0,0,0});

        for(unsigned i=0; i<m_cells.size(); i++){
            assert(m_migrate_in_by_dir[i][Centre]==0);
        }

        assert(m_migrate_in_by_dir.size()==m_cells.size());

    }

    void on_bead_migrate(vec3i_t src_cell, vec3i_t dst_cell, const Bead *b) override 
    {
        assert(m_migrate_in_by_dir.size()==m_cells.size());
        assert(src_cell!=dst_cell);

        for(unsigned i=0; i<m_cells.size(); i++){
            assert(m_migrate_in_by_dir[i][Centre]==0);
        }
        auto dir=calc_interaction_dir(src_cell, dst_cell);
        assert(dir!=Centre);

        m_migrate_out_by_dir[ cell_pos_to_index(src_cell) ][All] += 1;
        m_migrate_in_by_dir[ cell_pos_to_index(dst_cell) ][All] += 1;

        m_migrate_out_by_dir[ cell_pos_to_index(src_cell) ][dir] += 1;
        m_migrate_in_by_dir[ cell_pos_to_index(dst_cell) ][dir] += 1;

        for(unsigned i=0; i<m_cells.size(); i++){
            assert(m_migrate_in_by_dir[i][Centre]==0);
        }
    }

    void on_bead_interaction(const Bead &home, const Bead &other, vec3r_t dx, double dr, vec3r_t f) override
    {
        vec3i_t h_cell=vec3_floor(home.x);
        vec3i_t o_cell=vec3_floor(other.x);

        auto dir=calc_interaction_dir(h_cell, o_cell);

        if(home.bead_id==other.bead_id){
            return;
        }

        auto &a=m_curr_stats.interactions_by_dir[dir];
        a.n += 1;
        a.hit += (dr < 1);
        a.r.add(dr);

        auto &b=m_curr_stats.interactions_all;
        b.n += 1;
        b.hit += (dr<1);
        b.r.add(dr);
    }

    void on_end_step() override
    {
        vec3r_t gv{0,0,0};
        double sum_vv=0;
        for(unsigned i=0; i<m_cells.size(); i++){
            for(int dir=0; dir<DirCount; dir++){
                if(dir==Centre){
                    assert(m_migrate_in_by_dir[i][dir]==0);
                }
                m_curr_stats.migrations_in_per_cell_by_dir[dir].add( m_migrate_in_by_dir[i][dir] );
                m_curr_stats.migrations_out_per_cell_by_dir[dir].add( m_migrate_out_by_dir[i][dir] );
            }
        
            m_curr_stats.beads_per_cell.add(m_cells[i].size());
        
            for(auto b : m_cells[i]){
                double vv=b->v.l2_norm();
                m_curr_stats.bead_velocity.add( vv );
                m_curr_stats.bead_force.add( b->f.l2_norm() );
                gv += b->v;
                sum_vv += vv;
            }
        }
        m_curr_stats.global_velocity.add(gv.l2_norm() / m_state->beads.size());
        m_curr_stats.global_temperature.add(sum_vv / m_state->beads.size());

        m_total_stats.add(m_curr_stats);

        std::ofstream dst("stats_dpd_engine_out.txt", std::ios::app);
        JSONDoc doc;
        dst<<m_curr_stats.to_json(doc.get()).pretty_print()<<"\n";
        dst<<m_total_stats.to_json(doc.get()).pretty_print()<<"\n";
        dst.flush();
    }



};

#endif
