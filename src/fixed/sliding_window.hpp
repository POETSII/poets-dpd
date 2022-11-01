#include "fixed_math_v2.hpp"

const int MAX_BEADS_PER_CELL = 4;

struct cell
{
    ap_uint<4> n;
    fixed_maths_v2::bead_info_t beads[MAX_BEADS_PER_CELL];
};

template<int N>
struct interact_many
{
    static void execute(
        const fixed_maths_v2::bead_info_t &hb,
        const fixed_maths_v2::bead_info_t *obs, 
        ap_fixed_vec3<2,1> ob_cell_delta,
        const ap_uint<64> &round_seed,
        fixed_maths_v2::force_t &f
    ){
        fixed_maths_v2::force_t f1, f2;
        interact_many<N/2>::execute(
            hb, obs, ob_cell_delta, round_seed, f1
        );
        interact_many<N-N/2>::execute(
            hb, obs+N/2, ob_cell_delta, round_seed, f2
        );
        #pragma HLS UNROLL
        for(int i=0; i<3; i++){
            f[i] = f1[i] + f2[i];
        }
    }
};

template<>
struct interact_many<1>
{
    static void execute(
        const fixed_maths_v2::bead_info_t &hb,
        const fixed_maths_v2::bead_info_t *obs, 
        ap_fixed_vec3<2,1> ob_cell_delta,
        const ap_uint<64> &round_seed,
        fixed_maths_v2::force_t &f
    ){
        calc_force_fix_v2(
            hb, obs[0],
            ob_cell_delta,
            round_seed,
            f
        );
    }
};

struct slider_1d
{
    cell prev, curr, next;

    void push(
        const cell &in,
        const ap_uint<64> &round_seed,
        std::array<fixed_maths_v2::force_t,MAX_BEADS_PER_CELL> &f
    ){
#pragma HLS PIPELINE II=1

        prev=curr;
        curr=next;
        next=in;

        for(int i=0; i<MAX_BEADS_PER_CELL; i++){
            fixed_maths_v2::force_t f_prev, f_curr, f_next;
            interact_many<MAX_BEADS_PER_CELL>::execute(
                prev.beads[i], &prev.beads[0], {-1,0,0}, round_seed, f_prev
            );
            interact_many<MAX_BEADS_PER_CELL>::execute(
                curr.beads[i], &curr.beads[0], {0,0,0},  round_seed,f_curr
            );
            interact_many<MAX_BEADS_PER_CELL>::execute(
                next.beads[i], &next.beads[0], {+1,0,0},  round_seed,f_next
            );

            for(int d=0; d<3; d++){
                f[i][d]=f_prev[d] + f_curr[d] + f_next[d];
                
            }
        };
    }
};

void slider_1d_push(
    const cell &in,
    const ap_uint<64> &round_seed,
    fixed_maths_v2::force_t &f_0,
    fixed_maths_v2::force_t &f_1,
    fixed_maths_v2::force_t &f_2,
    fixed_maths_v2::force_t &f_3
){
#pragma HLS DATA_PACK variable=in
#pragma HLS INTERFACE ap_ctrl_none port=return
#pragma HLS PIPELINE

    static slider_1d inst;
    std::array<fixed_maths_v2::force_t,MAX_BEADS_PER_CELL> f;
    inst.push(in, round_seed, f);
    f_0=f[0];
    f_1=f[1];
    f_2=f[2];
    f_3=f[3];
}


struct slider_2d_point
{
	/*     | +-+  +-+
	 *    C1 | C2 | C3
	 *    C4 | C5 | C6
	 *    C7 | C8 | C9
	       +-+  +-+
	 */

	cell C1, C2, C3, C4, C5, C6, C7, C8, C9;
	cell Q1[64];
	cell Q2[64];
	ap_uint<6> q_index=0;

    void push(
        ap_uint<4> in_index,
    	cell &in,
        const ap_uint<64> &round_seed,
        fixed_maths_v2::force_t &f
    ){
#pragma HLS PIPELINE II=1

    	if(in_index==0){
			C9 = C6;
			C6 = C3;
			C3  = Q2[q_index];
			Q2[q_index] = C8;
			C8 = C5;
			C5 = C2;
			C2 = Q1[q_index];
			Q1[q_index] = C7;
			C7 = C4;
			C4 = C1;
			C1 = in;
			++q_index;
    	}

    	auto sel = C5.beads[in_index];

        fixed_maths_v2::force_t f1, f2, f3, f4, f5, f6, f7, f8, f9;
		interact_many<MAX_BEADS_PER_CELL>::execute(
			sel, &C1.beads[0], {-1,-1,0}, round_seed, f1
		);
		interact_many<MAX_BEADS_PER_CELL>::execute(
			sel, &C2.beads[0], {0,-1,0}, round_seed, f2
		);
		interact_many<MAX_BEADS_PER_CELL>::execute(
			sel, &C3.beads[0], {+1,-1,0}, round_seed, f3
		);
		interact_many<MAX_BEADS_PER_CELL>::execute(
			sel, &C4.beads[0], {-1,0,0}, round_seed, f4
		);
		interact_many<MAX_BEADS_PER_CELL>::execute(
			sel, &C5.beads[0], {0,0,0}, round_seed, f5
		);
		interact_many<MAX_BEADS_PER_CELL>::execute(
			sel, &C6.beads[0], {+1,0,0}, round_seed, f6
		);
		interact_many<MAX_BEADS_PER_CELL>::execute(
			sel, &C7.beads[0], {-1,+1,0}, round_seed, f7
		);
		interact_many<MAX_BEADS_PER_CELL>::execute(
			sel, &C8.beads[0], {0,+1,0}, round_seed, f8
		);
		interact_many<MAX_BEADS_PER_CELL>::execute(
			sel, &C9.beads[0], {+1,+1,0}, round_seed, f9
		);


		for(int d=0; d<3; d++){
			f[d]=f1[d]+f2[d]+f3[d]+f4[d]+f5[d]+f6[d]+f7[d]+f8[d]+f9[d];
		}
    }
};

void slider_2d_push_point(
		ap_uint<4> in_index,
		cell &in,
		const ap_uint<64> &round_seed,
		fixed_maths_v2::force_t &out
){
#pragma HLS DATA_PACK variable=in
#pragma HLS INTERFACE ap_ctrl_none port=return
#pragma HLS PIPELINE

    static slider_2d_point inst;
    inst.push(in_index, in, round_seed, out);
}


struct slider_3d_point
{
	/*     | +-+  +-+
	 *    C1 | C2 | C3
	 *    C4 | C5 | C6
	 *    C7 | C8 | C9
	       +-+  +-+
	 */


	cell C[27];
	cell Q2324[32], Q2021[32], Q1415[32], Q1112[32], Q67[32], Q23[32];
	cell Q89[32*32], Q1718[32*32];
	ap_uint<5> q_short;
	ap_uint<10> q_long;

    void push(
        ap_uint<4> in_index,
    	cell &in,
        const ap_uint<64> &round_seed,
        fixed_maths_v2::force_t &f
    ){
#pragma HLS PIPELINE II=1

    	if(in_index==0){
    		C[26]=C[25];
			C[25]=C[24];
			C[24]=Q2324[q_short];
			Q2324[q_short]=C[23];
			C[23]=C[22];
			C[22]=C[21];
			C[21]=Q2021[q_short];
			Q2021[q_short]=C[20];
			C[20]=C[19];
			C[19]=C[18];
			C[18]=Q1718[q_long];

    		Q1718[q_long]=C[17];
    		C[17]=C[16];
    		C[16]=C[15];
    		C[15]=Q1415[q_short];
    		Q1415[q_short]=C[14];
    		C[14]=C[12];
    		C[13]=C[12];
    		C[12]=Q1112[q_short];
    		Q1112[q_short]=C[11];
    		C[11]=C[10];
    		C[10]=C[9];
    		C[9]=Q89[q_long];

    		Q89[q_long]=C[8];
    		C[8]=C[7];
    		C[7]=C[6];
    		C[6]=Q67[q_short];
    		Q67[q_short]=C[5];
    		C[5]=C[4];
    		C[4]=C[3];
    		C[3]=Q23[q_short];
    		Q23[q_short]=C[2];
    		C[2]=C[1];
    		C[1]=C[0];
    		C[0]=in;

    		q_short++;
    		q_long++;
    	}

    	auto sel = C[13].beads[in_index];

        fixed_maths_v2::force_t ftot{0,0,0};
#pragma HLS UNROLL
        for(int i=0; i<27; i++){
        	const int d0=(i/9)-1;
        	const int d1=((i/3)%3)-1;
        	const int d2=(i%3)-1;

        	fixed_maths_v2::force_t floc;
    		interact_many<MAX_BEADS_PER_CELL>::execute(
    			sel, &C[i].beads[0], {d0,d1,d2}, round_seed, floc
    		);
#pragma HLS UNROLL
    		for(int d=0; d<3; d++){
    			ftot[d] += floc[d];
    		}
        }

#pragma HLS UNROLL
        for(int d=0; d<3; d++){
			f[d] = ftot[d];
		}
    }
};

void slider_3d_push_point(
		ap_uint<4> in_index,
		cell &in,
		const ap_uint<64> &round_seed,
		fixed_maths_v2::force_t &out
){
#pragma HLS DATA_PACK variable=in
#pragma HLS INTERFACE ap_ctrl_none port=return
#pragma HLS PIPELINE

    static slider_3d_point inst;
#pragma HLS ARRAY_PARTITION variable=inst.C complete dim=1
    inst.push(in_index, in, round_seed, out);
}



