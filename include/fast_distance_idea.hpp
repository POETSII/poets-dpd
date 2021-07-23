

#include <cstdint>
#include <vector>

/*
Assume we are simulating a 3d box with each side of
length "length", where length is a positive integer. So co-ordinates
within that space will be real numbers in [0,length)^3.

A big part of DPD is working out which beads are close
enough to interact - specifically, those which are within
unit distance of each other. For this version of DPD
we split the world into unit-length cells, and only
look for interactions between beads in adjacent
cells. There are 27 cells in each cell's neighbourhood,
and a lot of beads in the neighbourhood are more than
distance one from the beads in the central cell.

The distance check between beads in the current cell
and beads in the neighbours is relatively expensive,
due to wrapping, number of of memory accesses, and
floating-point. So if we can much more cheaply reject
beads more than 1+eps distance away, we only need
to do all the calculations for thos remaining.

For fast rejection we quantise the world again at a resolution
of 0.25, so each cell is split into a 4x4x4 sub-box. This
gives a 30-bit packed location. If we record each bead's
packed location then we can do fast rejection whenever
we can prove that they can't be within distance 1. This
is particularly effective for "corner" and "edge" neighbours,
where it is quite likely that beads are more than 1.25
away. However, the coarse nature of 0.25 means we must
be conservative.

It would be nice to use a resolution of 0.125, but that
means we can't do a 256^3 space, which Julian would like.

Another option is to do 11-bit in two axis, and 10-bit
in another axis.

*/

struct Bead
{
    uint32_t x_packed;  // Packed position
    float x[3];         // Full position
    // Other stuff
};

struct Cell
{
    std::vector<Bead> beads;  // Conceptually a vector
};


/* This is the hot-path of the whole system. */
void on_receive(Cell &c, Bead *b)
{
    float a_x[3];

    // We always have to calculate wrapped position
    for(int i=0; i<3; i++){
        a_x[i] = wrap_length(b->x[i], i);
    }

    for(Bead *a : c.beads){
        ////////////////////////////////////////////////////////////////
        // Do the full interaction test in floating-point.
        float dx[3];
        dx[0] = a->x[0] - b->x[0];
        dx[1] = a->x[1] - b->x[1];
        dx[2] = a->x[2] - b->x[2];
        float r_sqr = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];

        if(r_sqr >= 1.0f){
            continue; // Low probability of rejection
        }

        //////////////////////////////////////////////////
        // Full DPD interaction.
        float r = sqrt(f);
        //...
        // About 50-60 more instructions
    }
}


void on_receive_fast_reject(Cell &c, Bead *b)
{
    uint32_t b_pos_packed = b->pos_packed;
    bool a_x_wrapped = false;
    float a_x[3];
    for(Bead *a : c.beads){
        uint32_t a_pos_packed = a->pos_packed;
        if( fast_reject(b_pos_packed, a_pos_packed) ){
            continue; // Large probability of rejection
        }

        ////////////////////////////////////////////////////////////////
        // Do the full interaction test in floating-point.

        // We only wrap _if_ any of the points are close enough. For 
        // some cells (e.g. corner cells), it is entirely possible that
        // no beads touch. For others it is only the cost of one/two
        // instructions to skip.
        if(!a_x_wrapped){
            for(int i=0; i<3; i++){
                a_x[i] = wrap_length(b->x[i], i);
            }
        }

        // Calculate dx and r. It is very likely we accept, and so this
        // is useful compute.
        float dx[3];
        dx[0] = a->x[0] - b->x[0];
        dx[1] = a->x[1] - b->x[1];
        dx[2] = a->x[2] - b->x[2];
        float r_sqr = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];

        if(r_sqr >= 1.0f){
            continue; // Low probability of rejection
        }

        //////////////////////////////////////////////////
        // Full DPD interaction.
        float r = sqrt(f);
        //...
    }
}

// Proposed new instruction. It calculate the dot product
// between two absolute positions, which encode the
bool fast_reject(uint32_t aa, uint32_t bb)
{
    // Global constant. Could be put in CSR?
    const uint32_t length_cells=?;

    // This could be compiled into hardware. Needs to be
    // chosen to be conservative
    const int32_t thresh=?;

    // Bit extraction. Technically a_cell is just the MSBs of a.
    int16_t a[3], b[3];
    uint8_t a_cell[3], b_cell[3];
    for(int i=0; i<3; i++){
        a[i] = (aa>>(i*10) ) &0x3FF;
        a_cell[i] = (aa>>(i*10+2) ) &0xFF;
    }

    // calculate wrapping
    int16_t offset[3];
    for(int i=0; i<3; i++){
        if(a_cell[i]==0 && b_cell[i]==length_cells-1){
            offset[i] = length_cells<<2;
        }else if(a_cell[i]==length_cells-1 && b_cell[i]==0){
            offset[i] = (-(int16_t)length_cells)<<2;
        }else{
            offset[i] = 0;
        }
    }

    // Calculate distance and wrap
    int16_t dab[3];
    for(int i=0; i<3; i++){
        dab[i] = a[i] - b[i] + offset[i];
    }

    // Dot product
    int32_t dist_sqr=0;
    for(int i=0; i<3; i++){
        dist_sqr += dab[i] * int32_t(dab[i]);
    }

    return (dist_sqr <= thresh);
}

