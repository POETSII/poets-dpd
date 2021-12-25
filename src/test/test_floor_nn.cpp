#include <cstdio>
#include <cmath>
#include <cfloat>

const float delta=nextafterf(0.5f, -FLT_MAX);

int tinsel_floor_nn(float x)
{
    return roundf(x-delta);
}

int main()
{
    fprintf(stderr, "delta=%.12f\n", delta);
    float x=0;
    while(x<100){
        int got=tinsel_floor_nn(x);
        int ref=floorf(x);

        if(got != ref){
            fprintf(stderr, "x=%.10g, got=%d, ref=%d\n", x, got, ref);
            exit(1);
        }

        x=nextafterf(x, FLT_MAX);
    }
}