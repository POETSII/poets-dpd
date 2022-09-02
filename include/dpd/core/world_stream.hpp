#ifndef world_stream_hpp
#define world_stream_hpp

#include <cstdint>
#include <vector>
#include <algorithm>
#include <functional>

#include "dpd/core/dpd_state.hpp"
#include "dpd/core/vec3.hpp"


template<class T>
struct is_blittable
{
    static const bool value = T::is_blittable;
};

template<>
struct is_blittable<uint8_t>
{ static const bool value = true; };

template<>
struct is_blittable<uint16_t>
{ static const bool value = true; };

template<>
struct is_blittable<uint32_t>
{ static const bool value = true; };

template<>
struct is_blittable<uint64_t>
{ static const bool value = true; };

template<>
struct is_blittable<float>
{ static const bool value = true; };

template<>
struct is_blittable<double>
{ static const bool value = true; };

template<>
struct is_blittable<std::string>
{ static const bool value = false; };

template<class A,class B>
struct is_blittable<std::pair<A,B>>
{
    static const bool value = is_blittable<A>::value && is_blittable<B>::value;
};

class protocol_error
    : public std::runtime_error
{
public:
    protocol_error(const std::string &msg)
        : std::runtime_error(msg)
    {}
};

struct packet_write_ctxt_t
{
    uint32_t pos;
    uint32_t max_length;
    char *buffer;

    void write(unsigned n, const void *p)
    {
        if(pos+n > max_length){
            throw std::runtime_error("Not enough space in packet.");
        }
        memcpy(buffer+pos, p, n);
        pos += n;
    }
};

inline void write(packet_write_ctxt_t &ctxt, unsigned n, const void *p)
{
    ctxt.write(n,p);
}

template<class T>
inline size_t calc_max_size(const T &x)
{
    static_assert(is_blittable<T>::value);
    return sizeof(x);
}

template<class T>
inline void write(packet_write_ctxt_t &dst, const T &x)
{
    static_assert(is_blittable<T>::value);
    write(dst, sizeof(x), &x);
}

template<class T>
inline size_t calc_max_size(const std::vector<T> &x)
{
    if(x.empty()){
        return 4;
    }
    if(is_blittable<T>::value){
        return 4+x.size()*calc_max_size(x[0]);
    }else{
        size_t res=4;
        for(const auto &xx : x){
            res += calc_max_size(xx);
        }
        return res;
    }
}

template<class T>
inline void write(packet_write_ctxt_t &dst, const std::vector<T> &x)
{
    write<uint32_t>(dst, x.size());
    if constexpr (is_blittable<T>::value){
        T dual[2];
        static_assert(sizeof(dual) == 2*sizeof(T));
        write(dst, sizeof(T)*x.size(), &x[0]);
    }else{
        for(unsigned i=0; i<x.size(); i++){
            write(dst, x[i]);
        }
    }
}

template<>
inline size_t calc_max_size<std::string>(const std::string &x)
{
    return 4+x.size();
}

inline void write(packet_write_ctxt_t &dst, const std::string &x)
{
    uint32_t size=x.size();
    write(dst, size);
    if(size){
        write(dst, size, x.data());
    }
}

template<class A, class B>
inline size_t calc_max_size(const std::pair<A,B> &x)
{
    if constexpr (is_blittable<std::pair<A,B>>::value){
        static_assert( sizeof(A)+sizeof(B) == sizeof(std::pair<A,B>) );
        return sizeof(x);
    }else{
        return calc_max_size(x.first)+calc_max_size(x.second);
    }
}

template<class A, class B>
inline void write(packet_write_ctxt_t &dst, const std::pair<A,B> &x)
{
    if constexpr (is_blittable<std::pair<A,B>>::value){
        static_assert( sizeof(A)+sizeof(B) == sizeof(std::pair<A,B>) );
        write(dst, sizeof(x), &x);
    }else{
        write(dst, x.first);
        write(dst, x.second);
    }
}

struct packet_read_ctxt_t
{
    uint32_t pos;
    uint32_t length;
    const char *buffer;

    void read(unsigned n, void *p)
    {
        if(pos+n > length){
            throw std::runtime_error("Attempt to read beyond end of packet.");
        }
        memcpy(p, buffer+pos, n);
        pos+=n;
    }
};

inline void read(packet_read_ctxt_t &ctxt, unsigned n, void *p)
{
    ctxt.read(n,p);
}

template<class T>
inline void read(packet_read_ctxt_t &src, T &res)
{
    read(src, sizeof(T), &res);
}

template<class T>
inline T read(packet_read_ctxt_t &src)
{
    T res;
    read(src, sizeof(T), &res);
    return res;
}

inline void read(packet_read_ctxt_t &ctxt, std::string &x)
{
    uint32_t size=read<uint32_t>(ctxt);
    x.resize(size);
    if(size){
        read(ctxt, size, &x[0]);
    }
}

template<class T>
inline void read(packet_read_ctxt_t &src, std::vector<T> &x)
{
    unsigned size=read<uint32_t>(src);
    x.resize(size);
    if(is_blittable<T>::value){
        read(src, x);
    }else{
        for(unsigned i=0; i<size; i++){
            x[i] = read<T>(src);
        }
    }
}

template<class T>
std::vector<char> make_packet(uint32_t code, const T &x)
{
    uint32_t size=calc_max_size(x)+8;
    std::vector<char> buffer(size);
    packet_write_ctxt_t ctxt;
    ctxt.pos=0;
    ctxt.max_length=size;
    ctxt.buffer=&buffer[0];

    write(ctxt, size);
    write(ctxt, code);
    write(ctxt, x);

    if(ctxt.pos!=size){
        throw std::runtime_error("calc_max_size lied to us.");
    }

    return buffer;
}


struct bit_writer_ctxt_t
{
private:
    packet_write_ctxt_t &dst;
    unsigned n=0;
    uint64_t acc=0;

public:
    bit_writer_ctxt_t(packet_write_ctxt_t &_dst)
        : dst(_dst)
    {}

    ~bit_writer_ctxt_t()
    {
        if(n!=0){
            write(dst, uint32_t(acc));
        }
    }


    void write_bits(unsigned w, uint32_t bits)
    {
        assert(w<=32);
        assert( (bits==32) || (bits>>w) == 0);
        assert( n < 32);
        acc |= uint64_t(bits)<<n;
        n += w;
        if(n >= 32){
            write(dst, uint32_t(acc));
            acc >>= 32;
            n -= 32;
        }
    }
};

struct bit_reader_ctxt_t
{
private:
    packet_read_ctxt_t &src;
    unsigned n=0;
    uint64_t acc=0;

public:
    bit_reader_ctxt_t(packet_read_ctxt_t &_src)
        : src(_src)
    {}

    ~bit_reader_ctxt_t()
    {
    }


    uint32_t read_bits(unsigned w)
    {
        assert(w <= 32);
        
        if(n < w){
            acc |= read<uint32_t>(src);
            n += 32;
        }

        uint32_t res=uint32_t(acc);

        acc >>= w;
        n -= w;

        return res;
    }
};





enum PacketType : uint32_t
{
    PacketType_FullWorld     =0,           // Transfers the entire world over
    PacketType_FullPositions =1        // Transfers all positions over
};


struct bead_info_t
{
    uint8_t bead_type;
    uint8_t polymer_offset;
    uint16_t polymer_id;
    uint32_t position[3];

    static const bool is_blittable = true;
};

struct world_info_t
{
    uint32_t t;
    uint32_t log2_wire_dims[3];         // Logarithm of world dimensions. Gaurantees masks can be used for wrapping
    double world_dims[3];          // Dimensions in world coordinates
    std::vector<std::string> bead_type_names;
    std::vector<bead_info_t> beads;

    static const bool is_blittable = false;
};

size_t calc_max_size(const world_info_t &x)
{
    return calc_max_size(x.bead_type_names)+calc_max_size(x.beads);
}

void write(packet_write_ctxt_t &ctxt, const world_info_t &x)
{
    write(ctxt, x.t);
    write(ctxt, sizeof(x.log2_wire_dims), x.log2_wire_dims);
    write(ctxt, sizeof(x.world_dims), x.world_dims);
    write(ctxt, x.bead_type_names );
    write(ctxt, x.beads);
}

void read(packet_read_ctxt_t &ctxt, world_info_t &x)
{
    read(ctxt, x.t);
    read(ctxt, sizeof(x.log2_wire_dims), x.log2_wire_dims);
    read(ctxt, sizeof(x.world_dims), x.world_dims);
    read(ctxt, x.bead_type_names );
    read(ctxt, x.beads);
}

/*
    Sends bead co-ordinates as bit-packed absolute co-ordinates.
    Each co-ordinate is sent as log2_wire_dims[0]+log2_wire_dims[1]+log2_wire_dims[2] bits,
    packed end to end.
*/
struct world_absolute_positions_t
{
    uint32_t t;
    uint32_t batch_size;
    std::vector<uint32_t> bits;
};

void write(packet_write_ctxt_t &dst, const world_absolute_positions_t &positions)
{
    write(dst, positions.t);
    write(dst, positions.batch_size);
    write(dst, positions.bits);
}


struct world_state_to_stream_t
{
    world_info_t info;

    const WorldState *state;
    std::vector<uint32_t> wire_index_to_bead_id;
    vec3r_t world_to_wire;
    vec3i_t wire_mask;

    void choose_dims_mapping()
    {
        static const double min_resolution=128;
        
        for(int d=0; d<3; d++){
            info.log2_wire_dims[d]=7;
            while(1){
                double resolution=ldexp( info.world_dims[d], -info.log2_wire_dims[d] );
                if(resolution >= min_resolution){
                    break;
                }
                ++info.log2_wire_dims[d];
                if(info.log2_wire_dims[d]>30){
                    throw std::logic_error("Couldn't find wire to world mapping.");
                }
            }

            wire_mask[d] = (1ul<<info.log2_wire_dims[d])-1;
            world_to_wire[d] = wire_mask[d]  / info.world_dims[d];
        }
    }

    world_state_to_stream_t(const WorldState *_state, std::function<bool(const BeadType &b)> filter_include)
    {
        state=_state;

        state->box.extract(info.world_dims);
        choose_dims_mapping();

        uint64_t bead_type_active_mask=0;
        for(const auto &bt : state->bead_types){
            if(filter_include(bt)){
                info.bead_type_names.push_back(bt.name);
                bead_type_active_mask |= 1ull<<bt.id;
            }
        }

        for(const auto &b : state->beads){
            if( (1ull<<b.bead_type) & bead_type_active_mask ){
                bead_info_t bi;
                bi.bead_type=b.bead_type;
                bi.polymer_id=b.polymer_id;
                bi.polymer_offset=b.polymer_offset;
                b.x.extract(bi.position);
                info.beads.push_back(bi);
                
                
                wire_index_to_bead_id.push_back(b.bead_id);
            }
        }
    }

    void update_bead_from_state(unsigned wire_bead_id)
    {
        auto &bi=info.beads.at(wire_bead_id);
        const auto &b = state->beads.at(wire_index_to_bead_id[wire_bead_id]);
        for(int d=0; d<3; d++){
            bi.position[wire_bead_id] = unsigned(floor(b.x[d] * world_to_wire[d] )) & wire_mask[d];
        }
    }

    void update_from_state()
    {
        for(unsigned i=0; i<info.beads.size(); i++){
            update_bead_from_state(i);
        }
    }

    void create_absolute_positions_packet(world_absolute_positions_t &res)
    {
        static const unsigned batch_size=256;

        const unsigned bits_per_pos=info.log2_wire_dims[0]+info.log2_wire_dims[1]+info.log2_wire_dims[2];
        const unsigned words_per_batch=(bits_per_pos*batch_size+31)/32;
        const unsigned num_batches=(info.beads()+batch_size-1)/batch_size;

        res.resize( num_batches*words_per_batch );

        packet_write_ctxt_t dst{#erree}

        res.t=info.t;
        for(unsigned i=0; i<info.beads.size(); i+=batch_size){
            bit_writer_ctxt_t bdst{}
        }
    }

};

/*
struct stream_to_world_info_t
{
    bool initialised=false;
    world_info_t info;

    vec3i_t wire_mask;
    vec3r_t wire_to_world;

    void on_recv(const std::vector<char> &x)
    {
        packet_read_ctxt_t ctxt;
        ctxt.buffer=&x[0];
        ctxt.length=x.size();
        ctxt.pos=0;

        uint32_t size=read<uint32_t>(ctxt);
        uint32_t type=read<uint32_t>(ctxt);

        if(size!=x.size()){
            throw protocol_error("Mismatch between header and packet size.");
        }

        if(type==PacketType_FullWorld){
            if(initialised){
                throw protocol_error("received FullWorld twice");
            }
            read(ctxt, info);
        }
    }
};
*/




struct bead_delta_t
{
    uint32_t shift : 5;
    int32_t dx : 9;
    int32_t dy : 9;
    int32_t dz : 9;

    static unsigned get_msb(unsigned x)
    {
        return x==0 ? 0 : sizeof(unsigned)*8 - __builtin_clz(x);
    }

    bead_delta_t(vec3i_t delta)
    {
        int msb=0;
        for(int d=0; d<3; d++){
            int curr=std::abs(delta[d]);
            msb=std::max<int>(msb, get_msb(curr));
        }
        shift=std::max<int>(0, msb-8);
        int bias = shift==0 ? 0 : 1<<(shift-1);
        
        dx= (delta[0]) >> shift;
        dy= (delta[1]) >> shift;
        dz= (delta[2]) >> shift;
    }

    vec3i_t to_delta()
    {
        vec3i_t res;
        int bias=shift==0 ? 0 : 1<<shift;
        res[0] = (dx<<shift);
        res[1] = (dy<<shift);
        res[2] = (dz<<shift);
        return res;
    }
};


#endif