#ifndef struct_to_c_hpp
#define struct_to_c_hpp

#include <iostream>
#include <string>

template<class T>
struct primitive_info
{
    static constexpr bool is_primitive=false;
};

template<>
struct primitive_info<uint8_t>
{
    static constexpr bool is_primitive=true;
    static std::string name() { return "uint8_t"; }
};

template<>
struct primitive_info<uint16_t>
{
    static constexpr bool is_primitive=true;
    static std::string name() { return "uint16_t"; }
};

template<>
struct primitive_info<uint32_t>
{
    static constexpr bool is_primitive=true;
    static std::string name() { return "uint32_t"; }
};


template<>
struct primitive_info<uint64_t>
{
    static constexpr bool is_primitive=true;
    static std::string name() { return "uint64_t"; }
};

template<>
struct primitive_info<int8_t>
{
    static constexpr bool is_primitive=true;
    static std::string name() { return "int8_t"; }
};

template<>
struct primitive_info<int16_t>
{
    static constexpr bool is_primitive=true;
    static std::string name() { return "int16_t"; }
};

template<>
struct primitive_info<int32_t>
{
    static constexpr bool is_primitive=true;
    static std::string name() { return "int32_t"; }
};

template<>
struct primitive_info<int64_t>
{
    static constexpr bool is_primitive=true;
    static std::string name() { return "int64_t"; }
};


template<>
struct primitive_info<float>
{
    static constexpr bool is_primitive=true;
    static std::string name() { return "float"; }
};

template<class T>
constexpr bool is_primitive()
{ return primitive_info<T>::is_primitive; }

template<class T>
std::string primitive_name()
{ return primitive_info<T>::name();  }

template<class T>
struct is_array
{
    static constexpr bool value()
    { return false;  }
};

template<class T, size_t N>
struct is_array<T[N]>
{
    static constexpr bool value()
    { return true;  }
};

struct StructToC
{
    std::ostream &dst;
    std::string indent;

    template<class T>
    void operator()(const char *name, T x)
    {
        if constexpr(is_primitive<T>()){
            (void)x;
            dst<<indent<<primitive_name<T>()<<" "<<name<<";\n";
        }else{
            dst<<indent<<"struct {\n";
            StructToC rec{dst, indent+"    "};
            x.walk(rec);
            dst<<indent<<"} "<<name<<";\n";
        }
    }

    template<class T>
    void operator()(const char *name, T *x, unsigned n)
    {
        if constexpr(is_primitive<T>()){
            (void)x;
            dst<<indent<<primitive_name<T>()<<" "<<name<<"["<<n<<"];\n";
        }else{
            dst<<indent<<"struct {\n";
            StructToC rec{dst, indent+"    "};
            x[0].walk(rec);
            dst<<indent<<"} "<<name<<"["<<n<<"];\n";
        }
    }
};

template<class T>
std::string struct_to_c_body()
{
    T x;
    std::stringstream dst;
    StructToC s2c{dst, "  "};
    x.walk(s2c);
    return dst.str();
}

template<class T>
std::string struct_to_c(const std::string &name)
{
    T x;
    std::stringstream dst;
    dst<<"struct "<<name<<"{\n";
    StructToC s2c{dst, "  "};
    x.walk(s2c);
    dst<<"};\n";
    return dst.str();
}

#endif
