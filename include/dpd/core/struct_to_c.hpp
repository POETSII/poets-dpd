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

struct c_init
{
    std::string name;
    // value.empty() != parts.empty(). Either it is a branch or a leaf.
    std::string value;
    std::vector<c_init> parts;
};

void render_c_init(const c_init &x, std::ostream &dst, const std::string &indent="", bool newLine=false)
{
    if(!x.value.empty()){
        dst<<indent<<x.value;  // "+x.name<<"\n";
    }else{
        dst<<indent<<"{"; // "+x.name+"\n";
        for(unsigned i=0; i<x.parts.size(); i++){
            if(i!=0){
                dst<<indent<<",";
            }
            render_c_init(x.parts[i], dst, newLine ? indent+"  " : indent, newLine);
        }
        dst<<indent<<"}";
    }
}

struct StructToCInit
{
    c_init res;

    template<class T>
    std::string prim_to_string(const T &x)
    {
        if(x==0){
            return "0";
        }else if(std::is_same<T,uint64_t>::value){
            return std::to_string(x)+"ull";
        }else if(std::is_same<T,int64_t>::value){
            return std::to_string(x)+"ll";
        }else{
            return std::to_string(x);
        }
    }

    template<class T>
    void operator()(const char *name, T &x)
    {
        if constexpr(is_primitive<T>()){
            res.parts.push_back({name, prim_to_string(x), {}});
        }else{
            StructToCInit rec{ {name, {}, {}} };
            x.walk(rec);

            res.parts.push_back(std::move(rec.res));
        }
    }

    template<class T>
    void operator()(const char *name, T *x, unsigned n)
    {
        c_init init{name, {}, {}};
        
        if constexpr(is_primitive<T>()){
            for(unsigned i=0; i<n; i++){    
                init.parts.push_back({ std::string(name)+"["+std::to_string(i)+"]", prim_to_string(x[i]), {} });
            }
        }else{
            for(unsigned i=0; i<n; i++){
                StructToCInit rec{ {std::string(name)+"["+std::to_string(i)+"]", {}, {}} };
                x[i].walk(rec);
                init.parts.push_back(std::move(rec.res));
            }
        }

        res.parts.push_back(std::move(init));
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
std::string struct_to_c_init(const T &x, bool newLine=false)
{
    StructToCInit s2c{"", {}, {}};
    const_cast<T&>(x).walk(s2c);
    std::stringstream dst;
    assert(s2c.res.value.empty());
    dst<<"{";
    if(newLine){
        dst<<"\n";
    }
    for(unsigned i=0; i<s2c.res.parts.size(); i++){
        if(i!=0){
            dst<<",";
        }
        render_c_init(s2c.res.parts[i], dst, newLine ? "    " : "", newLine);
    }
    dst<<"}";
    if(newLine){
        dst<<"\n";
    }
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
    dst<<"}\n";
    return dst.str();
}


#endif
