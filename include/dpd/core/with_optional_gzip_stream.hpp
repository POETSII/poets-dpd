#ifndef with_zip_stream_hpp
#define with_zip_stream_hpp

#include <memory>

#ifdef __GLIBCXX__
#include <ext/stdio_filebuf.h>
#endif

void with_optional_gzip_istream(
    std::string src,
    std::function<void (istream &)> cb
){
    if(src.size()>=4 && src.substr(src.size()-3)==".gz"){
#ifdef __GLIBCXX__
        std::string cmd="gunzip -k -c "+src;
        FILE *f=popen(cmd.c_str(), "r");
        if(!f){
            throw std::runtime_error("Error when spawning gunzip command '"+cmd+"'");
        }

        std::unique_ptr<FILE*> fholder(f, pclose);
        
        __gnu_cxx::stdio_filebuf<char> buf(f, std::ios_base::in);
        std::istream fs(&buf);

        cb(fs);
#else
        throw std::runtime_error("Attempt to open gzip stream, but not using GNU libc++. File="+src);
#endif
    }else{
        std::ifstream ss(src, std::ios_base::in);
        if(!ss.is_open()){
            throw std::runtime_error("Couldn't open file "+src);
        }

        cb(ss);
    }
}

void with_optional_gzip_ostream(
    std::string dst,
    std::function<void (ostream &)> cb
){
    if(dst.size()>=4 && dst.substr(dst.size()-3)==".gz"){
#ifdef __GLIBCXX__
        std::string cmd="gzip -9 -c > "+src;
        FILE *f=popen(dst.c_str(), "w");
        if(!f){
            throw std::runtime_error("Error when spawning gzip command '"+cmd+"'");
        }

        std::unique_ptr<FILE*> fholder(f, pclose);
        
        __gnu_cxx::stdio_filebuf<char> buf(f, std::ios_base::out);
        std::ostream fs(&buf);

        cb(fs);
#else
        throw std::runtime_error("Attempt to open gzip stream, but not using GNU libc++. File="+src);
#endif
    }else{
        std::ofstream ss(dst, std::ios_base::out);
        if(!ss.is_open()){
            throw std::runtime_error("Couldn't open file "+dst);
        }

        cb(ss);
    }
}

#endif