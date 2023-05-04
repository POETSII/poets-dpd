#ifndef dpd_core_json_helper_hpp
#define dpd_core_json_helper_hpp

#include <rapidjson/document.h>
#include "rapidjson/prettywriter.h"

#include <string>

class JSON
{
private:
    rapidjson::Document &doc;
    rapidjson::Value val;
public:
    JSON() = delete;
    JSON(const JSON &) = delete;
    JSON & operator()(const JSON &x) = delete;


    JSON(rapidjson::Document &_doc)
        : doc(_doc)
    {
    }

    JSON(JSON &&_val)
        : doc(_val.doc)
        , val(std::move(_val.val))
    {
    }

    JSON(rapidjson::Value &&_val, rapidjson::Document &_doc)
        : doc(_doc)
        , val(std::move(_val))
    {
    }

    rapidjson::Value &&get()
    {
        return std::move(val);
    }

    rapidjson::Value Wrap(rapidjson::Value &&x)
    { return std::move(x); }

    rapidjson::Value Wrap(const std::string &x)
    { return rapidjson::Value(x.c_str(), doc.GetAllocator()); }

    rapidjson::Value Wrap(uint64_t x)
    { return rapidjson::Value(x); }

    rapidjson::Value Wrap(double x)
    { return rapidjson::Value(x); }

    template<class T>
    rapidjson::Value Wrap(std::vector<T> x)
    {
        rapidjson::Value res;
        res.SetArray();
        res.Reserve(x.size(), doc.GetAllocator());
        for(const auto &e : x){
            res.PushBack( Wrap(e), doc.GetAllocator() );
        }
        return res;
    }


    template<class T>
    void AddMember(const std::string &name, const T &value)
    {
        if(val.IsNull()){
            val.SetObject();
        }else{
            assert(val.IsObject());
        }
        val.AddMember(Wrap(name), Wrap(value), doc.GetAllocator());
    }

    void AddMember(const std::string &name, rapidjson::Value x)
    {
        if(val.IsNull()){
            val.SetObject();
        }else{
            assert(val.IsObject());
        }
        val.AddMember(Wrap(name), x, doc.GetAllocator());
    }

    void AddMember(const std::string &name, JSON &&x)
    {
        if(val.IsNull()){
            val.SetObject();
        }else{
            assert(val.IsObject());
        }
        val.AddMember(Wrap(name), std::move(x.val), doc.GetAllocator());
    }

    std::string pretty_print() const
    {
        rapidjson::StringBuffer sb;
        rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(sb);
        val.Accept(writer);    // Accept() traverses the DOM and generates Handler events.
        return sb.GetString();
    }
};

class JSONDoc
    : public JSON
{
private:
    rapidjson::Document thedoc;
public:
    JSONDoc()
        : JSON(thedoc)
    {
    }

    rapidjson::Document &get()
    {
        return thedoc;
    }
};

#endif
