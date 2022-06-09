#ifndef expression_hpp
#define expression_hpp

#include <unordered_map>
#include <unordered_set>
#include <string>
#include <iostream>

#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <memory>
#include <vector>
#include <variant>
#include <cassert>

class Expression
{
public:
    using leaf_value_t = double;

    virtual bool is_constant() const
    { return false; }

    virtual std::string as_string() const =0;

    virtual void collect_variables(std::unordered_set<std::string> &params) const=0;

    virtual double evaluate(
        const std::unordered_map<std::string,leaf_value_t> &bindings
    ) const = 0;
};

struct ExprConstant
    : public Expression
{
    leaf_value_t m_value;

public:
    ExprConstant(leaf_value_t value )
        : m_value(value)
    {}

    virtual bool is_constant() const
    { return true; }

    std::string as_string() const override
    { return std::to_string(m_value); }

    void collect_variables(std::unordered_set<std::string> &params) const override
    {}

    double evaluate(
        const std::unordered_map<std::string,leaf_value_t> &bindings
    ) const override {
        return m_value;
    }
};

struct ExprVariable
    : public Expression
{
    std::string m_name;

public:
    ExprVariable(std::string name)
        : m_name(name)
    {}


    std::string as_string() const override
    { return m_name; }

    void collect_variables(std::unordered_set<std::string> &params) const override
    {
        params.insert(m_name);
    }

    double evaluate(
        const std::unordered_map<std::string,leaf_value_t> &bindings
    ) const override {
        auto it=bindings.find(m_name);
        if(it==bindings.end()){
            throw std::runtime_error("Could not find variable binding for '"+m_name+"'");
        }
        return it->second;
    }
};

struct ExprAdd
    : public Expression
{
    std::shared_ptr<Expression> m_left, m_right;
public:
    ExprAdd(std::shared_ptr<Expression> left, std::shared_ptr<Expression> right)
        : m_left(left)
        , m_right(right)
    {}

    std::string as_string() const override
    { return "(" + m_left->as_string() + "+" + m_right->as_string() + ")"; }

    void collect_variables(std::unordered_set<std::string> &params) const override
    {
        m_left->collect_variables(params);
        m_right->collect_variables(params);
    }

    double evaluate(
        const std::unordered_map<std::string,leaf_value_t> &bindings
    ) const override {
        return m_left->evaluate(bindings) + m_right->evaluate(bindings);
    }
};

struct ExprSub
    : public Expression
{
    std::shared_ptr<Expression> m_left, m_right;
public:
    ExprSub(std::shared_ptr<Expression> left, std::shared_ptr<Expression> right)
        : m_left(left)
        , m_right(right)
    {}

    std::string as_string() const override
    { return "(" + m_left->as_string() + "-" + m_right->as_string() + ")"; }

    void collect_variables(std::unordered_set<std::string> &params) const override
    {
        m_left->collect_variables(params);
        m_right->collect_variables(params);
    }

    double evaluate(
        const std::unordered_map<std::string,leaf_value_t> &bindings
    ) const override {
        return m_left->evaluate(bindings) - m_right->evaluate(bindings);
    }
};

struct ExprMul
    : public Expression
{
    std::shared_ptr<Expression> m_left, m_right;
public:

    ExprMul(std::shared_ptr<Expression> left, std::shared_ptr<Expression> right)
        : m_left(left)
        , m_right(right)
    {}

    std::string as_string() const override
    { return "(" + m_left->as_string() + "*" + m_right->as_string() + ")"; }

    void collect_variables(std::unordered_set<std::string> &params) const override
    {
        m_left->collect_variables(params);
        m_right->collect_variables(params);
    }

    double evaluate(
        const std::unordered_map<std::string,leaf_value_t> &bindings
    ) const override {
        return m_left->evaluate(bindings) * m_right->evaluate(bindings);
    }
};

struct ExprDiv
    : public Expression
{
    std::shared_ptr<Expression> m_left, m_right;
public:
    ExprDiv(std::shared_ptr<Expression> left, std::shared_ptr<Expression> right)
        : m_left(left)
        , m_right(right)
    {}


    std::string as_string() const override
    { return "(" + m_left->as_string() + "/" + m_right->as_string() + ")"; }

    void collect_variables(std::unordered_set<std::string> &params) const override
    {
        m_left->collect_variables(params);
        m_right->collect_variables(params);
    }

    double evaluate(
        const std::unordered_map<std::string,leaf_value_t> &bindings
    ) const override {
        double right=m_right->evaluate(bindings);
        if(right==0){
            throw std::runtime_error("Divide by 0 in "+as_string());
        }
        return m_left->evaluate(bindings) / right;
    }
};

std::pair<std::shared_ptr<Expression>,int> parse_expression_terms(const std::vector<std::variant<std::string,double> > &tokens, int pos);
std::pair<std::shared_ptr<Expression>,int> parse_expression_factors(const std::vector<std::variant<std::string,double> > &tokens, int pos);
std::pair<std::shared_ptr<Expression>,int> parse_expression_base(const std::vector<std::variant<std::string,double> > &tokens, int pos);



inline std::pair<std::shared_ptr<Expression>,int> parse_expression_terms(const std::vector<std::variant<std::string,double> > &tokens, int pos)
{
    auto curr=parse_expression_factors(tokens, pos);
    while(curr.second < (int)tokens.size() ){
        auto next= std::get<std::string>(tokens[curr.second]);
        //std::cerr<<"terms "<<curr.second<<" "<<next<<"\n";

        if(next==")"){
            break;
        }else if(next=="+"){
            auto right=parse_expression_factors(tokens, curr.second+1);
            assert(curr.second < right.second);
            curr={ std::make_shared<ExprAdd>( curr.first, right.first ), right.second };
        }else if(next=="-"){
            auto right=parse_expression_factors(tokens, curr.second+1);
            assert(curr.second < right.second);
            curr={ std::make_shared<ExprSub>( curr.first, right.first ), right.second };
        }else{
            throw std::runtime_error("Unexpected token in parse_expression_terms.");
        }
    }
    return curr;
}

inline std::pair<std::shared_ptr<Expression>,int> parse_expression_factors(const std::vector<std::variant<std::string,double> > &tokens, int pos)
{
    auto curr=parse_expression_base(tokens, pos);
    while(curr.second < (int)tokens.size() ){
        auto next=std::get<std::string>(tokens[curr.second]);

        //std::cerr<<"factors "<<curr.second<<" "<<next<<"\n";

        if(next==")" || next=="+" || next=="-"){
            break;
        }else if(next=="*"){
            auto right=parse_expression_base(tokens, curr.second+1);
            assert(curr.second < right.second);
            curr={ std::make_shared<ExprMul>( curr.first, right.first ), right.second };
        }else if(next=="/"){
            auto right=parse_expression_base(tokens, curr.second+1);
            assert(curr.second < right.second);
            curr={ std::make_shared<ExprDiv>( curr.first, right.first ), right.second };
        }else{
            throw std::runtime_error("Unexpected token in parse_expression_factors.");
        }
    }
    return curr;
}


inline std::pair<std::shared_ptr<Expression>,int> parse_expression_base(const std::vector<std::variant<std::string,double> > &tokens, int pos)
{
    if(std::holds_alternative<double>(tokens[pos])){
        return { std::make_shared<ExprConstant>( std::get<double>(tokens[pos]) ), pos+1 };
    }

    std::string token=std::get<std::string>(tokens[pos]);
    if(isalpha(token[0]) || token[0]=='_'){
        return { std::make_shared<ExprVariable>( token ), pos+1 };
    }

    if(token=="("){
        auto cont=parse_expression_terms(tokens, pos+1);
        auto npos=cont.second;
        if(npos >= (int)tokens.size() || std::get<std::string>(tokens[npos])!=")"){
            throw std::runtime_error("Un-matched '(' bracket.");
        }
        return {cont.first, npos+1};
    }

    throw std::runtime_error("Unexpected token in parse_expression_base, token = "+token);
}

inline std::shared_ptr<Expression> parse_expression(const std::string &src)
{
    std::vector<std::variant<std::string,double> > tokens;
    const char *curr=src.c_str(); // Gaurantees null terminated
    while( *curr ){
        auto ch=*curr;
        if(isspace(ch)){
            ++curr;
        }else if( ch=='(' || ch==')' || ch=='+' || ch=='-' || ch=='*' || ch=='/' ){
            char cc[2]={ch,0};
            tokens.push_back(std::string(cc));
            ++curr;
        }else if( isalpha(ch) || ch=='_' ){
            auto end=curr+1;
            while( isalnum(*end) || *end=='_' ){
                ++end;
            }
            tokens.push_back(std::string{curr, end});
            assert(curr < end);
            curr=end;
        }else{
            char *end;
            double v=strtod( curr, &end );
            if(curr==end){
                throw std::runtime_error("Couldn't parse expression from "+std::string(curr));
            }
            tokens.push_back( v );
            assert(curr<end);
            curr=end;
        }
    }



    auto res=parse_expression_terms(tokens, 0);
    if(res.second!=(int)tokens.size()){
        throw std::runtime_error("Unconsumed tokens.");
    }
    return res.first;
}

/*
Well-known parameters:
- ux, uy, uz : Normalised positions within grid, i.e. in [0,1)^3
- x, y, z : Absolute positions in the grid, in [0,w), [0,h), [0,d)
*/
struct Parameter
{
private:
    double m_value;
    std::shared_ptr<Expression> m_expr;

    void set_expr(std::shared_ptr<Expression> expr)
    {
        if(expr->is_constant()){
            auto p=std::dynamic_pointer_cast<ExprConstant>(expr);
            *this = p->m_value;
        }else{
            m_value=std::nan("");
            m_expr=expr;
        }
    }
public:
    Parameter()
        : m_value(0.0)
    {}

    Parameter(double x)
        : m_value(x)
    {}

    Parameter(std::shared_ptr<Expression> x)
    {
        set_expr(x);
    }

    Parameter(const std::string &str)
    {
        set_expr(parse_expression(str));
    }

    Parameter(const char *str)
    {
        set_expr(parse_expression(str));
    }

    Parameter &operator=(double x)
    {
        m_expr.reset();
        m_value=x;
        return *this;
    }

    Parameter &operator=(std::shared_ptr<Expression> x)
    {
        set_expr(x);
        return *this;
    }

    bool operator==(const Parameter &o)
    {
        if(m_expr){
            if(!m_expr){
                return false;
            }
            return m_expr->as_string()==o.m_expr->as_string();
        }else{
            if(m_expr){
                return false;
            }
            return m_value==o.m_value;
        }
    }

    operator double() const
    {
        if(__builtin_expect(!!m_expr,0)){
            throw std::runtime_error("Attempt to use expression with variables as a constant.");
        }
        return m_value;
    }

    double evaluate(const std::unordered_map<std::string,double> &x)
    {
        if(!m_expr){
            return m_value;
        }else{
            return m_expr->evaluate(x);
        }
    }

    std::string as_string() const
    {
        if(m_expr){
            return m_expr->as_string();
        }else{
            return std::to_string(m_value);
        }
    }

    bool is_constant() const
    { return !m_expr; }

    void collect_variables(std::unordered_set<std::string> &params) const
    {
        if(m_expr){
            m_expr->collect_variables(params);
        }
    }
};

inline std::ostream &operator<<(std::ostream &dst, const Parameter &p)
{
    return dst<<p.as_string();
}

#endif
