#include "dpd/core/expression.hpp"

void check(std::string src, double value, const std::unordered_map<std::string,double> &bindings)
{
    std::cerr<<"expr = "<<src<<"\n";
    auto e=parse_expression(src);
    std::cerr<<"     = "<<e->as_string()<<"\n";
    auto got=e->evaluate(bindings);
    std::cerr<<"     = "<<got<<"\n";
    if( std::abs(value-got) > 0.0000001 ){
        throw std::runtime_error("Value is wrong.");
    }
}

int main()
{
    check("0", 0.0, {});
    check("0.0", 0.0, {});
    check("10", 10.0, {});
    check("10.0", 10.0, {});
    check("10000.0", 10000.0, {});

    check("x", 10, {{"x",10}});
    check("y", -10, {{"y",-10}});

    check("3+4", 7, {});
    check("3-4", -1, {});
    check("3*4", 3*4, {});
    check("3/4", 3/4.0, {});

    check("3+3+3", 9, {});
    check("3*3*3", 27, {});
    check("3-3-3", -3, {});
    check("3/3/3", 1.0/3, {});

    check("3+4*5", 23, {});
    check("3*4+5", 17, {});

    check("3*4/5", 3*4.0/5, {});
    check("3/4*5", 3/4.0*5, {});

    check("3+x*5", 23, {{"x",4}});
    check("3*x+5", 17, {{"x",4}});

    check("(3+x)*5", (3+4)*5, {{"x",4}});
    check("3+(x*5)", 3+(4*5), {{"x",4}});
    check("(3*x+5)", (3*4+5), {{"x",4}});
    check("(3*(x)+5)", (3*(4)+5), {{"x",4}});
    check("(3.1*(x)+5)", (3.1*(4)+5), {{"x",4}});
    check("(3.1*(x)+(y))", (3.1*(4)+11.2), {{"x",4},{"y",11.2}});
    check("(3.1*((x)+(y)))", (3.1*((4)+11.2)), {{"x",4},{"y",11.2}});

    fprintf(stdout, "Success\n");
}