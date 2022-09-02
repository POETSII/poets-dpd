

#include <variant>

class Optimisable
{
public:
    virtual ~Optimisable()
    {}

    enum parameter_type
    {
        EnumType,
        IntegerRange,
        RealRange
    };

    struct parameter_info
    {
        unsigned index;
        std::string name;
        parameter_type type;
        std::vector<std::string> enumValues; 
        double minValue;  // Used both forinteger and real
        double maxValue;
    };

    void GetOptimisationParameters(
        std::vector<parameter_info> &parameters
    ) const =0;

    /*
        Returns true if the parameters were successfully set.
        If it returns false, then no parameters should be changed.
        Returning false may indicate an invalid combination, but
        where possible the parameters should be orthogonal.
    */
    bool SetOptimisationParameterValues(
        const std::vector<double> &bindings
    ) =0;

    void GetOptimisationParameterValues(
        std::vector<double> &bindings
    ) const =0;
};


class OptimisableBase
    : public Optimisable
{
private:
    struct parameter_binding
    {
        parameter_info info;
        void *pvalue;
    };
    std::vector<std::string,parameter_binding> m_bindings;
    bool m_bindings_locked=false;

    void add_binding(parameter_info &&info, void *pvalue)
    {
        if(m_bindings_locked){
            throw std::runtime_error("")
        }
        for(unsigned i=0; i<m_bindings.size(); i++){
            if(m_bindings[i].name==info.name){
                throw std::runtime_error("Duplicate parameter name");
            }
        }
        if(info.index!=m_bindings.size()){
            throw std::runtime_error("Invalid parameter index.");
        }
        m_bindings.push_back({std::move(info), pvalue});
    }
protected:
    // Adapt to new parameters. return false if they could not be set
    virtual bool on_parameters_changed(const std::vector<double> &values)=0;


    void add_parameter(const std::string &name, int &enumValue, const std::vector<std::string> &options)
    {
        if(enumValue < 0 || options.size() <= intValue){
            throw std::runtime_error("Initial valud is out of range.");
        }
        parameter_info{
            m_bindings.size(),
            name,
            EnumType,
            options,
            nanf(), nanf()
        };
        add_binding(std::move(info), &enumValue);
    }

    void add_parameter(const std::string &name, int &intValue, int minVal, int maxVal)
    {
        if(intValue < minVal || maxVal < intValue){
            throw std::runtime_error("Initial valud is out of range.");
        }
        parameter_info{
            m_bindings.size(),
            name,
            IntType,
            {},
            minVal, maxVal
        };
        add_binding(std::move(info), &intValue);
    }

    void add_parameter(const std::string &name, double &realValue, double minVal, double maxVal)
    {
        if(realValue < minVal || maxVal < realValue){
            throw std::runtime_error("Initial valud is out of range.");
        }
        parameter_info{
            m_bindings.size(),
            name,
            RealType,
            {},
            minVal, maxVal
        };
        add_binding(std::move(info), &realValue);
    }

public:
    void GetOptimisationParameters(
        std::vector<parameter_info> &parameters
    ) const
    {
        m_bindings_locked=true;
        parameters.resize(m_bindings.size());
        for(unsigned i=0; i<parameters.size(); i++){
            parameters[i]=m_bindings[i].info;
        }
    }

    bool SetOptimisationParameterValues(
        const std::vector<double> &bindings
    ){
        if(!m_bindings_locked){
            throw std::runtime_error("Attempt to set optimisation parameters without getting info first.");
        }

        std::vector<double> current;
        GetOptimisationParameterValues(current);

        if(on_parameters_changed(bindings)){
            return true;
        }

        bool success=on_parameters_changed(current);
        if(!success){
            throw std::logic_error("Restoring valid parameters failed.");
        }
        return false;
    }

    void GetOptimisationParameterValues(
        std::vector<double> &bindings
    ) const
    {
        if(!m_bindings_locked){
            throw std::runtime_error("Attempt to set optimisation parameters without getting info first.");
        }
        bindings.resize(m_bindings.size());
        for(unsigned i=0; i<m_bindings.size(); i++){
            const auto &info=m_bindings[i].info;
            switch(info.type){
            case EnumType:
            {
                int val=*(int*)m_bindings[i].pvalue;
                if(val<0 || val>= info.options.size()){
                    throw std::runtime_error("Current value is out of range.");
                }
                bindings[i]=val;
            }
            case IntType:
            {
                int val=*(int*)m_bindings[i].pvalue;
                if(val<info.minVal || val> info.maxVal){
                    throw std::runtime_error("Current value is out of range.");
                }
                bindings[i]=val;
            }
            case RealType:
            {
                double val=*(double*)m_bindings[i].pvalue;
                if(val<info.minVal || val> info.maxVal){
                    throw std::runtime_error("Current value is out of range.");
                }
                bindings[i]=val;
            }
            default:
                assert(0);
                throw std::logic_error("Invalid parameter type.");
            }
        }
    }


};