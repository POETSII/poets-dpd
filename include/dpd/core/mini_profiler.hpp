

class MiniProfiler
{
private:
    std::unordered_map<std::string,intptr_t> name_to_token;
    std::vector<std::string> token_to_name;


    struct thread_info
    {
        std::unordered_map<std::vector<uint64_t>>
    };
    
    struct section_info
    {
        std::vector<uint16_t> path;
        uint64_t count;
        uint64_t time_ns;
    };
public:
    intptr_t DeclareSection(const std::string &name);

    void EnterSection(intptr_t token);
    void ExitSection(intptr_t token);


};