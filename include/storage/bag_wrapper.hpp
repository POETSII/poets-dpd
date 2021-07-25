#ifndef bag_wrapper_hpp
#define bag_wrapper_hpp

#include <iterator>


template<class T, unsigned MAX_N>
struct bag_storage_concept
{
    T elements[MAX_N];
    uint16_t n=0;
    uint16_t lost=0;
};


template<class TStorage>
struct bag_wrapper
{
    using T = typename std::remove_reference<decltype(TStorage::elements[0])>::type;
    static constexpr size_t MAX_N = sizeof(TStorage::elements) / sizeof(TStorage::elements[0]);

    using iterator = T *;
    using const_iterator = const T *;

    TStorage &storage;

    size_t size() const
    { return storage.n; }

    void clear()
    { storage.n=0; }

    bool empty() const
    { return storage.n==0; }

    const T *begin() const
    { return &storage.elements[0]; }

    T *begin()
    { return &storage.elements[0]; }

    const T *end() const
    { return &storage.elements[storage.n]; }

    T *end()
    { return &storage.elements[storage.n]; }

    const T &back() const
    { assert(!empty()); return storage.elements[storage.n-1]; }

    T &back()
    { assert(!empty()); return storage.elements[storage.n-1]; }

    const T &front() const
    { assert(!empty()); return storage.elements[0]; }

    T &front()
    { assert(!empty()); return storage.elements[0]; }

    const T &operator[](unsigned index) const
    { assert(index<storage.n); return storage.elements[index]; }

    T &operator[](unsigned index)
    { assert(index<storage.n); return storage.elements[index]; }

    void resize(unsigned new_n)
    {
        assert(new_n <= MAX_N);
        if(new_n < MAX_N){
            storage.n=new_n;
        }else{
            assert(false);
            storage.lost += new_n-storage.n;
        }
    }

    void push_back(const T &x)
    {
        if(storage.n<MAX_N){
            // We are defensive here. Prefer to lose entries rather than overflow and corrupt memory
            storage.elements[storage.n]=x;
            storage.n=storage.n+1;
        }else{
            assert(false);
            storage.lost++;
        }
    }

    void pop_back()
    {
        assert(!empty());
        --storage.n;
    }

    void erase(iterator it)
    {
        assert( (begin() <= it) && (it < end()) );
        std::move(it+1, end(), it);
    }
};

template<class TStorage>
auto make_bag_wrapper(TStorage &stg)
{
    return bag_wrapper<TStorage>{stg};
}

#endif
