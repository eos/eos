/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_UTILS_MEMOISE_HH
#define WFITTER_GUARD_SRC_UTILS_MEMOISE_HH 1

#include <src/utils/instantiation_policy.hh>
#include <src/utils/instantiation_policy-impl.hh>
#include <src/utils/lock.hh>
#include <src/utils/mutex.hh>

#include <tuple>
#include <unordered_map>

/* We are using a custom std::hash for tuple<FunctinPtr, double, ..., double> */
namespace std
{
    template <typename U_> static unsigned long __hash_one(const U_ & u)
    {
        static_assert(sizeof(U_) == sizeof(double), "Need to specialize __hash_one for sizeof(U_) != sizeof(double)");
        return *reinterpret_cast<const unsigned long *>(&u);
    }

    template <unsigned n_, typename ... T_>
    struct __TupleHasher
    {
        static unsigned long hash(const std::tuple<T_ ...> & t)
        {
            return __TupleHasher<n_ - 1, T_ ...>::hash(t) ^ __hash_one<decltype(std::get<n_>(t))>(std::get<n_>(t));
        }
    };

    template <typename ... T_>
    struct __TupleHasher<0, T_ ...>
    {
        static unsigned long hash(const std::tuple<T_ ...> & t)
        {
            return __hash_one<decltype(std::get<0>(t))>(std::get<0>(t));
        }
    };

    template <typename ... T_>
    struct hash<std::tuple<T_ ...>>
    {
        size_t operator() (const std::tuple<T_ ...> & t) const
        {
            return __TupleHasher<sizeof...(T_) - 1, T_ ...>::hash(t);
        }
    };
}

namespace wf
{
    namespace implementation
    {
        template <typename T_> struct ResultOf;

        template <typename Result_, typename Class_, typename ... Args_>
        struct ResultOf<Result_ (Class_::*) (Args_ ...)>
        {
            typedef Result_ Type;
        };

        template <typename Result_, typename ... Args_>
        struct ResultOf<Result_ (*) (Args_ ...)>
        {
            typedef Result_ Type;
        };
    }

    template <typename Result_, typename ... Params_>
    class Memoiser :
        public InstantiationPolicy<Memoiser<Result_, Params_ ...>, Singleton>
    {
        public:
            Mutex * const mutex;
            typedef Result_ (*FunctionType)(const Params_ & ...);
            typedef std::tuple<FunctionType, Params_...> KeyType;

        private:
            std::unordered_map<KeyType, Result_> memoisations;

        public:
            Memoiser() :
                mutex(new Mutex)
            {
            }

            ~Memoiser()
            {
                delete mutex;
            }

            Result_ operator() (const FunctionType & f, const Params_ & ... p)
            {
                Lock l(*mutex);

                KeyType key(f, p ...);
                auto i = memoisations.find(key);

                if (memoisations.end() != i)
                    return i->second;

                Result_ result = f(p ...);
                memoisations.insert(std::pair<KeyType, Result_>(key, result));

                return result;
            }
    };

    template <typename FunctionType_, typename ... Params>
    typename implementation::ResultOf<FunctionType_>::Type memoise(FunctionType_ f, const Params & ... p)
    {
        return (*Memoiser<typename implementation::ResultOf<FunctionType_>::Type, Params ...>::instance())(f, p ...);
    }
}

#endif
