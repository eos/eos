/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/random_number_engine.hh>

#include <tr1/random>

namespace eos
{
    template <>
    struct Implementation<RandomNumberEngine>
    {
        std::tr1::mt19937 engine;
    };

    RandomNumberEngine::RandomNumberEngine() :
        PrivateImplementationPattern<RandomNumberEngine>(new Implementation<RandomNumberEngine>)
    {
    }

    RandomNumberEngine::~RandomNumberEngine()
    {
    }

    unsigned int
    RandomNumberEngine::operator() (void)
    {
        return _imp->engine();
    }
}

