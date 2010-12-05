/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/random_number_engine.hh>

#include <random>

#include <iostream>

namespace eos
{
    template <>
    struct Implementation<RandomNumberEngine>
    {
        std::mt19937 engine;
    };

    RandomNumberEngine::RandomNumberEngine() :
        PrivateImplementationPattern<RandomNumberEngine>(new Implementation<RandomNumberEngine>)
    {
    }

    RandomNumberEngine::~RandomNumberEngine()
    {
    }

    double
    RandomNumberEngine::operator() (void)
    {
        return _imp->engine() / static_cast<double>((1L << 32) - 1);
    }
}

