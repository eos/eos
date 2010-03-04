/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/accumulator.hh>
#include <src/utils/lock.hh>
#include <src/utils/mutex.hh>
#include <src/utils/private_implementation_pattern-impl.hh>

namespace wf
{
    template <>
    struct Implementation<Accumulator>
    {
        Mutex mutex;

        double value;

        Implementation(const double & value) :
            value(value)
        {
        }
    };

    Accumulator::Accumulator(const double & value) :
        PrivateImplementationPattern<Accumulator>(new Implementation<Accumulator>(value))
    {
    }

    Accumulator::~Accumulator()
    {
    }

    double
    Accumulator::operator+= (const double & value)
    {
        Lock(_imp->mutex);

        return _imp->value += value;
    }

    double
    Accumulator::value() const
    {
        Lock(_imp->mutex);

        return _imp->value;
    }
}
