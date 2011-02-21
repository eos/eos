/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/likelihood.hh>
#include <src/utils/private_implementation_pattern-impl.hh>

#include <cmath>
#include <vector>

namespace eos
{
    template <>
    struct Implementation<Likelihood>
    {
        std::vector<std::tuple<ObservablePtr, double, double, double>> observables;

        Parameters parameters;

        Implementation(const Parameters & parameters) :
            parameters(parameters)
        {
        }
    };

    Likelihood::Likelihood(const Parameters & parameters) :
        PrivateImplementationPattern<Likelihood>(new Implementation<Likelihood>(parameters))
    {
    }

    Likelihood::~Likelihood()
    {
    }

    void
    Likelihood::add(const ObservablePtr & observable, const double & min, const double & central, const double & max)
    {
        if (observable->parameters() != _imp->parameters)
        {
            throw InternalError("Likelihood::add(): Encountered observable whose parameters doesn't fit ours");
        }

        _imp->observables.push_back(std::make_tuple(observable, min, central, max));
    }

    Likelihood
    Likelihood::clone() const
    {
        Parameters parameters(_imp->parameters.clone());
        Likelihood result(parameters);
        for (auto i = _imp->observables.cbegin(), i_end = _imp->observables.cend() ; i != i_end ; ++i)
        {
            result._imp->observables.push_back(std::make_tuple(std::get<0>(*i)->clone(parameters), std::get<1>(*i), std::get<2>(*i), std::get<3>(*i)));
        }

        return result;
    }

    Parameters
    Likelihood::parameters() const
    {
        return _imp->parameters;
    }

    double
    Likelihood::operator() () const
    {
        double chi_squared = 0.0;

        for (auto i = _imp->observables.cbegin(), i_end = _imp->observables.cend() ; i != i_end ; ++i)
        {
            double value = std::get<0>(*i)->evaluate();
            double chi = (value - std::get<2>(*i)) / (std::get<3>(*i) - std::get<1>(*i));

            chi_squared += chi * chi;
        }

        return std::exp(-0.5 * chi_squared);
    }
}
