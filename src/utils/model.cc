/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/model.hh>
#include <src/utils/private_implementation_pattern-impl.hh>

#include <cmath>

namespace wf
{
    using std::sqrt;

    Model::~Model()
    {
    }

    template <> struct Implementation<StandardModel>
    {
        /* CKM Wolfenstein parameters */
        Parameter A;
        Parameter lambda;
        Parameter rhobar;
        Parameter etabar;

        /* Measurement of V_cb is quite good, so use that */
        Parameter ckm_cb;

        Implementation(const Parameters & p) :
            A(p["CKM::A"]),
            lambda(p["CKM::lambda"]),
            rhobar(p["CKM::rhobar"]),
            etabar(p["CKM::etabar"]),
            ckm_cb(p["CKM::|V_cb|"])
        {
        }
    };

    StandardModel::StandardModel(const Parameters & p) :
        PrivateImplementationPattern<StandardModel>(new Implementation<StandardModel>(p))
    {
    }

    StandardModel::~StandardModel()
    {
    }

    complex<double>
    StandardModel::ckm_cb() const
    {
        return complex<double>(_imp->ckm_cb, 0.0);
    }

    complex<double>
    StandardModel::ckm_us() const
    {
        return complex<double>(_imp->lambda, 0.0);
    }

    complex<double>
    StandardModel::ckm_ub() const
    {
        double lambda2 = _imp->lambda * _imp->lambda, lambda3 = _imp->lambda * lambda2, lambda4 = lambda2 * lambda3;
        double A2 = _imp->A * _imp->A;
        complex<double> num = complex<double>(_imp->rhobar, -1.0 * _imp->etabar);
        complex<double> denom = 1.0 - num * A2 * lambda4;
        double rest =  _imp->A * lambda3 * sqrt((1.0 - A2 * lambda4) / (1.0 - lambda2));

        return (num / denom) * rest;
    }

    complex<double>
    StandardModel::ckm_ts() const
    {
        return complex<double>(-1.0 * _imp->A * _imp->lambda * _imp->lambda, 0.0);
    }

    complex<double>
    StandardModel::ckm_tb() const
    {
        return complex<double>(1.0, 0.0);
    }
}
