/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/exclusive-b-to-s-gamma.hh>
#include <src/utils/private_implementation_pattern-impl.hh>

#include <cmath>

namespace wf
{
    template <>
    struct Implementation<BToKstarGamma>
    {
        Parameter c7;

        Parameter c7prime;

        Parameter ckm_A;

        Parameter ckm_lambda;

        Parameter ckm_etabar;

        Parameter ckm_rhobar;

        Implementation(const Parameters & p) :
            c7(p["c7"]),
            c7prime(p["c7prime"]),
            ckm_A(p["CKM::A"]),
            ckm_lambda(p["CKM::lambda"]),
            ckm_etabar(p["CKM::etabar"]),
            ckm_rhobar(p["CKM::rhobar"])
        {
        }

        double beta() const
        {
            double A = ckm_A(), lambda = ckm_lambda(), etabar = ckm_etabar(), rhobar = ckm_rhobar();
            double A2 = A * A, lambda2 = lambda * lambda, lambda4 = lambda2 * lambda2;

            Complex<double> a = Complex<double>::Cartesian(rhobar, etabar);
            Complex<double> b = 1.0 - A2 * lambda4 * a.conjugate();
            double c = std::sqrt((1.0 - A2 * lambda4) / (1.0 - lambda2));
            double d = std::pow(1.0 - A2 * lambda4 * rhobar, 2.0) + std::pow(A2 * lambda4 * etabar, 2.0);

            return -1.0 * (1.0 - a * b * (c / d)).phase();
        }

        double s_kstar_gamma() const
        {
            double r = std::abs(c7prime / c7);

            return -2.0 * r / (1.0 + r * r) * std::sin(2.0 * beta());
        }
    };

    BToKstarGamma::BToKstarGamma(const Parameters & parameters, const ObservableOptions &) :
        PrivateImplementationPattern<BToKstarGamma>(new Implementation<BToKstarGamma>(parameters))
    {
    }

    BToKstarGamma::~BToKstarGamma()
    {
    }

    double
    BToKstarGamma::s_kstar_gamma() const
    {
        return _imp->s_kstar_gamma();
    }
}
