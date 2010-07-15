/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/rare-b-decays/inclusive-b-to-s-gamma.hh>
#include <src/utils/concrete_observable.hh>
#include <src/utils/integrate.hh>
#include <src/utils/kinematic.hh>
#include <src/utils/model.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/qcd.hh>

namespace wf
{
    /* Minimal */

    template <>
    struct Implementation<BToXsGamma<Minimal>>
    {
        std::shared_ptr<Model> model;

        Parameter c7;

        Parameter m_b_MSbar;

        Parameter br_bcsl;

        Parameter uncertainty;

        Implementation(const Parameters & p) :
            model(new StandardModel(p)),
            c7(p["c7"]),
            m_b_MSbar(p["mass::b(MSbar)"]),
            br_bcsl(p["exp::BR(B->X_clnu)"]),
            uncertainty(p["B->X_sgamma::uncertainty"])
        {
        }

        double m_b_pole() const
        {
            // Return the m_b pole mass
            return QCD::mb_pole(m_b_MSbar);
        }

        double m_c_pole() const
        {
            return 1.6;
        }

        double branching_ratio() const
        {
            static const double sm = 3.15e-4;
            static const double sm_delta = 0.23e-4;
            static const double c7sm = -0.3;
            static const double alpha_e = 1.0 / 133.0; // cf. [BHP2008]

            double m_c_hat = m_c_pole() / m_b_pole();
            double z = pow(m_c_hat, 2);
            double z2 = pow(z, 2), z3 = z * z2, z4 = z3 * z, lnz = log(z);

            // cf. [BMU1999], Eq. (46), p. 16
            double g = 1.0 - 8.0 * z + 8.0 * z3 - z4 - 12.0 * z2 * lnz;
            double kappa = 1.0 - 2.0/3.0 * QCD::alpha_s(m_b_pole()) / M_PI * (1.5 + (M_PI * M_PI - 31.0 / 4.0) * pow(1.0 - m_c_hat, 2));

            double ckm = norm(model->ckm_tb() * conj(model->ckm_ts()) / model->ckm_cb());
            double c7np = c7 - c7sm;

            return (sm + sm_delta * uncertainty)
                + 6.0 * alpha_e / M_PI * br_bcsl * ckm / g / kappa * (pow(c7np, 2) + 2.0 * c7sm * c7np);
        }
    };

    BToXsGamma<Minimal>::BToXsGamma(const Parameters & parameters, const ObservableOptions &) :
        PrivateImplementationPattern<BToXsGamma<Minimal>>(new Implementation<BToXsGamma<Minimal>>(parameters))
    {
    }

    BToXsGamma<Minimal>::~BToXsGamma()
    {
    }

    double
    BToXsGamma<Minimal>::integrated_branching_ratio() const
    {
        return _imp->branching_ratio();
    }
}
