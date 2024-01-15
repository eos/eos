/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Danny van Dyk
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <eos/c-decays/dstarq-to-l-nu.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <string>

namespace eos
{
    using std::norm;

    /*
     * Decay: D^*_q -> lbar nu, based on [PS:2023A]
     */
    template <>
    struct Implementation<DstarqToLeptonNeutrino>
    {
        SpecifiedOption opt_model;

        std::shared_ptr<Model> model;

        QuarkFlavorOption opt_q;

        UsedParameter hbar;

        UsedParameter g_fermi;

        UsedParameter m_Dstarq;

        UsedParameter f_Dstarq;

        UsedParameter f_perp_Dstarq;

        UsedParameter tau_Dstarq;

        LeptonFlavorOption opt_l;

        UsedParameter m_l;

        SpecifiedOption opt_cp_conjugate;

        bool cp_conjugate;

        UsedParameter mu;

        std::function<double (const double &)> m_D_msbar;
        std::function<complex<double> ()> v_cD;
        std::function<WilsonCoefficients<ChargedCurrent> (LeptonFlavor, bool)> wc;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            opt_model(o, options, "model"),
            model(Model::make(opt_model.value(), p, o)),
            opt_q(o, options, "q"),
            hbar(p["QM::hbar"], u),
            g_fermi(p["WET::G_Fermi"], u),
            m_Dstarq(p["mass::D_" + opt_q.str() + "^*"], u),
            f_Dstarq(p["decay-constant::D_" + opt_q.str() + "^*"], u),
            f_perp_Dstarq(p["decay-constant::D_" + opt_q.str() + "^*,T"], u),
            tau_Dstarq(p["life_time::D_" + opt_q.str() + "^*"], u),
            opt_l(o, options, "l"),
            m_l(p["mass::" + opt_l.str()], u),
            opt_cp_conjugate(o, options, "cp-conjugate"),
            cp_conjugate(destringify<bool>(opt_cp_conjugate.value())),
            mu(p[opt_q.str() + "c" + "nu" + opt_l.str() + opt_l.str() + "::mu"], u)
        {
            Context ctx("When constructing D_q^*+->l^+nu observable");

            switch (opt_q.value())
            {
                case QuarkFlavor::strange:
                    m_D_msbar = [this](const double & mu) -> double { return model->m_s_msbar(mu); };
                    v_cD      = [this]() -> complex<double> { return model->ckm_cs(); };
                    wc        = [this](LeptonFlavor l, bool cp) -> WilsonCoefficients<ChargedCurrent> { return model->wet_scnul(l, cp); };
                    break;
                default:
                    throw InternalError("Invalid quark flavor: " + stringify(opt_q.value()));
            }
            u.uses(*model);
        }

        inline double beta_l() const
        {
            return 1.0 - power_of<2>(m_l() / m_Dstarq());
        }

        double decay_width() const
        {
            const WilsonCoefficients<ChargedCurrent> wc = this->wc(opt_l.value(), cp_conjugate);

            // masses
            const double m_Dstarq = this->m_Dstarq(), m_Dstarq2 = m_Dstarq * m_Dstarq;
            const double m_l      = this->m_l(),      m_l2      = m_l      * m_l;

            // decay constants
            const double f_Dstarq      = this->f_Dstarq(),      f_Dstarq2      = f_Dstarq      * f_Dstarq;
            const double f_perp_Dstarq = this->f_perp_Dstarq(), f_perp_Dstarq2 = f_perp_Dstarq * f_perp_Dstarq;

            // cf. [PS:2023A], eq. (54), p. 17, using that all operators with right-handed neutrinos do not contribute
            return power_of<2>(g_fermi * std::abs(this->v_cD()) * beta_l())
                * m_Dstarq / (24.0 * M_PI)
                * (
                    f_Dstarq2 * (m_l2 + 2.0 * m_Dstarq2) * norm(wc.cvl() + wc.cvr())
                    + 16.0 * f_perp_Dstarq2 * (2.0 * m_l2 + m_Dstarq2) * norm(wc.ct())
                );
        }

        double branching_ratio() const
        {
            return decay_width() * tau_Dstarq / hbar;
        }
    };

    const std::vector<OptionSpecification>
    Implementation<DstarqToLeptonNeutrino>::options
    {
        Model::option_specification(),
        { "cp-conjugate", { "true", "false" },  "false" },
        { "l",            { "e", "mu", "tau" }, "mu"    },
        { "q",            { "s" },              "s"     }
    };

    DstarqToLeptonNeutrino::DstarqToLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<DstarqToLeptonNeutrino>(new Implementation<DstarqToLeptonNeutrino>(parameters, options, *this))
    {
    }

    DstarqToLeptonNeutrino::~DstarqToLeptonNeutrino()
    {
    }

    double
    DstarqToLeptonNeutrino::branching_ratio() const
    {
        return _imp->branching_ratio();
    }

    double
    DstarqToLeptonNeutrino::decay_width() const
    {
        return _imp->decay_width();
    }

    const std::set<ReferenceName>
    DstarqToLeptonNeutrino::references
    {
        "DBG:2013A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    DstarqToLeptonNeutrino::begin_options()
    {
        return Implementation<DstarqToLeptonNeutrino>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    DstarqToLeptonNeutrino::end_options()
    {
        return Implementation<DstarqToLeptonNeutrino>::options.cend();
    }
}
