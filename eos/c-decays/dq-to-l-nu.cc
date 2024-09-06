/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015-2023 Danny van Dyk
 * Copyright (c) 2018 Ahmet Kokulu
 * Copyright (c) 2018 Christoph Bobeth
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

#include <eos/c-decays/dq-to-l-nu.hh>
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
     * Decay: D_q -> lbar nu, based on B_q -> l nubar, cf. [DBG:2013A]
     */
    template <>
    struct Implementation<DqToLeptonNeutrino>
    {
        SpecifiedOption opt_model;

        std::shared_ptr<Model> model;

        QuarkFlavorOption opt_q;

        UsedParameter hbar;

        UsedParameter g_fermi;

        UsedParameter m_Dq;

        UsedParameter f_Dq;

        UsedParameter tau_Dq;

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
            m_Dq(p["mass::D_" + opt_q.str()], u),
            f_Dq(p["decay-constant::D_" + opt_q.str()], u),
            tau_Dq(p["life_time::D_" + opt_q.str()], u),
            opt_l(o, options, "l"),
            m_l(p["mass::" + opt_l.str()], u),
            opt_cp_conjugate(o, options, "cp-conjugate"),
            cp_conjugate(destringify<bool>(opt_cp_conjugate.value())),
            mu(p[opt_q.str() + "c" + "nu" + opt_l.str() + opt_l.str() + "::mu"], u)
        {
            Context ctx("When constructing D_q^+->l^+nu observable");

            switch (opt_q.value())
            {
                case QuarkFlavor::down:
                    m_D_msbar = [this](const double & mu) -> double { return model->m_d_msbar(mu); };
                    v_cD      = [this]() -> complex<double> { return model->ckm_cd(); };
                    wc        = [this](LeptonFlavor l, bool cp) -> WilsonCoefficients<ChargedCurrent> { return model->wet_dcnul(l, cp); };
                    break;
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
            return 1.0 - power_of<2>(m_l() / m_Dq());
        }

        double decay_width() const
        {
            const WilsonCoefficients<ChargedCurrent> wc = this->wc(opt_l.value(), cp_conjugate);

            // cf. [DBG2013], eq. (5), p. 5
            const complex<double> ga = wc.cvl() - wc.cvr();
            const complex<double> gp = wc.csl() - wc.csr();

            // masses
            const double m_Dq   = this->m_Dq(), m_Dq2 = m_Dq * m_Dq;
            const double m_l    = this->m_l();
            const double mcatmu = model->m_c_msbar(mu);
            const double mDatmu = this->m_D_msbar(mu);

            return power_of<2>(g_fermi * std::abs(this->v_cD()) * f_Dq * beta_l())
                * m_Dq / (8.0 * M_PI)
                * norm(ga * m_l - gp * m_Dq2 / (mcatmu + mDatmu));
        }

        double branching_ratio() const
        {
            return decay_width() * tau_Dq / hbar;
        }
    };

    const std::vector<OptionSpecification>
    Implementation<DqToLeptonNeutrino>::options
    {
        Model::option_specification(),
        { "cp-conjugate", { "true", "false" },  "false" },
        { "l",            { "e", "mu", "tau" }, "mu"    },
        { "q",            { "d", "s" },         ""      }
    };

    DqToLeptonNeutrino::DqToLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<DqToLeptonNeutrino>(new Implementation<DqToLeptonNeutrino>(parameters, options, *this))
    {
    }

    DqToLeptonNeutrino::~DqToLeptonNeutrino()
    {
    }

    double
    DqToLeptonNeutrino::branching_ratio() const
    {
        return _imp->branching_ratio();
    }

    double
    DqToLeptonNeutrino::decay_width() const
    {
        return _imp->decay_width();
    }

    const std::set<ReferenceName>
    DqToLeptonNeutrino::references
    {
        "DBG:2013A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    DqToLeptonNeutrino::begin_options()
    {
        return Implementation<DqToLeptonNeutrino>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    DqToLeptonNeutrino::end_options()
    {
        return Implementation<DqToLeptonNeutrino>::options.cend();
    }
}
