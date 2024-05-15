/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015-2024 Danny van Dyk
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

#include <eos/b-decays/b-to-l-nu.hh>
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
     * Decay: B_q -> l nubar, cf. [DBG:2013A]
     */
    template <>
    struct Implementation<BToLeptonNeutrino>
    {
        SpecifiedOption opt_model;

        std::shared_ptr<Model> model;

        QuarkFlavorOption opt_q;

        UsedParameter hbar;

        UsedParameter g_fermi;

        UsedParameter m_B;

        UsedParameter f_B;

        UsedParameter tau_B;

        LeptonFlavorOption opt_l;

        UsedParameter m_l;

        BooleanOption opt_cp_conjugate;

        UsedParameter mu;

        std::function<double (const double &)> m_U_msbar;
        std::function<complex<double> ()> v_Ub;
        std::function<WilsonCoefficients<ChargedCurrent> (LeptonFlavor)> wc;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            opt_model(o, options, "model"),
            model(Model::make(opt_model.value(), p, o)),
            opt_q(o, options, "q"),
            hbar(p["QM::hbar"], u),
            g_fermi(p["WET::G_Fermi"], u),
            m_B(p["mass::B_" + opt_q.str()], u),
            f_B(p["decay-constant::B_" + opt_q.str()], u),
            tau_B(p["life_time::B_" + opt_q.str()], u),
            opt_l(o, options, "l"),
            m_l(p["mass::" + opt_l.str()], u),
            opt_cp_conjugate(o, options, "cp-conjugate"),
            mu(p[opt_q.str() + "b" + opt_l.str() + "nu" + opt_l.str() + "::mu"], u)
        {
            Context ctx("When constructing B_q->lnu observable");

            using std::placeholders::_1;
            using std::placeholders::_2;
            switch (opt_q.value())
            {
                case QuarkFlavor::up:
                    m_U_msbar = std::bind(&ModelComponent<components::QCD>::m_u_msbar, model.get(), _1);
                    v_Ub      = std::bind(&ModelComponent<components::CKM>::ckm_ub, model.get());
                    wc        = std::bind(&ModelComponent<components::WET::UBLNu>::wet_ublnu, model.get(), _1, opt_cp_conjugate.value());
                    break;
                case QuarkFlavor::charm:
                    m_U_msbar = std::bind(&ModelComponent<components::QCD>::m_c_msbar, model.get(), _1);
                    v_Ub      = std::bind(&ModelComponent<components::CKM>::ckm_cb, model.get());
                    wc        = std::bind(&ModelComponent<components::WET::CBLNu>::wet_cblnu, model.get(), _1, opt_cp_conjugate.value());
                    break;
                default:
                    throw InternalError("Invalid quark flavor: " + stringify(opt_q.value()));
            }
            u.uses(*model);
        }

        inline double beta_l() const
        {
            return (1.0 - power_of<2>(m_l() / m_B()));
        }

        double decay_width() const
        {
            const WilsonCoefficients<ChargedCurrent> wc = this->wc(opt_l.value());

            // cf. [DBG2013], eq. (5), p. 5
            const complex<double> ga = wc.cvl() - wc.cvr();
            const complex<double> gp = wc.csl() - wc.csr();

            // masses
            const double m_B = this-> m_B(), m_B2 = m_B * m_B;
            const double m_l = this->m_l();
            const double mbatmu = model->m_b_msbar(mu);
            const double mUatmu = this->m_U_msbar(mu);

            return power_of<2>(g_fermi * std::abs(this->v_Ub()) * f_B * beta_l())
                * m_B / (8.0 * M_PI)
                * norm(ga * m_l - gp * m_B2 / (mbatmu + mUatmu));
        }

        double branching_ratio() const
        {
            return decay_width() * tau_B / hbar;
        }
    };

    const std::vector<OptionSpecification>
    Implementation<BToLeptonNeutrino>::options
    {
        Model::option_specification(),
        { "cp-conjugate", { "true", "false" },  "false" },
        { "l",            { "e", "mu", "tau" }, "mu"    },
        { "q",            { "c", "u" },         "c"     }
    };

    BToLeptonNeutrino::BToLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToLeptonNeutrino>(new Implementation<BToLeptonNeutrino>(parameters, options, *this))
    {
    }

    BToLeptonNeutrino::~BToLeptonNeutrino()
    {
    }

    double
    BToLeptonNeutrino::branching_ratio() const
    {
        return _imp->branching_ratio();
    }

    double
    BToLeptonNeutrino::decay_width() const
    {
        return _imp->decay_width();
    }

    const std::set<ReferenceName>
    BToLeptonNeutrino::references
    {
        "DBG:2013A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    BToLeptonNeutrino::begin_options()
    {
        return Implementation<BToLeptonNeutrino>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToLeptonNeutrino::end_options()
    {
        return Implementation<BToLeptonNeutrino>::options.cend();
    }
}
