/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016-2025 Danny van Dyk
 * Copyright (c) 2021      Méril Reboud
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

#include <eos/form-factors/mesonic.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/rare-b-decays/b-to-kstar-charmonium.hh>
#include <eos/nonlocal-form-factors/nonlocal-formfactors.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <cmath>
#include <complex>

namespace eos
{
    using std::abs;
    using std::arg;
    using std::conj;
    using std::norm;
    using std::real;
    using std::sqrt;
    using namespace std::literals::string_literals;

    /*!
     * Implementation for the decay @f$\bar{B} \to \bar{K}^* \psi@f$.
     */
    template <>
    struct Implementation<BToKstarCharmonium>
    {
        UsedParameter g_fermi;

        UsedParameter hbar;

        std::shared_ptr<Model> model;

        QuarkFlavorOption opt_q;

        UsedParameter m_B;

        UsedParameter tau_B;

        UsedParameter m_Kstar;

        SwitchOption opt_nonlocal_formfactor;

        NonlocalFormFactorPtr<PToV> nonlocal_formfactor;

        SwitchOption opt_psi;

        UsedParameter m_psi;

        UsedParameter f_psi;

        std::function<complex<double> ()> residue_H_long;
        std::function<complex<double> ()> residue_H_perp;
        std::function<complex<double> ()> residue_H_para;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            g_fermi(p["WET::G_Fermi"], u),
            hbar(p["QM::hbar"], u),
            model(Model::make(o.get("model"_ok, "SM"), p, o)),
            opt_q(o, options, "q"_ok),
            m_B(p["mass::B_" + opt_q.str()], u),
            tau_B(p["life_time::B_" + opt_q.str()], u),
            m_Kstar(p["mass::K_" + opt_q.str() + "^*"], u),
            opt_nonlocal_formfactor(o, "nonlocal-formfactor"_ok, { "GvDV2020", "naive", "GRvDV2022order5" }, "GvDV2020"),
            nonlocal_formfactor(NonlocalFormFactor<PToV>::make("B->K^*::" + opt_nonlocal_formfactor.value(), p, o)),
            opt_psi(o, "psi"_ok, { "J/psi", "psi(2S)" }, "J/psi"),
            m_psi(p["mass::" + opt_psi.value()], u),
            f_psi(p["decay-constant::" + opt_psi.value()], u)
        {
            Context ctx("When constructing B->K^*psi observables");

            if (! nonlocal_formfactor.get())
                throw InternalError("Cannot construct the nonlocal formfactor");

            if ("J/psi" == opt_psi.value())
            {
                residue_H_long = std::bind(&NonlocalFormFactor<PToV>::H_long_residue_jpsi, nonlocal_formfactor);
                residue_H_perp = std::bind(&NonlocalFormFactor<PToV>::H_perp_residue_jpsi, nonlocal_formfactor);
                residue_H_para = std::bind(&NonlocalFormFactor<PToV>::H_para_residue_jpsi, nonlocal_formfactor);
            }
            else
            {
                residue_H_long = std::bind(&NonlocalFormFactor<PToV>::H_long_residue_psi2s, nonlocal_formfactor);
                residue_H_perp = std::bind(&NonlocalFormFactor<PToV>::H_perp_residue_psi2s, nonlocal_formfactor);
                residue_H_para = std::bind(&NonlocalFormFactor<PToV>::H_para_residue_psi2s, nonlocal_formfactor);
            }

            u.uses(*model);
            u.uses(*nonlocal_formfactor);
        }

        ~Implementation() = default;

        // The amplitudes in the conventions of [BCvDV2016], eq. (B14)
        struct AmplitudesBCvDV2016
        {
            complex<double> A_perp, A_para, A_long;
        };

        // The amplitudes in the conventions of [T2002], eq. (2.38)
        struct AmplitudesExperimental
        {
            complex<double> A_perp, A_para, A_long;
        };

        AmplitudesBCvDV2016 amplitudes_bcvdv2016() const
        {
            const complex<double> res_H_long = this->residue_H_long();
            const complex<double> res_H_perp = this->residue_H_perp();
            const complex<double> res_H_para = this->residue_H_para();

            const double m_B = this->m_B(), m_B2 = power_of<2>(m_B);
            const double m_psi = this->m_psi();

            complex<double> A_perp = m_B2 / (f_psi * m_psi) * res_H_perp;
            complex<double> A_para = m_B2 / (f_psi * m_psi) * res_H_para;
            complex<double> A_long = m_B2 / (f_psi * m_psi) * res_H_long;

            return { A_perp, A_para, A_long };
        }

        // Amplitudes are CP invariant according to [BRY:2006A].
        AmplitudesExperimental amplitudes_experimental() const
        {
            static complex<double> I(0.0, 1.0);

            const auto amps = this->amplitudes_bcvdv2016();

            return {
                    -I * amps.A_perp,
                    -I * amps.A_para,
                    +I * (m_B / m_psi) * amps.A_long
                };
        }

        double branching_ratio() const
        {
            const auto amps = amplitudes_bcvdv2016();
            const auto lambda = eos::lambda(power_of<2>(m_B), power_of<2>(m_Kstar), power_of<2>(m_psi));
            const auto prefactor = pow(g_fermi * abs(model->ckm_cb() * conj(model->ckm_cs())), 2)
                    * tau_B() / hbar() * sqrt(lambda) / (2.0 * M_PI * m_B);

            return prefactor * (norm(amps.A_perp) + norm(amps.A_para) + power_of<2>(m_B/m_psi) * norm(amps.A_long));
        }

    };

    BToKstarCharmonium::BToKstarCharmonium(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<BToKstarCharmonium>(new Implementation<BToKstarCharmonium>(p, o, *this))
    {
    }

    BToKstarCharmonium::~BToKstarCharmonium() = default;

    const std::vector<OptionSpecification>
    Implementation<BToKstarCharmonium>::options
    {
        {"q"_ok, { "d"s, "u"s }, "d"s},
        {"psi"_ok, { "J/psi"s, "psi(2S)"s }, "J/psi"s}
    };

    double
    BToKstarCharmonium::branching_ratio() const
    {
        return _imp->branching_ratio();
    }

    double
    BToKstarCharmonium::perp_polarization() const
    {
        const auto amps = _imp->amplitudes_experimental();

        return norm(amps.A_perp) / (norm(amps.A_perp) + norm(amps.A_para) + norm(amps.A_long));
    }

    double
    BToKstarCharmonium::para_polarization() const
    {
        const auto amps = _imp->amplitudes_experimental();

        return norm(amps.A_para) / (norm(amps.A_perp) + norm(amps.A_para) + norm(amps.A_long));
    }

    double
    BToKstarCharmonium::long_polarization() const
    {
        const auto amps = _imp->amplitudes_experimental();

        return norm(amps.A_long) / (norm(amps.A_perp) + norm(amps.A_para) + norm(amps.A_long));
    }

    double
    BToKstarCharmonium::long_phase() const
    {
        const auto amps = _imp->amplitudes_experimental();

        return arg(amps.A_long);
    }

    double
    BToKstarCharmonium::delta_perp_long() const
    {
        const auto amps = _imp->amplitudes_experimental();

        const auto result = arg(amps.A_perp / amps.A_long);

        // clamp this result between 0 and 2 pi
        if (result < 0)
        {
            return result + 2.0 * M_PI;
        }
        else
        {
            return result;
        }
    }

    double
    BToKstarCharmonium::delta_para_long() const
    {
        const auto amps = _imp->amplitudes_experimental();

        const auto result = arg(amps.A_para / amps.A_long);

        // clamp this result between -2 pi and 0
        if (result > 0)
        {
            return result - 2.0 * M_PI;
        }
        else
        {
            return result;
        }
    }

    double
    BToKstarCharmonium::S_1c_LHCb() const
    {
        return this->long_polarization();
    }

    double
    BToKstarCharmonium::S_1s_LHCb() const
    {
        return 3.0 / 4.0 * (this->perp_polarization() + this->para_polarization());
    }

    double
    BToKstarCharmonium::S_3_LHCb() const
    {
        return 1.0 / 2.0 * (this->perp_polarization() - this->para_polarization());
    }

    double
    BToKstarCharmonium::S_4_LHCb() const
    {
        return +sqrt(this->long_polarization() * this->para_polarization() / 2.0) * cos(this->delta_para_long());
    }

    double
    BToKstarCharmonium::S_8_LHCb() const
    {
        return +sqrt(this->long_polarization() * this->perp_polarization() / 2.0) * sin(-this->delta_perp_long());
    }

    double
    BToKstarCharmonium::S_9_LHCb() const
    {
        return +sqrt(this->para_polarization() * this->perp_polarization()) * sin(this->delta_perp_long() - this->delta_para_long());
    }

    const std::set<ReferenceName>
    BToKstarCharmonium::references
    {
        "KMPW:2010A"_rn,
        "GvDV:2020A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    BToKstarCharmonium::begin_options()
    {
        return Implementation<BToKstarCharmonium>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToKstarCharmonium::end_options()
    {
        return Implementation<BToKstarCharmonium>::options.cend();
    }
}
