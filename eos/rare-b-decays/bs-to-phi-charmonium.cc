/*
 * Copyright (c) 2021 Méril Reboud
 * Copyright (c) 2025 Danny van Dyk
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
#include <eos/rare-b-decays/bs-to-phi-charmonium.hh>
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

    /*!
     * Implementation for the decay @f$\bar{B_s} \to \phi \psi@f$.
     */
    template <>
    struct Implementation<BsToPhiCharmonium>
    {
        UsedParameter g_fermi;

        UsedParameter hbar;

        std::shared_ptr<Model> model;

        UsedParameter m_Bs;

        UsedParameter tau_Bs;

        UsedParameter m_phi;

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
            m_Bs(p["mass::B_s"], u),
            tau_Bs(p["life_time::B_s"], u),
            m_phi(p["mass::phi"], u),
            opt_nonlocal_formfactor(o, "nonlocal-formfactor"_ok, { "GvDV2020", "naive", "GRvDV2022order5" }, "GvDV2020"),
            nonlocal_formfactor(NonlocalFormFactor<PToV>::make("B_s->phi::" + opt_nonlocal_formfactor.value(), p, o)),
            opt_psi(o, "psi"_ok, { "J/psi", "psi(2S)" }, "J/psi"),
            m_psi(p["mass::" + opt_psi.value()], u),
            f_psi(p["decay-constant::" + opt_psi.value()], u)
        {
            Context ctx("When constructing Bs->Phipsi observables");

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

        // The amplitudes in the conventions of ??
        struct AmplitudesExperimental
        {
            complex<double> A_perp, A_para, A_long;
        };

        AmplitudesBCvDV2016 amplitudes_bcvdv2016() const
        {
            const complex<double> res_H_long = this->residue_H_long();
            const complex<double> res_H_perp = this->residue_H_perp();
            const complex<double> res_H_para = this->residue_H_para();

            const double m_Bs = this->m_Bs(), m_Bs2 = power_of<2>(m_Bs);
            const double m_psi = this->m_psi();

            complex<double> A_perp = m_Bs2 / (f_psi * m_psi) * res_H_perp;
            complex<double> A_para = m_Bs2 / (f_psi * m_psi) * res_H_para;
            complex<double> A_long = m_Bs2 / (f_psi * m_psi) * res_H_long;

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
                    +I * (m_Bs / m_psi) * amps.A_long
                };
        }

        double branching_ratio() const
        {
            const auto amps = amplitudes_bcvdv2016();
            const auto lambda = eos::lambda(power_of<2>(m_Bs), power_of<2>(m_phi), power_of<2>(m_psi));
            const auto prefactor = power_of<2>(g_fermi * abs(model->ckm_cb() * conj(model->ckm_cs())))
                    * tau_Bs() / hbar() * sqrt(lambda) / (2.0 * M_PI * m_Bs);

            return prefactor * (norm(amps.A_perp) + norm(amps.A_para) + pow(m_Bs/m_psi, 2) * norm(amps.A_long));
        }

    };

    BsToPhiCharmonium::BsToPhiCharmonium(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<BsToPhiCharmonium>(new Implementation<BsToPhiCharmonium>(p, o, *this))
    {
    }

    BsToPhiCharmonium::~BsToPhiCharmonium() = default;

    const std::vector<OptionSpecification>
    Implementation<BsToPhiCharmonium>::options
    {
        {"psi"_ok, { "J/psi", "psi(2S)" }, "J/psi"}
    };

    double
    BsToPhiCharmonium::branching_ratio() const
    {
        return _imp->branching_ratio();
    }

    double
    BsToPhiCharmonium::perp_polarization() const
    {
        const auto amps = _imp->amplitudes_experimental();

        return norm(amps.A_perp) / (norm(amps.A_perp) + norm(amps.A_para) + norm(amps.A_long));
    }

    double
    BsToPhiCharmonium::para_polarization() const
    {
        const auto amps = _imp->amplitudes_experimental();

        return norm(amps.A_para) / (norm(amps.A_perp) + norm(amps.A_para) + norm(amps.A_long));
    }

    double
    BsToPhiCharmonium::long_polarization() const
    {
        const auto amps = _imp->amplitudes_experimental();

        return norm(amps.A_long) / (norm(amps.A_perp) + norm(amps.A_para) + norm(amps.A_long));
    }

    double
    BsToPhiCharmonium::long_phase() const
    {
        const auto amps = _imp->amplitudes_experimental();

        return arg(amps.A_long);
    }

    double
    BsToPhiCharmonium::delta_perp_long() const
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
    BsToPhiCharmonium::delta_para_long() const
    {
        const auto amps = _imp->amplitudes_experimental();

        const auto result = arg(amps.A_para / amps.A_long);

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
    BsToPhiCharmonium::S_1c_LHCb() const
    {
        return this->long_polarization();
    }

    double
    BsToPhiCharmonium::S_1s_LHCb() const
    {
        return 3.0 / 4.0 * (this->perp_polarization() + this->para_polarization());
    }

    double
    BsToPhiCharmonium::S_3_LHCb() const
    {
        return 1.0 / 2.0 * (this->perp_polarization() - this->para_polarization());
    }

    double
    BsToPhiCharmonium::S_4_LHCb() const
    {
        return +sqrt(this->long_polarization() * this->para_polarization() / 2.0) * cos(this->delta_para_long());
    }

    double
    BsToPhiCharmonium::S_8_LHCb() const
    {
        return +sqrt(this->long_polarization() * this->perp_polarization() / 2.0) * sin(-this->delta_perp_long());
    }

    double
    BsToPhiCharmonium::S_9_LHCb() const
    {
        return +sqrt(this->para_polarization() * this->perp_polarization()) * sin(this->delta_perp_long() - this->delta_para_long());
    }

    const std::set<ReferenceName>
    BsToPhiCharmonium::references
    {
        "KMPW:2010A"_rn,
        "GvDV:2020A"_rn
    };

    std::vector<OptionSpecification>::const_iterator
    BsToPhiCharmonium::begin_options()
    {
        return Implementation<BsToPhiCharmonium>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BsToPhiCharmonium::end_options()
    {
        return Implementation<BsToPhiCharmonium>::options.cend();
    }
}
