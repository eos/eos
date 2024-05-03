/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Méril Reboud
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

#include <eos/b-decays/b-to-psd-psd.hh>
#include <eos/models/model.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/maths/power-of.hh>

using namespace std::literals::complex_literals;

namespace eos
{

    template <>
    struct Implementation<BToPseudoscalarPseudoscalar>
    {
        QuarkFlavorOption opt_q;
        LightMesonOption opt_p1;
        LightMesonOption opt_p2;
        UsedParameter hbar;
        UsedParameter tau;
        UsedParameter mB;
        UsedParameter mP1;
        UsedParameter mP2;
        SpecifiedOption opt_rep;
        std::shared_ptr<NonleptonicAmplitudes<PToPP>> nl_amplitudes;
        std::shared_ptr<NonleptonicAmplitudes<PToPP>> cp_nl_amplitudes;
        std::shared_ptr<NonleptonicAmplitudes<PToPP>> Bbar_nl_amplitudes;
        std::shared_ptr<Model> model;

        // B mixing
        std::shared_ptr<UsedParameter> life_time_difference;
        std::function<double ()> yq;
        std::function<double ()> phiB;
        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            opt_q(o, options, "q"),
            opt_p1(o, options, "P1"),
            opt_p2(o, options, "P2"),
            hbar(p["QM::hbar"], u),
            tau(p["life_time::B_" + opt_q.str()], u),
            mB(p["mass::B_" + opt_q.str()], u),
            mP1(p["mass::" + opt_p1.str()], u),
            mP2(p["mass::" + opt_p2.str()], u),
            opt_rep(o, options, "representation"),
            nl_amplitudes(NonleptonicAmplitudeFactory<PToPP>::create("B->PP::" + opt_rep.value(), p, o + Options{{"cp-conjugate", "false"}})),
            cp_nl_amplitudes(NonleptonicAmplitudeFactory<PToPP>::create("B->PP::" + opt_rep.value(), p, o + Options{{"cp-conjugate", "true"}})),
            Bbar_nl_amplitudes(NonleptonicAmplitudeFactory<PToPP>::create("B->PP::" + opt_rep.value(), p, o + Options{{"cp-conjugate", "false"}} + Options{{"B_bar", "true"}}))
        {
            Context ctx("When constructing B->PP observable");

            switch (opt_q.value())
            {
                case QuarkFlavor::up:
                    life_time_difference = nullptr;
                    yq = []() { return 0.0; };
                    phiB = []() { return 0.0; };
                    break;

                case QuarkFlavor::down:
                    life_time_difference = std::make_shared<UsedParameter>(p["life_time::Delta_B_d"], u);
                    yq = [this]() { return this->life_time_difference->evaluate() / (2.0 / tau); };
                    phiB = [this]() { return 2.0 * arg(model->ckm_tb() * conj(model->ckm_td())); };
                    break;

                case QuarkFlavor::strange:
                    life_time_difference = std::make_shared<UsedParameter>(p["life_time::Delta_B_s"], u);
                    yq = [this]() { return this->life_time_difference->evaluate() / (2.0 / tau); };
                    phiB = [this]() { return 2.0 * arg(model->ckm_tb() * conj(model->ckm_ts())); };
                    break;

                default:
                    throw InternalError("Invalid quark flavor: " + stringify(opt_q.value()));
            }

            u.uses(*model);
            u.uses(*nl_amplitudes);
            u.uses(*cp_nl_amplitudes);
            u.uses(*Bbar_nl_amplitudes);
        }

        double decay_width() const
        {
            const double mB = this->mB();
            const double mP1 = this->mP1();
            const double mP2 = this->mP2();
            double prefactor = 1.0 / 16.0 / M_PI / power_of<3>(mB) * std::sqrt(lambda(mB * mB, mP1 * mP1, mP2 * mP2));

            // The decay width is a physical observable, the inputs are symmetrized
            // BR(B -> P1 P2) = prefactor * S * |A(B -> P1 P2) + A(B -> P2 P1)|^2
            if (opt_p1.value() == opt_p2.value())
            {
                prefactor *= 0.5;
            }

            return  prefactor * std::norm(nl_amplitudes->amplitude());
        }

        double cp_decay_width() const
        {
            const double mB = this->mB();
            const double mP1 = this->mP1();
            const double mP2 = this->mP2();
            double prefactor = 1.0 / 16.0 / M_PI / power_of<3>(mB) * std::sqrt(lambda(mB * mB, mP1 * mP1, mP2 * mP2));

            // The decay width is a physical observable, the inputs are symmetrized
            // BR(B -> P1 P2) = prefactor * S * |A(B -> P1 P2) + A(B -> P2 P1)|^2
            if (opt_p1.value() == opt_p2.value())
            {
                prefactor *= 0.5;
            }

            return  prefactor * std::norm(cp_nl_amplitudes->amplitude());
        }

        double mixing_induced_cp_asymmetry() const
        {
            const complex<double> amp = nl_amplitudes->amplitude();
            const complex<double> Bbar_amp = Bbar_nl_amplitudes->amplitude();

            // Assumes the mixing parameter ratio q / p to be a pure phase
            const complex<double> xif = - std::exp(-1.0i * phiB()) * Bbar_amp / amp;


            return 2 * std::imag(xif) / (1 + std::norm(xif)) ;
        }

        double a_Delta_Gamma() const
        {
            const complex<double> amp = nl_amplitudes->amplitude();
            const complex<double> Bbar_amp = Bbar_nl_amplitudes->amplitude();

            // Assumes the mixing parameter ratio q / p to be a pure phase
            const complex<double> xif = - std::exp(-1.0i * phiB()) * Bbar_amp / amp;


            return 2 * std::real(xif) / (1 + std::norm(xif)) ;
        }
    };

    const std::vector<OptionSpecification>
    Implementation<BToPseudoscalarPseudoscalar>::options
    {
        Model::option_specification(),
        NonleptonicAmplitudeFactory<PToPP>::option_specification(),
        { "q",              { "u", "d", "s" },                                                                       "" },
        { "P1",             { "pi^0", "pi^+", "pi^-", "K_d", "Kbar_d", "K_S", "K_u", "Kbar_u", "eta", "eta_prime" }, "" },
        { "P2",             { "pi^0", "pi^+", "pi^-", "K_d", "Kbar_d", "K_S", "K_u", "Kbar_u", "eta", "eta_prime" }, "" },
        { "representation", { "topological", "SU3F", "QCDF" },                                                       "" }
    };

    BToPseudoscalarPseudoscalar::BToPseudoscalarPseudoscalar(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToPseudoscalarPseudoscalar>(new Implementation<BToPseudoscalarPseudoscalar>(parameters, options, *this))
    {

    }

    BToPseudoscalarPseudoscalar::~BToPseudoscalarPseudoscalar()
    {
    }

    double
    BToPseudoscalarPseudoscalar::decay_width() const
    {
        return _imp->decay_width();
    }

    double
    BToPseudoscalarPseudoscalar::branching_ratio() const
    {
        return _imp->decay_width() * _imp->tau() / _imp->hbar();
    }

    double
    BToPseudoscalarPseudoscalar::cp_branching_ratio() const
    {
        return _imp->cp_decay_width() * _imp->tau() / _imp->hbar();
    }

    double
    BToPseudoscalarPseudoscalar::avg_branching_ratio() const
    {
        return 0.5 * (branching_ratio() + cp_branching_ratio());
    }

    double
    BToPseudoscalarPseudoscalar::exp_branching_ratio() const
    {
        const double yq = _imp->yq();
        return avg_branching_ratio() * (1 + _imp->a_Delta_Gamma() * yq) / (1 - power_of<2>(yq));
    }

    double
    BToPseudoscalarPseudoscalar::cp_asymmetry() const
    {
        return (branching_ratio() - cp_branching_ratio()) / (branching_ratio() + cp_branching_ratio());
    }

    double
    BToPseudoscalarPseudoscalar::mixing_induced_cp_asymmetry() const
    {
        return _imp->mixing_induced_cp_asymmetry();
    }

    double
    BToPseudoscalarPseudoscalar::a_Delta_Gamma() const
    {
        return _imp->a_Delta_Gamma();
    }


    const std::string
    BToPseudoscalarPseudoscalar::description = "\
    The decay B->PP, where all states are pseudoscalars.";

    const std::set<ReferenceName>
    BToPseudoscalarPseudoscalar::references
    {
        "HTX:2021A"_rn,
    };

    std::vector<OptionSpecification>::const_iterator
    BToPseudoscalarPseudoscalar::begin_options()
    {
        return Implementation<BToPseudoscalarPseudoscalar>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BToPseudoscalarPseudoscalar::end_options()
    {
        return Implementation<BToPseudoscalarPseudoscalar>::options.cend();
    }
}
