/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021-2025 Danny van Dyk
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

#include <eos/maths/power-of.hh>
#include <eos/meson-mixing/bq-mixing.hh>
#include <eos/models/model.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

namespace eos
{
    using std::norm;

    template <>
    struct Implementation<BMixing>
    {
        std::shared_ptr<Model> model;

        UsedParameter mu;

        UsedParameter hbar;

        UsedParameter g_fermi;

        SwitchOption opt_q;

        UsedParameter m_B;

        UsedParameter f_B;

        UsedParameter tau_B;

        UsedParameter R_1;
        UsedParameter R_2;
        UsedParameter R_3;
        UsedParameter R_4;
        UsedParameter R_5;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "SM"), p, o)),
            mu(p["sbsb::mu"], u),
            hbar(p["QM::hbar"], u),
            g_fermi(p["WET::G_Fermi"], u),
            opt_q(o, "q"_ok, {"s"}),
            m_B(p["mass::B_" + opt_q.value()], u),
            f_B(p["decay-constant::B_" + opt_q.value()], u),
            tau_B(p["life_time::B_" + opt_q.value()], u),
            R_1(p["B_s<->Bbar_s::R^1"], u),
            R_2(p["B_s<->Bbar_s::R^2"], u),
            R_3(p["B_s<->Bbar_s::R^3"], u),
            R_4(p["B_s<->Bbar_s::R^4"], u),
            R_5(p["B_s<->Bbar_s::R^5"], u)
        {
            u.uses(*model);
        }

        complex<double> M_12() const
        {
            const auto wc = model->wet_sbsb();

            // cf. [DDHLMSW:2019A]
            // TODO: still needs to be evolved to scale mu from reference scale 4.2 GeV.
            const std::array<complex<double>, 8> contributions{{
                wc.c1()  * R_1(),
                wc.c2()  * R_2(),
                wc.c3()  * R_3(),
                wc.c4()  * R_4(),
                wc.c5()  * R_5(),
                wc.c1p() * R_1(), // primed operators share the hadronic matrix elements of their unprimed partners
                wc.c2p() * R_2(),
                wc.c3p() * R_3(),
            }};

            complex<double> result = 0.0;
            for (auto & c : contributions)
            {
                result += c;
            }

            // cf. [BBL:1995A], eq. (XVIII.17), p. 153
            return 4.0 * g_fermi / std::sqrt(2.0) * power_of<2>(model->ckm_tb() * std::conj(model->ckm_ts()))
                * f_B() * f_B() * m_B() / 2.0 * result;

        }
    };

    const std::vector<OptionSpecification>
    Implementation<BMixing>::options
    {
    };

    BMixing::BMixing(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BMixing>(new Implementation<BMixing>(parameters, options, *this))
    {
    }

    BMixing::~BMixing()
    {
    }

    double
    BMixing::delta_m() const
    {
        // cf. [BBL:1995A], eq. (XVIII.16), p. 153
        return 2.0 * std::abs(_imp->M_12()) / _imp->hbar() * 1.0e-12; // return value in ps^-1
    }

    const std::set<ReferenceName>
    BMixing::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    BMixing::begin_options()
    {
        return Implementation<BMixing>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BMixing::end_options()
    {
        return Implementation<BMixing>::options.cend();
    }
}
