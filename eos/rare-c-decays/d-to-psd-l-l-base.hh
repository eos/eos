/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026    Carolina Bolognani
 * Copyright (c) 2026    Dominik Suelmann
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

#ifndef MASTER_GUARD_EOS_RARE_C_DECAYS_D_TO_PSDEUDOSCALAR_LEPTON_LEPTON_BASE_HH
#define MASTER_GUARD_EOS_RARE_C_DECAYS_D_TO_PSDEUDOSCALAR_LEPTON_LEPTON_BASE_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/models/model.hh>
#include <eos/rare-c-decays/d-to-psd-l-l.hh>
#include <eos/utils/options-impl.hh>

namespace eos
{
    class DToPseudoscalarLeptonLepton::AmplitudeGenerator : public ParameterUser
    {
        public:
            std::shared_ptr<Model>             model;
            LeptonFlavorOption                 opt_l;
            BooleanOption                      opt_cp_conjugate;
            QuarkFlavorOption                  opt_q;
            SpecifiedOption                    opt_P;
            bool                               cp_conjugate;
            LeptonFlavor                       lepton_flavor;
            QuarkFlavor                        q;
            std::shared_ptr<FormFactors<PToP>> form_factors;

            UsedParameter mu;
            UsedParameter alpha_e;
            UsedParameter g_fermi;
            UsedParameter hbar;

            const double isospin_factor;

            UsedParameter m_D;
            UsedParameter m_P;
            UsedParameter m_l;

            static const std::vector<OptionSpecification> options;

            // { q, P } -> { process, D_name, P_name }
            // q: u, d, s: the spectator quark flavor
            // P: K, pi, eta: the type of daughter meson
            // process: string that can be used to obtain the form factor
            // D_name: name of the D meson
            // P_name: name of the daughter meson
            // c_I: isospin factor by which the amplitudes are multiplied
            static const std::map<std::tuple<QuarkFlavor, std::string>, std::tuple<std::string, std::string, std::string, std::string, double>> process_map;

            complex<double> ckm_factor_ds;

            inline std::string
            _process_q() const
            {
                const std::string P = opt_P.value();
                const auto        p = process_map.find(std::make_tuple(q, P));

                if (p == process_map.end())
                {
                    throw InternalError("Unsupported combination of q = " + stringify(q) + ", P = " + P);
                }

                return std::get<0>(p->second);
            }

            inline std::string
            _process() const
            {
                const std::string P = opt_P.value();
                const auto        p = process_map.find(std::make_tuple(q, P));

                if (p == process_map.end())
                {
                    throw InternalError("Unsupported combination of q = " + stringify(q) + ", P = " + P);
                }

                return std::get<1>(p->second);
            }

            inline std::string
            _D() const
            {
                const std::string P = opt_P.value();
                const auto        p = process_map.find(std::make_tuple(q, P));

                if (p == process_map.end())
                {
                    throw InternalError("Unsupported combination of q = " + stringify(q) + ", P = " + P);
                }

                return std::get<2>(p->second);
            }

            inline std::string
            _P() const
            {
                const std::string P = opt_P.value();
                const auto        p = process_map.find(std::make_tuple(q, P));

                if (p == process_map.end())
                {
                    throw InternalError("Unsupported combination of q = " + stringify(q) + ", P = " + P);
                }

                return std::get<3>(p->second);
            }

            inline double
            _isospin_factor() const
            {
                const std::string P = opt_P.value();
                const auto        p = process_map.find(std::make_tuple(q, P));

                if (p == process_map.end())
                {
                    throw InternalError("Unsupported combination of q = " + stringify(q) + ", P = " + P);
                }

                return std::get<4>(p->second);
            }

            AmplitudeGenerator(const Parameters &, const Options &);

            double mqatmu() const;
            double beta_l(const double & q2) const;
            double energy(const double & q2) const;
            double lambda(const double & q2) const;
            double xi_pseudo(const double & q2) const;
            double normalisation(const double & q2) const;

            virtual ~AmplitudeGenerator();
            virtual DToPseudoscalarLeptonLepton::Amplitudes amplitudes(const double & q2) const = 0;
    };

    template <typename Tag_> class DToPseudoscalarLeptonLeptonAmplitudes;

    namespace tag
    {
        /*
         * Amplitudes according to [BGHT:2019A].
         */
        struct BGHT2019;
    } // namespace tag
} // namespace eos

#endif
