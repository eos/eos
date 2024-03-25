/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Marta Burgos
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

#ifndef EOS_GUARD_EOS_NONLEPTONIC_AMPLITUDES_SU3_AMPLITUDES_HH
#define EOS_GUARD_EOS_NONLEPTONIC_AMPLITUDES_SU3_AMPLITUDES_HH 1

#include <eos/nonleptonic-amplitudes/nonleptonic-amplitudes.hh>
#include <eos/nonleptonic-amplitudes/nonleptonic-amplitudes-fwd.hh>
#include <eos/maths/complex.hh>
#include <eos/models/model.hh>
#include <eos/utils/parameters.hh>

#include <array>
#include <map>

namespace eos
{
    template <typename Transition_> class SU3FRepresentation;

    template <>
    class SU3FRepresentation<PToPP> :
        public NonleptonicAmplitudes<PToPP>
    {
        private:
            std::shared_ptr<Model> model;
            QuarkFlavorOption opt_q;
            LightMesonOption opt_p1;
            LightMesonOption opt_p2;
            BooleanOption opt_cp_conjugate;

            UsedParameter theta_18;

            su3f::rank1 B, H3bar, H3tilde;
            su3f::rank2 P1, P2;
            su3f::rank3 H6bar, H6tilde, H15bar, H15tilde;

            UsedParameter Gfermi;

            UsedParameter re_AT3, im_AT3;
            UsedParameter re_CT3, im_CT3;
            UsedParameter re_AT6, im_AT6;
            UsedParameter re_CT6, im_CT6;
            UsedParameter re_AT15, im_AT15;
            UsedParameter re_CT15, im_CT15;
            UsedParameter re_BT3, im_BT3;
            UsedParameter re_BT6, im_BT6;
            UsedParameter re_BT15, im_BT15;
            UsedParameter re_DT3, im_DT3;

            UsedParameter re_AP3, im_AP3;
            UsedParameter re_CP3, im_CP3;
            UsedParameter re_AP6, im_AP6;
            UsedParameter re_CP6, im_CP6;
            UsedParameter re_AP15, im_AP15;
            UsedParameter re_CP15, im_CP15;
            UsedParameter re_BP3, im_BP3;
            UsedParameter re_BP6, im_BP6;
            UsedParameter re_BP15, im_BP15;
            UsedParameter re_DP3, im_DP3;

            static const std::vector<OptionSpecification> options;

            complex<double> lamdu;
            complex<double> lamsu;
            complex<double> lamdt;
            complex<double> lamst;

        public:
            SU3FRepresentation(const Parameters & p, const Options & o) :
                model(Model::make(o.get("model", "SM"), p, o)),
                opt_q(o, options, "q"),
                opt_p1(o, options, "P1"),
                opt_p2(o, options, "P2"),
                opt_cp_conjugate(o, options, "cp-conjugate"),
                theta_18(p["eta::theta_18"], *this),
                B(su3f::psd_b_triplet.find(opt_q.value())->second),
                H3bar({}),
                H3tilde({}),
                P1(su3f::psd_octet(theta_18()).find(opt_p1.value())->second),
                P2(su3f::psd_octet(theta_18()).find(opt_p2.value())->second),
                H6bar({}),
                H6tilde({}),
                H15bar({}),
                H15tilde({}),
                Gfermi(p["WET::G_Fermi"], *this),
                re_AT3( p["nonleptonic::Re{AT3}@SU3F"] , *this),
                im_AT3( p["nonleptonic::Im{AT3}@SU3F"] , *this),
                re_CT3( p["nonleptonic::Re{CT3}@SU3F"] , *this),
                im_CT3( p["nonleptonic::Im{CT3}@SU3F"] , *this),
                re_AT6( p["nonleptonic::Re{AT6}@SU3F"] , *this),
                im_AT6( p["nonleptonic::Im{AT6}@SU3F"] , *this),
                re_CT6( p["nonleptonic::Re{CT6}@SU3F"] , *this),
                im_CT6( p["nonleptonic::Im{CT6}@SU3F"] , *this),
                re_AT15(p["nonleptonic::Re{AT15}@SU3F"], *this),
                im_AT15(p["nonleptonic::Im{AT15}@SU3F"], *this),
                re_CT15(p["nonleptonic::Re{CT15}@SU3F"], *this),
                im_CT15(p["nonleptonic::Im{CT15}@SU3F"], *this),
                re_BT3( p["nonleptonic::Re{BT3}@SU3F"] , *this),
                im_BT3( p["nonleptonic::Im{BT3}@SU3F"] , *this),
                re_BT6( p["nonleptonic::Re{BT6}@SU3F"] , *this),
                im_BT6( p["nonleptonic::Im{BT6}@SU3F"] , *this),
                re_BT15(p["nonleptonic::Re{BT15}@SU3F"], *this),
                im_BT15(p["nonleptonic::Im{BT15}@SU3F"], *this),
                re_DT3( p["nonleptonic::Re{DT3}@SU3F"] , *this),
                im_DT3( p["nonleptonic::Im{DT3}@SU3F"] , *this),
                re_AP3( p["nonleptonic::Re{AP3}@SU3F"] , *this),
                im_AP3( p["nonleptonic::Im{AP3}@SU3F"] , *this),
                re_CP3( p["nonleptonic::Re{CP3}@SU3F"] , *this),
                im_CP3( p["nonleptonic::Im{CP3}@SU3F"] , *this),
                re_AP6( p["nonleptonic::Re{AP6}@SU3F"] , *this),
                im_AP6( p["nonleptonic::Im{AP6}@SU3F"] , *this),
                re_CP6( p["nonleptonic::Re{CP6}@SU3F"] , *this),
                im_CP6( p["nonleptonic::Im{CP6}@SU3F"] , *this),
                re_AP15(p["nonleptonic::Re{AP15}@SU3F"], *this),
                im_AP15(p["nonleptonic::Im{AP15}@SU3F"], *this),
                re_CP15(p["nonleptonic::Re{CP15}@SU3F"], *this),
                im_CP15(p["nonleptonic::Im{CP15}@SU3F"], *this),
                re_BP3( p["nonleptonic::Re{BP3}@SU3F"] , *this),
                im_BP3( p["nonleptonic::Im{BP3}@SU3F"] , *this),
                re_BP6( p["nonleptonic::Re{BP6}@SU3F"] , *this),
                im_BP6( p["nonleptonic::Im{BP6}@SU3F"] , *this),
                re_BP15(p["nonleptonic::Re{BP15}@SU3F"], *this),
                im_BP15(p["nonleptonic::Im{BP15}@SU3F"], *this),
                re_DP3( p["nonleptonic::Re{DP3}@SU3F"] , *this),
                im_DP3( p["nonleptonic::Im{DP3}@SU3F"] , *this)
            {
                Context ctx("When constructing B->PP SU3 amplitudes");

                if (opt_cp_conjugate.value())
                {
                    lamdu = conj(model->ckm_ub()) * model->ckm_ud();
                    lamsu = conj(model->ckm_ub()) * model->ckm_us();
                    lamdt = conj(model->ckm_tb()) * model->ckm_td();
                    lamst = conj(model->ckm_tb()) * model->ckm_ts();
                }
                else
                {
                    lamdu = model->ckm_ub() * conj(model->ckm_ud());
                    lamsu = model->ckm_ub() * conj(model->ckm_us());
                    lamdt = model->ckm_tb() * conj(model->ckm_td());
                    lamst = model->ckm_tb() * conj(model->ckm_ts());
                }

                H3bar[1] = +lamdu;
                H3bar[2] = +lamsu;

                H6bar[0][1][0] = +lamdu;
                H6bar[1][0][0] = -lamdu;
                H6bar[1][2][2] = +lamdu;
                H6bar[2][1][2] = -lamdu;

                H6bar[0][2][0] = +lamsu;
                H6bar[2][0][0] = -lamsu;
                H6bar[2][1][1] = +lamsu; // Corrected with respect to typo in [HTX:2021A]
                H6bar[1][2][1] = -lamsu; // Corrected with respect to typo in [HTX:2021A]

                H15bar[0][1][0] = +3.0 * lamdu;
                H15bar[1][0][0] = +3.0 * lamdu;
                H15bar[1][1][1] = -2.0 * lamdu;
                H15bar[1][2][2] = -lamdu;
                H15bar[2][1][2] = -lamdu;

                H15bar[0][2][0] = +3.0 * lamsu;
                H15bar[2][0][0] = +3.0 * lamsu;
                H15bar[2][2][2] = -2.0 * lamsu;
                H15bar[2][1][1] = -lamsu; // Corrected with respect to typo in [HTX:2021A]
                H15bar[1][2][1] = -lamsu; // Corrected with respect to typo in [HTX:2021A]


                H3tilde[1] = +lamdt;
                H3tilde[2] = +lamst;

                H6tilde[0][1][0] = +lamdt;
                H6tilde[1][0][0] = -lamdt;
                H6tilde[1][2][2] = +lamdt;
                H6tilde[2][1][2] = -lamdt;

                H6tilde[0][2][0] = +lamst;
                H6tilde[2][0][0] = -lamst;
                H6tilde[2][1][1] = +lamst; // Corrected with respect to typo in [HTX:2021A]
                H6tilde[1][2][1] = -lamst; // Corrected with respect to typo in [HTX:2021A]

                H15tilde[0][1][0] = +3.0 * lamdt;
                H15tilde[1][0][0] = +3.0 * lamdt;
                H15tilde[1][1][1] = -2.0 * lamdt;
                H15tilde[1][2][2] = -lamdt;
                H15tilde[2][1][2] = -lamdt;

                H15tilde[0][2][0] = +3.0 * lamst;
                H15tilde[2][0][0] = +3.0 * lamst;
                H15tilde[2][2][2] = -2.0 * lamst;
                H15tilde[2][1][1] = -lamst; // Corrected with respect to typo in [HTX:2021A]
                H15tilde[1][2][1] = -lamst; // Corrected with respect to typo in [HTX:2021A]
            };

            ~SU3FRepresentation() {};

            static NonleptonicAmplitudes<PToPP> * make(const Parameters &, const Options &);

            // Helper functions
            complex<double> tree_amplitude(su3f::rank2 p1, su3f::rank2 p2) const;
            complex<double> penguin_amplitude(su3f::rank2 p1, su3f::rank2 p2) const;

            // Diagnostic functions
            complex<double> tree_amplitude() const { return tree_amplitude(P1, P2); }
            complex<double> penguin_amplitude() const { return penguin_amplitude(P1, P2); };

            // Amplitude for B -> P1 P2
            complex<double> ordered_amplitude() const;
            // Amplitude for B -> P2 P1
            complex<double> inverse_amplitude() const;
    };
}
#endif
