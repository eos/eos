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

#ifndef EOS_GUARD_EOS_NONLEPTONIC_AMPLITUDES_QCDF_AMPLITUDES_HH
#define EOS_GUARD_EOS_NONLEPTONIC_AMPLITUDES_QCDF_AMPLITUDES_HH 1

#include <eos/maths/complex.hh>
#include <eos/models/model.hh>
#include <eos/nonleptonic-amplitudes/nonleptonic-amplitudes-fwd.hh>
#include <eos/nonleptonic-amplitudes/nonleptonic-amplitudes.hh>
#include <eos/utils/parameters.hh>

#include <array>
#include <map>

namespace eos
{
    template <typename Transition_> class QCDFRepresentation;

    template <> class QCDFRepresentation<PToPP> : public NonleptonicAmplitudes<PToPP>
    {
        private:
            std::shared_ptr<Model> model;
            QuarkFlavorOption      opt_q;
            LightMesonOption       opt_p1;
            LightMesonOption       opt_p2;
            BooleanOption          opt_cp_conjugate;
            BooleanOption          opt_B_bar;

            UsedParameter theta_18;

            su3f::rank1                  B;
            std::function<su3f::rank1()> Lambda_u, Lambda_c;
            mutable su3f::rank2          P1, P2, U, I;

            UsedParameter Gfermi;
            UsedParameter mB;
            UsedParameter mB_q_0;
            UsedParameter mP1;
            UsedParameter mP2;
            UsedParameter FP1;
            UsedParameter FP2;
            UsedParameter fB;
            UsedParameter fP1;
            UsedParameter fP2;

            UsedParameter re_alpha1, im_alpha1;
            UsedParameter re_alpha2, im_alpha2;
            UsedParameter re_b1, im_b1;
            UsedParameter re_b2, im_b2;
            UsedParameter re_bS1, im_bS1;
            UsedParameter re_bS2, im_bS2;

            UsedParameter re_alpha3_u, im_alpha3_u;
            UsedParameter re_alpha3_c, im_alpha3_c;
            UsedParameter re_alpha4_u, im_alpha4_u;
            UsedParameter re_alpha4_c, im_alpha4_c;
            UsedParameter re_b4_u, im_b4_u;
            UsedParameter re_b4_c, im_b4_c;
            UsedParameter re_bS4_u, im_bS4_u;
            UsedParameter re_bS4_c, im_bS4_c;

            UsedParameter re_alpha3EW_c, im_alpha3EW_c;
            UsedParameter re_alpha4EW_c, im_alpha4EW_c;
            UsedParameter re_b3EW_c, im_b3EW_c;
            UsedParameter re_bS3EW_c, im_bS3EW_c;
            UsedParameter re_b4EW_c, im_b4EW_c;
            UsedParameter re_bS4EW_c, im_bS4EW_c;

            static const std::vector<OptionSpecification> options;

            std::function<complex<double>()> lamdu;
            std::function<complex<double>()> lamsu;
            std::function<complex<double>()> lamdc;
            std::function<complex<double>()> lamsc;

        public:
            QCDFRepresentation(const Parameters & p, const Options & o);

            ~QCDFRepresentation() {}

            inline void
            update() const
            {
                const double theta_18 = this->theta_18.evaluate();
                su3f::psd_octet.find(opt_p1.value())->second(theta_18, P1);
                su3f::psd_octet.find(opt_p2.value())->second(theta_18, P2);
            }

            static NonleptonicAmplitudes<PToPP> * make(const Parameters &, const Options &);

            // Helper functions
            complex<double> alpha_amplitude(su3f::rank2 & p1, su3f::rank2 & p2) const;
            complex<double> b_amplitude(su3f::rank2 & p1, su3f::rank2 & p2) const;

            // Diagnostic functions
            complex<double>
            alpha_amplitude() const
            {
                update();
                return alpha_amplitude(P1, P2);
            }

            complex<double>
            b_amplitude() const
            {
                update();
                return b_amplitude(P1, P2);
            }

            // Amplitude for B -> P1 P2
            complex<double> ordered_amplitude() const;
            // Amplitude for B -> P2 P1
            complex<double> inverse_amplitude() const;
    };
} // namespace eos
#endif
