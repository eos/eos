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

#include <eos/nonleptonic-amplitudes/nonleptonic-amplitudes.hh>
#include <eos/nonleptonic-amplitudes/nonleptonic-amplitudes-fwd.hh>
#include <eos/maths/complex.hh>
#include <eos/models/model.hh>
#include <eos/utils/parameters.hh>

#include <array>
#include <map>

namespace eos
{
    template <typename Transition_> class QCDFRepresentation;

    template <>
    class QCDFRepresentation<PToPP> :
        public NonleptonicAmplitudes<PToPP>
    {
        private:
            std::shared_ptr<Model> model;
            QuarkFlavorOption opt_q;
            LightMesonOption opt_p1;
            LightMesonOption opt_p2;
            BooleanOption opt_cp_conjugate;

            UsedParameter theta_18;

            su3f::rank1 B; 
            std::array<complex<double>, 6> T, P1_c, P1_u, P2_c, Lambda_u, Lambda_c;
            mutable su3f::rank2 P1, P2, U, I;
            std::array<std::array<std::array<complex<double>, 6>, 6>, 6> C_u, C_c;


            UsedParameter Gfermi;
            UsedParameter mB;
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

            complex<double> lamdu;
            complex<double> lamsu;
            complex<double> lamdc;
            complex<double> lamsc;
            

        public:
            QCDFRepresentation(const Parameters & p, const Options & o) :
                model(Model::make(o.get("model", "SM"), p, o)),
                opt_q(o, options, "q"),
                opt_p1(o, options, "P1"),
                opt_p2(o, options, "P2"),
                opt_cp_conjugate(o, options, "cp-conjugate"),
                theta_18(p["eta::theta_18"], *this),
                B(su3f::psd_b_triplet.find(opt_q.value())->second),
                T({}),
                P1_c({}),
                P1_u({}),
                P2_c({}),
                P1{{}},
                P2{{}},
                C_u({}),
                C_c({}),
                Gfermi(p["WET::G_Fermi"], *this),
                mB(p["mass::B_" + opt_q.str()], *this),
                FP1(p["B_" + opt_q.str() + "->" + opt_p1.str() + "::f_+(0)"], *this),
                FP2(p["B_" + opt_q.str() + "->" + opt_p2.str() + "::f_+(0)"], *this),
                fB(p["decay-constant::B_" + opt_q.str()], *this),
                fP1(p["decay-constant::" + opt_p1.str()], *this),
                fP2(p["decay-constant::" + opt_p2.str()], *this),
                
                re_alpha1(    p["nonleptonic::Re{alpha1}@QCDF"] , *this),
                im_alpha1(    p["nonleptonic::Im{alpha1}@QCDF"] , *this),
                re_alpha2(    p["nonleptonic::Re{alpha2}@QCDF"] , *this),
                im_alpha2(    p["nonleptonic::Im{alpha2}@QCDF"] , *this),
                re_b1(     p["nonleptonic::Re{b1}@QCDF"] , *this),
                im_b1(     p["nonleptonic::Im{b1}@QCDF"] , *this),
                re_b2(     p["nonleptonic::Re{b2}@QCDF"] , *this),
                im_b2(     p["nonleptonic::Im{b2}@QCDF"] , *this),
                re_bS1(    p["nonleptonic::Re{bS1}@QCDF"], *this),
                im_bS1(    p["nonleptonic::Im{bS1}@QCDF"], *this),
                re_bS2(    p["nonleptonic::Re{bS2}@QCDF"], *this),
                im_bS2(    p["nonleptonic::Im{bS2}@QCDF"], *this),

                re_alpha3_u(  p["nonleptonic::Re{alpha3_u}@QCDF"] , *this),
                im_alpha3_u(  p["nonleptonic::Im{alpha3_u}@QCDF"] , *this),
                re_alpha3_c(  p["nonleptonic::Re{alpha3_c}@QCDF"] , *this),
                im_alpha3_c(  p["nonleptonic::Im{alpha3_c}@QCDF"] , *this),
                re_alpha4_u(  p["nonleptonic::Re{alpha4_u}@QCDF"] , *this),
                im_alpha4_u(  p["nonleptonic::Im{alpha4_u}@QCDF"] , *this),
                re_alpha4_c(  p["nonleptonic::Re{alpha4_c}@QCDF"] , *this),
                im_alpha4_c(  p["nonleptonic::Im{alpha4_c}@QCDF"] , *this),
                re_b4_u(   p["nonleptonic::Re{b4_u}@QCDF"] , *this),
                im_b4_u(   p["nonleptonic::Im{b4_u}@QCDF"] , *this),
                re_b4_c(   p["nonleptonic::Re{b4_c}@QCDF"] , *this),
                im_b4_c(   p["nonleptonic::Im{b4_c}@QCDF"] , *this),
                re_bS4_u(  p["nonleptonic::Re{bS4_u}@QCDF"] , *this),
                im_bS4_u(  p["nonleptonic::Im{bS4_u}@QCDF"] , *this),
                re_bS4_c(  p["nonleptonic::Re{bS4_c}@QCDF"] , *this),
                im_bS4_c(  p["nonleptonic::Im{bS4_c}@QCDF"] , *this),

                re_alpha3EW_c(p["nonleptonic::Re{alpha3EW_c}@QCDF"] , *this),
                im_alpha3EW_c(p["nonleptonic::Im{alpha3EW_c}@QCDF"] , *this),
                re_alpha4EW_c(p["nonleptonic::Re{alpha4EW_c}@QCDF"] , *this),
                im_alpha4EW_c(p["nonleptonic::Im{alpha4EW_c}@QCDF"] , *this),
                re_b3EW_c( p["nonleptonic::Re{b3EW_c}@QCDF"], *this),
                im_b3EW_c( p["nonleptonic::Im{b3EW_c}@QCDF"], *this),
                re_bS3EW_c(p["nonleptonic::Re{bS3EW_c}@QCDF"], *this),
                im_bS3EW_c(p["nonleptonic::Im{bS3EW_c}@QCDF"], *this),
                re_b4EW_c(    p["nonleptonic::Re{b4EW_c}@QCDF"] , *this),
                im_b4EW_c(    p["nonleptonic::Im{b4EW_c}@QCDF"] , *this),
                re_bS4EW_c(   p["nonleptonic::Re{bS4EW_c}@QCDF"] , *this),
                im_bS4EW_c(   p["nonleptonic::Im{bS4EW_c}@QCDF"] , *this)
            {
                Context ctx("When constructing B->PP QCD amplitudes");

                if (opt_cp_conjugate.value())
                {
                    lamdu = model->ckm_ub() * conj(model->ckm_ud());
                    lamsu = model->ckm_ub() * conj(model->ckm_us());
                    lamdc = model->ckm_cb() * conj(model->ckm_cd());
                    lamsc = model->ckm_cb() * conj(model->ckm_cs());
                }
                else
                {
                    lamdu = conj(model->ckm_ub()) * model->ckm_ud();
                    lamsu = conj(model->ckm_ub()) * model->ckm_us();
                    lamdc = conj(model->ckm_cb()) * model->ckm_cd();
                    lamsc = conj(model->ckm_cb()) * model->ckm_cs();
                }



                U[0][0] = 1;
                
                I[0][0] = 1;
                I[1][1] = 1;
                I[2][2] = 1;

                Lambda_u[1] = lamdu;
                Lambda_u[2] = lamsu;

                Lambda_c[1] = lamdc;
                Lambda_c[2] = lamsc;

    
            };

            void update() const
            {
                const double theta_18 = this->theta_18.evaluate();
                su3f::psd_octet.find(opt_p1.value())->second(theta_18, P1);
                su3f::psd_octet.find(opt_p2.value())->second(theta_18, P2);
            }

            ~QCDFRepresentation() {};

            static NonleptonicAmplitudes<PToPP> * make(const Parameters &, const Options &);

            // Helper functions
            complex<double> alpha_amplitude(su3f::rank2 p1, su3f::rank2 p2) const;
            complex<double> b_amplitude(su3f::rank2 p1, su3f::rank2 p2) const;
            
            // Diagnostic functions
            complex<double> alpha_amplitude() const { update(); return alpha_amplitude(P1, P2); }
            complex<double> b_amplitude() const { update(); return b_amplitude(P1, P2); };

            // Amplitude for B -> P1 P2
            complex<double> ordered_amplitude() const;
            // Amplitude for B -> P2 P1
            complex<double> inverse_amplitude() const;
    };
}
#endif
