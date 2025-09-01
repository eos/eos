/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018, 2019 Ahmet Kokulu
 * Copyright (c) 2019, 2021 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_B_DECAYS_B_TO_VEC_L_NU_IMPL_HH
#define EOS_GUARD_EOS_B_DECAYS_B_TO_VEC_L_NU_IMPL_HH 1

#include <eos/observable.hh>
#include <eos/b-decays/b-to-vec-l-nu.hh>
#include <eos/maths/complex.hh>

#include <array>

namespace eos
{
    namespace b_to_vec_l_nu
    {
        struct Amplitudes
        {
            complex<double> a_0;
            complex<double> a_0_T;
            complex<double> a_plus;
            complex<double> a_plus_T;
            complex<double> a_minus;
            complex<double> a_minus_T;
            complex<double> a_P;
            complex<double> a_t;
            complex<double> a_para;
            complex<double> a_para_T;
            complex<double> a_perp;
            complex<double> a_perp_T;
            double mlH;
            double NF;
        };

        class AngularObservables
        {
            private:
                std::array<double, 12> _vv;

            public:
                friend class BToVectorLeptonNeutrino;
                friend class Implementation<BToVectorLeptonNeutrino>;

                // angular observables V's. cf. from [DDS:2014A], p. 16, redifined V's in order to include NF
                AngularObservables(const Amplitudes & a)
                {
                    // charged lepton velocity in the dilepton rest frame
                    const double mlH = a.mlH;
                    const double mlH2 = mlH * mlH;
                    const double NF = a.NF;

                    _vv[0] = NF * 2.0 * (
                                (1.0 + mlH2) * (std::norm(a.a_0) + 16.0 * std::norm(a.a_0_T))
                              + 2.0 * mlH2 * std::norm(a.a_t)
                              + 2.0 * std::norm(a.a_P)
                              + 4.0 * mlH * std::real(a.a_t * std::conj(a.a_P))
                              - 16.0 * mlH * std::real(a.a_0_T * std::conj(a.a_0))
                             );

                    _vv[1] = NF * 2.0 * (1.0 - mlH2) * ( - std::norm(a.a_0) + 16.0 * std::norm(a.a_0_T) );

                    _vv[2] = - NF * 8.0 * std::real(
                                mlH * (mlH * a.a_t + a.a_P) * std::conj(a.a_0)
                              - 4.0 * (mlH * a.a_t + a.a_P) * std::conj(a.a_0_T)
                            );

                    _vv[3] = NF * (
                                (3.0 + mlH2) * (std::norm(a.a_para) + std::norm(a.a_perp)) / 2.0
                              + 8.0 * (1.0 + 3.0 * mlH2) * (std::norm(a.a_para_T) + std::norm(a.a_perp_T))
                              - 16.0 * mlH * std::real(a.a_para_T * std::conj(a.a_para) + a.a_perp_T * std::conj(a.a_perp))
                            );

                    _vv[4] = NF * (1.0 - mlH2) * (
                                (std::norm(a.a_para) + std::norm(a.a_perp)) / 2.0
                              - 8.0 * (std::norm(a.a_para_T) + std::norm(a.a_perp_T))
                            );

                    _vv[5] = NF * 4.0 * std::real(
                              - a.a_para * std::conj(a.a_perp)
                              - 16.0 * mlH2 * a.a_para_T * std::conj(a.a_perp_T)
                              + 4.0 * mlH * (a.a_perp_T * std::conj(a.a_para) + a.a_para_T * std::conj(a.a_perp))
                            );

                    _vv[6] = NF * (1.0 - mlH2) * (
                              - (std::norm(a.a_para) - std::norm(a.a_perp))
                              + 16.0 * (std::norm(a.a_para_T) - std::norm(a.a_perp_T))
                            );

                    _vv[7] = NF * 2.0 * (1.0 - mlH2) * std::imag( a.a_para * std::conj(a.a_perp));

                    _vv[8] = NF * std::sqrt(2.0) * (1.0 - mlH2) * std::real(
                                a.a_para * std::conj(a.a_0)
                              - 16.0 * a.a_para_T * std::conj(a.a_0_T)
                            );

                    _vv[9] = NF * 2.0 * std::sqrt(2.0) * std::real(
                              - a.a_perp * std::conj(a.a_0)
                              + a.a_para * mlH * std::conj(mlH * a.a_t + a.a_P)
                              - 16.0 * mlH2 * a.a_perp_T * std::conj(a.a_0_T)
                              + 4.0 * mlH * (a.a_0_T * std::conj(a.a_perp) + a.a_perp_T * std::conj(a.a_0))
                              - 4.0 * a.a_para_T * std::conj(mlH * a.a_t + a.a_P)
                            );

                    _vv[10] = NF * 2.0 * std::sqrt(2.0) * std::imag(
                              - a.a_para * std::conj(a.a_0)
                              + mlH * a.a_perp * std::conj(mlH * a.a_t + a.a_P)
                              + 4.0 * mlH * ( a.a_0_T * std::conj(a.a_para) - a.a_para_T * std::conj(a.a_0))
                              + 4.0 * a.a_perp_T * std::conj(mlH * a.a_t + a.a_P)
                            );

                    _vv[11] = NF * std::sqrt(2.0) * (1.0 - mlH2) * std::imag( a.a_perp * std::conj(a.a_0));
                }

                AngularObservables()
                {
                    _vv.fill(0.0);
                }

                AngularObservables(const std::array<double, 12> & vv) :
                    _vv(vv)
                {
                }

                inline double vv10()    const  { return _vv[0]; }  // J_1c in B->K* ll literature
                inline double vv20()    const  { return _vv[1]; }  // J_2c
                inline double vv30()    const  { return _vv[2]; }  // J_6c
                inline double vv1T()    const  { return _vv[3]; }  // J_1s
                inline double vv2T()    const  { return _vv[4]; }  // J_2s
                inline double vv3T()    const  { return _vv[5]; }  // J_6s
                inline double vv4T()    const  { return _vv[6]; }  // J_3
                inline double vv5T()    const  { return _vv[7]; }  // J_9
                inline double vv10T()   const  { return _vv[8]; }  // J_4
                inline double vv20T()   const  { return _vv[9]; }  // J_5
                inline double vv30T()   const  { return _vv[10];}  // J_7
                inline double vv40T()   const  { return _vv[11];}  // J_8

                // longitudinal polarization amplitude
                inline double normalized_amplitude_polarization_L() const
                {
                    return  vv10() - vv20() / 3.0;
                }

                // transverse polarization amplitude
                inline double normalized_amplitude_polarization_T() const
                {
                    return  2.0 * (vv1T() - vv2T() / 3.0);
                }

                // redefined decay width
                inline double normalized_decay_width() const
                {
                    return 3.0 / 4.0 * ( normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T() );
                }

                // polarization fraction
                inline double f_L() const
                {
                    return  normalized_amplitude_polarization_L() / (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T());
                }

                // polarization fraction from cos(theta_l) distribution; identical to F_L in the SM and the limit m_l -> 0.
                inline double ftilde_L() const
                {
                    // (1 - 3 Ftilde_L)  == 16/3 (S2s + S2c/2)
                    return 1.0 / 3.0 - 16.0 / 9.0 * (vv2T() + vv20() / 2.0) / (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T());
                }

                // a_fb leptonic
                inline double a_fb_leptonic() const
                {
                    return  (vv3T() + vv30() / 2.0) / (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T());
                }

                // transverse azimuthal asymmetries
                inline double a_c_1() const
                {
                    return  4.0 * vv4T() / (3.0 * (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T()));
                }

                inline double a_c_2() const
                {
                    return  vv20T() / (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T());
                }

                inline double a_c_3() const
                {
                    return  vv10T() / (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T());
                }

                // T-odd CP asymmetries
                inline double a_t_1() const
                {
                    return  4.0 * vv5T() / (3.0 * (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T()));
                }

                inline double a_t_2() const
                {
                    return  vv30T() / (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T());
                }

                inline double a_t_3() const
                {
                    return  vv40T() / (normalized_amplitude_polarization_L() + normalized_amplitude_polarization_T());
                }
        };
    }

    class BToVectorLeptonNeutrino::IntermediateResult :
        public CacheableObservable::IntermediateResult
    {
        public:
            b_to_vec_l_nu::AngularObservables ao;

            IntermediateResult()
            {
            }

            ~IntermediateResult() = default;
    };
}

#endif
