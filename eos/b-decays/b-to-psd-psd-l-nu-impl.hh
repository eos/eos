/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Florian Herren
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

#ifndef EOS_GUARD_EOS_B_DECAYS_B_TO_PSD_PSD_L_NU_IMPL_HH
#define EOS_GUARD_EOS_B_DECAYS_B_TO_PSD_PSD_L_NU_IMPL_HH 1

#include <eos/observable.hh>
#include <eos/b-decays/b-to-psd-psd-l-nu.hh>
#include <eos/maths/angular-integrals.hh>
#include <eos/maths/complex.hh>
#include <eos/maths/power-of.hh>

#include <array>
#include <iostream>

namespace eos
{
    namespace b_to_psd_psd_l_nu
    {
        struct Amplitudes
        {
            std::array<complex<double>, 4> f_perp;
            std::array<complex<double>, 4> f_para;
            std::array<complex<double>, 4> f_long;
            std::array<complex<double>, 4> f_time;
            double q2;
            double betaL;
            double betaPi;
            double pref;
        };

        class AngularObservables
        {
            private:
                std::array<std::array<double, 5>, 9> _m;
                std::array<std::array<std::array<double, 5>, 4>, 4> _ints00, _ints11, _ints01;

            public:
                friend class BToPPLeptonNeutrino;
                friend class Implementation<BToPPLeptonNeutrino>;

                AngularObservables(const Amplitudes & a)
                {
                    _m.fill({0.0, 0.0, 0.0, 0.0, 0.0});
                    init_mats();

                    for (unsigned l = 0; l < 4; l++)
                    {
                        for (unsigned m = 0; m < 4; m++)
                        {
                            double amppm = (1.0 - a.betaL / 4.0) * power_of<2>(a.betaPi) * (std::real(a.f_perp[l] * std::conj(a.f_perp[m])) + std::real(a.f_para[l] * std::conj(a.f_para[m])));
                            double amp00 = ((1.0 - a.betaL) * std::real(a.q2 * a.f_time[l] * std::conj(a.q2 * a.f_time[m])) + (1.0 - a.betaL / 2.0) * std::real(a.f_long[l] * std::conj(a.f_long[m])));

                            for (unsigned i = 0; i < 5; i++)
                            {
                                _m[0][i] += amppm * _ints11[l][m][i];
                                _m[0][i] += amp00 * _ints00[l][m][i];
                            }
                        }
                    }

                    for (unsigned l = 0; l < 4; l++)
                    {
                        for (unsigned m = 0; m < 4; m++)
                        {
                            double amppm = a.betaL / 4.0 * power_of<2>(a.betaPi) * (std::real(a.f_perp[l] * std::conj(a.f_perp[m])) + std::real(a.f_para[l] * std::conj(a.f_para[m])));
                            double amp00 = -(a.betaL / 2.0 *  std::real(a.f_long[l] * std::conj(a.f_long[m])));

                            for (unsigned i = 0; i < 5; i++)
                            {
                                _m[1][i] += amppm * _ints11[l][m][i];
                                _m[1][i] += amp00 * _ints00[l][m][i];
                            }
                        }
                    }

                    for (unsigned l = 1; l < 4; l++)
                    {
                        for (unsigned m = 1; m < 4; m++)
                        {
                            double amppm = a.betaL / 2.0 * power_of<2>(a.betaPi) * (std::real(a.f_perp[l] * std::conj(a.f_perp[m])) - std::real(a.f_para[l] * std::conj(a.f_para[m])));

                            for (unsigned i = 0; i < 5; i++)
                            {
                                _m[2][i] += amppm * _ints11[l][m][i];
                            }
                        }
                    }

                    // Additional minus signs in front of ampmix are due to P_l^1 vs d P_l / dx
                    for (unsigned l = 0; l < 4; l++)
                    {
                        for (unsigned m = 1; m < 4; m++)
                        {
                            double ampmix = a.betaL * a.betaPi * std::real(a.f_long[l] * std::conj(a.f_para[m]));

                            for (unsigned i = 0; i < 5; i++)
                            {
                                _m[3][i] -= ampmix * _ints01[l][m][i];
                            }
                        }
                    }

                    for (unsigned l = 0; l < 4; l++)
                    {
                        for (unsigned m = 1; m < 4; m++)
                        {
                            double ampmix = 2.0 * a.betaPi * (std::real(a.f_long[l] * std::conj(a.f_perp[m])) + (1.0 - a.betaL) * std::real(a.q2 * a.f_time[l] * std::conj(a.f_para[m])));

                            for (unsigned i = 0; i < 5; i++)
                            {
                                _m[4][i] -= ampmix * _ints01[l][m][i];
                            }
                        }
                    }

                    for (unsigned l = 0; l < 4; l++)
                    {
                        for (unsigned m = 0; m < 4; m++)
                        {
                            double amppm = 2.0 * power_of<2>(a.betaPi) * std::real(a.f_perp[l] * std::conj(a.f_para[m]));
                            double amp00 = -2.0 * (1.0 - a.betaL) *  std::real(a.q2 * a.f_time[l] * std::conj(a.f_long[m]));

                            for (unsigned i = 0; i < 5; i++)
                            {
                                _m[5][i] += amppm * _ints11[l][m][i];
                                _m[5][i] += amp00 * _ints00[l][m][i];
                            }
                        }
                    }

                    for (unsigned l = 0; l < 4; l++)
                    {
                        for (unsigned m = 1; m < 4; m++)
                        {
                            double ampmix = -2.0 * a.betaPi * (std::imag(a.f_long[l] * std::conj(a.f_para[m])) - (1.0 - a.betaL) * std::imag(a.q2 * a.f_time[l] * std::conj(a.f_perp[m])));

                            for (unsigned i = 0; i < 5; i++)
                            {
                                _m[6][i] -= ampmix * _ints01[l][m][i];
                            }
                        }
                    }

                    for (unsigned l = 0; l < 4; l++)
                    {
                        for (unsigned m = 1; m < 4; m++)
                        {
                            double ampmix = a.betaL * a.betaPi * std::imag(a.f_long[l] * std::conj(a.f_perp[m]));

                            for (unsigned i = 0; i < 5; i++)
                            {
                                _m[7][i] -= ampmix * _ints01[l][m][i];
                            }
                        }
                    }

                    for (unsigned l = 1; l < 4; l++)
                    {
                        for (unsigned m = 1; m < 4; m++)
                        {
                            double amppm = -a.betaL * power_of<2>(a.betaPi) * std::imag(a.f_perp[l] * std::conj(a.f_para[m]));

                            for (unsigned i = 0; i < 5; i++)
                            {
                                _m[8][i] += amppm * _ints11[l][m][i];
                            }
                        }
                    }

                    // Multiply by prefactor
                    for (unsigned m = 0; m < 9; m++)
                    {
                        for (unsigned i = 0; i < 5; i++) _m[m][i] *= a.pref;
                    }
                }

                AngularObservables()
                {
                    _m.fill({0.0, 0.0, 0.0, 0.0, 0.0});
                    init_mats();
                }

                AngularObservables(const std::array<std::array<double, 5>, 9> & m) :
                    _m(m)
                {
                    init_mats();
                }

                inline void init_mats()
                {
                    for (unsigned l = 0; l < 4; l++)
                    {
                        for (unsigned m = 0; m < 4; m++)
                        {
                            double pref = std::sqrt((2 * l + 1) * (2 * m + 1));
                            for (unsigned i = 0; i < 5; i++)
                            {
                                _ints00[l][m][i] = pref * three_legendre_integral(l, 0, m, 0, i, 0);
                                _ints11[l][m][i] = pref * three_legendre_integral(l, 1, m, 1, i, 0);
                                _ints01[l][m][i] = pref * three_legendre_integral(l, 0, m, 1, i, 0);
                            }
                        }
                    }
                }

                inline double M1(unsigned i) const  { return _m[0][i]; }
                inline double M2(unsigned i) const  { return _m[1][i]; }
                inline double M3(unsigned i) const  { return _m[2][i]; }
                inline double M4(unsigned i) const  { return _m[3][i]; }
                inline double M5(unsigned i) const  { return _m[4][i]; }
                inline double M6(unsigned i) const  { return _m[5][i]; }
                inline double M7(unsigned i) const  { return _m[6][i]; }
                inline double M8(unsigned i) const  { return _m[7][i]; }
                inline double M9(unsigned i) const  { return _m[8][i]; }

                inline double double_differential_decay_width() const
                {
                    return M1(0) - M2(0) / 3.0;
                }

                inline double double_differential_mesonic_afb() const
                {
                    return M1(1) - M2(1) / 3.0;
                }
        };
    }
}

#endif
