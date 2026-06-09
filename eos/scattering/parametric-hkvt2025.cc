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

#include <eos/scattering/parametric-hkvt2025.hh>
#include <eos/maths/omnes-factor-impl.hh>
#include <eos/maths/outer-function.hh>

#include <functional>
#include <iostream>

namespace eos
{
    double HKVT2025ScatteringAmplitudes::_calc_w(const double & s, const double & s0) const
    {
        if (s > s0)
            throw InternalError("The real conformal mapping is used above threshold: " + stringify(s) + " > " + stringify(s0));

        return (std::sqrt(s) - std::sqrt(s0 - s)) / (std::sqrt(s) + std::sqrt(s0 - s));
    }

    double HKVT2025ScatteringAmplitudes::_calc_s(const complex<double> & z, const double & sp, const double & s0) const
    {
        if (s0 > sp)
            throw InternalError("The inverse conformal mapping is used with s_+ < s_0: " + stringify(sp) + " < " + stringify(s0));

        return std::abs((-4.0 * sp * z + s0 * power_of<2>(1.0 + z)) / power_of<2>(z - 1.0));
    }

    complex<double> HKVT2025ScatteringAmplitudes::_calc_z(const double & s, const double & sp, const double & s0) const
    {
        if (s0 > sp)
            throw InternalError("The conformal mapping is used with s_+ < s_0: " + stringify(sp) + " < " + stringify(s0));

        return (std::sqrt(complex<double>(sp - s, 0.0)) - std::sqrt(sp - s0)) / (std::sqrt(complex<double>(sp - s, 0.0)) + std::sqrt(sp - s0));
    }

    // S0 wave Omnes factor from [DHK:2015A] and [RHK:2018A]
    complex<double> HKVT2025ScatteringAmplitudes::_omnes_S0(const double & s) const
    {
        const static thread_local OmnesInterpolation Om11_interpolation = {
                {
                    #include "svalues-om-dhk.hh"
                },
                {
                    #include "reom11-dhk.hh"
                },
                {
                    #include "imom11-dhk.hh"
                }
        };

        const static thread_local OmnesInterpolation Om12_interpolation = {
                {
                    #include "svalues-om-dhk.hh"
                },
                {
                    #include "reom12-dhk.hh"
                },
                {
                    #include "imom12-dhk.hh"
                }
        };

        if (std::sqrt(s) <= 9.88873)
        {
            return complex<double>(_Gamma_pi_0) * Om11_interpolation(std::sqrt(s)) + complex<double>(_Gamma_K_0) * 2.0 / std::sqrt(3.0) * Om12_interpolation(std::sqrt(s));
        }
        else
        {
            return (complex<double>(_Gamma_pi_0) * Om11_interpolation(9.88873) + complex<double>(_Gamma_K_0) * 2.0 / std::sqrt(3.0) * Om12_interpolation(9.88873)) * 9.88873 * 9.88873 / s;
        }
    }

    // P1 wave phase-shift from [CHS:2018A]
    double HKVT2025ScatteringAmplitudes::_phase_P1(const double & s) const
    {
        const static thread_local PhaseInterpolation CHS_interpolation = {
                {
                    #include "svalues-chs.hh"
                },
                {
                    #include "delvalues-chs.hh"
                }
        };

        if (s <= 1.69)
        {
            return CHS_interpolation(s);
        }
        else
        {
            return M_PI + (CHS_interpolation(1.69) - M_PI) * 2.0 / (1.0 + std::pow(s / 1.69, _cont_pow_P1));
        }

    }

    // D0 wave phase-shift from [GMKPRDEY:2011A]
    double HKVT2025ScatteringAmplitudes::_phase_D0(const double & s) const
    {
        double mpi2 = _mPi * _mPi;
        double mf22 = _mF2 * _mF2;
        double sqrts2 = std::sqrt(s) / 2.0;
        double k = std::sqrt(s / 4.0 - mpi2);

        if (s <= _s0_D0)
        {
            return std::atan( power_of<5>(k) / sqrts2 / (mf22 - s) / mpi2 / (_params_D0[0] + _params_D0[1] * _calc_w(s, _s0_D0)) );
        }
        else if (s <= 2.0164)
        {
            double Bh0 = _params_D0[0] + _params_D0[1] - _params_D0[2] * _calc_w(_s0_D0, _sh_D0);
            return M_PI/2 - std::atan( sqrts2 / power_of<5>(k) * (mf22 - s) * mpi2 * (Bh0 + _params_D0[2] * _calc_w(s, _sh_D0)) );
        }
        else
        {
            return M_PI + (_phase_D0(2.0164) - M_PI) * 2.0 / (1.0 + std::pow(s / 2.0164, _cont_pow_D0));
        }
    }

    HKVT2025ScatteringAmplitudes::HKVT2025ScatteringAmplitudes(const Parameters & p, const Options &) :
        _params_D0 {    UsedParameter(p[_par_name("D0", "B", 0)],   *this),
                        UsedParameter(p[_par_name("D0", "B", 1)],   *this),
                        UsedParameter(p[_par_name("D0", "Bh", 1)],  *this)
        },
        _mPi(UsedParameter(p["mass::pi^+@GMKPRDEY2011"], *this)),
        _mF2(UsedParameter(p["mass::f_2@GMKPRDEY2011"], *this)),
        _mOmega(UsedParameter(p["mass::omega@GMKPRDEY2011"], *this)),
        _GammaOmega(UsedParameter(p["width::omega@GMKPRDEY2011"], *this)),
        _kappa(UsedParameter(p["mixing::kappaEM@GMKPRDEY2011"], *this)),
        _s0_D0(UsedParameter(p[_par_name("D0", "s0")],  *this)),
        _sh_D0(UsedParameter(p[_par_name("D0", "sh")],  *this)),
        _cont_pow_P1(UsedParameter(p[_par_name("P1", "n")],  *this)),
        _cont_pow_D0(UsedParameter(p[_par_name("D0", "n")],  *this)),
        _Gamma_pi_0(UsedParameter(p["pipi->pipi::Gamman0_pi@HKvT2025"], *this)),
        _Gamma_K_0(UsedParameter(p["pipi->pipi::Gamman0_K@HKvT2025"], *this)),
        _intervals_P1({4.0 * _mPi * _mPi , 0.5 , 1.0 , 2.0 }),
        _intervals_D0({4.0 * _mPi * _mPi , 0.7, 1.1, 1.45, 2.0 }),
        _f_phase_P1(std::bind(&HKVT2025ScatteringAmplitudes::_phase_P1, this, std::placeholders::_1)),
        _f_phase_D0(std::bind(&HKVT2025ScatteringAmplitudes::_phase_D0, this, std::placeholders::_1)),
        _omnes_P1(_intervals_P1, _f_phase_P1, 0.0),
        _omnes_D0(_intervals_D0, _f_phase_D0, 0.0)
    {
    }

    HKVT2025ScatteringAmplitudes::~HKVT2025ScatteringAmplitudes() = default;

    ScatteringAmplitudes<PPToPP> *
    HKVT2025ScatteringAmplitudes::make(const Parameters & parameters, const Options & options)
    {
        return new HKVT2025ScatteringAmplitudes(parameters, options);
    }

    QualifiedName
    HKVT2025ScatteringAmplitudes::_par_name(const std::string & partial_wave, const std::string & par_name, unsigned idx) const
    {
        return QualifiedName(stringify(PiPiToPiPi::label) + "::" + partial_wave + "_" + par_name + "_" + stringify(idx) + "@GMKPRDEY2011");
    }

    QualifiedName
    HKVT2025ScatteringAmplitudes::_par_name(const std::string & partial_wave, const std::string & par_name) const
    {
        return QualifiedName(stringify(PiPiToPiPi::label) + "::" + partial_wave + "_" + par_name + "@GMKPRDEY2011");
    }

    complex<double> HKVT2025ScatteringAmplitudes::scattering_amplitude(const double & s, const unsigned & l, const IsospinRepresentation & i) const
    {
        if (s <= 4 * _mPi * _mPi)
        {
            return 0.0;
        }

        const double rho = std::sqrt(1 - 4 * _mPi * _mPi / s);

        if ((l == 0) && (i == IsospinRepresentation::zero))
        {
            double del = std::arg(_omnes_S0(s));
            return std::exp(complex<double>(0, del)) * std::sin(del) / rho;
        }
        else if ((l == 1) && (i == IsospinRepresentation::one))
        {
            double del = _phase_P1(s);
            return std::exp(complex<double>(0, del)) * std::sin(del) / power_of<3>(rho);
        }
        else if ((l == 2) && (i == IsospinRepresentation::zero))
        {
            double del = _phase_D0(s);
            return std::exp(complex<double>(0, del)) * std::sin(del) / power_of<5>(rho);
        }
        else
        {
            return 0.0;
        }
    }

    complex<double> HKVT2025ScatteringAmplitudes::omnes_factor(const double & s, const unsigned & l, const IsospinRepresentation & i) const
    {
        if ((l == 0) && (i == IsospinRepresentation::zero))
        {
            return _omnes_S0(s);
        }
        else if ((l == 1) && (i == IsospinRepresentation::one))
        {
            return _omnes_P1(s);
        }
        else if ((l == 2) && (i == IsospinRepresentation::zero))
        {
            return _omnes_D0(s);
        }
        else
        {
            return 0.0;
        }
    }

    // Simplified isospin-breaking correction following [CHS:2018A]
    complex<double> HKVT2025ScatteringAmplitudes::isospin_breaking(const double & s, const unsigned & l, const IsospinRepresentation & i) const
    {
        if ( (l == 1) && (i == IsospinRepresentation::one) )
        {
            return 1.0 + s * _kappa / (_mOmega * _mOmega - s - complex<double>(0, _mOmega * _GammaOmega));
        }
        else
        {
            return 1.0;
        }
    }

    // Note: all our omnes factors go like 1/s for large s. Thus we need to take out a factor of (1 - z)^2 which would cause issues with the integration
    complex<double> HKVT2025ScatteringAmplitudes::omnes_outer_function(const double & s, const double & sp, const double & s0, const double & prec, const unsigned & l, const IsospinRepresentation & i) const
    {
        // Point to extract asymptotic behaviour at
        const double sM = 1000000.0;

        if ((l == 0) && (i == IsospinRepresentation::zero))
        {
            complex<double> zeval = _calc_z(s, sp, s0);

            std::function<complex<double>(const complex<double> &)> integrand = [&] (const complex<double> & z)
            {
                // Ensure we provide a real s to Omnes factor
                double sArg = _calc_s(z, sp, s0);

                if (!isfinite(sArg) || (sArg > sM))
                {
                    // Here we fix the constant Omega(s) ~ c/s for large s and take into account the rest of the conformal mapping
                    return sM * _omnes_S0(sM) / (-4.0 * sp * z + s0 * power_of<2>(z + 1.0));
                }
                else
                {
                    return _omnes_S0(sArg) / power_of<2>(z - 1.0);
                }
            };

            return power_of<2>(zeval - 1.0) * outer(integrand, zeval, prec);
        }
        else if ((s < sp) && (s0 < sp) && (l == 1) && (i == IsospinRepresentation::one))
        {
            complex<double> zeval = _calc_z(s, sp, s0);

            std::function<complex<double>(const complex<double> &)> integrand = [&] (const complex<double> & z)
            {
                // Ensure we provide a real s to Omnes factor
                double sArg = _calc_s(z, sp, s0);

                if (!isfinite(sArg) || (sArg > sM))
                {
                    // Here we fix the constant Omega(s) ~ c/s for large s and take into account the rest of the conformal mapping
                    return sM * _omnes_P1(sM) / (-4.0 * sp * z + s0 * power_of<2>(z + 1.0));
                }
                else
                {
                    return _omnes_P1(sArg) / power_of<2>(z - 1.0);
                }
            };

            return power_of<2>(zeval - 1.0) * outer(integrand, zeval, prec);
        }
        else if ((s < sp) && (s0 < sp) && (l == 2) && (i == IsospinRepresentation::zero))
        {
            complex<double> zeval = _calc_z(s, sp, s0);

            std::function<complex<double>(const complex<double> &)> integrand = [&] (const complex<double> & z)
            {
                // Ensure we provide a real s to Omnes factor
                double sArg = _calc_s(z, sp, s0);

                if (!isfinite(sArg) || (sArg > sM))
                {
                    // Here we fix the constant Omega(s) ~ c/s for large s and take into account the rest of the conformal mapping
                    return sM * _omnes_D0(sM) / (-4.0 * sp * z + s0 * power_of<2>(z + 1.0));
                }
                else
                {
                    return _omnes_D0(sArg) / power_of<2>(z - 1.0);
                }
            };

            return power_of<2>(zeval - 1.0) * outer(integrand, zeval, prec);
        }
        else
        {
            return complex<double>(1.0, 0.0);
        }
    }

    Diagnostics
    HKVT2025ScatteringAmplitudes::diagnostics() const
    {
        Diagnostics results;

        results.add({ _calc_w(0.0,  _s0_D0), "w_D0(s =  0.0)" });
        results.add({ _calc_w(1.0,  _s0_D0), "w_D0(s =  1.0)" });

        results.add({ _phase_P1(0.25),  "del_P1(s =  0.25)" });
        results.add({ _phase_P1(0.9),   "del_P1(s =  0.9)"  });
        results.add({ _phase_P1(1.44),  "del_P1(s =  1.44)" });
        results.add({ _phase_P1(4.0),   "del_P1(s =  4.0)" });

        results.add({ _phase_D0(0.25),  "del_D0(s =  0.25)" });
        results.add({ _phase_D0(0.9),   "del_D0(s =  0.9)"  });
        results.add({ _phase_D0(1.44),  "del_D0(s =  1.44)" });
        results.add({ _phase_D0(4.0),   "del_D0(s =  4.0)"  });

        return results;
    }

    const std::set<ReferenceName> HKVT2025ScatteringAmplitudes::references
    {
        "HKvT:2025A"_rn,
        "DHK:2015A"_rn,
        "RHK:2018A"_rn,
        "CHS:2018A"_rn,
        "GMKPRDEY:2011A"_rn
    };

    const std::vector<OptionSpecification> HKVT2025ScatteringAmplitudes::options
    {
    };

    std::vector<OptionSpecification>::const_iterator
    HKVT2025ScatteringAmplitudes::begin_options()
    {
        return options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    HKVT2025ScatteringAmplitudes::end_options()
    {
        return options.cend();
    }
}
