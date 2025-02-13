/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Florian Herren
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

#include <eos/scattering/parametric-gmkprdey2011.hh>
#include <eos/maths/omnes-factor-impl.hh>
#include <eos/maths/outer-function.hh>


#include <functional>

namespace eos
{

    double GMKPRDEY2011ScatteringAmplitudes::_calc_w(const double & s, const double & s0) const
    {
        if (s > s0)
            throw InternalError("The real conformal mapping is used above threshold: " + stringify(s) + " > " + stringify(s0));

        return (std::sqrt(s) - std::sqrt(s0 - s)) / (std::sqrt(s) + std::sqrt(s0 - s));
    }

    double GMKPRDEY2011ScatteringAmplitudes::_calc_s(const complex<double> & z, const double & sp, const double & s0) const
    {
        if (s0 > sp)
            throw InternalError("The inverse conformal mapping is used with s_+ < s_0: " + stringify(sp) + " < " + stringify(s0));

        return std::abs((-4.0 * sp * z + s0 * power_of<2>(1.0 + z)) / power_of<2>(z - 1.0));
    }

    complex<double> GMKPRDEY2011ScatteringAmplitudes::_calc_z(const double & s, const double & sp, const double & s0) const
    {
        if (s0 > sp)
            throw InternalError("The conformal mapping is used with s_+ < s_0: " + stringify(sp) + " < " + stringify(s0));

        return (std::sqrt(complex<double>(sp - s, 0.0)) - std::sqrt(complex<double>(sp - s0, 0.0))) / (std::sqrt(complex<double>(sp - s, 0.0)) + std::sqrt(complex<double>(sp - s0, 0.0)));
    }

    // S0 wave
    double GMKPRDEY2011ScatteringAmplitudes::_phase_S0(const double & s) const
    {
        double mpi2 = _mPi * _mPi;
        double mk2 = _mK * _mK;
        double meta2 = _mEta * _mEta;
        double sqrts = std::sqrt(s);

        if (s <= _sM_S0)
        {
            double k = std::sqrt(s / 4.0 - mpi2);
            double ws = _calc_w(s, 4.0 * mk2);
            if (s <= 0.6) // Avoid division by 0 for s = 4*mpi2 and s approx 0.7
            {
                 return std::atan(2.0 * k / sqrts * (s - mpi2 / 2.0) / mpi2 / (_mPi / sqrts + _params_S0[0] + _params_S0[1] * ws + _params_S0[2] * power_of<2>(ws) + _params_S0[3] * power_of<3>(ws)));
            }
            else
            {
                 return M_PI / 2.0 - std::atan(sqrts * mpi2 / 2.0 / k / (s - mpi2 / 2.0) * (_mPi / sqrts + _params_S0[0] + _params_S0[1] * ws + _params_S0[2] * power_of<2>(ws) + _params_S0[3] * power_of<3>(ws)));
            }
        }
        else if (s <= 4.0 * mk2)
        {
            // Quantities entering the derivative of the phase at s = _sM_s0
            double mpi3 = mpi2 * _mPi;
            double sqrtsM = std::sqrt(_sM_S0);
            double kMpi = std::sqrt(_sM_S0 / 4.0 - mpi2);
            double sqrts0sM = std::sqrt(4.0 * mk2 - _sM_S0);
            double wM = _calc_w(_sM_S0, 4.0 * mk2);
            double Bsum = std::sqrt(mpi2 / _sM_S0) + _params_S0[0] + _params_S0[1] * wM + _params_S0[2] * power_of<2>(wM) + _params_S0[3] * power_of<3>(wM);

            // Individual pieces entering the derivative
            double x1 = mpi2 * ( (mpi2 - 2.0 * _sM_S0) * _sM_S0 - 4.0 * power_of<2>(kMpi) * (mpi2 + 2.0 * _sM_S0)) * Bsum / (8.0 * sqrtsM * power_of<2>(mpi2 - 2.0 * _sM_S0) * power_of<3>(kMpi));
            double x2 = mpi3 / kMpi / 2.0 / _sM_S0 / (mpi2 - 2.0 * _sM_S0);
            double x3 = -mpi2 / 2.0 / (mpi2 - 2.0 * _sM_S0) / kMpi / sqrts0sM * ((1.0 - wM) * sqrts0sM + (1.0 + wM) * sqrtsM) / (sqrtsM + sqrts0sM);
            double x4 = _sM_S0 * power_of<2>(mpi2 * Bsum / (mpi2 - 2.0 * _sM_S0) / kMpi);

            // Derivative of the phase at s = _sM_s0
            double delpM = -(x1 + x2 + x3 * (_params_S0[1] + 2.0 * wM * _params_S0[2] + 3.0 * power_of<2>(wM) * _params_S0[3])) / (1 + x4);

            // Remaining quantities entering the parametrization of the phase
            double absk2 = std::sqrt(mk2 - s / 4.0);
            double kM = std::sqrt(mk2 - _sM_S0 / 4.0);
            double mk3 = mk2 * _mK;

            return _params_S0[4] * power_of<2>( 1.0 - absk2 / kM ) + _phase_S0(_sM_S0) * absk2 / kM * (2.0 - absk2 / kM) + absk2 * (kM - absk2) * ( 8.0 * delpM + _params_S0[5] * (kM - absk2) / mk3);
        }
        else if (s <= 4.0 * meta2)
        {
            double k22norm = (s / 4.0 / mk2 - 1.0);

            return _params_S0[4] + _params_S0[6] * k22norm + _params_S0[7] * power_of<2>(k22norm);
        }
        else if (s <= 2.0164)
        {
            double k22norm = (s / 4.0 / mk2 - 1.0);
            double k32norm = (s / 4.0 / meta2 - 1.0);

            return _params_S0[4] + _params_S0[6] * k22norm + _params_S0[7] * power_of<2>(k22norm) + _params_S0[8] * k32norm;
        }
        else
        {
            return 2 * M_PI + (_phase_S0(2.0164) - 2.0 * M_PI) * 2.0 / (1.0 + std::pow(s / 2.0164, _cont_pow_S0));
        }
    }

    // P wave
    double GMKPRDEY2011ScatteringAmplitudes::_phase_P1(const double & s) const
    {
        double mpi2 = _mPi * _mPi;
        double mk2 = _mK * _mK;
        double sqrts = std::sqrt(s);

        if (s <= 4.0 * mk2)
        {
            double mrho2 = _mRho * _mRho;
            double mpi3 = _mPi * mpi2;
            double k = std::sqrt(s / 4.0 - mpi2);
            if ( s <= 0.5 ) // Avoid divide by 0 for s = 4*mpi2 and s = mrho2
            {
                return std::atan(2.0 * power_of<3>(k) / sqrts / (mrho2 - s) / (2.0 * mpi3 / mrho2 / sqrts + _params_P1[0] + _params_P1[1] * _calc_w(s, _s0_P1)));
            }
            else
            {
                return M_PI / 2.0 - std::atan(sqrts / 2.0 / power_of<3>(k) * (mrho2 - s) * (2.0 * mpi3 / mrho2 / sqrts + _params_P1[0] + _params_P1[1] * _calc_w(s, _s0_P1)));
            }
        }
        else if (s <= 2.0164)
        {
            return _phase_P1(4.0 * mk2) + _params_P1[2] * (sqrts / 2.0 / _mK - 1.0) + _params_P1[3] * power_of<2>(sqrts / 2.0 / _mK - 1.0);
        }
        else
        {
            return M_PI + (_phase_P1(2.0164) - M_PI) * 2.0 / (1.0 + std::pow(s / 2.0164, _cont_pow_P1));
        }
    }

    // D0 wave
    double GMKPRDEY2011ScatteringAmplitudes::_phase_D0(const double & s) const
    {
        double mpi2 = _mPi * _mPi;
        double mf22 = _mF2 * _mF2;
        double sqrts = std::sqrt(s);
        double k = std::sqrt(s / 4.0 - mpi2);

        if (s <= _s0_D0)
        {
            return std::atan(2.0 * power_of<5>(k) / sqrts / (mf22 - s) / mpi2 / (_params_D0[0] + _params_D0[1] * _calc_w(s, _s0_D0)));
        }
        else if (s <= 2.0164)
        {
            double Bh0 = _params_D0[0] + _params_D0[1] - _params_D0[2] * _calc_w(_s0_D0, _sh_D0);
            return M_PI / 2.0 - std::atan(sqrts / 2.0 / power_of<5>(k) * (mf22 - s) * mpi2 * (Bh0 + _params_D0[2] * _calc_w(s, _sh_D0)));
        }
        else
        {
            return M_PI + (_phase_D0(2.0164) - M_PI) * 2.0 / (1.0 + std::pow(s / 2.0164, _cont_pow_D0));
        }
    }

    GMKPRDEY2011ScatteringAmplitudes::GMKPRDEY2011ScatteringAmplitudes(const Parameters & p, const Options &) :
        _params_S0 {    UsedParameter(p[_par_name("S0", "B", 0)],   *this),
                        UsedParameter(p[_par_name("S0", "B", 1)],   *this),
                        UsedParameter(p[_par_name("S0", "B", 2)],   *this),
                        UsedParameter(p[_par_name("S0", "B", 3)],   *this),
                        UsedParameter(p[_par_name("S0", "d", 0)],   *this),
                        UsedParameter(p[_par_name("S0", "c")],      *this),
                        UsedParameter(p[_par_name("S0", "B")],      *this),
                        UsedParameter(p[_par_name("S0", "C")],      *this),
                        UsedParameter(p[_par_name("S0", "D")],      *this)
        },
        _params_P1 {    UsedParameter(p[_par_name("P1", "B", 0)],   *this),
                        UsedParameter(p[_par_name("P1", "B", 1)],   *this),
                        UsedParameter(p[_par_name("P1", "lam", 1)], *this),
                        UsedParameter(p[_par_name("P1", "lam", 2)], *this)
        },
        _params_D0 {    UsedParameter(p[_par_name("D0", "B", 0)],   *this),
                        UsedParameter(p[_par_name("D0", "B", 1)],   *this),
                        UsedParameter(p[_par_name("D0", "Bh", 1)],  *this)
        },
        _mPi(UsedParameter(p["mass::pi^+@GMKPRDEY2011"], *this)),
        _mK(UsedParameter(p["mass::K_u@GMKPRDEY2011"], *this)),
        _mEta(UsedParameter(p["mass::eta@GMKPRDEY2011"], *this)),
        _mRho(UsedParameter(p["mass::rho^0@GMKPRDEY2011"], *this)),
        _mF2(UsedParameter(p["mass::f_2@GMKPRDEY2011"], *this)),
        _sM_S0(UsedParameter(p[_par_name("S0", "sM")],  *this)),
        _s0_P1(UsedParameter(p[_par_name("P1", "s0")],  *this)),
        _s0_D0(UsedParameter(p[_par_name("D0", "s0")],  *this)),
        _sh_D0(UsedParameter(p[_par_name("D0", "sh")],  *this)),
        _cont_pow_S0(UsedParameter(p[_par_name("S0", "n")],  *this)),
        _cont_pow_P1(UsedParameter(p[_par_name("P1", "n")],  *this)),
        _cont_pow_D0(UsedParameter(p[_par_name("D0", "n")],  *this)),
        _intervals_P1({4 * _mPi * _mPi , 0.5 , 1.0 , 2.0 }),
        _intervals_D0({4 * _mPi * _mPi , 0.7, 1.1, 1.45, 2.0 }),
        _f_phase_P1(std::bind(&GMKPRDEY2011ScatteringAmplitudes::_phase_P1, this, std::placeholders::_1)),
        _f_phase_D0(std::bind(&GMKPRDEY2011ScatteringAmplitudes::_phase_D0, this, std::placeholders::_1)),
        _omnes_P1(_intervals_P1, _f_phase_P1, 0.0),
        _omnes_D0(_intervals_D0, _f_phase_D0, 0.0)
    {
    }

    GMKPRDEY2011ScatteringAmplitudes::~GMKPRDEY2011ScatteringAmplitudes() = default;

    ScatteringAmplitudes<PPToPP> *
    GMKPRDEY2011ScatteringAmplitudes::make(const Parameters & parameters, const Options & options)
    {
        return new GMKPRDEY2011ScatteringAmplitudes(parameters, options);
    }

    QualifiedName
    GMKPRDEY2011ScatteringAmplitudes::_par_name(const std::string & partial_wave, const std::string & par_name, unsigned idx) const
    {
        return QualifiedName(stringify(PiPiToPiPi::label) + "::" + partial_wave + "_" + par_name + "_" + stringify(idx) + "@GMKPRDEY2011");
    }

    QualifiedName
    GMKPRDEY2011ScatteringAmplitudes::_par_name(const std::string & partial_wave, const std::string & par_name) const
    {
        return QualifiedName(stringify(PiPiToPiPi::label) + "::" + partial_wave + "_" + par_name + "@GMKPRDEY2011");
    }

    complex<double> GMKPRDEY2011ScatteringAmplitudes::scattering_amplitude(const double & s, const unsigned & l, const IsospinRepresentation & i) const
    {
        if (s <= 4 * _mPi * _mPi)
        {
            return 0.0;
        }

        double rho = std::sqrt(1 - 4 * _mPi * _mPi / s);

        if ((l == 0) && (i == IsospinRepresentation::zero))
        {
            double del = _phase_S0(s);
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

    complex<double> GMKPRDEY2011ScatteringAmplitudes::omnes_factor(const double & s, const unsigned & l, const IsospinRepresentation & i) const
    {
        if ((l == 0) && (i == IsospinRepresentation::zero))
        {
            throw InternalError("Current Omnes factor solution strategy does not allow for phases exceeding 2 Pi! Consider implementing coupled-channel treatment!");
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
            return 1.0;
        }
    }

    // Note: all our omnes factors go like 1/s for large s. Thus we need to take out a factor of (1 - z)^2 which would cause issues with the integration
    complex<double> GMKPRDEY2011ScatteringAmplitudes::omnes_outer_function(const double & s, const double & sp, const double & s0, const unsigned & npoints, const unsigned & l, const IsospinRepresentation & i) const
    {
        // Point to extract asymptotic behaviour at.
        const double sM = 1000000.0;

        if ((l == 0) && (i == IsospinRepresentation::zero))
        {
            throw InternalError("Current Omnes factor solution strategy does not allow for phases exceeding 2 Pi! Consider implementing coupled-channel treatment!");
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

            return power_of<2>(zeval - 1.0) * outer(integrand, zeval, npoints);
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

            return power_of<2>(zeval - 1.0) * outer(integrand, zeval, npoints);
        }
        else
        {
            return complex<double>(1.0, 0.0);
        }
    }

    Diagnostics
    GMKPRDEY2011ScatteringAmplitudes::diagnostics() const
    {
        Diagnostics results;

        results.add({ _calc_w(0.0,  _s0_P1), "w_P1(s =  0.0)" });
        results.add({ _calc_w(1.0,  _s0_P1), "w_P1(s =  1.0)" });

        results.add({ _phase_S0(0.25),  "del_S0(s =  0.25)" });
        results.add({ _phase_S0(0.72),  "del_S0(s =  0.72)" });
        results.add({ _phase_S0(0.9),   "del_S0(s =  0.9)"  });
        results.add({ _phase_S0(1.44),  "del_S0(s =  1.44)" });
        results.add({ _phase_S0(4.0),   "del_S0(s =  4.0)"  });

        results.add({ _phase_P1(0.25),  "del_P1(s =  0.25)" });
        results.add({ _phase_P1(0.9),   "del_P1(s =  0.9)"  });
        results.add({ _phase_P1(1.0),   "del_P1(s =  1.0)"  });
        results.add({ _phase_P1(4.0),   "del_P1(s =  4.0)"  });

        results.add({ _phase_D0(0.25),  "del_D0(s =  0.25)" });
        results.add({ _phase_D0(0.9),   "del_D0(s =  0.9)"  });
        results.add({ _phase_D0(1.44),  "del_D0(s =  1.44)" });
        results.add({ _phase_D0(4.0),   "del_D0(s =  4.0)"  });

        return results;
    }

    const std::set<ReferenceName> GMKPRDEY2011ScatteringAmplitudes::references
    {
        "GMKPRDEY:2011A"_rn
    };

    const std::vector<OptionSpecification> GMKPRDEY2011ScatteringAmplitudes::options
    {
    };

    std::vector<OptionSpecification>::const_iterator
    GMKPRDEY2011ScatteringAmplitudes::begin_options()
    {
        return options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    GMKPRDEY2011ScatteringAmplitudes::end_options()
    {
        return options.cend();
    }
}
