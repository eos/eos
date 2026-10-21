/*
 * Copyright (c) 2026 Danny van Dyk
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

#include <eos/maths/cauchy-transform.hh>
#include <eos/maths/faddeeva.hh>
#include <eos/maths/integrate.hh>
#include <eos/utils/exception.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <functional>
#include <limits>
#include <utility>
#include <vector>

namespace eos
{
    CauchyTransform::~CauchyTransform() = default;

    namespace
    {
        // Kernel::Gaussian(u, mu, sigma): S(zeta) = -i sqrt(pi / 2) / sigma w((zeta - mu) / (sqrt(2) sigma))
        class GaussianCauchyTransform : public CauchyTransform
        {
            private:
                double _mu    = 0.0;
                double _sigma = 1.0;

            public:
                virtual void
                prepare(std::span<const double> parameters) override
                {
                    _mu    = parameters[0];
                    _sigma = parameters[1];
                }

                virtual complex<double>
                evaluate(const complex<double> & zeta) const override
                {
                    return -1.0i * std::sqrt(M_PI / 2.0) / _sigma * faddeeva_w((zeta - _mu) / (std::sqrt(2.0) * _sigma));
                }
        };

        /*
         * Kernel::CrystalBall(u, mu, sigma, alpha_L, n_L, alpha_R, n_R): S = S_G + S_Delta, where Delta = R_CB - R_G
         * is smooth on the scale sigma. With zeta - mu = sigma (s + i g), and to first order in g,
         *   sigma S_Delta = pi H[delta](s) - i pi delta(s) + i g (pi H[delta'](s) - i pi delta'(s)),
         * where delta(t) = sigma Delta(sigma t) and pi H[f](s) = p.v. int dt f(t) / (s - t). The Hilbert
         * transforms depend only on the tail parameters, and are tabulated once per set of them.
         */
        class CrystalBallCauchyTransform : public CauchyTransform
        {
            private:
                double                  _extent;
                GaussianCauchyTransform _gaussian;
                double                  _mu    = 0.0;
                double                  _sigma = 1.0;

                // the tail parameters alpha_L, n_L, alpha_R, n_R of the table, and the derived constants
                std::array<double, 4> _tails;
                double                _norm, _E_L, _E_R;

                // pi H[delta] and pi H[delta'] on s = sinh(u), for u uniform in [-asinh(s_max), +asinh(s_max)]
                double              _s_max = 0.0;
                double              _u_max, _du;
                std::vector<double> _hilbert, _hilbert_prime;

                // the tails as exp(-alpha^2 / 2) (1 + alpha (|t| - alpha) / n)^-n, since (n / alpha)^n overflows for large n
                double
                shape(const double & t) const
                {
                    const auto [alpha_L, n_L, alpha_R, n_R] = _tails;
                    if (t < -alpha_L)
                    {
                        return _E_L * std::pow(1.0 + alpha_L * (-t - alpha_L) / n_L, -n_L);
                    }
                    if (t > alpha_R)
                    {
                        return _E_R * std::pow(1.0 + alpha_R * (t - alpha_R) / n_R, -n_R);
                    }

                    return std::exp(-0.5 * t * t);
                }

                double
                shape_prime(const double & t) const
                {
                    const auto [alpha_L, n_L, alpha_R, n_R] = _tails;
                    if (t < -alpha_L)
                    {
                        return _E_L * alpha_L * std::pow(1.0 + alpha_L * (-t - alpha_L) / n_L, -n_L - 1.0);
                    }
                    if (t > alpha_R)
                    {
                        return -_E_R * alpha_R * std::pow(1.0 + alpha_R * (t - alpha_R) / n_R, -n_R - 1.0);
                    }

                    return -t * std::exp(-0.5 * t * t);
                }

                double
                delta(const double & t) const
                {
                    return shape(t) / _norm - std::exp(-0.5 * t * t) / std::sqrt(2.0 * M_PI);
                }

                double
                delta_prime(const double & t) const
                {
                    return shape_prime(t) / _norm + t * std::exp(-0.5 * t * t) / std::sqrt(2.0 * M_PI);
                }

                // p.v. int dt f(t) / (s - t), with f(s) subtracted and integrated analytically within [-L, L]; f decays
                // like a power, so the tails beyond L, which matter for n close to 1, are integrated in u = L / |t|
                double
                principal_value(const std::function<double(const double &)> & f, const double & s) const
                {
                    const double L   = 10.0 * std::max(_s_max, std::abs(s)) + 100.0;
                    const double f_s = f(s);

                    std::vector<double> points{ -L, L, -_tails[0], _tails[2], s - 1.0, s, s + 1.0 };
                    for (double p = 10.0; p < L; p *= 10.0)
                    {
                        points.push_back(-p);
                        points.push_back(+p);
                    }
                    std::sort(points.begin(), points.end());

                    // merge breakpoints that are too close to resolve, keeping s: next to s, the nodes of an interval of a
                    // few ulps would round to s, where the integrand below is 0 / 0
                    std::vector<double> merged;
                    for (const double p : points)
                    {
                        if (merged.empty() || (p - merged.back() > 1.0e-9 * (1.0 + std::abs(p))))
                        {
                            merged.push_back(p);
                        }
                        else if (p == s)
                        {
                            merged.back() = s;
                        }
                    }
                    points = std::move(merged);

                    const auto integrand = [&](const double & t) -> double { return (f(t) - f_s) / (s - t); };
                    const auto config    = GSL::QAGS::Config().epsrel(1.0e-10).epsabs(1.0e-13);

                    // a failed integration reports which part of the principal value it belongs to
                    const auto checked_integrate = [&](const std::function<double(const double &)> & g, const double & a, const double & b, const char * part) -> double
                    {
                        try
                        {
                            return integrate<GSL::QAGS>(g, a, b, config);
                        }
                        catch (const Exception & e)
                        {
                            throw InternalError(std::format("CauchyTransform: the principal value at s = {} failed for the {} on [{}, {}], with the tails ({}, {}, {}, {}): {}",
                                                            s,
                                                            part,
                                                            a,
                                                            b,
                                                            _tails[0],
                                                            _tails[1],
                                                            _tails[2],
                                                            _tails[3],
                                                            e.what()));
                        }
                    };

                    double result = f_s * std::log((L + s) / (L - s));
                    for (std::size_t i = 1; i < points.size(); ++i)
                    {
                        result += checked_integrate(integrand, points[i - 1], points[i], "range within L");
                    }

                    // t = +-L u^-a, with dt = a |t| du / u; for a tail ~ |t|^-n, a = max(1, 1 / (n - 1)) makes the integrand
                    // vanish at least linearly at u = 0, also for n close to 1
                    for (const auto & [sign, n] : {
                             std::pair{ -1.0, _tails[1] },
                             std::pair{ +1.0, _tails[3] }
                    })
                    {
                        const double a    = std::max(1.0, 1.0 / (n - 1.0));
                        const auto   tail = [&](const double & u) -> double
                        {
                            const double t = sign * L * std::pow(u, -a);
                            if (! std::isfinite(t))
                            {
                                return 0.0;
                            }

                            const double f_t = f(t);
                            return (0.0 == f_t) ? 0.0 : a * f_t * (std::abs(t) / (s - t)) / u;
                        };
                        result += checked_integrate(tail, 0.0, 1.0, sign < 0.0 ? "left tail in u" : "right tail in u");
                    }

                    return result;
                }

                void
                tabulate(const double & s_max)
                {
                    _s_max = s_max;
                    _u_max = std::asinh(s_max);

                    const std::size_t points = 2 * static_cast<std::size_t>(std::ceil(_u_max / 0.005)) + 1;
                    _du                      = 2.0 * _u_max / static_cast<double>(points - 1);
                    _hilbert.resize(points);
                    _hilbert_prime.resize(points);

                    const auto f       = [this](const double & t) { return delta(t); };
                    const auto f_prime = [this](const double & t) { return delta_prime(t); };
                    for (std::size_t i = 0; i < points; ++i)
                    {
                        const double s    = std::sinh(-_u_max + static_cast<double>(i) * _du);
                        _hilbert[i]       = principal_value(f, s);
                        _hilbert_prime[i] = principal_value(f_prime, s);
                    }
                }

                // four-point Lagrange interpolation in u = asinh(s)
                double
                interpolate(const std::vector<double> & table, const double & s) const
                {
                    const double      x = (std::asinh(s) + _u_max) / _du;
                    const std::size_t i = std::clamp<std::size_t>(static_cast<std::size_t>(x), 1, table.size() - 3);
                    const double      w = x - static_cast<double>(i);

                    return -w * (w - 1.0) * (w - 2.0) / 6.0 * table[i - 1] + (w + 1.0) * (w - 1.0) * (w - 2.0) / 2.0 * table[i] - (w + 1.0) * w * (w - 2.0) / 2.0 * table[i + 1]
                           + (w + 1.0) * w * (w - 1.0) / 6.0 * table[i + 2];
                }

            public:
                CrystalBallCauchyTransform(const double & extent) :
                    _extent(extent)
                {
                    _tails.fill(std::numeric_limits<double>::quiet_NaN());
                }

                virtual void
                prepare(std::span<const double> parameters) override
                {
                    _mu    = parameters[0];
                    _sigma = parameters[1];
                    _gaussian.prepare(parameters.subspan(0, 2));

                    const std::array<double, 4> tails{ parameters[2], parameters[3], parameters[4], parameters[5] };
                    const double                s_max = _extent / _sigma;
                    if ((tails == _tails) && (s_max <= _s_max))
                    {
                        return;
                    }

                    const auto [alpha_L, n_L, alpha_R, n_R] = tails;

                    _tails = tails;
                    _E_L   = std::exp(-0.5 * alpha_L * alpha_L);
                    _E_R   = std::exp(-0.5 * alpha_R * alpha_R);
                    _norm  = std::sqrt(M_PI / 2.0) * (std::erf(alpha_L / std::sqrt(2.0)) + std::erf(alpha_R / std::sqrt(2.0))) + n_L * _E_L / (alpha_L * (n_L - 1.0))
                            + n_R * _E_R / (alpha_R * (n_R - 1.0));

                    tabulate(s_max);
                }

                virtual complex<double>
                evaluate(const complex<double> & zeta) const override
                {
                    const double s = (zeta.real() - _mu) / _sigma;
                    const double g = zeta.imag() / _sigma;

                    double hilbert, hilbert_prime;
                    if (std::abs(s) <= _s_max)
                    {
                        hilbert       = interpolate(_hilbert, s);
                        hilbert_prime = interpolate(_hilbert_prime, s);
                    }
                    else
                    {
                        hilbert       = principal_value([this](const double & t) { return delta(t); }, s);
                        hilbert_prime = principal_value([this](const double & t) { return delta_prime(t); }, s);
                    }

                    const complex<double> S_delta = (hilbert - 1.0i * M_PI * delta(s) + 1.0i * g * (hilbert_prime - 1.0i * M_PI * delta_prime(s))) / _sigma;

                    return _gaussian.evaluate(zeta) + S_delta;
                }
        };
    } // namespace

    std::unique_ptr<CauchyTransform>
    CauchyTransform::make(const std::string & kernel, const double & extent)
    {
        if ("Kernel::Gaussian" == kernel)
        {
            return std::make_unique<GaussianCauchyTransform>();
        }

        if ("Kernel::CrystalBall" == kernel)
        {
            return std::make_unique<CrystalBallCauchyTransform>(extent);
        }

        throw InternalError("CauchyTransform: no Cauchy transform is known for the kernel '" + kernel + "'");
    }
} // namespace eos
