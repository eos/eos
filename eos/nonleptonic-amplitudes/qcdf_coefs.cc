/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Marta Burgos
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

#include <eos/models/model.hh>
#include <eos/models/wilson-coefficients.hh>
#include <eos/nonleptonic-amplitudes/qcdf_coefs.hh>
#include <eos/maths/matrix.hh>
#include <eos/maths/power-of.hh>
#include <gsl/gsl_sf_dilog.h>
#include <eos/maths/integrate.hh>
#include <eos/utils/options-impl.hh>

#include <array>
#include <cmath>
#include <vector>
#include <map>


namespace eos
{
    using std::sqrt;

    NonleptonicAmplitudes<PToPP> *
    QCDFCoefficients<PToPP>::make(const Parameters & p, const Options & o)
    {
        return new QCDFCoefficients<PToPP>(p, o);
    }

        QCDFCoefficients<PToPP>::QCDFCoefficients(const Parameters & p, const Options & o):
                model(Model::make(o.get("model"_ok, "SM"), p, o)),
                opt_q(o, options, "q"_ok),
                opt_p1(o, options, "P1"_ok),
                opt_p2(o, options, "P2"_ok),
                mB(p["mass::B_" + opt_q.str()], *this),
                mB_q_0(p["mass::B_" + opt_q.str() + ",0@BSZ2015"], *this),
                mP1(p["mass::" + opt_p1.str()], *this),
                mP2(p["mass::" + opt_p2.str()], *this),
                FP1(p["B_" + opt_q.str() + "->" + opt_p1.str() + "::f_+(0)"], *this),
                FP2(p["B_" + opt_q.str() + "->" + opt_p2.str() + "::f_+(0)"], *this),
                fB(p["decay-constant::B_" + opt_q.str()], *this),
                fP1(p["decay-constant::" + opt_p1.str()], *this),
                mb(p["mass::b(MSbar)"], *this),
                mus(p["mass::b(MSbar)"], *this),
                //mu(p[stringify(opt_q.value() == QuarkFlavor::down ? "s" : "d") + "bcu::mu"], u), //what is this? do I need u for this?
                fP2(p["decay-constant::" + opt_p2.str()], *this)
        {
                if (opt_p1.str() == "pi^0" || opt_p1.str() == "pi^+" || opt_p1.str() == "pi^-")
                {
                    lcdasP1 = PseudoscalarLCDAs::make("pi", p, o);
                }
                else if (opt_p1.str() == "K_d" || opt_p1.str() == "K_u")
                {
                    lcdasP1 = PseudoscalarLCDAs::make("K", p, o);
                }
                else if (opt_p1.str() == "Kbar_d" ||  opt_p1.str() == "Kbar_u")
                {
                    lcdasP1 = PseudoscalarLCDAs::make("Kbar", p, o);
                }


                if (opt_p2.str() == "pi^0" || opt_p2.str() == "pi^+" || opt_p2.str() == "pi^-")
                {
                    lcdasP2 = PseudoscalarLCDAs::make("pi", p, o);
                }
                else if (opt_p2.str() == "K_d" || opt_p2.str() == "K_u")
                {
                    lcdasP2 = PseudoscalarLCDAs::make("K", p, o);
                }
                else if (opt_p2.str() == "Kbar_d" ||  opt_p2.str() == "Kbar_u")
                {
                    lcdasP2 = PseudoscalarLCDAs::make("Kbar", p, o);
                }


        }

        const std::vector<OptionSpecification>  QCDFCoefficients<PToPP>::options{
            Model::option_specification(),
            { "q"_ok, { "u", "d", "s" } },
            { "P1"_ok, { "pi^0", "pi^+", "pi^-", "K_d", "Kbar_d", "K_s", "K_u", "Kbar_u", "eta", "eta_prime", "eta_q", "eta_s" } },
            { "P2"_ok, { "pi^0", "pi^+", "pi^-", "K_d", "Kbar_d", "K_s", "K_u", "Kbar_u", "eta", "eta_prime", "eta_q", "eta_s" } },
        };



        complex<double>
        QCDFCoefficients<PToPP>::gvertex(const double & x) const
            {
                return 3.0 * ((1.0 - 2.0 * x) / (1.0 - x) * log(x) - complex<double>(0.0, 1.0) * M_PI) + 2.0 * gsl_sf_dilog(x) - power_of<2>(log(x)) + 2.0 * log(x) / (1.0 - x)
                - log(x) * (3.0 + 2.0 * complex<double>(0.0, 1.0) * M_PI) - (2.0 * gsl_sf_dilog(1.0 - x) - power_of<2>(log(1.0 - x)) + 2.0 * log(1.0 - x) / x
                - log(1.0 - x) * (3.0 + 2.0 * complex<double>(0.0, 1.0) * M_PI));
            }

        std::array<complex<double>, 11>
        QCDFCoefficients<PToPP>::Vertex(const double & x) const
        {
            double eps = 1e-5;
            std::array<complex<double>, 11> Vi;

            for (unsigned i = 0; i < 11; i++)
                {
                        if (i == 1 || i == 2 || i == 3 || i == 4 || i == 9 || i == 10)
                        {
                            Vi[i] = integrate1D(
                                std::function<complex<double>(const double&)>(
                                    [=, this](const double &x)
                                    {
                                        return lcdasP2->phi(x, mu()) * (12.0 * std::log(mb / mus) - 18.0 + this->gvertex(x));
                                    }
                                ),64, 0.0, 1.0);
                        }

                        else if (i == 5 || i == 7)
                        {

                            Vi[i] = integrate1D(
                                std::function<complex<double>(const double&)>(
                                    [=, this](const double &x)
                                    {
                                        return lcdasP2->phi(x, mu()) * (-12.0 * std::log(mb / mus) + 6.0 - this->gvertex(1.0 - x));
                                    }
                                ),64, 0.0, 1.0);

                        }

                        else if (i == 6 || i == 8)
                        {
                            Vi[i] = -6.0;
                        }


                }

            return Vi;
        }

        std::array<complex<double>, 11>
        QCDFCoefficients<PToPP>::HardSpec(const double & x, const double & y) const{

            double eps = 1e-5;

            std::array<complex<double>, 11> Haux;

            double rchi = 1.0 / (1.0 - power_of<2>(mP2() / mB_q_0()));
            double XH = 1.0 / (1.0 - power_of<2>(mP2() / mB_q_0())) * (1.0 - power_of<2>(mP1() / mB_q_0()));
            double prefac = power_of<2>(mB()) * FP1()  / ((1 - power_of<2>(mP2() / mB_q_0())) * fB() * fP1());

            for (unsigned i = 0; i < 11; i++)
                {
                    std::function<double(const double &)> integrand;

                    if (i == 1 || i == 2 || i == 3 || i == 4 || i == 9 || i == 10)
                    {

                            std::function<complex<double>(const double&)> outer_integral = [=, this](const double &x)
                            {
                                auto inner_integral = std::function<complex<double>(const double&)>(
                                    [=, this, x](const double &y)
                                    {
                                        return lcdasP1->phi(y,mu()) * lcdasP2->phi(x,mu())/ ((1 - x) * (1 - y)) + rchi * lcdasP2->phi(x,mu()) * XH / x;
                                    }
                                );
                                return integrate1D(inner_integral, 32, 0.0, 1.0);
                            };

                            Haux[i] = prefac * integrate1D(outer_integral, 32, 0.0, 1.0);
                    }

                    else if (i == 5 || i == 7)
                    {

                        std::function<complex<double>(const double&)> outer_integral = [=, this](const double &x)
                            {
                                auto inner_integral = std::function<complex<double>(const double&)>(
                                    [=, this, x](const double &y)
                                    {
                                        return lcdasP1->phi(y,mu()) * lcdasP2->phi(x,mu())/ (x * (1 - y)) + rchi * lcdasP2->phi(x,mu()) * XH / x;
                                    }
                                );
                                return integrate1D(inner_integral, 32, 0.0, 1.0);
                            };

                            Haux[i] = -1.0 * prefac * integrate1D(outer_integral, 32, 0.0, 1.0);
                    }

                    else if (i == 6 || i == 8)
                    {

                        Haux[i] = 0.0;
                    }


                };

            return Haux;

        }

        complex<double>
        QCDFCoefficients<PToPP>::Gsx(const double & s, const double & x) const
        {
            return 2.0 * (12.0 * s + 5.0 * x - 3.0 * x * log(s))/(9.0 * x) - 4.0 * sqrt(4.0 * s - x) * (2.0 * s + x) / (3.0 * sqrt(power_of<3>(x))) * atan(sqrt(x/(4.0 * s - x)));
        }

        complex<double>
        QCDFCoefficients<PToPP>::GM2(const double & s) const
        {
            return integrate1D(
                                std::function<complex<double>(const double&)>(
                                    [=, this](const double &x)
                                    {
                                        return lcdasP2->phi(x, mu()) * this->Gsx(s, x);
                                    }
                                ),64, 0.0, 1.0);
        }

        complex<double>
        QCDFCoefficients<PToPP>::GM2hat(const double & s) const
        {
            double eps = 1e-5;
            return integrate1D(
                                std::function<complex<double>(const double&)>(
                                    [=, this](const double &x)
                                    {
                                        return this->Gsx(s - eps, 1.0 - x);
                                    }
                                ),64, 0.0, 1.0);
        }

        std::array<complex<double>, 11>
        QCDFCoefficients<PToPP>::Penguin(const double & x) const
        {
            std::array<complex<double>, 11> peng;
            complex<double> CF = 4.0/3.0;
            double alphas = 0.25; //change
            double Nc = 3.0;
            double nf = 5.0;
            double sp = 0.0; //change
            double sc = 0.0; //change
            double alpha = 1.0 / 137.0; //change

            std::array<complex<double>, 11> wc;

            for( unsigned i = 0; i < 11; i++)
            {
                if (i == 4)
                {
                    peng[i] = CF * alphas / (4.0 * M_PI * Nc) * (wc2.c1() * (4.0 / 3.0 * log(mb / mus) + 2.0 / 3.0 - this->GM2(sp)) +
                    wc[3] * (8.0 / 3.0 * log(mb / mus) + 4.0 / 3.0 - this->GM2(0.0) - this->GM2(1.0)) +
                    (wc[4] + wc[6]) * (4.0 * nf / 3.0 * log(mb / mus) - (nf - 2.0) * this->GM2(0.0) - this->GM2(sc) - this->GM2(1.0)) -
                    wc[8] * integrate1D(std::function<complex<double>(const double&)>([=, this](const double &x)
                                    {
                                        return lcdasP2->phi(x, mu()) / (1.0 - x);
                                    }
                                ),64, 0.0, 1.0));
                }

                else if (i == 6)
                {
                    peng[i] = CF * alphas / (4.0 * M_PI * Nc) * (wc[1] * (4.0 / 3.0 * log(mb / mus) + 2.0 / 3.0 - this->GM2hat(sp)) +
                    wc[3] * (8.0 / 3.0 * log(mb / mus) + 4.0 / 3.0 - this->GM2hat(0.0) - this->GM2hat(1.0)) +
                    (wc[4] + wc[6]) * (4.0 * nf / 3.0 * log(mb / mus) - (nf - 2.0) * this->GM2hat(0.0) - this->GM2hat(sc) - this->GM2hat(1.0)) -
                    2.0 * wc[8]);
                }

                else if (i == 8)
                {
                    peng[i] = alpha / (9.0 * M_PI * Nc) * ((wc[1] + Nc * wc[2]) * (4.0 / 3.0 * log(mb/mus) + 2.0 / 3.0 - this->GM2hat(sp)) - 3.0 * wc[7]);
                }

                else if (i == 10)
                {
                    peng[i] = alpha / (9.0 * M_PI * Nc) * ((wc[1] + Nc * wc[2]) * (4.0 / 3.0 * log(mb/mus) + 2.0 / 3.0 - this->GM2(sp)) -
                    3.0 * wc[7] * integrate1D(std::function<complex<double>(const double&)>([=, this](const double &x)
                                    {
                                        return lcdasP2->phi(x, mu()) / (1.0 - x);
                                    }
                                ),64, 0.0, 1.0));
                }

                else
                {
                    peng[i] = 0.0;
                }

            }
            return peng;
        }


}
