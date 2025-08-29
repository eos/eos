/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017 Danny van Dyk
 * Copyright (c) 2017 Elena Graverini
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

#include <eos/form-factors/hqet-b-to-c.hh>
#include <eos/models/model.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/maths/polylog.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/stringify.hh>

namespace eos
{
    template <>
    struct Implementation<HQETBToC>
    {
        std::shared_ptr<Model> model;

        UsedParameter m_b_msbar;
        UsedParameter m_c_msbar;


        Implementation(const Parameters & p, const Options &, ParameterUser & u) :
            model(Model::make("SM", p, Options{ })),
            m_b_msbar(p["mass::b(MSbar)"], u),
            m_c_msbar(p["mass::c"], u)
        {
            u.uses(*model);
        }

        ~Implementation() = default;

        /* auxilliary functions from [N:1993] */

        // r(omega) is defined in [N:1993] eq. (3.104), p. 63
        inline static double r(const double & omega)
        {
            if (omega < 1.0)
                throw InternalError("HQETBToC::r omega = '" + stringify(omega) + "' outside its domain of validity");

            // for small omega - 1, taylor expand r up to second order
            if ((omega - 1.0) < 1.0e-5)
            {
                const double c0 = +1.0;
                const double c1 = -1.0 / 3.0;
                const double c2 = +4.0 / 30.0;

                return c0 + (omega - 1.0) * c1 + power_of<2>(omega - 1.0) * c2;
            }

            return log(omega + sqrt(omega * omega - 1.0)) / sqrt(omega * omega - 1.0);
        }

        // f(omega) is defined in [N:1993] eq. (3.117), p. 65
        inline static double f(const double & omega)
        {
            if (omega < 1.0)
                throw InternalError("HQETBToC::f omega = '" + stringify(omega) + "' outside its domain of validity");

            // for small omega - 1, taylor expand f up to second order
            if ((omega - 1.0) < 1.0e-5)
            {
                const double c0 = -3.0;
                const double c1 = -10.0 / 9.0;
                const double c2 = +2.0 / 150.0;

                return c0 + (omega - 1.0) * c1 + power_of<2>(omega - 1.0) * c2;
            }

            const double omega_m = omega - sqrt(omega * omega - 1.0);
            const double L2      = real(dilog(complex<double>(1.0 - omega_m * omega_m, 0.0)));
            const double r_omega = r(omega);

            return omega * r_omega - 2.0 - omega / sqrt(omega * omega - 1.0) * (L2 + (omega * omega - 1.0) * r_omega * r_omega);
        }

        // g(z, omega) is defined in [N:1993] eq. (3.129), p. 70
        inline static double g(const double & z, const double & omega)
        {
            if (omega < 1.0)
                throw InternalError("HQETBToC::g omega = '" + stringify(omega) + "' outside its domain of validity");

            // for small omega - 1, taylor expand g up to second order
            if ((omega - 1.0) < 1.0e-5)
            {
                const double c0 = z * log(z) / (z - 1.0);
                const double c1 = z * (2.0 - 2.0 * z + (1.0 + z) * log(z)) / power_of<3>(z - 1.0);
                const double c2 = z * (1.0 + 9.0 * z - 9.0 * z * z - z * z + z + 6.0 * z * (1.0 + z) * log(z)) / (3.0 * power_of<5>(z - 1.0));

                return c0 + (omega - 1.0) * c1 + power_of<2>(omega - 1.0) * c2;
            }

            const double omega_m = omega - sqrt(omega * omega - 1.0);
            const double omega_p = omega + sqrt(omega * omega - 1.0);
            const double L2_p    = real(dilog(complex<double>(1.0 - z * omega_p, 0.0)));
            const double L2_m    = real(dilog(complex<double>(1.0 - z * omega_m, 0.0)));
            const double r_omega = r(omega);

            return omega / sqrt(omega * omega - 1.0) * (L2_m - L2_p)
                - z / (1.0 - 2.0 * omega * z + z * z) * ((omega * omega - 1.0) * r_omega + (omega - z) * log(z));
        }

        /* anomalous dimensions and auxilliaries for next-to-leading log terms */

        // a_hh(omega) is defined in [N:1993] eq. (3.119), p. 66
        inline static double a_hh(const double & omega)
        {
            return 8.0 / 27.0 * (omega * r(omega) - 1.0);
        }

        // Z_hh(omega) is defined in [N:1993] eq. (3.119), p. 66
        // Note that we use only the Taylor expansion in (omega - 1) up to second order.
        inline static double Z_hh(const double & omega)
        {
            return 8.0 / (81.0) * (94.0 / 9.0 - M_PI * M_PI) * (omega - 1.0)
                - 4.0 / 135.0 * (92.0 / 9.0 - M_PI * M_PI) * power_of<2>(omega - 1.0);
        }

        // S_{1,2,3}^{(5)} for the two currents are defined in [N:1993] eq. (3.145)
        inline static double S_1(const double & x, const double & omega)
        {
            return omega * (17.0 / 27.0 - 5.0 / 9.0 * pow(x, -6.0 / 25.0) - 2.0 / 27.0 * pow(x, -9.0 / 25.0) + 8.0 / 25.0 * log(x))
                + (1.0 / 6.0 - 5.0 / 9.0 * pow(x, -6.0 / 25.0) + 4.0 / 9.0 * pow(x, -9.0 / 25.0) - 1.0 / 18.0 * pow(x, -12.0 / 25.0));
        }
        inline static double S_1_5(const double & x, const double & omega)
        {
            return omega * (17.0 / 27.0 - 5.0 / 9.0 * pow(x, -6.0 / 25.0) - 2.0 / 27.0 * pow(x, -9.0 / 25.0) + 8.0 / 25.0 * log(x))
                - (1.0 / 6.0 - 5.0 / 9.0 * pow(x, -6.0 / 25.0) + 4.0 / 9.0 * pow(x, -9.0 / 25.0) - 1.0 / 18.0 * pow(x, -12.0 / 25.0));
        }
        inline static double S_2(const double & x, const double & omega)
        {
            return -omega * (14.0 / 27.0 + 10.0 / 9.0 * pow(x, -6.0 / 25.0) - 44.0 / 27.0 * pow(x, -9.0 / 25.0))
                + (2.0 / 3.0 + 5.0 / 9.0 * pow(x, -6.0 / 25.0) + 2.0 / 9.0 * pow(x, -9.0 / 25.0) - 13.0 / 9.0 * pow(x, -12.0 / 25.0));
        }
        inline static double S_2_5(const double & x, const double & omega)
        {
            return +omega * (14.0 / 27.0 + 10.0 / 9.0 * pow(x, -6.0 / 25.0) - 44.0 / 27.0 * pow(x, -9.0 / 25.0))
                + (2.0 / 3.0 + 5.0 / 9.0 * pow(x, -6.0 / 25.0) + 2.0 / 9.0 * pow(x, -9.0 / 25.0) - 13.0 / 9.0 * pow(x, -12.0 / 25.0));
        }
        inline static double S_3(const double & x)
        {
            return 1.0 - 5.0 / 3.0 * pow(x, -6.0 / 25.0) + 2.0 / 3.0 * pow(x, -9.0 / 25.0);
        }
        inline static double S_3_5(const double & x) { return S_3(x); }

        // we use the form factors at a fixed scale mu = m_c
        inline double mu() const
        {
            return m_c_msbar();
        }

        // we use a fixed matching scale mu_match = sqrt(m_b * m_c)
        inline double mu_match() const
        {
            return sqrt(m_b_msbar() * m_c_msbar());
        }

        // universal mu-dependence of the Wilson coefficients
        inline double K_hh(const double & omega) const
        {
            const double alpha_s_mu = model->alpha_s(mu());

            return pow(alpha_s_mu, -a_hh(omega)) * (1.0 - alpha_s_mu / M_PI * Z_hh(omega));
        }

        // universal prefactor A
        inline double A(const double & omega) const
        {
            const double alpha_s_mb = model->alpha_s(m_b_msbar());
            const double alpha_s_mc = model->alpha_s(m_c_msbar());

            return pow(alpha_s_mc / alpha_s_mb, 6.0 / 25.0) * pow(alpha_s_mc, a_hh(omega));
        }

        // h_2(z, omega) is defined in [N:1993] eq. (3.129), p. 70
        inline static double h_2(const double & z, const double & omega)
        {
            const double denom = 1.0 - 2.0 * omega * z + z * z;

            return z / power_of<2>(denom) * (
                    2.0 * (omega - 1.0) * z * (1.0 + z) * log(z)
                    - (
                        (omega + 1.0) - 2.0 * omega * (2.0 * omega + 1.0) * z
                        + (5.0 * omega + 2.0 * omega * omega - 1.0) * z * z - 2.0 * z * z * z
                    ) * r(omega)
                ) - z / denom * (log(z) - 1.0 + z);
        }
        // h_2_5(z, omega) is defined in [N:1993] eq. (3.129), p. 70
        inline static double h_2_5(const double & z, const double & omega)
        {
            const double denom = 1.0 - 2.0 * omega * z + z * z;

            return z / power_of<2>(denom) * (
                    2.0 * (omega + 1.0) * z * (1.0 - z) * log(z)
                    - (
                        (omega - 1.0) - 2.0 * omega * (2.0 * omega - 1.0) * z
                        + (5.0 * omega - 2.0 * omega * omega + 1.0) * z * z - 2.0 * z * z * z
                    ) * r(omega)
                ) - z / denom * (log(z) - 1.0 - z);
        }

        // hatted Wilson coefficients without the universal mu-dependence,
        // cf. [N:1993] eq. (3.142), p. 74
        double Chat_1_v(const double & omega) const
        {
            static const double Z_4 = -9403.0 / 7500.0 - 7.0 * M_PI * M_PI / 225.0;

            const double alpha_s_mb = model->alpha_s(m_b_msbar());
            const double alpha_s_mc = model->alpha_s(m_c_msbar());
            const double alpha_s_m  = model->alpha_s(mu_match());
            const double x          = alpha_s_mc / alpha_s_mb;
            const double z          = m_c_msbar() / m_b_msbar();
            const double G          = g(z, omega) + 3.0 * omega * z * log(z);

            return A(omega) * (
                1.0
                + (alpha_s_mb - alpha_s_mc) / M_PI * Z_4 + alpha_s_mc / M_PI * (Z_hh(omega) + 2.0 / 3.0 * (f(omega) + r(omega)))
                + z * S_1(x, omega) + 2.0 * alpha_s_m / (3.0 * M_PI) * G
            );
        }
        double Chat_2_v(const double & omega) const
        {
            const double alpha_s_mb = model->alpha_s(m_b_msbar());
            const double alpha_s_mc = model->alpha_s(m_c_msbar());
            const double alpha_s_m  = model->alpha_s(mu_match());
            const double x          = alpha_s_mc / alpha_s_mb;
            const double z          = m_c_msbar() / m_b_msbar();
            const double h1         = h_2(1.0 / z, omega) - 2.0 * r(omega) + 1.0;
            const double H1         = h1 - (3.0 - 2.0 * omega) * z * log(z);

            return A(omega) * (
                + 2.0 * alpha_s_mb / (3.0 * M_PI)
                - 4.0 * alpha_s_mc / (3.0 * M_PI) * r(omega)
                + z * S_2(x, omega)
                - 2.0 * alpha_s_m / (3.0 * M_PI) * H1
            );
        }
        double Chat_3_v(const double & omega) const
        {
            const double alpha_s_mb = model->alpha_s(m_b_msbar());
            const double alpha_s_mc = model->alpha_s(m_c_msbar());
            const double alpha_s_m  = model->alpha_s(mu_match());
            const double x          = alpha_s_mc / alpha_s_mb;
            const double z          = m_c_msbar() / m_b_msbar();
            const double h2         = h_2(z, omega);
            const double H2         = h2 + z * log(z);
            const double S_3        = this->S_3(x);

            return -A(omega) * (
                + z * S_3
                + 2.0 * alpha_s_m / (3.0 * M_PI) * H2
            );
        }
        double Chat_1_a(const double & omega) const
        {
            static const double Z_4 = -9403.0 / 7500.0 - 7.0 * M_PI * M_PI / 225.0;

            const double alpha_s_mb = model->alpha_s(m_b_msbar());
            const double alpha_s_mc = model->alpha_s(m_c_msbar());
            const double alpha_s_m  = model->alpha_s(mu_match());
            const double x          = alpha_s_mc / alpha_s_mb;
            const double z          = m_c_msbar() / m_b_msbar();
            const double G          = g(z, omega) + 3.0 * omega * z * log(z);

            return A(omega) * (
                1.0
                + (alpha_s_mb - alpha_s_mc) / M_PI * Z_4 + alpha_s_mc / M_PI * (Z_hh(omega) + 2.0 / 3.0 * (f(omega) - r(omega)))
                + z * S_1_5(x, omega) + 2.0 * alpha_s_m / (3.0 * M_PI) * G
            );
        }
        double Chat_2_a(const double & omega) const
        {
            const double alpha_s_mb = model->alpha_s(m_b_msbar());
            const double alpha_s_mc = model->alpha_s(m_c_msbar());
            const double alpha_s_m  = model->alpha_s(mu_match());
            const double x          = alpha_s_mc / alpha_s_mb;
            const double z          = m_c_msbar() / m_b_msbar();
            const double h1         = h_2_5(1.0 / z, omega) - 2.0 * r(omega) - 1.0;
            const double H1         = h1 - (3.0 + 2.0 * omega) * z * log(z);

            return A(omega) * (
                - 2.0 * alpha_s_mb / (3.0 * M_PI)
                - 4.0 * alpha_s_mc / (3.0 * M_PI) * r(omega)
                + z * S_2_5(x, omega)
                - 2.0 * alpha_s_m / (3.0 * M_PI) * H1
            );
        }
        double Chat_3_a(const double & omega) const
        {
            const double alpha_s_mb = model->alpha_s(m_b_msbar());
            const double alpha_s_mc = model->alpha_s(m_c_msbar());
            const double alpha_s_m  = model->alpha_s(mu_match());
            const double x          = alpha_s_mc / alpha_s_mb;
            const double z          = m_c_msbar() / m_b_msbar();
            const double h2         = h_2_5(z, omega);
            const double H2         = h2 + z * log(z);

            return +A(omega) * (
                + z * S_3_5(x)
                + 2.0 * alpha_s_m / (3.0 * M_PI) * H2
            );
        }

        Diagnostics diagnostics() const
        {
            Diagnostics results;

            // Inputs
            {
                const double alpha_s_mb = model->alpha_s(m_b_msbar());
                const double alpha_s_mc = model->alpha_s(m_c_msbar());
                const double alpha_s_m  = model->alpha_s(mu_match());
                const double x          = alpha_s_mc / alpha_s_mb;
                const double z          = m_c_msbar() / m_b_msbar();

                results.add(Diagnostics::Entry{ alpha_s_mb, "alpha_s(m_b^MSbar)" });
                results.add(Diagnostics::Entry{ alpha_s_mc, "alpha_s(m_c^MSbar)" });
                results.add(Diagnostics::Entry{ alpha_s_m,  "alpha_s(mu_match)"  });

                results.add(Diagnostics::Entry{ x, "x = alpha_s_mc / alpha_s_mb" });
                results.add(Diagnostics::Entry{ z, "z = m_c / m_b"               });
            }

            // Universal mu dependence
            {
                results.add(Diagnostics::Entry{ K_hh(1.0), "K_hh(1.0)" });
                results.add(Diagnostics::Entry{ K_hh(1.1), "K_hh(1.1)" });
                results.add(Diagnostics::Entry{ K_hh(1.2), "K_hh(1.2)" });
            }

            // Chat_1_v
            {
                static const double z = 0.305024;
                static const double Z_4 = -9403.0 / 7500.0 - 7.0 * M_PI * M_PI / 225.0;
                results.add(Diagnostics::Entry{ Z_4, "Z_4" });

                results.add(Diagnostics::Entry{ A(1.0), "A(1.0)" });
                results.add(Diagnostics::Entry{ A(1.1), "A(1.1)" });
                results.add(Diagnostics::Entry{ A(1.2), "A(1.2)" });

                results.add(Diagnostics::Entry{ f(1.0), "f(1.0)" });
                results.add(Diagnostics::Entry{ f(1.1), "f(1.1)" });
                results.add(Diagnostics::Entry{ f(1.2), "f(1.2)" });

                results.add(Diagnostics::Entry{ g(z, 1.0), "g(z = 0.305024, omega = 1.0)" });
                results.add(Diagnostics::Entry{ g(z, 1.1), "g(z = 0.305024, omega = 1.1)" });
                results.add(Diagnostics::Entry{ g(z, 1.2), "g(z = 0.305024, omega = 1.2)" });

                results.add(Diagnostics::Entry{ g(z, 1.0) + 3.0 * 1.0 * z * log(z), "G(z = 0.305024, omega = 1.0)" });
                results.add(Diagnostics::Entry{ g(z, 1.1) + 3.0 * 1.1 * z * log(z), "G(z = 0.305024, omega = 1.1)" });
                results.add(Diagnostics::Entry{ g(z, 1.2) + 3.0 * 1.2 * z * log(z), "G(z = 0.305024, omega = 1.2)" });

                results.add(Diagnostics::Entry{ r(1.0), "r(1.0)" });
                results.add(Diagnostics::Entry{ r(1.1), "r(1.1)" });
                results.add(Diagnostics::Entry{ r(1.2), "r(1.2)" });

                results.add(Diagnostics::Entry{ S_1(1.7589, 1.0), "S_1(x = 1.7589, omega = 1.0)" });
                results.add(Diagnostics::Entry{ S_1(1.7589, 1.1), "S_1(x = 1.7589, omega = 1.1)" });
                results.add(Diagnostics::Entry{ S_1(1.7589, 1.2), "S_1(x = 1.7589, omega = 1.2)" });

                results.add(Diagnostics::Entry{ Z_hh(1.0), "Z_hh(1.0)" });
                results.add(Diagnostics::Entry{ Z_hh(1.1), "Z_hh(1.1)" });
                results.add(Diagnostics::Entry{ Z_hh(1.2), "Z_hh(1.2)" });

            }

            // Hatted Wilson coefficients
            {
                results.add(Diagnostics::Entry{ Chat_1_v(1.0), "Chat_1_v(1.0)" });
                results.add(Diagnostics::Entry{ Chat_1_v(1.1), "Chat_1_v(1.1)" });
                results.add(Diagnostics::Entry{ Chat_1_v(1.2), "Chat_1_v(1.2)" });

                results.add(Diagnostics::Entry{ Chat_2_v(1.0), "Chat_2_v(1.0)" });
                results.add(Diagnostics::Entry{ Chat_2_v(1.1), "Chat_2_v(1.1)" });
                results.add(Diagnostics::Entry{ Chat_2_v(1.2), "Chat_2_v(1.2)" });

                results.add(Diagnostics::Entry{ Chat_3_v(1.0), "Chat_3_v(1.0)" });
                results.add(Diagnostics::Entry{ Chat_3_v(1.1), "Chat_3_v(1.1)" });
                results.add(Diagnostics::Entry{ Chat_3_v(1.2), "Chat_3_v(1.2)" });

                results.add(Diagnostics::Entry{ Chat_1_a(1.0), "Chat_1_a(1.0)" });
                results.add(Diagnostics::Entry{ Chat_1_a(1.1), "Chat_1_a(1.1)" });
                results.add(Diagnostics::Entry{ Chat_1_a(1.2), "Chat_1_a(1.2)" });

                results.add(Diagnostics::Entry{ Chat_2_a(1.0), "Chat_2_a(1.0)" });
                results.add(Diagnostics::Entry{ Chat_2_a(1.1), "Chat_2_a(1.1)" });
                results.add(Diagnostics::Entry{ Chat_2_a(1.2), "Chat_2_a(1.2)" });

                results.add(Diagnostics::Entry{ Chat_3_a(1.0), "Chat_3_a(1.0)" });
                results.add(Diagnostics::Entry{ Chat_3_a(1.1), "Chat_3_a(1.1)" });
                results.add(Diagnostics::Entry{ Chat_3_a(1.2), "Chat_3_a(1.2)" });
            }

            return results;
        }
    };

    HQETBToC::HQETBToC(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<HQETBToC>(new Implementation<HQETBToC>(p, o, *this))
    {
    }

    HQETBToC::~HQETBToC() = default;

    double
    HQETBToC::c_1_vector(const double & omega) const
    {
        const double Khh = _imp->K_hh(omega);

        return Khh * _imp->Chat_1_v(omega);
    }

    double
    HQETBToC::c_2_vector(const double & omega) const
    {
        const double Khh = _imp->K_hh(omega);

        return Khh * _imp->Chat_2_v(omega);
    }

    double
    HQETBToC::c_3_vector(const double & omega) const
    {
        const double Khh = _imp->K_hh(omega);

        return Khh * _imp->Chat_3_v(omega);
    }

    double
    HQETBToC::c_1_axialvector(const double & omega) const
    {
        const double Khh = _imp->K_hh(omega);

        return Khh * _imp->Chat_1_a(omega);
    }

    double
    HQETBToC::c_2_axialvector(const double & omega) const
    {
        const double Khh = _imp->K_hh(omega);

        return Khh * _imp->Chat_2_a(omega);
    }

    double
    HQETBToC::c_3_axialvector(const double & omega) const
    {
        const double Khh = _imp->K_hh(omega);

        return Khh * _imp->Chat_3_a(omega);
    }

    Diagnostics
    HQETBToC::diagnostics() const
    {
        return _imp->diagnostics();
    }

    const std::set<ReferenceName>
    HQETBToC::references
    {
    };
}
