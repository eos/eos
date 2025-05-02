/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017-2025 Danny van Dyk
 * Copyright (c) 2019 Nico Gubernari
 * Copyright (c) 2021 Méril Reboud
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

#include <eos/maths/complex.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
#include <eos/maths/complex.hh>
#include <eos/nonlocal-form-factors/nonlocal-formfactors.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/memoise.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/stringify.hh>


#include <map>
#include <numeric>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

namespace eos
{
    using std::abs;

    namespace nff
    {
        struct BToKstar
        {
            static constexpr const char * label = "B->K^*";
        };
        constexpr const char * BToKstar::label;

        struct BsToPhi
        {
            static constexpr const char * label = "B_s->phi";
        };
        constexpr const char * BsToPhi::label;
    }

    namespace nff_p_to_v
    {
        class Naive :
            public NonlocalFormFactor<PToV>
        {
            public:
                Naive(const Parameters &, const Options &)
                {
                }

                ~Naive() = default;

                virtual complex<double> H_perp(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> H_para(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> H_long(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> H_perp(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> H_para(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> H_long(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> Hhat_perp(const double &) const
                {
                    return 0.0;
                }
                virtual complex<double> Hhat_para(const double &) const
                {
                    return 0.0;
                }
                virtual complex<double> Hhat_long(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_perp(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_para(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_long(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_perp(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_para(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_long(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> F_ratio_perp(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> F_ratio_para(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> F_ratio_long(const complex<double> &) const
                {
                    return 0.0;
                }

                virtual complex<double> get_orthonormal_perp_coefficients(const unsigned &) const
                {
                    return 0.0;
                }

                virtual complex<double> get_orthonormal_para_coefficients(const unsigned &) const
                {
                    return 0.0;
                }

                virtual complex<double> get_orthonormal_long_coefficients(const unsigned &) const
                {
                    return 0.0;
                }

                virtual double weak_bound() const
                {
                    return 0.0;
                }

                virtual double strong_bound() const
                {
                    return 0.0;
                }

                virtual double weak_bound_log_likelihood() const
                {
                    return 0.0;
                }

                virtual double strong_bound_log_likelihood() const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V1(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V2(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V23(const double &) const
                {
                    return 0.0;
                }

                static NonlocalFormFactorPtr<PToV> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<PToV>(new Naive(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    return {};
                }
        };

        /*
         * Parametrize the entire formfactor, i.e., both leading and all sub-leading powers as described in [GvDV:2020].
         */
        template <typename Process_>
        class GvDV2020 :
            public NonlocalFormFactor<PToV>
        {
            private:
                std::shared_ptr<FormFactors<PToV>> form_factors;

                // spectator quark option
                SwitchOption opt_q;

                //Polynomial expansion parameters
                UsedParameter re_alpha_0_perp;
                UsedParameter im_alpha_0_perp;
                UsedParameter re_alpha_1_perp;
                UsedParameter im_alpha_1_perp;
                UsedParameter re_alpha_2_perp;
                UsedParameter im_alpha_2_perp;
                UsedParameter re_alpha_3_perp;
                UsedParameter im_alpha_3_perp;
                UsedParameter re_alpha_4_perp;
                UsedParameter im_alpha_4_perp;
                UsedParameter re_alpha_5_perp;
                UsedParameter im_alpha_5_perp;

                UsedParameter re_alpha_0_para;
                UsedParameter im_alpha_0_para;
                UsedParameter re_alpha_1_para;
                UsedParameter im_alpha_1_para;
                UsedParameter re_alpha_2_para;
                UsedParameter im_alpha_2_para;
                UsedParameter re_alpha_3_para;
                UsedParameter im_alpha_3_para;
                UsedParameter re_alpha_4_para;
                UsedParameter im_alpha_4_para;
                UsedParameter re_alpha_5_para;
                UsedParameter im_alpha_5_para;

                UsedParameter re_alpha_0_long;
                UsedParameter im_alpha_0_long;
                UsedParameter re_alpha_1_long;
                UsedParameter im_alpha_1_long;
                UsedParameter re_alpha_2_long;
                UsedParameter im_alpha_2_long;
                UsedParameter re_alpha_3_long;
                UsedParameter im_alpha_3_long;
                UsedParameter re_alpha_4_long;
                UsedParameter im_alpha_4_long;
                UsedParameter re_alpha_5_long;
                UsedParameter im_alpha_5_long;

                //Charmonium masses
                UsedParameter m_Jpsi;
                UsedParameter m_psi2S;

                // B-meson parameters
                UsedParameter m_B;

                // final state meson parameters
                UsedParameter m_V;

                UsedParameter m_D0;
                UsedParameter t_0;

                // Subtraction point for the dispersion relation...
                UsedParameter t_s;
                // ...and value of the dispersion bound at that point in the OPE
                UsedParameter chiOPE;

                // Orthogonal polynomials on an arc of the unit circle
                std::shared_ptr<SzegoPolynomial<5u>> polynomials;

                std::string _final_state() const
                {
                    switch (opt_q.value()[0])
                    {
                        case 's':
                            return "phi";
                            break;

                        default:
                            return "K_d^*";
                    }
                }

            public:
                GvDV2020(const Parameters & p, const Options & o) :
                    form_factors(FormFactorFactory<PToV>::create(stringify(Process_::label) + "::" + o.get("form-factors"_ok, "BSZ2015"), p)),
                    opt_q(o, "q"_ok, { "u", "d", "s" }),

                    re_alpha_0_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_0^perp}@GvDV2020"], *this),
                    im_alpha_0_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_0^perp}@GvDV2020"], *this),
                    re_alpha_1_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_1^perp}@GvDV2020"], *this),
                    im_alpha_1_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_1^perp}@GvDV2020"], *this),
                    re_alpha_2_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_2^perp}@GvDV2020"], *this),
                    im_alpha_2_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_2^perp}@GvDV2020"], *this),
                    re_alpha_3_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_3^perp}@GvDV2020"], *this),
                    im_alpha_3_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_3^perp}@GvDV2020"], *this),
                    re_alpha_4_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_4^perp}@GvDV2020"], *this),
                    im_alpha_4_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_4^perp}@GvDV2020"], *this),
                    re_alpha_5_perp(p[stringify(Process_::label) + "ccbar::Re{alpha_5^perp}@GvDV2020"], *this),
                    im_alpha_5_perp(p[stringify(Process_::label) + "ccbar::Im{alpha_5^perp}@GvDV2020"], *this),

                    re_alpha_0_para(p[stringify(Process_::label) + "ccbar::Re{alpha_0^para}@GvDV2020"], *this),
                    im_alpha_0_para(p[stringify(Process_::label) + "ccbar::Im{alpha_0^para}@GvDV2020"], *this),
                    re_alpha_1_para(p[stringify(Process_::label) + "ccbar::Re{alpha_1^para}@GvDV2020"], *this),
                    im_alpha_1_para(p[stringify(Process_::label) + "ccbar::Im{alpha_1^para}@GvDV2020"], *this),
                    re_alpha_2_para(p[stringify(Process_::label) + "ccbar::Re{alpha_2^para}@GvDV2020"], *this),
                    im_alpha_2_para(p[stringify(Process_::label) + "ccbar::Im{alpha_2^para}@GvDV2020"], *this),
                    re_alpha_3_para(p[stringify(Process_::label) + "ccbar::Re{alpha_3^para}@GvDV2020"], *this),
                    im_alpha_3_para(p[stringify(Process_::label) + "ccbar::Im{alpha_3^para}@GvDV2020"], *this),
                    re_alpha_4_para(p[stringify(Process_::label) + "ccbar::Re{alpha_4^para}@GvDV2020"], *this),
                    im_alpha_4_para(p[stringify(Process_::label) + "ccbar::Im{alpha_4^para}@GvDV2020"], *this),
                    re_alpha_5_para(p[stringify(Process_::label) + "ccbar::Re{alpha_5^para}@GvDV2020"], *this),
                    im_alpha_5_para(p[stringify(Process_::label) + "ccbar::Im{alpha_5^para}@GvDV2020"], *this),

                    re_alpha_0_long(p[stringify(Process_::label) + "ccbar::Re{alpha_0^long}@GvDV2020"], *this),
                    im_alpha_0_long(p[stringify(Process_::label) + "ccbar::Im{alpha_0^long}@GvDV2020"], *this),
                    re_alpha_1_long(p[stringify(Process_::label) + "ccbar::Re{alpha_1^long}@GvDV2020"], *this),
                    im_alpha_1_long(p[stringify(Process_::label) + "ccbar::Im{alpha_1^long}@GvDV2020"], *this),
                    re_alpha_2_long(p[stringify(Process_::label) + "ccbar::Re{alpha_2^long}@GvDV2020"], *this),
                    im_alpha_2_long(p[stringify(Process_::label) + "ccbar::Im{alpha_2^long}@GvDV2020"], *this),
                    re_alpha_3_long(p[stringify(Process_::label) + "ccbar::Re{alpha_3^long}@GvDV2020"], *this),
                    im_alpha_3_long(p[stringify(Process_::label) + "ccbar::Im{alpha_3^long}@GvDV2020"], *this),
                    re_alpha_4_long(p[stringify(Process_::label) + "ccbar::Re{alpha_4^long}@GvDV2020"], *this),
                    im_alpha_4_long(p[stringify(Process_::label) + "ccbar::Im{alpha_4^long}@GvDV2020"], *this),
                    re_alpha_5_long(p[stringify(Process_::label) + "ccbar::Re{alpha_5^long}@GvDV2020"], *this),
                    im_alpha_5_long(p[stringify(Process_::label) + "ccbar::Im{alpha_5^long}@GvDV2020"], *this),

                    m_Jpsi(p["mass::J/psi"], *this),
                    m_psi2S(p["mass::psi(2S)"], *this),

                    m_B(p["mass::B_" + opt_q.value()], *this),

                    m_V(p["mass::" + _final_state()], *this),

                    m_D0(p["mass::D^0"], *this),
                    t_0(p["b->sccbar::t_0"], *this),

                    t_s(p["b->sccbar::t_s"], *this),

                    chiOPE(p["b->sccbar::chiOPE@GvDV2020"], *this),

                    // The parameters of the polynomial expension are computed using t0 = 4.0 and
                    // the masses are set to the same values as for local form-factors
                    polynomials(PolynomialsFactory::create(opt_q.value()))
                {
                    this->uses(*form_factors);
                }

                ~GvDV2020() = default;

                inline complex<double> phi(const complex<double> & q2, const std::array<unsigned, 4> & phi_parameters) const
                {
                    // Values of a, b, c and d depends on the form factor:
                    // FF                        a    b    c    d
                    // 0(P->P) aka plus          3    3    2    2
                    // perp(P->V) = par(P->V)    3    1    3    0
                    // 0(P->V) aka long          3    1    2    2

                    const double m_V2  = power_of<2>(m_V);
                    const double m_B2  = power_of<2>(m_B),  m_B4 =  power_of<4>(m_B);
                    const double m_D02 = power_of<2>(m_D0), m_D04 = power_of<4>(m_D0);
                    const double s_0   = this->t_0();
                    const auto   z     = eos::nff_utils::z(q2, 4.0 * power_of<2>(m_D0), s_0);
                    const double Q2   = this->t_s();
                    const double chi = this->chiOPE();

                    const double a = phi_parameters[0], b = phi_parameters[1], c = phi_parameters[2], d = phi_parameters[3];

                    const double Nlambda = 4 * M_PI * pow(m_B2, 0.5 * (a - b + c + d) - 1.) * pow(2 * (4 * m_D02 - s_0) / 3 / chi, 0.5); //(C6)
                    const complex<double> phi1 = -pow(2 * pow((4 * m_D02 - Q2) * (4 * m_D02 - s_0), 0.5) + 8 * m_D02 - Q2 - s_0, 0.5) /
                                                (2 * pow((4 * m_D02 - Q2) * (4 * m_D02 - s_0), 0.5) + 8 * m_D02 + Q2 * (z - 1.) - s_0 * (z + 1.)); //(C7)
                    const complex<double> phi2 = pow(m_B4 * power_of<4>(z - 1.) - 2 * m_B2 * power_of<2>(z - 1.) * (-16 * m_D02 * z + m_V2 * power_of<2>(z - 1.) + s_0 * power_of<2>(z + 1.)) +
						power_of<2>(16 * m_D02 * z + m_V2 * power_of<2>(z - 1.) - s_0 * power_of<2>(z + 1.)), 0.5);//(C8)
                    const complex<double> phi3 = pow(8 * m_D02 + 4 * pow(4 * m_D04 - s_0 * m_D02, 0.5) - s_0, 0.5) /
                                                (-8 * m_D02 - 4 * pow(4 * m_D04 - s_0 * m_D02, 0.5) + s_0*(z+1.)); //(C9)
                    const complex<double> phi4 = pow(s_0 * power_of<2>(z + 1.) - 16. * z * m_D02, -0.5); //(C10)

                    return Nlambda * pow(1.+z, 0.5) * pow(1.-z, a-b+c+d-1.5) * pow(phi1, a) * pow(phi2, 0.5*b) * pow(phi3, c) * pow(phi4, d); //(C5)
                }

                inline complex<double> phi(const double & q2, const std::array<unsigned, 4> & phi_parameters) const
                {
                    return phi(complex<double>(q2, 0.0), phi_parameters);
                }

                // Residue of H at s = m_Jpsi2 computed as the residue wrt z -z_Jpsi divided by dz/ds evaluated at s = m_Jpsi2
                inline complex<double> H_residue_jpsi(const std::array<unsigned, 4> & phi_parameters, const std::array<complex<double>, 6> & alpha) const
                {
                    const double m_Jpsi2  = power_of<2>(m_Jpsi);
                    const double m_psi2S2 = power_of<2>(m_psi2S);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2, s_p, s_0);

                    const auto & polynomials_at_z = (*polynomials)(z_Jpsi);
                    const complex<double> p_at_z = std::inner_product(alpha.begin(), alpha.end(), polynomials_at_z.begin(), complex<double>(0, 0));

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_Jpsi2, -0.5) * pow(pow(s_p - m_Jpsi2, 0.5) + pow(s_p - s_0, 0.5), -2);

                    return p_at_z / phi(m_Jpsi2, phi_parameters) * (1 - norm(z_Jpsi)) * (1. - z_Jpsi * std::conj(z_psi2S)) / (z_Jpsi - z_psi2S) / dzds;
                }

                // Residue of H at s = m_psi2S2 computed as the residue wrt z -z_psi2S divided by dz/ds evaluated at s = m_psi2S2
                inline complex<double> H_residue_psi2s(const std::array<unsigned, 4> & phi_parameters, const std::array<complex<double>, 6> & alpha) const
                {
                    const double m_Jpsi2  = power_of<2>(m_Jpsi);
                    const double m_psi2S2 = power_of<2>(m_psi2S);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2, s_p, s_0);

                    const auto & polynomials_at_z = (*polynomials)(z_psi2S);
                    const complex<double> p_at_z = std::inner_product(alpha.begin(), alpha.end(), polynomials_at_z.begin(), complex<double>(0, 0));

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_psi2S2, -0.5) * pow(pow(s_p - m_psi2S2, 0.5) + pow(s_p - s_0, 0.5), -2);

                    return p_at_z / phi(m_psi2S2, phi_parameters) *(1 - norm(z_psi2S)) * (1. - z_psi2S * std::conj(z_Jpsi)) / (z_psi2S - z_Jpsi) / dzds;
                }

                virtual complex<double> H_perp(const complex<double> & q2) const
                {
                    const std::array<complex<double>, 6> alpha_perp{
                        complex<double>(re_alpha_0_perp, im_alpha_0_perp),
                        complex<double>(re_alpha_1_perp, im_alpha_1_perp),
                        complex<double>(re_alpha_2_perp, im_alpha_2_perp),
                        complex<double>(re_alpha_3_perp, im_alpha_3_perp),
                        complex<double>(re_alpha_4_perp, im_alpha_4_perp),
                        complex<double>(re_alpha_5_perp, im_alpha_5_perp),
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2,                     s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(power_of<2>(m_Jpsi),    s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(power_of<2>(m_psi2S),   s_p, s_0);

                    const auto & polynomials_at_z = (*polynomials)(z);
                    const complex<double> p_at_z = std::inner_product(alpha_perp.begin(), alpha_perp.end(), polynomials_at_z.begin(), complex<double>(0, 0));

                    const complex<double> blaschke_factor = eos::nff_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 3, 0};

                    return p_at_z / phi(q2, phi_parameters) / blaschke_factor;
                }

                virtual complex<double> H_perp(const double & q2) const
                {
                    return H_perp(complex<double>(q2, 0.0));
                }

                virtual complex<double> Hhat_perp(const double & q2) const
                {
                    const std::array<complex<double>, 6> alpha_perp{
                        complex<double>(re_alpha_0_perp, im_alpha_0_perp),
                        complex<double>(re_alpha_1_perp, im_alpha_1_perp),
                        complex<double>(re_alpha_2_perp, im_alpha_2_perp),
                        complex<double>(re_alpha_3_perp, im_alpha_3_perp),
                        complex<double>(re_alpha_4_perp, im_alpha_4_perp),
                        complex<double>(re_alpha_5_perp, im_alpha_5_perp),
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);

                    const auto & polynomials_at_z = (*polynomials)(z);

                    return std::inner_product(alpha_perp.begin(), alpha_perp.end(), polynomials_at_z.begin(), complex<double>(0, 0));
                }

                virtual complex<double> H_para(const complex<double> & q2) const
                {
                    const std::array<complex<double>, 6> alpha_para{
                        complex<double>(re_alpha_0_para, im_alpha_0_para),
                        complex<double>(re_alpha_1_para, im_alpha_1_para),
                        complex<double>(re_alpha_2_para, im_alpha_2_para),
                        complex<double>(re_alpha_3_para, im_alpha_3_para),
                        complex<double>(re_alpha_4_para, im_alpha_4_para),
                        complex<double>(re_alpha_5_para, im_alpha_5_para),
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2,                     s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(power_of<2>(m_Jpsi),    s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(power_of<2>(m_psi2S),   s_p, s_0);

                    const complex<double> blaschke_factor = eos::nff_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 3, 0};

                    const auto & polynomials_at_z = (*polynomials)(z);
                    const complex<double> p_at_z = std::inner_product(alpha_para.begin(), alpha_para.end(), polynomials_at_z.begin(), complex<double>(0, 0));

                    return p_at_z / phi(q2, phi_parameters) / blaschke_factor;
                }

                virtual complex<double> H_para(const double & q2) const
                {
                    return H_para(complex<double>(q2, 0.0));
                }

                virtual complex<double> Hhat_para(const double & q2) const
                {
                    const std::array<complex<double>, 6> alpha_para{
                        complex<double>(re_alpha_0_para, im_alpha_0_para),
                        complex<double>(re_alpha_1_para, im_alpha_1_para),
                        complex<double>(re_alpha_2_para, im_alpha_2_para),
                        complex<double>(re_alpha_3_para, im_alpha_3_para),
                        complex<double>(re_alpha_4_para, im_alpha_4_para),
                        complex<double>(re_alpha_5_para, im_alpha_5_para),
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);

                    const auto & polynomials_at_z = (*polynomials)(z);

                    return std::inner_product(alpha_para.begin(), alpha_para.end(), polynomials_at_z.begin(), complex<double>(0, 0));
                }

                virtual complex<double> H_long(const complex<double> & q2) const
                {
                    const std::array<complex<double>, 6> alpha_long{
                        complex<double>(re_alpha_0_long, im_alpha_0_long),
                        complex<double>(re_alpha_1_long, im_alpha_1_long),
                        complex<double>(re_alpha_2_long, im_alpha_2_long),
                        complex<double>(re_alpha_3_long, im_alpha_3_long),
                        complex<double>(re_alpha_4_long, im_alpha_4_long),
                        complex<double>(re_alpha_5_long, im_alpha_5_long),
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2,                     s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(power_of<2>(m_Jpsi),    s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(power_of<2>(m_psi2S),   s_p, s_0);

                    const complex<double> blaschke_factor = eos::nff_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 2, 2};

                    const auto & polynomials_at_z = (*polynomials)(z);
                    const complex<double> p_at_z = std::inner_product(alpha_long.begin(), alpha_long.end(), polynomials_at_z.begin(), complex<double>(0, 0));

                    return p_at_z / phi(q2, phi_parameters) / blaschke_factor;
                }

                virtual complex<double> H_long(const double & q2) const
                {
                    return H_long(complex<double>(q2, 0.0));
                }

                virtual complex<double> Hhat_long(const double & q2) const
                {
                    const std::array<complex<double>, 6> alpha_long{
                        complex<double>(re_alpha_0_long, im_alpha_0_long),
                        complex<double>(re_alpha_1_long, im_alpha_1_long),
                        complex<double>(re_alpha_2_long, im_alpha_2_long),
                        complex<double>(re_alpha_3_long, im_alpha_3_long),
                        complex<double>(re_alpha_4_long, im_alpha_4_long),
                        complex<double>(re_alpha_5_long, im_alpha_5_long),
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);

                    const auto & polynomials_at_z = (*polynomials)(z);

                    return std::inner_product(alpha_long.begin(), alpha_long.end(), polynomials_at_z.begin(), complex<double>(0, 0));
                }


                virtual complex<double> H_perp_residue_jpsi() const
                {
                    const std::array<complex<double>, 6> alpha_perp{
                        complex<double>(re_alpha_0_perp, im_alpha_0_perp),
                        complex<double>(re_alpha_1_perp, im_alpha_1_perp),
                        complex<double>(re_alpha_2_perp, im_alpha_2_perp),
                        complex<double>(re_alpha_3_perp, im_alpha_3_perp),
                        complex<double>(re_alpha_4_perp, im_alpha_4_perp),
                        complex<double>(re_alpha_5_perp, im_alpha_5_perp),
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 3, 0};

                    return H_residue_jpsi(phi_parameters, alpha_perp);
                }

                virtual complex<double> H_perp_residue_psi2s() const
                {
                    const std::array<complex<double>, 6> alpha_perp{
                        complex<double>(re_alpha_0_perp, im_alpha_0_perp),
                        complex<double>(re_alpha_1_perp, im_alpha_1_perp),
                        complex<double>(re_alpha_2_perp, im_alpha_2_perp),
                        complex<double>(re_alpha_3_perp, im_alpha_3_perp),
                        complex<double>(re_alpha_4_perp, im_alpha_4_perp),
                        complex<double>(re_alpha_5_perp, im_alpha_5_perp),
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 3, 0};

                    return H_residue_psi2s(phi_parameters, alpha_perp);
                }

                virtual complex<double> H_para_residue_jpsi() const
                {
                    const std::array<complex<double>, 6> alpha_para{
                        complex<double>(re_alpha_0_para, im_alpha_0_para),
                        complex<double>(re_alpha_1_para, im_alpha_1_para),
                        complex<double>(re_alpha_2_para, im_alpha_2_para),
                        complex<double>(re_alpha_3_para, im_alpha_3_para),
                        complex<double>(re_alpha_4_para, im_alpha_4_para),
                        complex<double>(re_alpha_5_para, im_alpha_5_para),
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 3, 0};

                    return H_residue_jpsi(phi_parameters, alpha_para);
                }

                virtual complex<double> H_para_residue_psi2s() const
                {
                    const std::array<complex<double>, 6> alpha_para{
                        complex<double>(re_alpha_0_para, im_alpha_0_para),
                        complex<double>(re_alpha_1_para, im_alpha_1_para),
                        complex<double>(re_alpha_2_para, im_alpha_2_para),
                        complex<double>(re_alpha_3_para, im_alpha_3_para),
                        complex<double>(re_alpha_4_para, im_alpha_4_para),
                        complex<double>(re_alpha_5_para, im_alpha_5_para),
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 3, 0};

                    return H_residue_psi2s(phi_parameters, alpha_para);
                }

                virtual complex<double> H_long_residue_jpsi() const
                {
                    const std::array<complex<double>, 6> alpha_long{
                        complex<double>(re_alpha_0_long, im_alpha_0_long),
                        complex<double>(re_alpha_1_long, im_alpha_1_long),
                        complex<double>(re_alpha_2_long, im_alpha_2_long),
                        complex<double>(re_alpha_3_long, im_alpha_3_long),
                        complex<double>(re_alpha_4_long, im_alpha_4_long),
                        complex<double>(re_alpha_5_long, im_alpha_5_long),
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 2, 2};

                    return H_residue_jpsi(phi_parameters, alpha_long);
                }

                virtual complex<double> H_long_residue_psi2s() const
                {
                    const std::array<complex<double>, 6> alpha_long{
                        complex<double>(re_alpha_0_long, im_alpha_0_long),
                        complex<double>(re_alpha_1_long, im_alpha_1_long),
                        complex<double>(re_alpha_2_long, im_alpha_2_long),
                        complex<double>(re_alpha_3_long, im_alpha_3_long),
                        complex<double>(re_alpha_4_long, im_alpha_4_long),
                        complex<double>(re_alpha_5_long, im_alpha_5_long),
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 2, 2};

                    return H_residue_psi2s(phi_parameters, alpha_long);
                }

                virtual complex<double> ratio_perp(const complex<double> & q2) const
                {
                    const double m_V2  = power_of<2>(m_V);
                    const double m_B2  = power_of<2>(m_B);
                    const complex<double> lambda = eos::lambda(complex<double>(m_B2, 0.0), complex<double>(m_V2, 0.0), q2);
                    const complex<double> F_perp = pow(2.0 * lambda, 0.5) / (m_B + m_V) / m_B() * form_factors->v(q2);

                    return H_perp(q2) / F_perp;
                }

                virtual complex<double> ratio_perp(const double & q2) const
                {
                    return ratio_perp(complex<double>(q2, 0.0));
                }

                virtual complex<double> F_ratio_perp(const complex<double> & q2) const
                {
                    return (m_B + m_V) / m_B * form_factors->t_1(q2) / form_factors->v(q2);
                }

                virtual complex<double> ratio_para(const complex<double> & q2) const
                {
                    const complex<double> F_para = sqrt(2) * (m_B + m_V) / m_B * form_factors->a_1(q2);

                    return H_para(q2) / F_para;
                }

                virtual complex<double> ratio_para(const double & q2) const
                {
                    return ratio_para(complex<double>(q2, 0.0));
                }

                virtual complex<double> F_ratio_para(const complex<double> & q2) const
                {
                    return (m_B - m_V) / m_B * form_factors->t_2(q2) / form_factors->a_1(q2);
                }

                virtual complex<double> ratio_long(const complex<double> & q2) const
                {
                    const double m_V2  = power_of<2>(m_V);
                    const double m_B2  = power_of<2>(m_B);
                    const complex<double> lambda = eos::lambda(complex<double>(m_B2, 0.0), complex<double>(m_V2, 0.0), q2);
                    const complex<double> F_long = ((m_B2 - m_V2 - q2) * power_of<2>(m_B + m_V) * form_factors->a_1(q2)
                            - lambda * form_factors->a_2(q2)) / (2 * m_V * m_B2 * (m_B + m_V));

                    return H_long(q2) / F_long;
                }

                virtual complex<double> ratio_long(const double & q2) const
                {
                    return ratio_long(complex<double>(q2, 0.0));
                }

                virtual complex<double> F_ratio_long(const complex<double> & q2) const
                {
                    const double m_V2  = power_of<2>(m_V);
                    const double m_B2  = power_of<2>(m_B);
                    const complex<double> lambda = eos::lambda(complex<double>(m_B2, 0.0), complex<double>(m_V2, 0.0), q2);
                    const complex<double> F_long = ((m_B2 - m_V2 - q2) * power_of<2>(m_B + m_V) * form_factors->a_1(q2)
                            - lambda * form_factors->a_2(q2)) / (2 * m_V * m_B2 * (m_B + m_V));
                    const complex<double> F_T_long = q2 * ((m_B2 + 3 * m_V2 - q2) * (m_B2 - m_V2) * form_factors->t_2(q2)
                            - lambda * form_factors->a_2(q2)) / (2 * m_V * m_B2 * m_B * (m_B2 - m_V2));

                    return F_T_long / F_long;
                }

                virtual complex<double> normalized_moment_V1(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V2(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V23(const double &) const
                {
                    return 0.0;
                }                virtual complex<double> get_orthonormal_perp_coefficients(const unsigned & i) const
                {
                    const std::array<complex<double>, 6> alpha_perp{
                        complex<double>(re_alpha_0_perp, im_alpha_0_perp),
                        complex<double>(re_alpha_1_perp, im_alpha_1_perp),
                        complex<double>(re_alpha_2_perp, im_alpha_2_perp),
                        complex<double>(re_alpha_3_perp, im_alpha_3_perp),
                        complex<double>(re_alpha_4_perp, im_alpha_4_perp),
                        complex<double>(re_alpha_5_perp, im_alpha_5_perp),
                    };

                    return alpha_perp[i];
                }

                virtual complex<double> get_orthonormal_para_coefficients(const unsigned & i) const
                {
                    const std::array<complex<double>, 6> alpha_para{
                        complex<double>(re_alpha_0_para, im_alpha_0_para),
                        complex<double>(re_alpha_1_para, im_alpha_1_para),
                        complex<double>(re_alpha_2_para, im_alpha_2_para),
                        complex<double>(re_alpha_3_para, im_alpha_3_para),
                        complex<double>(re_alpha_4_para, im_alpha_4_para),
                        complex<double>(re_alpha_5_para, im_alpha_5_para),
                    };

                    return alpha_para[i];
                }

                virtual complex<double> get_orthonormal_long_coefficients(const unsigned & i) const
                {
                    const std::array<complex<double>, 6> alpha_long{
                        complex<double>(re_alpha_0_long, im_alpha_0_long),
                        complex<double>(re_alpha_1_long, im_alpha_1_long),
                        complex<double>(re_alpha_2_long, im_alpha_2_long),
                        complex<double>(re_alpha_3_long, im_alpha_3_long),
                        complex<double>(re_alpha_4_long, im_alpha_4_long),
                        complex<double>(re_alpha_5_long, im_alpha_5_long),
                    };

                    return alpha_long[i];
                }

                virtual double weak_bound() const
                {
                    return 0.0;
                }

                virtual double strong_bound() const
                {
                    return 0.0;
                }

                virtual double weak_bound_log_likelihood() const
                {
                    return 0.0;
                }

                virtual double strong_bound_log_likelihood() const
                {
                    return 0.0;
                }

                static NonlocalFormFactorPtr<PToV> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<PToV>(new GvDV2020<Process_>(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    Diagnostics results;

                    const std::array<unsigned, 4> phi_parameters_long = {3, 1, 2, 2}; //long polarization
                    results.add({ real(1./this->phi(0.0, phi_parameters_long)), "Re{1/phi_long(q2 = 0.0)}" });
                    results.add({ imag(1./this->phi(0.0, phi_parameters_long)), "Im{1/phi_long(q2 = 0.0)}" });
                    results.add({ real(this->phi(16.0, phi_parameters_long)), "Re{phi_long(q2 = 16.0)}" });
                    results.add({ imag(this->phi(16.0, phi_parameters_long)), "Im{phi_long(q2 = 16.0)}" });

                    const std::array<unsigned, 4> phi_parameters_perp = {3, 1, 3, 0}; //perp or para polarization
                    results.add({ real(this->phi(16.0, phi_parameters_perp)), "Re{phi_perp(q2 = 16.0)}" });
                    results.add({ imag(this->phi(16.0, phi_parameters_perp)), "Im{phi_perp(q2 = 16.0)}" });

                    return results;
                }
        };


        /*
         * Parametrize the entire formfactor, i.e., both leading and all sub-leading powers as described in [GvDV:2020],
         * but using an ad hoc Lagrange polynomial.
         */
        template <typename Process_>
        class GRvDV2022order5 :
            public NonlocalFormFactor<PToV>
        {
            private:
                std::shared_ptr<FormFactors<PToV>> form_factors;

                // spectator quark option
                SwitchOption opt_q;

                //Polynomial expansion parameters
                UsedParameter re_at_m7_perp;
                UsedParameter im_at_m7_perp;
                UsedParameter re_at_m5_perp;
                UsedParameter im_at_m5_perp;
                UsedParameter re_at_m3_perp;
                UsedParameter im_at_m3_perp;
                UsedParameter re_at_m1_perp;
                UsedParameter im_at_m1_perp;
                UsedParameter abs_at_Jpsi_perp;
                UsedParameter arg_at_Jpsi_perp_minus_long;
                UsedParameter abs_at_psi2S_perp;
                UsedParameter arg_at_psi2S_perp_minus_long;

                UsedParameter re_at_m7_para;
                UsedParameter im_at_m7_para;
                UsedParameter re_at_m5_para;
                UsedParameter im_at_m5_para;
                UsedParameter re_at_m3_para;
                UsedParameter im_at_m3_para;
                UsedParameter re_at_m1_para;
                UsedParameter im_at_m1_para;
                UsedParameter abs_at_Jpsi_para;
                UsedParameter arg_at_Jpsi_para_minus_long;
                UsedParameter abs_at_psi2S_para;
                UsedParameter arg_at_psi2S_para_minus_long;

                UsedParameter re_at_m7_long;
                UsedParameter im_at_m7_long;
                UsedParameter re_at_m5_long;
                UsedParameter im_at_m5_long;
                UsedParameter re_at_m3_long;
                UsedParameter im_at_m3_long;
                UsedParameter re_at_m1_long;
                UsedParameter im_at_m1_long;
                UsedParameter abs_at_Jpsi_long;
                UsedParameter arg_at_Jpsi_long;
                UsedParameter abs_at_psi2S_long;
                UsedParameter arg_at_psi2S_long;

                //Charmonium masses
                UsedParameter m_Jpsi;
                UsedParameter m_psi2S;

                // B-meson parameters
                UsedParameter m_B;

                // final state meson parameters
                UsedParameter m_V;

                UsedParameter m_D0;
                UsedParameter t_0;

                // Subtraction point for the dispersion relation...
                UsedParameter t_s;
                // ...and value of the dispersion bound at that point in the OPE
                UsedParameter chiOPE;
                UsedParameter bound;
                UsedParameter bound_uncertainty;

                // Lagrange interpolating polynomial
                const static unsigned interpolation_order = 5;
                const LagrangePolynomial<interpolation_order> lagrange;

                // Orthogonal polynomials on an arc of the unit circle used for the computation of dispersive bounds
                std::shared_ptr<SzegoPolynomial<5u>> orthonormal_polynomials;

                std::string _final_state() const
                {
                    switch (opt_q.value()[0])
                    {
                        case 's':
                            return "phi";
                            break;

                        default:
                            return "K_d^*";
                    }
                }

            public:
                GRvDV2022order5(const Parameters & p, const Options & o) :
                    form_factors(FormFactorFactory<PToV>::create(stringify(Process_::label) + "::" + o.get("form-factors"_ok, "BSZ2015"), p)),
                    opt_q(o, "q"_ok, { "u", "d", "s" }),

                    re_at_m7_perp(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_m7_perp@GRvDV2022"], *this),
                    im_at_m7_perp(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_m7_perp@GRvDV2022"], *this),
                    re_at_m5_perp(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_m5_perp@GRvDV2022"], *this),
                    im_at_m5_perp(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_m5_perp@GRvDV2022"], *this),
                    re_at_m3_perp(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_m3_perp@GRvDV2022"], *this),
                    im_at_m3_perp(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_m3_perp@GRvDV2022"], *this),
                    re_at_m1_perp(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_m1_perp@GRvDV2022"], *this),
                    im_at_m1_perp(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_m1_perp@GRvDV2022"], *this),
                    abs_at_Jpsi_perp(p[stringify(Process_::label) + "ccbar::Abs_Hhat_at_Jpsi_perp@GRvDV2022"], *this),
                    arg_at_Jpsi_perp_minus_long(p[stringify(Process_::label) + "ccbar::Arg_Hhat_at_Jpsi_perp_minus_long@GRvDV2022"], *this),
                    abs_at_psi2S_perp(p[stringify(Process_::label) + "ccbar::Abs_Hhat_at_psi2S_perp@GRvDV2022"], *this),
                    arg_at_psi2S_perp_minus_long(p[stringify(Process_::label) + "ccbar::Arg_Hhat_at_psi2S_perp_minus_long@GRvDV2022"], *this),

                    re_at_m7_para(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_m7_para@GRvDV2022"], *this),
                    im_at_m7_para(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_m7_para@GRvDV2022"], *this),
                    re_at_m5_para(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_m5_para@GRvDV2022"], *this),
                    im_at_m5_para(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_m5_para@GRvDV2022"], *this),
                    re_at_m3_para(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_m3_para@GRvDV2022"], *this),
                    im_at_m3_para(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_m3_para@GRvDV2022"], *this),
                    re_at_m1_para(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_m1_para@GRvDV2022"], *this),
                    im_at_m1_para(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_m1_para@GRvDV2022"], *this),
                    abs_at_Jpsi_para(p[stringify(Process_::label) + "ccbar::Abs_Hhat_at_Jpsi_para@GRvDV2022"], *this),
                    arg_at_Jpsi_para_minus_long(p[stringify(Process_::label) + "ccbar::Arg_Hhat_at_Jpsi_para_minus_long@GRvDV2022"], *this),
                    abs_at_psi2S_para(p[stringify(Process_::label) + "ccbar::Abs_Hhat_at_psi2S_para@GRvDV2022"], *this),
                    arg_at_psi2S_para_minus_long(p[stringify(Process_::label) + "ccbar::Arg_Hhat_at_psi2S_para_minus_long@GRvDV2022"], *this),

                    re_at_m7_long(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_m7_long@GRvDV2022"], *this),
                    im_at_m7_long(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_m7_long@GRvDV2022"], *this),
                    re_at_m5_long(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_m5_long@GRvDV2022"], *this),
                    im_at_m5_long(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_m5_long@GRvDV2022"], *this),
                    re_at_m3_long(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_m3_long@GRvDV2022"], *this),
                    im_at_m3_long(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_m3_long@GRvDV2022"], *this),
                    re_at_m1_long(p[stringify(Process_::label) + "ccbar::Re_Hhat_at_m1_long@GRvDV2022"], *this),
                    im_at_m1_long(p[stringify(Process_::label) + "ccbar::Im_Hhat_at_m1_long@GRvDV2022"], *this),
                    abs_at_Jpsi_long(p[stringify(Process_::label) + "ccbar::Abs_Hhat_at_Jpsi_long@GRvDV2022"], *this),
                    arg_at_Jpsi_long(p[stringify(Process_::label) + "ccbar::Arg_Hhat_at_Jpsi_long@GRvDV2022"], *this),
                    abs_at_psi2S_long(p[stringify(Process_::label) + "ccbar::Abs_Hhat_at_psi2S_long@GRvDV2022"], *this),
                    arg_at_psi2S_long(p[stringify(Process_::label) + "ccbar::Arg_Hhat_at_psi2S_long@GRvDV2022"], *this),

                    m_Jpsi(p["mass::J/psi"], *this),
                    m_psi2S(p["mass::psi(2S)"], *this),

                    m_B(p["mass::B_" + opt_q.value()], *this),

                    m_V(p["mass::" + _final_state()], *this),

                    m_D0(p["mass::D^0"], *this),
                    t_0(p["b->sccbar::t_0"], *this),

                    t_s(p["b->sccbar::t_s"], *this),

                    chiOPE(p["b->sccbar::chiOPE@GvDV2020"], *this),
                    bound(p["b->sccbar::bound@GvDV2020"], *this),
                    bound_uncertainty(p["b->sccbar::bound_uncertainty@GvDV2020"], *this),

                    lagrange({eos::nff_utils::z(-7.0, 4.0 * power_of<2>(m_D0), t_0()),
                              eos::nff_utils::z(-5.0, 4.0 * power_of<2>(m_D0), t_0()),
                              eos::nff_utils::z(-3.0, 4.0 * power_of<2>(m_D0), t_0()),
                              eos::nff_utils::z(-1.0, 4.0 * power_of<2>(m_D0), t_0()),
                              eos::nff_utils::z(power_of<2>(m_Jpsi),  4.0 * power_of<2>(m_D0), t_0()),
                              eos::nff_utils::z(power_of<2>(m_psi2S), 4.0 * power_of<2>(m_D0), t_0())}),

                    // The parameters of the polynomial expension are computed using t0 = 4.0 and
                    // the masses are set to mB(s) = 5.279 (5.366) and mKst(phi) = 0.896 (1.02) (same values as for local form-factors)
                    orthonormal_polynomials(PolynomialsFactory::create(opt_q.value()))
                {
                    this->uses(*form_factors);
                }

                ~GRvDV2022order5() = default;

                inline complex<double> phi(const complex<double> & q2, const std::array<unsigned, 4> & phi_parameters) const
                {
                    // Values of a, b, c and d depends on the form factor:
                    // FF                        a    b    c    d
                    // 0(P->P) aka plus          3    3    2    2
                    // perp(P->V) = par(P->V)    3    1    3    0
                    // 0(P->V) aka long          3    1    2    2

                    const double m_V2  = power_of<2>(m_V);
                    const double m_B2  = power_of<2>(m_B),  m_B4 =  power_of<4>(m_B);
                    const double m_D02 = power_of<2>(m_D0), m_D04 = power_of<4>(m_D0);
                    const double s_0   = this->t_0();
                    const auto   z     = eos::nff_utils::z(q2, 4.0 * m_D02, s_0);
                    const double Q2   = this->t_s();
                    const double chi = this->chiOPE();

                    const double a = phi_parameters[0], b = phi_parameters[1], c = phi_parameters[2], d = phi_parameters[3];

                    const double Nlambda = 4 * M_PI * pow(m_B2, 0.5 * (a - b + c + d) - 1.) * pow(2 * (4 * m_D02 - s_0) / 3 / chi, 0.5); //(C6)
                    const complex<double> phi1 = -pow(2 * pow((4 * m_D02 - Q2) * (4 * m_D02 - s_0), 0.5) + 8 * m_D02 - Q2 - s_0, 0.5) /
                                                (2 * pow((4 * m_D02 - Q2) * (4 * m_D02 - s_0), 0.5) + 8 * m_D02 + Q2 * (z - 1.) - s_0 * (z + 1.)); //(C7)
                    const complex<double> phi2 = pow(m_B4 * power_of<4>(z - 1.) - 2 * m_B2 * power_of<2>(z - 1.) * (-16 * m_D02 * z + m_V2 * power_of<2>(z - 1.) + s_0 * power_of<2>(z + 1.)) +
						power_of<2>(16 * m_D02 * z + m_V2 * power_of<2>(z - 1.) - s_0 * power_of<2>(z + 1.)), 0.5);//(C8)
                    const complex<double> phi3 = pow(8 * m_D02 + 4 * pow(4 * m_D04 - s_0 * m_D02, 0.5) - s_0, 0.5) /
                                                (-8 * m_D02 - 4 * pow(4 * m_D04 - s_0 * m_D02, 0.5) + s_0*(z+1.)); //(C9)
                    const complex<double> phi4 = pow(s_0 * power_of<2>(z + 1.) - 16. * z * m_D02, -0.5); //(C10)

                    return Nlambda * pow(1.+z, 0.5) * pow(1.-z, a-b+c+d-1.5) * pow(phi1, a) * pow(phi2, 0.5*b) * pow(phi3, c) * pow(phi4, d); //(C5)
                }

                inline complex<double> phi(const double & q2, const std::array<unsigned, 4> & phi_parameters) const
                {
                    return phi(complex<double>(q2, 0.0), phi_parameters);
                }

                // Residue of H at s = m_Jpsi2 computed as the residue wrt z -z_Jpsi divided by dz/ds evaluated at s = m_Jpsi2
                inline complex<double> H_residue_jpsi(const std::array<unsigned, 4> & phi_parameters, const std::array<complex<double>, interpolation_order + 1> & interpolation_values) const
                {
                    const double m_Jpsi2  = power_of<2>(m_Jpsi);
                    const double m_psi2S2 = power_of<2>(m_psi2S);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2, s_p, s_0);

                    const complex<double> p_at_z = lagrange(interpolation_values, z_Jpsi);

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_Jpsi2, -0.5) * pow(pow(s_p - m_Jpsi2, 0.5) + pow(s_p - s_0, 0.5), -2);

                    return p_at_z / phi(m_Jpsi2, phi_parameters) * (1 - norm(z_Jpsi)) * (1. - z_Jpsi * std::conj(z_psi2S)) / (z_Jpsi - z_psi2S) / dzds;
                }

                // Residue of H at s = m_psi2S2 computed as the residue wrt z -z_psi2S divided by dz/ds evaluated at s = m_psi2S2
                inline complex<double> H_residue_psi2s(const std::array<unsigned, 4> & phi_parameters, const std::array<complex<double>, interpolation_order + 1> & interpolation_values) const
                {
                    const double m_Jpsi2  = power_of<2>(m_Jpsi);
                    const double m_psi2S2 = power_of<2>(m_psi2S);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2, s_p, s_0);

                    const complex<double> p_at_z = lagrange(interpolation_values, z_psi2S);

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_psi2S2, -0.5) * pow(pow(s_p - m_psi2S2, 0.5) + pow(s_p - s_0, 0.5), -2);

                    return p_at_z / phi(m_psi2S2, phi_parameters) *(1 - norm(z_psi2S)) * (1. - z_psi2S * std::conj(z_Jpsi)) / (z_psi2S - z_Jpsi) / dzds;
                }

                virtual complex<double> H_perp(const complex<double> & q2) const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_perp, im_at_m7_perp),
                        complex<double>(re_at_m5_perp, im_at_m5_perp),
                        complex<double>(re_at_m3_perp, im_at_m3_perp),
                        complex<double>(re_at_m1_perp, im_at_m1_perp),
                        polar<double>(abs_at_Jpsi_perp, arg_at_Jpsi_perp_minus_long + arg_at_Jpsi_long),
                        polar<double>(abs_at_psi2S_perp, arg_at_psi2S_perp_minus_long + arg_at_psi2S_long)
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2,                     s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(power_of<2>(m_Jpsi),    s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(power_of<2>(m_psi2S),   s_p, s_0);

                    const complex<double> blaschke_factor = eos::nff_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 3, 0};

                    const complex<double> p_at_z = lagrange(interpolation_values, z);

                    return p_at_z / phi(q2, phi_parameters) / blaschke_factor;
                }

                virtual complex<double> H_perp(const double & q2) const
                {
                    return H_perp(complex<double>(q2, 0.0));
                }

                virtual complex<double> Hhat_perp(const double & q2) const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_perp, im_at_m7_perp),
                        complex<double>(re_at_m5_perp, im_at_m5_perp),
                        complex<double>(re_at_m3_perp, im_at_m3_perp),
                        complex<double>(re_at_m1_perp, im_at_m1_perp),
                        polar<double>(abs_at_Jpsi_perp, arg_at_Jpsi_perp_minus_long + arg_at_Jpsi_long),
                        polar<double>(abs_at_psi2S_perp, arg_at_psi2S_perp_minus_long + arg_at_psi2S_long)
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);

                    return lagrange(interpolation_values, z);
                }

                virtual complex<double> H_para(const complex<double> & q2) const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_para, im_at_m7_para),
                        complex<double>(re_at_m5_para, im_at_m5_para),
                        complex<double>(re_at_m3_para, im_at_m3_para),
                        complex<double>(re_at_m1_para, im_at_m1_para),
                        polar<double>(abs_at_Jpsi_para, arg_at_Jpsi_para_minus_long + arg_at_Jpsi_long),
                        polar<double>(abs_at_psi2S_para, arg_at_psi2S_para_minus_long + arg_at_psi2S_long)
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2,                     s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(power_of<2>(m_Jpsi),    s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(power_of<2>(m_psi2S),   s_p, s_0);

                    const complex<double> blaschke_factor = eos::nff_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 3, 0};

                    const complex<double> p_at_z = lagrange(interpolation_values, z);

                    return p_at_z / phi(q2, phi_parameters) / blaschke_factor;
                }

                virtual complex<double> H_para(const double & q2) const
                {
                    return H_para(complex<double>(q2, 0.0));
                }

                virtual complex<double> Hhat_para(const double & q2) const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_para, im_at_m7_para),
                        complex<double>(re_at_m5_para, im_at_m5_para),
                        complex<double>(re_at_m3_para, im_at_m3_para),
                        complex<double>(re_at_m1_para, im_at_m1_para),
                        polar<double>(abs_at_Jpsi_para, arg_at_Jpsi_para_minus_long + arg_at_Jpsi_long),
                        polar<double>(abs_at_psi2S_para, arg_at_psi2S_para_minus_long + arg_at_psi2S_long)
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);

                    return lagrange(interpolation_values, z);
                }

                virtual complex<double> H_long(const complex<double> & q2) const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_long, im_at_m7_long),
                        complex<double>(re_at_m5_long, im_at_m5_long),
                        complex<double>(re_at_m3_long, im_at_m3_long),
                        complex<double>(re_at_m1_long, im_at_m1_long),
                        polar<double>(abs_at_Jpsi_long, arg_at_Jpsi_long),
                        polar<double>(abs_at_psi2S_long, arg_at_psi2S_long)
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2,                     s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(power_of<2>(m_Jpsi),    s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(power_of<2>(m_psi2S),   s_p, s_0);

                    const complex<double> blaschke_factor = eos::nff_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 2, 2};

                    const complex<double> p_at_z = lagrange(interpolation_values, z);

                    return p_at_z / phi(q2, phi_parameters) / blaschke_factor;
                }

                virtual complex<double> H_long(const double & q2) const
                {
                    return H_long(complex<double>(q2, 0.0));
                }

                virtual complex<double> Hhat_long(const double & q2) const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_long, im_at_m7_long),
                        complex<double>(re_at_m5_long, im_at_m5_long),
                        complex<double>(re_at_m3_long, im_at_m3_long),
                        complex<double>(re_at_m1_long, im_at_m1_long),
                        polar<double>(abs_at_Jpsi_long, arg_at_Jpsi_long),
                        polar<double>(abs_at_psi2S_long, arg_at_psi2S_long)
                    };

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * power_of<2>(m_D0);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);

                    return lagrange(interpolation_values, z);
                }

                virtual complex<double> H_perp_residue_jpsi() const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_perp, im_at_m7_perp),
                        complex<double>(re_at_m5_perp, im_at_m5_perp),
                        complex<double>(re_at_m3_perp, im_at_m3_perp),
                        complex<double>(re_at_m1_perp, im_at_m1_perp),
                        polar<double>(abs_at_Jpsi_perp, arg_at_Jpsi_perp_minus_long + arg_at_Jpsi_long),
                        polar<double>(abs_at_psi2S_perp, arg_at_psi2S_perp_minus_long + arg_at_psi2S_long)
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 3, 0};

                    return H_residue_jpsi(phi_parameters, interpolation_values);
                }

                virtual complex<double> H_perp_residue_psi2s() const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_perp, im_at_m7_perp),
                        complex<double>(re_at_m5_perp, im_at_m5_perp),
                        complex<double>(re_at_m3_perp, im_at_m3_perp),
                        complex<double>(re_at_m1_perp, im_at_m1_perp),
                        polar<double>(abs_at_Jpsi_perp, arg_at_Jpsi_perp_minus_long + arg_at_Jpsi_long),
                        polar<double>(abs_at_psi2S_perp, arg_at_psi2S_perp_minus_long + arg_at_psi2S_long)
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 3, 0};

                    return H_residue_psi2s(phi_parameters, interpolation_values);
                }

                virtual complex<double> H_para_residue_jpsi() const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_para, im_at_m7_para),
                        complex<double>(re_at_m5_para, im_at_m5_para),
                        complex<double>(re_at_m3_para, im_at_m3_para),
                        complex<double>(re_at_m1_para, im_at_m1_para),
                        polar<double>(abs_at_Jpsi_para, arg_at_Jpsi_para_minus_long + arg_at_Jpsi_long),
                        polar<double>(abs_at_psi2S_para, arg_at_psi2S_para_minus_long + arg_at_psi2S_long)
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 3, 0};

                    return H_residue_jpsi(phi_parameters, interpolation_values);
                }

                virtual complex<double> H_para_residue_psi2s() const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_para, im_at_m7_para),
                        complex<double>(re_at_m5_para, im_at_m5_para),
                        complex<double>(re_at_m3_para, im_at_m3_para),
                        complex<double>(re_at_m1_para, im_at_m1_para),
                        polar<double>(abs_at_Jpsi_para, arg_at_Jpsi_para_minus_long + arg_at_Jpsi_long),
                        polar<double>(abs_at_psi2S_para, arg_at_psi2S_para_minus_long + arg_at_psi2S_long)
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 3, 0};

                    return H_residue_psi2s(phi_parameters, interpolation_values);
                }

                virtual complex<double> H_long_residue_jpsi() const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_long, im_at_m7_long),
                        complex<double>(re_at_m5_long, im_at_m5_long),
                        complex<double>(re_at_m3_long, im_at_m3_long),
                        complex<double>(re_at_m1_long, im_at_m1_long),
                        polar<double>(abs_at_Jpsi_long, arg_at_Jpsi_long),
                        polar<double>(abs_at_psi2S_long, arg_at_psi2S_long)
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 2, 2};

                    return H_residue_jpsi(phi_parameters, interpolation_values);
                }

                virtual complex<double> H_long_residue_psi2s() const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_long, im_at_m7_long),
                        complex<double>(re_at_m5_long, im_at_m5_long),
                        complex<double>(re_at_m3_long, im_at_m3_long),
                        complex<double>(re_at_m1_long, im_at_m1_long),
                        polar<double>(abs_at_Jpsi_long, arg_at_Jpsi_long),
                        polar<double>(abs_at_psi2S_long, arg_at_psi2S_long)
                    };

                    const std::array<unsigned, 4> phi_parameters = {3, 1, 2, 2};

                    return H_residue_psi2s(phi_parameters, interpolation_values);
                }

                virtual complex<double> ratio_perp(const complex<double> & q2) const
                {
                    const double m_V2  = power_of<2>(m_V);
                    const double m_B2  = power_of<2>(m_B);
                    const complex<double> lambda = eos::lambda(complex<double>(m_B2, 0.0), complex<double>(m_V2, 0.0), q2);
                    const complex<double> F_perp = pow(2.0 * lambda, 0.5) / (m_B + m_V) / m_B() * form_factors->v(q2);

                    return H_perp(q2) / F_perp;
                }

                virtual complex<double> ratio_perp(const double & q2) const
                {
                    return ratio_perp(complex<double>(q2, 0.0));
                }

                virtual complex<double> F_ratio_perp(const complex<double> & q2) const
                {
                    return (m_B + m_V) / m_B * form_factors->t_1(q2) / form_factors->v(q2);
                }

                virtual complex<double> ratio_para(const complex<double> & q2) const
                {
                    const complex<double> F_para = sqrt(2) * (m_B + m_V) / m_B * form_factors->a_1(q2);

                    return H_para(q2) / F_para;
                }

                virtual complex<double> ratio_para(const double & q2) const
                {
                    return ratio_para(complex<double>(q2, 0.0));
                }

                virtual complex<double> F_ratio_para(const complex<double> & q2) const
                {
                    return (m_B - m_V) / m_B * form_factors->t_2(q2) / form_factors->a_1(q2);
                }

                virtual complex<double> ratio_long(const complex<double> & q2) const
                {
                    const double m_V2  = power_of<2>(m_V);
                    const double m_B2  = power_of<2>(m_B);
                    const complex<double> lambda = eos::lambda(complex<double>(m_B2, 0.0), complex<double>(m_V2, 0.0), q2);
                    const complex<double> F_long = ((m_B2 - m_V2 - q2) * power_of<2>(m_B + m_V) * form_factors->a_1(q2)
                            - lambda * form_factors->a_2(q2)) / (2 * m_V * m_B2 * (m_B + m_V));

                    return H_long(q2) / F_long;
                }

                virtual complex<double> ratio_long(const double & q2) const
                {
                    return ratio_long(complex<double>(q2, 0.0));
                }

                virtual complex<double> F_ratio_long(const complex<double> & q2) const
                {
                    const double m_V2  = power_of<2>(m_V);
                    const double m_B2  = power_of<2>(m_B);
                    const complex<double> lambda = eos::lambda(complex<double>(m_B2, 0.0), complex<double>(m_V2, 0.0), q2);
                    const complex<double> F_long = ((m_B2 - m_V2 - q2) * power_of<2>(m_B + m_V) * form_factors->a_1(q2)
                            - lambda * form_factors->a_2(q2)) / (2 * m_V * m_B2 * (m_B + m_V));
                    const complex<double> F_T_long = q2 * ((m_B2 + 3 * m_V2 - q2) * (m_B2 - m_V2) * form_factors->t_2(q2)
                            - lambda * form_factors->a_2(q2)) / (2 * m_V * m_B2 * m_B * (m_B2 - m_V2));

                    return F_T_long / F_long;
                }

                inline std::pair<gsl_vector *, gsl_vector *> orthonormal_perp_coefficients() const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_perp, im_at_m7_perp),
                        complex<double>(re_at_m5_perp, im_at_m5_perp),
                        complex<double>(re_at_m3_perp, im_at_m3_perp),
                        complex<double>(re_at_m1_perp, im_at_m1_perp),
                        polar<double>(abs_at_Jpsi_perp, arg_at_Jpsi_perp_minus_long + arg_at_Jpsi_long),
                        polar<double>(abs_at_psi2S_perp, arg_at_psi2S_perp_minus_long + arg_at_psi2S_long)
                    };

                    std::array<complex<double>, interpolation_order + 1> L_coeffs = lagrange.get_coefficients(interpolation_values);

                    // Split array of coefficients to real and imaginary parts
                    gsl_vector * L_coeffs_real_part = gsl_vector_calloc(interpolation_order + 1);
                    gsl_vector * L_coeffs_imag_part = gsl_vector_calloc(interpolation_order + 1);

                    for (unsigned i = 0; i <= interpolation_order; ++i)
                    {
                        gsl_vector_set(L_coeffs_real_part, i, real(L_coeffs[i]));
                        gsl_vector_set(L_coeffs_imag_part, i, imag(L_coeffs[i]));
                    }

                    const gsl_matrix * coefficient_matrix = orthonormal_polynomials->coefficient_matrix();

                    // Solve the system by computing (coefficient_matrix)^(-1) . L_coeffs_real_part and idem for imag
                    gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, coefficient_matrix, L_coeffs_real_part);
                    gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, coefficient_matrix, L_coeffs_imag_part);

                    return std::make_pair(L_coeffs_real_part, L_coeffs_imag_part);
                }

                inline std::pair<gsl_vector *, gsl_vector *> orthonormal_para_coefficients() const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_para, im_at_m7_para),
                        complex<double>(re_at_m5_para, im_at_m5_para),
                        complex<double>(re_at_m3_para, im_at_m3_para),
                        complex<double>(re_at_m1_para, im_at_m1_para),
                        polar<double>(abs_at_Jpsi_para, arg_at_Jpsi_para_minus_long + arg_at_Jpsi_long),
                        polar<double>(abs_at_psi2S_para, arg_at_psi2S_para_minus_long + arg_at_psi2S_long)
                    };

                    std::array<complex<double>, interpolation_order + 1> L_coeffs = lagrange.get_coefficients(interpolation_values);

                    // Split array of coefficients to real and imaginary parts
                    gsl_vector * L_coeffs_real_part = gsl_vector_calloc(interpolation_order + 1);
                    gsl_vector * L_coeffs_imag_part = gsl_vector_calloc(interpolation_order + 1);

                    for (unsigned i = 0; i <= interpolation_order; ++i)
                    {
                        gsl_vector_set(L_coeffs_real_part, i, real(L_coeffs[i]));
                        gsl_vector_set(L_coeffs_imag_part, i, imag(L_coeffs[i]));
                    }

                    const gsl_matrix * coefficient_matrix = orthonormal_polynomials->coefficient_matrix();

                    // Solve the system by computing (coefficient_matrix)^(-1) . L_coeffs_real_part and idem for imag
                    gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, coefficient_matrix, L_coeffs_real_part);
                    gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, coefficient_matrix, L_coeffs_imag_part);

                    return std::make_pair(L_coeffs_real_part, L_coeffs_imag_part);
                }

                inline std::pair<gsl_vector *, gsl_vector *> orthonormal_long_coefficients() const
                {
                    const std::array<complex<double>, interpolation_order + 1> interpolation_values{
                        complex<double>(re_at_m7_long, im_at_m7_long),
                        complex<double>(re_at_m5_long, im_at_m5_long),
                        complex<double>(re_at_m3_long, im_at_m3_long),
                        complex<double>(re_at_m1_long, im_at_m1_long),
                        polar<double>(abs_at_Jpsi_long, arg_at_Jpsi_long),
                        polar<double>(abs_at_psi2S_long, arg_at_psi2S_long)
                    };

                    std::array<complex<double>, interpolation_order + 1> L_coeffs = lagrange.get_coefficients(interpolation_values);

                    // Split array of coefficients to real and imaginary parts
                    gsl_vector * L_coeffs_real_part = gsl_vector_calloc(interpolation_order + 1);
                    gsl_vector * L_coeffs_imag_part = gsl_vector_calloc(interpolation_order + 1);

                    for (unsigned i = 0; i <= interpolation_order; ++i)
                    {
                        gsl_vector_set(L_coeffs_real_part, i, real(L_coeffs[i]));
                        gsl_vector_set(L_coeffs_imag_part, i, imag(L_coeffs[i]));
                    }

                    const gsl_matrix * coefficient_matrix = orthonormal_polynomials->coefficient_matrix();

                    // Solve the system by computing (coefficient_matrix)^(-1) . L_coeffs_real_part and idem for imag
                    gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, coefficient_matrix, L_coeffs_real_part);
                    gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, coefficient_matrix, L_coeffs_imag_part);

                    return std::make_pair(L_coeffs_real_part, L_coeffs_imag_part);
                }

                virtual complex<double> get_orthonormal_perp_coefficients(const unsigned & i) const
                {
                    auto coefficients = orthonormal_perp_coefficients();

                    return complex<double>(gsl_vector_get(coefficients.first,  i),
                                           gsl_vector_get(coefficients.second, i));
                }

                virtual complex<double> get_orthonormal_para_coefficients(const unsigned & i) const
                {
                    auto coefficients = orthonormal_para_coefficients();

                    return complex<double>(gsl_vector_get(coefficients.first,  i),
                                           gsl_vector_get(coefficients.second, i));
                }

                virtual complex<double> get_orthonormal_long_coefficients(const unsigned & i) const
                {
                    auto coefficients = orthonormal_long_coefficients();

                    return complex<double>(gsl_vector_get(coefficients.first,  i),
                                           gsl_vector_get(coefficients.second, i));
                }

                virtual double weak_bound() const
                {
                    auto perp_coefficients = orthonormal_perp_coefficients();
                    auto para_coefficients = orthonormal_para_coefficients();
                    auto long_coefficients = orthonormal_long_coefficients();

                    double largest_absolute_coeff = 0.0, coeff;

                    for (unsigned i = 0; i <= interpolation_order; ++i)
                    {
                        coeff = power_of<2>(gsl_vector_get(perp_coefficients.first,  i))
                              + power_of<2>(gsl_vector_get(perp_coefficients.second, i));
                        if (coeff > largest_absolute_coeff)
                        {
                            largest_absolute_coeff = coeff;
                        }
                        coeff = power_of<2>(gsl_vector_get(para_coefficients.first,  i))
                              + power_of<2>(gsl_vector_get(para_coefficients.second, i));
                        if (coeff > largest_absolute_coeff)
                        {
                            largest_absolute_coeff = coeff;
                        }
                        coeff = power_of<2>(gsl_vector_get(long_coefficients.first,  i))
                              + power_of<2>(gsl_vector_get(long_coefficients.second, i));
                        if (coeff > largest_absolute_coeff)
                        {
                            largest_absolute_coeff = coeff;
                        }
                    }

                    return largest_absolute_coeff;
                }

                virtual double strong_bound() const
                {
                    auto perp_coefficients = orthonormal_perp_coefficients();
                    auto para_coefficients = orthonormal_para_coefficients();
                    auto long_coefficients = orthonormal_long_coefficients();

                    double coefficient_sum = 0.0;

                    for (unsigned i = 0; i <= interpolation_order; ++i)
                    {
                        coefficient_sum += power_of<2>(gsl_vector_get(perp_coefficients.first,  i))
                                          + power_of<2>(gsl_vector_get(perp_coefficients.second, i));
                        coefficient_sum += power_of<2>(gsl_vector_get(para_coefficients.first,  i))
                                          + power_of<2>(gsl_vector_get(para_coefficients.second, i));
                        coefficient_sum += power_of<2>(gsl_vector_get(long_coefficients.first,  i))
                                          + power_of<2>(gsl_vector_get(long_coefficients.second, i));
                    }

                    return coefficient_sum;
                }

                virtual double weak_bound_log_likelihood() const
                {
                    const double saturation = weak_bound();
                    if (saturation < this->bound)
                    {
                        return 0.;
                    }
                    else
                    {
                        // Halfnormal constraint
                        return -0.5 * power_of<2>( (saturation - this->bound) / this->bound_uncertainty );
                    }
                }

                virtual double strong_bound_log_likelihood() const
                {
                    const double saturation = strong_bound();
                    if (saturation < this->bound)
                    {
                        return 0.;
                    }
                    else
                    {
                        // Halfnormal constraint
                        return -0.5 * power_of<2>( (saturation - this->bound) / this->bound_uncertainty );
                    }
                }

                virtual complex<double> normalized_moment_V1(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V2(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_V23(const double &) const
                {
                    return 0.0;
                }


                static NonlocalFormFactorPtr<PToV> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<PToV>(new GRvDV2022order5<Process_>(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    Diagnostics results;

                    return results;
                }
        };
    }

    NonlocalFormFactorPtr<PToV>
    NonlocalFormFactor<PToV>::make(const QualifiedName & name, const Parameters & p, const Options & o)
    {
        typedef QualifiedName KeyType;
        typedef std::function<NonlocalFormFactorPtr<PToV> (const Parameters &, const Options &)> ValueType;
        std::map<KeyType, ValueType> entries
        {
            // trivial
            std::make_pair("B->K^*::naive",             &nff_p_to_v::Naive::make),
            // parametrizations
            std::make_pair("B->K^*::GvDV2020",          &nff_p_to_v::GvDV2020<nff::BToKstar>::make),
            std::make_pair("B->K^*::GRvDV2022order5",   &nff_p_to_v::GRvDV2022order5<nff::BToKstar>::make),
            std::make_pair("B_s->phi::GvDV2020",        &nff_p_to_v::GvDV2020<nff::BsToPhi>::make),
            std::make_pair("B_s->phi::GRvDV2022order5", &nff_p_to_v::GRvDV2022order5<nff::BsToPhi>::make),
        };

        auto i = entries.find(name);
        if (entries.end() == i)
        {
            return NonlocalFormFactorPtr<PToV>(nullptr);
        }

        return i->second(p, o);
    }

    template <typename Process_>
    struct Implementation<NonlocalFormFactorObservable<Process_, PToV>>
    {
        NameOption opt_formfactor;
        NonlocalFormFactorPtr<PToV> nff;
        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            opt_formfactor(o, "nonlocal-formfactor"_ok, qnp::Name("GvDV2020")),
            nff(NonlocalFormFactor<PToV>::make(QualifiedName(qnp::Prefix(Process_::label), opt_formfactor.value()), p, o))
        {
            u.uses(*nff);
        }

        ~Implementation() = default;
    };

    template <typename Process_>
    NonlocalFormFactorObservable<Process_, PToV>::NonlocalFormFactorObservable(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<NonlocalFormFactorObservable<Process_, PToV>>(new Implementation<NonlocalFormFactorObservable<Process_, PToV>>(p, o, *this))
    {
    }

    template <typename Process_>
    NonlocalFormFactorObservable<Process_, PToV>::~NonlocalFormFactorObservable() = default;

    template <typename Process_>
    const std::vector<OptionSpecification>
    Implementation<NonlocalFormFactorObservable<Process_, PToV>>::options
    {
    };

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::re_H_perp(const double & q2) const
    {
        return real(this->_imp->nff->H_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::im_H_perp(const double & q2) const
    {
        return imag(this->_imp->nff->H_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::abs_H_perp(const double & q2) const
    {
        return abs(this->_imp->nff->H_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::re_Hhat_perp(const double & q2) const
    {
        return real(this->_imp->nff->Hhat_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::im_Hhat_perp(const double & q2) const
    {
        return imag(this->_imp->nff->Hhat_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::abs_Hhat_perp(const double & q2) const
    {
        return abs(this->_imp->nff->Hhat_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::re_H_para(const double & q2) const
    {
        return real(this->_imp->nff->H_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::im_H_para(const double & q2) const
    {
        return imag(this->_imp->nff->H_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::abs_H_para(const double & q2) const
    {
        return abs(this->_imp->nff->H_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::re_Hhat_para(const double & q2) const
    {
        return real(this->_imp->nff->Hhat_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::im_Hhat_para(const double & q2) const
    {
        return imag(this->_imp->nff->Hhat_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::abs_Hhat_para(const double & q2) const
    {
        return abs(this->_imp->nff->Hhat_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::re_H_long(const double & q2) const
    {
        return real(this->_imp->nff->H_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::im_H_long(const double & q2) const
    {
        return imag(this->_imp->nff->H_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::abs_H_long(const double & q2) const
    {
        return abs(this->_imp->nff->H_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::re_Hhat_long(const double & q2) const
    {
        return real(this->_imp->nff->Hhat_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::im_Hhat_long(const double & q2) const
    {
        return imag(this->_imp->nff->Hhat_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::abs_Hhat_long(const double & q2) const
    {
        return abs(this->_imp->nff->Hhat_long(q2));
    }


    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::re_ratio_perp(const double & q2) const
    {
        return real(this->_imp->nff->ratio_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::im_ratio_perp(const double & q2) const
    {
        return imag(this->_imp->nff->ratio_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::abs_ratio_perp(const double & q2) const
    {
        return abs(this->_imp->nff->ratio_perp(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::re_ratio_para(const double & q2) const
    {
        return real(this->_imp->nff->ratio_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::im_ratio_para(const double & q2) const
    {
        return imag(this->_imp->nff->ratio_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::abs_ratio_para(const double & q2) const
    {
        return abs(this->_imp->nff->ratio_para(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::re_ratio_long(const double & q2) const
    {
        return real(this->_imp->nff->ratio_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::im_ratio_long(const double & q2) const
    {
        return imag(this->_imp->nff->ratio_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::abs_ratio_long(const double & q2) const
    {
        return abs(this->_imp->nff->ratio_long(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::re_ratio_perp_complex(const double & re_q2, const double & im_q2) const
    {
        return real(this->_imp->nff->ratio_perp(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::im_ratio_perp_complex(const double & re_q2, const double & im_q2) const
    {
        return imag(this->_imp->nff->ratio_perp(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::re_ratio_para_complex(const double & re_q2, const double & im_q2) const
    {
        return real(this->_imp->nff->ratio_para(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::im_ratio_para_complex(const double & re_q2, const double & im_q2) const
    {
        return imag(this->_imp->nff->ratio_para(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::re_ratio_long_complex(const double & re_q2, const double & im_q2) const
    {
        return real(this->_imp->nff->ratio_long(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::im_ratio_long_complex(const double & re_q2, const double & im_q2) const
    {
        return imag(this->_imp->nff->ratio_long(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::re_F_ratio_perp_complex(const double & re_q2, const double & im_q2) const
    {
        return real(this->_imp->nff->F_ratio_perp(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::im_F_ratio_perp_complex(const double & re_q2, const double & im_q2) const
    {
        return imag(this->_imp->nff->F_ratio_perp(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::re_F_ratio_para_complex(const double & re_q2, const double & im_q2) const
    {
        return real(this->_imp->nff->F_ratio_para(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::im_F_ratio_para_complex(const double & re_q2, const double & im_q2) const
    {
        return imag(this->_imp->nff->F_ratio_para(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::re_F_ratio_long_complex(const double & re_q2, const double & im_q2) const
    {
        return real(this->_imp->nff->F_ratio_long(complex<double>(re_q2, im_q2)));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::im_F_ratio_long_complex(const double & re_q2, const double & im_q2) const
    {
        return imag(this->_imp->nff->F_ratio_long(complex<double>(re_q2, im_q2)));
    }


    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::re_normalized_moment_V1(const double & q2) const
    {
        return real(this->_imp->nff->normalized_moment_V1(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::re_normalized_moment_V2(const double & q2) const
    {
        return real(this->_imp->nff->normalized_moment_V2(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::re_normalized_moment_V23(const double & q2) const
    {
        return real(this->_imp->nff->normalized_moment_V23(q2));
    }



    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::get_real_perp_alpha(const unsigned & i) const
    {
        return real(this->_imp->nff->get_orthonormal_perp_coefficients(i));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::get_imag_perp_alpha(const unsigned & i) const
    {
        return imag(this->_imp->nff->get_orthonormal_perp_coefficients(i));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::get_real_para_alpha(const unsigned & i) const
    {
        return real(this->_imp->nff->get_orthonormal_para_coefficients(i));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::get_imag_para_alpha(const unsigned & i) const
    {
        return imag(this->_imp->nff->get_orthonormal_para_coefficients(i));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::get_real_long_alpha(const unsigned & i) const
    {
        return real(this->_imp->nff->get_orthonormal_long_coefficients(i));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::get_imag_long_alpha(const unsigned & i) const
    {
        return imag(this->_imp->nff->get_orthonormal_long_coefficients(i));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::weak_bound() const
    {
        return this->_imp->nff->weak_bound();
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::weak_bound_log_likelihood() const
    {
        return this->_imp->nff->weak_bound_log_likelihood();
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::strong_bound() const
    {
        return this->_imp->nff->strong_bound();
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, PToV>::strong_bound_log_likelihood() const
    {
        return this->_imp->nff->strong_bound_log_likelihood();
    }


    template <typename Process_>
    const std::set<ReferenceName>
    NonlocalFormFactorObservable<Process_, PToV>::references
    {
        "GvDV:2020A"_rn,
        "GRvDV:2022A"_rn,
    };

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    NonlocalFormFactorObservable<Process_, PToV>::begin_options()
    {
        return Implementation<NonlocalFormFactorObservable<Process_, PToV>>::options.cbegin();
    }

    template <typename Process_>
    std::vector<OptionSpecification>::const_iterator
    NonlocalFormFactorObservable<Process_, PToV>::end_options()
    {
        return Implementation<NonlocalFormFactorObservable<Process_, PToV>>::options.cend();
    }

    template class NonlocalFormFactorObservable<nff::BToKstar, PToV>;

    template class NonlocalFormFactorObservable<nff::BsToPhi, PToV>;
}
