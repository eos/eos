/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017-2020 Danny van Dyk
 * Copyright (c) 2020 Nico Gubernari
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

#include <eos/rare-b-decays/nonlocal-formfactors.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/model.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/stringify.hh>

#include <map>

namespace eos
{
    using std::abs;

    namespace nff
    {
        struct BToK
        {
            constexpr static const char * label = "B->K";
        };
        constexpr const char * BToK::label;
    }

    namespace nff_p_to_p
    {
        class Naive :
            public NonlocalFormFactor<nff::PToP>
        {
            public:
                Naive(const Parameters &, const Options &)
                {
                }

                ~Naive() = default;

                virtual complex<double> H_plus(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> Hhat_plus(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_plus(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> normalized_moment_A(const double &) const
                {
                    return 0.0;
                }

                static NonlocalFormFactorPtr<nff::PToP> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nff::PToP>(new Naive(p, o));
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
            public NonlocalFormFactor<nff::PToP>
        {
            public:
                std::shared_ptr<FormFactors<PToP>> form_factors;

                //Polynomial expansion parameters
                UsedParameter re_alpha_0_plus;
                UsedParameter im_alpha_0_plus;
                UsedParameter re_alpha_1_plus;
                UsedParameter im_alpha_1_plus;
                UsedParameter re_alpha_2_plus;
                UsedParameter im_alpha_2_plus;

                //Charmonium masses
                UsedParameter m_Jpsi;
                UsedParameter m_psi2S;

                // B-meson parameters
                UsedParameter m_B;

                // final state meson parameters
                UsedParameter m_P;

                UsedParameter m_D0;
                UsedParameter t_0;

                // Subtraction point for the dispersion relation...
                UsedParameter t_s;
                // ...and value of the dispersion bound at that point in the OPE
                UsedParameter chiOPE;

                GvDV2020(const Parameters & p, const Options & o) :
                    form_factors(FormFactorFactory<PToP>::create(stringify(Process_::label) + "::" + o.get("form-factors", "BSZ2015"), p)),

                    re_alpha_0_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_0^plus}@GvDV2020"], *this),
                    im_alpha_0_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_0^plus}@GvDV2020"], *this),
                    re_alpha_1_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_1^plus}@GvDV2020"], *this),
                    im_alpha_1_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_1^plus}@GvDV2020"], *this),
                    re_alpha_2_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_2^plus}@GvDV2020"], *this),
                    im_alpha_2_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_2^plus}@GvDV2020"], *this),

                    m_Jpsi(p["mass::J/psi"], *this),
                    m_psi2S(p["mass::psi(2S)"], *this),

                    m_B(p["mass::B_d"], *this),

                    m_P(p["mass::K_d"], *this),

                    m_D0(p["mass::D^0"], *this),
                    t_0(p["b->sccbar::t_0"], *this),

                    t_s(p["b->sccbar::t_s"], *this),
                    chiOPE(p["b->sccbar::chiOPE@GvDV2020"], *this)
                {
                    this->uses(*form_factors);
                }

                ~GvDV2020() = default;

                inline complex<double> phi(const double & q2, const unsigned phiParam[4]) const
                {
                    // Values of a, b, c and d depends on the form factor:
                    // FF                        a    b    c    d
                    // 0(P->P) aka plus          3    3    2    2
                    // perp(P->V) = par(P->V)    3    1    3    0
                    // 0(P->V) aka long          3    1    2    2

                    const double m_P2  = pow(m_P, 2);
                    const double m_B2  = pow(m_B, 2),  m_B4 =  pow(m_B, 4);
                    const double m_D02 = pow(m_D0, 2), m_D04 = pow(m_D0, 4);
                    const double s_0   = this->t_0();
                    const auto   z     = eos::nff_utils::z(q2, 4.0 * pow(m_D0, 2), s_0);
                    const double Q2   = this->t_s();
                    const double chi = this->chiOPE();

                    const double a = phiParam[0], b = phiParam[1], c = phiParam[2], d = phiParam[3];

                    const double Nlambda = 4 * M_PI * pow(m_B2, 0.5 * (a - b + c + d) - 1.) * pow(2 * (4 * m_D02 - s_0) / 3 / chi, 0.5); //(C6)
                    const complex<double> phi1 = -pow(2 * pow((4 * m_D02 - Q2) * (4 * m_D02 - s_0), 0.5) + 8 * m_D02 - Q2 - s_0, 0.5) /
                                                (2 * pow((4 * m_D02 - Q2) * (4 * m_D02 - s_0), 0.5) + 8 * m_D02 + Q2 * (z - 1.) - s_0 * (z + 1.)); //(C7)
                    const complex<double> phi2 = pow(m_B4 * pow(z - 1., 4) - 2 * m_B2 * pow(z - 1., 2) * (-16 * m_D02 * z + m_P2 * pow(z - 1., 2) + s_0 * pow(z + 1., 2)) +
                                                pow(16 * m_D02 * z + m_P2 * pow(z - 1., 2) - s_0 * pow(z + 1., 2), 2), 0.5);//(C8)
                    const complex<double> phi3 = pow(8 * m_D02 + 4 * pow(4 * m_D04 - s_0 * m_D02, 0.5) - s_0, 0.5) /
                                                (-8 * m_D02 - 4 * pow(4 * m_D04 - s_0 * m_D02, 0.5) + s_0*(z+1.)); //(C9)
                    const complex<double> phi4 = pow(s_0 * pow(z + 1., 2.) - 16. * z * m_D02, -0.5); //(C10)

                    return Nlambda * pow(1.+z, 0.5) * pow(1.-z, a-b+c+d-1.5) * pow(phi1, a) * pow(phi2, 0.5*b) * pow(phi3, c) * pow(phi4, d); //(C5)
                }

                // Residue of H at s = m_Jpsi2 computed as the residue wrt z -z_Jpsi divided by dz/ds evaluated at s = m_Jpsi2
                inline complex<double> H_residue_jpsi(const unsigned phiParam[4], const complex<double> & alpha_0, const complex<double> & alpha_1,
                                                      const complex<double> & alpha_2) const
                {
                    const double m_Jpsi2  = pow(m_Jpsi, 2);
                    const double m_psi2S2 = pow(m_psi2S, 2);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto zBP     = eos::nff_utils::z(pow(m_B + m_P, 2), s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,           s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2,          s_p, s_0);

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_Jpsi2, -0.5) * pow(pow(s_p - m_Jpsi2, 0.5) + pow(s_p - s_0, 0.5), -2);

                    return eos::nff_utils::PGvDV2020(z_Jpsi, zBP, alpha_0, alpha_1, alpha_2) / phi(m_Jpsi2, phiParam) *
                            (1 - norm(z_Jpsi)) * (1. - z_Jpsi * std::conj(z_psi2S)) / (z_Jpsi - z_psi2S) / dzds;
                }

                // Residue of H at s = m_psi2S2 computed as the residue wrt z -z_psi2S divided by dz/ds evaluated at s = m_psi2S2
                inline complex<double> H_residue_psi2s(const unsigned phiParam[4], const complex<double> & alpha_0, const complex<double> & alpha_1,
                                                       const complex<double> & alpha_2) const
                {
                    const double m_Jpsi2  = pow(m_Jpsi, 2);
                    const double m_psi2S2 = pow(m_psi2S, 2);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto zBP     = eos::nff_utils::z(pow(m_B + m_P, 2), s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,           s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2,          s_p, s_0);

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_psi2S2, -0.5) * pow(pow(s_p - m_psi2S2, 0.5) + pow(s_p - s_0, 0.5), -2);

                    return eos::nff_utils::PGvDV2020(z_psi2S, zBP, alpha_0, alpha_1, alpha_2) / phi(m_psi2S2, phiParam) *
                            (1 - norm(z_psi2S)) * (1. - z_psi2S * std::conj(z_Jpsi)) / (z_psi2S - z_Jpsi) / dzds;
                }


                virtual complex<double> H_plus(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_plus, im_alpha_0_plus);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_plus, im_alpha_1_plus);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_plus, im_alpha_2_plus);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nff_utils::z(q2,              s_p, s_0);
                    const auto zBP     = eos::nff_utils::z(pow(m_B+m_P, 2), s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(pow(m_Jpsi, 2),  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(pow(m_psi2S, 2), s_p, s_0);

                    const complex<double> blaschke_factor = eos::nff_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const unsigned phiParam[4] = {3, 3, 2, 2};

                    return eos::nff_utils::PGvDV2020(z, zBP, alpha_0, alpha_1, alpha_2) / phi(q2, phiParam) / blaschke_factor;
                }


                virtual complex<double> Hhat_plus(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_plus, im_alpha_0_plus);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_plus, im_alpha_1_plus);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_plus, im_alpha_2_plus);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);
                    const auto zBP     = eos::nff_utils::z(pow(m_B + m_P, 2), s_p, s_0);

                    return eos::nff_utils::PGvDV2020(z, zBP, alpha_0, alpha_1, alpha_2);
                }


                virtual complex<double> H_plus_residue_jpsi() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_plus, im_alpha_0_plus);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_plus, im_alpha_1_plus);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_plus, im_alpha_2_plus);

                    const unsigned phiParam[4] = {3, 3, 2, 2};

                    return H_residue_jpsi(phiParam, alpha_0, alpha_1, alpha_2);
                }

                virtual complex<double> H_plus_residue_psi2s() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_plus, im_alpha_0_plus);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_plus, im_alpha_1_plus);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_plus, im_alpha_2_plus);

                    const unsigned phiParam[4] = {3, 3, 2, 2};

                    return H_residue_psi2s(phiParam, alpha_0, alpha_1, alpha_2);
                }

                virtual complex<double> normalized_moment_A(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_plus(const double & q2) const
                {
                    const double F_plus = form_factors->f_p(q2);

                    return H_plus(q2) / F_plus;
                }

                static NonlocalFormFactorPtr<nff::PToP> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nff::PToP>(new GvDV2020<Process_>(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    Diagnostics results;

                    const unsigned phiParam[4] = {3, 3, 2, 2}; //plus polarization

                    results.add({ real(1./this->phi(0.0, phiParam)), "Re{1/phi_+(q2 = 0.0)}" });
                    results.add({ imag(1./this->phi(0.0, phiParam)), "Im{1/phi_+(q2 = 0.0)}" });
                    results.add({ real(this->phi(16.0, phiParam)), "Re{phi_+(q2 = 16.0)}" });
                    results.add({ imag(this->phi(16.0, phiParam)), "Im{phi_+(q2 = 16.0)}" });

                    const double s_0   = this->t_0();
                    const auto z1 = eos::nff_utils::z(1.0, 4.0 * pow(m_D0, 2), s_0);

                    results.add({ real(eos::nff_utils::PGvDV2020(z1, complex<double>(0.6,0.8), 2.0, 3.0, 4.0)),
                                "Re{PGvDV2020(q2 = 1.0, sXY = 0.6+0.8i, 2.0, 3.0, 4.0)}" });
                    results.add({ imag(eos::nff_utils::PGvDV2020(z1, complex<double>(0.6,0.8), 2.0, 3.0, 4.0)),
                                "Im{PGvDV2020(q2 = 1.0, sXY = 0.6+0.8i, 2.0, 3.0, 4.0)}" });

                    return results;
                }
        };



        /*
         * Parametrize the entire formfactor, i.e., both leading and all sub-leading powers as described in [GRvDV:2021].
         */
        template <typename Process_>
        class GRvDV2021 :
            public NonlocalFormFactor<nff::PToP>
        {
            public:
                std::shared_ptr<FormFactors<PToP>> form_factors;

                //Polynomial expansion parameters
                UsedParameter re_alpha_0_plus;
                UsedParameter im_alpha_0_plus;
                UsedParameter re_alpha_1_plus;
                UsedParameter im_alpha_1_plus;
                UsedParameter re_alpha_2_plus;
                UsedParameter im_alpha_2_plus;

                //Charmonium masses
                UsedParameter m_Jpsi;
                UsedParameter m_psi2S;

                // B-meson parameters
                UsedParameter m_B;
                UsedParameter m_Bsst;

                // final state meson parameters
                UsedParameter m_P;

                UsedParameter m_D0;
                UsedParameter t_0;

                // Subtraction point for the dispersion relation...
                UsedParameter t_s;
                // ...and value of the dispersion bound at that point in the OPE
                UsedParameter chiOPE;

                GRvDV2021(const Parameters & p, const Options & o) :
                    form_factors(FormFactorFactory<PToP>::create(stringify(Process_::label) + "::" + o.get("form-factors", "BSZ2015"), p)),

                    re_alpha_0_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_0^plus}@GRvDV2021"], *this),
                    im_alpha_0_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_0^plus}@GRvDV2021"], *this),
                    re_alpha_1_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_1^plus}@GRvDV2021"], *this),
                    im_alpha_1_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_1^plus}@GRvDV2021"], *this),
                    re_alpha_2_plus(p[stringify(Process_::label) + "ccbar::Re{alpha_2^plus}@GRvDV2021"], *this),
                    im_alpha_2_plus(p[stringify(Process_::label) + "ccbar::Im{alpha_2^plus}@GRvDV2021"], *this),

                    m_Jpsi(p["mass::J/psi"], *this),
                    m_psi2S(p["mass::psi(2S)"], *this),

                    m_B(p["mass::B_d"], *this),
                    m_Bsst(p["mass::B_s^*"], *this),

                    m_P(p["mass::K_d"], *this),

                    m_D0(p["mass::D^0"], *this),
                    t_0(p["b->sccbar::t_0"], *this),

                    t_s(p["b->sccbar::t_s"], *this),
                    chiOPE(p["b->sccbar::chiOPE@GRvDV2021"], *this)
                {
                    this->uses(*form_factors);
                }

                ~GRvDV2021() = default;

                inline complex<double> phi(const double & q2, const unsigned phiParam[5]) const
                {
                    // Values of a, b, c and d depends on the form factor:
                    // FF                        a    b    c    d    e
                    // 0(P->P) aka plus          5    3    2    2    2
                    // perp(P->V) = par(P->V)    5    1    3    0    2
                    // 0(P->V) aka long          5    1    2    2    2

                    const double m_P2  = pow(m_P, 2);
                    const double m_Bsst2 = pow(m_Bsst, 2);
                    const double m_B2  = pow(m_B, 2),  m_B4 =  pow(m_B, 4);
                    const double m_D02 = pow(m_D0, 2), m_D04 = pow(m_D0, 4);
                    const double s_0   = this->t_0();
                    const auto   z     = eos::nff_utils::z(q2, 4.0 * pow(m_D0, 2), s_0);
                    const double Q2   = this->t_s();
                    const double chi = this->chiOPE();

                    const double a = phiParam[0], b = phiParam[1], c = phiParam[2], d = phiParam[3], e = phiParam[4];

                    const complex<double> Nlambda = 4. * M_PI * pow(m_B2, 0.5 * (a - b + c + d - e) - 1.) * pow(2. * (4. * m_D02 - s_0) / 3. / chi, 0.5);
                    const complex<double> phi1 = -pow(2. * pow((4. * m_D02 - Q2) * (4. * m_D02 - s_0), 0.5) + 8. * m_D02 - Q2 - s_0, 0.5) /
                                                (2. * pow((4. * m_D02 - Q2) * (4. * m_D02 - s_0), 0.5) + 8. * m_D02 + Q2 * (z - 1.) - s_0*(z + 1.));
                    const complex<double> phi2 = pow(m_B4 * pow(z - 1., 4.) - 2. * m_B2 * pow(z - 1., 2) * (-16 * m_D02 * z + m_P2 * pow(z - 1., 2) +
                                                s_0 * pow(z + 1., 2)) + pow(16 * m_D02 * z + m_P2 * pow(z - 1., 2) - s_0 * pow(z + 1., 2), 2), 0.5);
                    const complex<double> phi3 = pow(8 * m_D02 + 4 * pow(4 * m_D04 - s_0 * m_D02, 0.5) - s_0, 0.5) /
                                                (-8 * m_D02 - 4 * pow(4 * m_D04 - s_0 * m_D02, 0.5) + s_0 * (z + 1.));
                    const complex<double> phi4 = pow(s_0 * pow(z + 1., 2.) - 16. * z * m_D02, -0.5);
                    const complex<double> phi5 = pow(s_0 * pow(z + 1., 2.) - 16. * z * m_D02 - m_Bsst2 * pow(-z + 1., 2.), 0.5);

                    return Nlambda * pow(1. + z, 0.5) * pow(1. - z, a - b + c + d - e - 1.5) * pow(phi1, a) * pow(phi2, 0.5*b) * pow(phi3, c) * pow(phi4, d) * pow(phi5, e);
                }

                // Residue of H at s = m_Jpsi2 computed as the residue wrt z -z_Jpsi divided by dz/ds evaluated at s = m_Jpsi2
                inline complex<double> H_residue_jpsi(const unsigned phiParam[5], const complex<double> & alpha_0, const complex<double> & alpha_1,
                                                      const complex<double> & alpha_2) const
                {
                    const double m_Jpsi2  = pow(m_Jpsi, 2);
                    const double m_psi2S2 = pow(m_psi2S, 2);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2, s_p, s_0);

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_Jpsi2, -0.5) * pow(pow(s_p - m_Jpsi2, 0.5) + pow(s_p - s_0, 0.5), -2);

                    return eos::nff_utils::P(z_Jpsi, alpha_0, alpha_1, alpha_2) / phi(m_Jpsi2, phiParam) *
                            (1 - norm(z_Jpsi)) * (1. - z_Jpsi * std::conj(z_psi2S)) / (z_Jpsi - z_psi2S) / dzds;
                }

                // Residue of H at s = m_psi2S2 computed as the residue wrt z -z_psi2S divided by dz/ds evaluated at s = m_psi2S2
                inline complex<double> H_residue_psi2s(const unsigned phiParam[5], const complex<double> & alpha_0, const complex<double> & alpha_1,
                                                       const complex<double> & alpha_2) const
                {
                    const double m_Jpsi2  = pow(m_Jpsi, 2);
                    const double m_psi2S2 = pow(m_psi2S, 2);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z_Jpsi  = eos::nff_utils::z(m_Jpsi2,  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(m_psi2S2, s_p, s_0);

                    const complex<double> dzds = -pow(s_p - s_0, 0.5) * pow(s_p - m_psi2S2, -0.5) * pow(pow(s_p - m_psi2S2, 0.5) + pow(s_p - s_0, 0.5), -2);

                    return eos::nff_utils::P(z_psi2S, alpha_0, alpha_1, alpha_2) / phi(m_psi2S2, phiParam) *
                            (1 - norm(z_psi2S)) * (1. - z_psi2S * std::conj(z_Jpsi)) / (z_psi2S - z_Jpsi) / dzds;
                }


                virtual complex<double> H_plus(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_plus, im_alpha_0_plus);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_plus, im_alpha_1_plus);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_plus, im_alpha_2_plus);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nff_utils::z(q2,              s_p, s_0);
                    const auto z_Jpsi  = eos::nff_utils::z(pow(m_Jpsi, 2),  s_p, s_0);
                    const auto z_psi2S = eos::nff_utils::z(pow(m_psi2S, 2), s_p, s_0);

                    const complex<double> blaschke_factor = eos::nff_utils::blaschke_cc(z, z_Jpsi, z_psi2S);

                    const unsigned phiParam[5] = {5, 3, 2, 2, 2};

                    return eos::nff_utils::P(z, alpha_0, alpha_1, alpha_2) / phi(q2, phiParam) / blaschke_factor;
                }


                virtual complex<double> Hhat_plus(const double & q2) const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_plus, im_alpha_0_plus);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_plus, im_alpha_1_plus);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_plus, im_alpha_2_plus);

                    const double s_0   = this->t_0();
                    const double s_p   = 4.0 * pow(m_D0, 2);
                    const auto z       = eos::nff_utils::z(q2, s_p, s_0);

                    return eos::nff_utils::P(z, alpha_0, alpha_1, alpha_2);
                }


                virtual complex<double> H_plus_residue_jpsi() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_plus, im_alpha_0_plus);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_plus, im_alpha_1_plus);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_plus, im_alpha_2_plus);

                    const unsigned phiParam[5] = {5, 3, 2, 2, 2};

                    return H_residue_jpsi(phiParam, alpha_0, alpha_1, alpha_2);
                }

                virtual complex<double> H_plus_residue_psi2s() const
                {
                    const complex<double> alpha_0 = complex<double>(re_alpha_0_plus, im_alpha_0_plus);
                    const complex<double> alpha_1 = complex<double>(re_alpha_1_plus, im_alpha_1_plus);
                    const complex<double> alpha_2 = complex<double>(re_alpha_2_plus, im_alpha_2_plus);

                    const unsigned phiParam[5] = {5, 3, 2, 2, 2};

                    return H_residue_psi2s(phiParam, alpha_0, alpha_1, alpha_2);
                };

                virtual complex<double> normalized_moment_A(const double &) const
                {
                    return 0.0;
                }

                virtual complex<double> ratio_plus(const double & q2) const
                {
                    const double F_plus = form_factors->f_p(q2);

                    return H_plus(q2) / F_plus;
                }

                static NonlocalFormFactorPtr<nff::PToP> make(const Parameters & p, const Options & o)
                {
                    return NonlocalFormFactorPtr<nff::PToP>(new GRvDV2021<Process_>(p, o));
                }

                virtual Diagnostics diagnostics() const
                {
                    Diagnostics results;

                    const unsigned phiParam[5] = {5, 3, 2, 2, 2}; //plus polarization

                    results.add({ real(this->phi(16.0, phiParam)), "Re{phi_+(q2 = 16.0)}" });
                    results.add({ imag(this->phi(16.0, phiParam)), "Im{phi_+(q2 = 16.0)}" });

                    const double s_0   = this->t_0();
                    const auto z1 = eos::nff_utils::z(1.0, 4.0 * pow(m_D0, 2), s_0);

                    results.add({ real(eos::nff_utils::P(z1, 2.0, 3.0, 4.0)), "Re{P(q2 = 1.0, 2.0, 3.0, 4.0)}" });
                    results.add({ imag(eos::nff_utils::P(z1, 2.0, 3.0, 4.0)), "Im{P(q2 = 1.0, 2.0, 3.0, 4.0)}" });
                    results.add({ real(eos::nff_utils::P(z1, complex<double>(2.0,5.0), complex<double>(3.0,6.0), complex<double>(4.0,7.0))),
                                "Re{P(q2 = 1.0, (2.0,5.0), (3.0,6.0), (4.0,7.0))}" });
                    results.add({ imag(eos::nff_utils::P(z1, complex<double>(2.0,5.0), complex<double>(3.0,6.0), complex<double>(4.0,7.0))),
                                "Im{P(q2 = 1.0, (2.0,5.0), (3.0,6.0), (4.0,7.0))}" });

                    return results;
                }
        };

    }

    NonlocalFormFactorPtr<nff::PToP>
    NonlocalFormFactor<nff::PToP>::make(const QualifiedName & name, const Parameters & p, const Options & o)
    {
        typedef QualifiedName KeyType;
        typedef std::function<NonlocalFormFactorPtr<nff::PToP> (const Parameters &, const Options &)> ValueType;
        std::map<KeyType, ValueType> entries
        {
            // trivial
            std::make_pair("B->K::naive",         &nff_p_to_p::Naive::make),
            // parametrizations
            std::make_pair("B->K::GvDV2020",      &nff_p_to_p::GvDV2020<nff::BToK>::make),
            std::make_pair("B->K::GRvDV2021",     &nff_p_to_p::GRvDV2021<nff::BToK>::make),
        };

        auto i = entries.find(name);
        if (entries.end() == i)
        {
            return NonlocalFormFactorPtr<nff::PToP>(nullptr);
        }

        return i->second(p, o);
    }

    template <typename Process_>
    struct Implementation<NonlocalFormFactorObservable<Process_, nff::PToP>>
    {
        NameOption opt_formfactor;
        NonlocalFormFactorPtr<nff::PToP> nff;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            opt_formfactor(o, "formfactor", qnp::Name("GvDV2020")),
            nff(NonlocalFormFactor<nff::PToP>::make(QualifiedName(qnp::Prefix(Process_::label), opt_formfactor.value()), p, o))
        {
            u.uses(*nff);
        }

        ~Implementation() = default;
    };

    template <typename Process_>
    NonlocalFormFactorObservable<Process_, nff::PToP>::NonlocalFormFactorObservable(const Parameters & p, const Options & o) :
        PrivateImplementationPattern<NonlocalFormFactorObservable<Process_, nff::PToP>>(new Implementation<NonlocalFormFactorObservable<Process_, nff::PToP>>(p, o, *this))
    {
    }

    template <typename Process_>
    NonlocalFormFactorObservable<Process_, nff::PToP>::~NonlocalFormFactorObservable() = default;

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::re_H_plus(const double & q2) const
    {
        return real(this->_imp->nff->H_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::im_H_plus(const double & q2) const
    {
        return imag(this->_imp->nff->H_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::abs_H_plus(const double & q2) const
    {
        return abs(this->_imp->nff->H_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::re_Hhat_plus(const double & q2) const
    {
        return real(this->_imp->nff->Hhat_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::im_Hhat_plus(const double & q2) const
    {
        return imag(this->_imp->nff->Hhat_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::abs_Hhat_plus(const double & q2) const
    {
        return abs(this->_imp->nff->Hhat_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::re_ratio_plus(const double & q2) const
    {
        return real(this->_imp->nff->ratio_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::im_ratio_plus(const double & q2) const
    {
        return imag(this->_imp->nff->ratio_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::abs_ratio_plus(const double & q2) const
    {
        return abs(this->_imp->nff->ratio_plus(q2));
    }

    template <typename Process_>
    double
    NonlocalFormFactorObservable<Process_, nff::PToP>::re_normalized_moment_A(const double & q2) const
    {
        return real(this->_imp->nff->normalized_moment_A(q2));
    }

    template class NonlocalFormFactorObservable<nff::BToK, nff::PToP>;
}
