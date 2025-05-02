/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015-2025 Danny van Dyk
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

#include <eos/form-factors/baryonic.hh>
#include <eos/form-factors/zero-recoil-sum-rule.hh>
#include <eos/models/model.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/maths/power-of.hh>

#include <iostream>

namespace eos
{
    template <>
    struct Implementation<ZeroRecoilSumRule<LambdaBToC>>
    {
        std::shared_ptr<Model> model;

        /* inclusive bounds */

        // renormalization scale (kinetic scheme)
        UsedParameter mu;

        // excitation energy cut off
        UsedParameter wM;

        // matrix elements of dim 5
        UsedParameter mu2_pi;

        // matrix elements of dim 6
        UsedParameter rho3_D;

        /* exclusive inelastic contributions */

        // Lambda_b -> Lambda_c(2595) form factors
        std::shared_ptr<FormFactors<OneHalfPlusToOneHalfMinus>> ff_2595;

        // Lambda_b -> Lambda_c(2625) form factors
        std::shared_ptr<FormFactors<OneHalfPlusToThreeHalfMinus>> ff_2625;

        // masses
        UsedParameter m_Lambda_b;
        UsedParameter m_Lambda_c_2595;
        UsedParameter m_Lambda_c_2625;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "SM"), p, o)),
            mu(p["Lambda_b->Lambda_c::mu@ZRSR"], u),
            wM(p["Lambda_b->Lambda_c::wM@ZRSR"], u),
            mu2_pi(p["Lambda_b->Lambda_b::mu_pi^2@1GeV"], u),
            rho3_D(p["Lambda_b->Lambda_b::rho_D^3@1GeV"], u),
            ff_2595(FormFactorFactory<OneHalfPlusToOneHalfMinus>::create("Lambda_b->Lambda_c(2595)@" + o.get("form-factors"_ok, "BBGIOvD2017"), p)),
            ff_2625(FormFactorFactory<OneHalfPlusToThreeHalfMinus>::create("Lambda_b->Lambda_c(2625)@" + o.get("form-factors"_ok, "BBGIOvD2017"), p)),
            m_Lambda_b(p["mass::Lambda_b"], u),
            m_Lambda_c_2595(p["mass::Lambda_c(2595)"], u),
            m_Lambda_c_2625(p["mass::Lambda_c(2625)"], u)
        {
            u.uses(*model);
            u.uses(*ff_2625);
        }

        inline double xiA() const
        {
            const double mb = model->m_b_kin(mu), mb2 = mb * mb;
            const double mc = model->m_c_kin(mu), mc2 = mc * mc;

            const double a_s = model->alpha_s(sqrt(mb * mc));
            const double mu2 = mu() * mu();
            const double mu_p = mu() - mc + std::sqrt(mc2 + mu2), mu_p2 = mu_p * mu_p;

            //std::cout << "alpha_s = " << a_s << std::endl;
            //std::cout << "m_b     = " << mb << std::endl;
            //std::cout << "m_c     = " << mc << std::endl;
            //std::cout << "mu'     = " << mu_p << std::endl;

            const double etaApert = 1.0 + a_s / M_PI * ((mb + mc) / (mb - mc) * log(mb / mc) - 8.0 / 3.0);
            const double etaAsoft = - a_s * mu2 / (3.0 * M_PI) * (1.0 / mc2 + 2.0 / (3.0 * mb * mc) + 1.0 / mb2);
            const double etaAspec = 4.0 * a_s / (3.0 * M_PI) * (
                        (wM - mu_p) * (wM + 2.0 * mc + mu_p) / (24.0 * mb2 * power_of<2>(wM + mc)) *
                        (2.0 * wM * (wM + 2.0 * mc) + mc2 * (mc2 - 3.0 * mb2 - 2.0 * mb * mc + 4.0 * mc * mu_p + 2.0 * mu_p2) / power_of<2>(mc + mu_p))
                        - (3.0 * mb - mc) * (mb + mc) / (12.0 * mb2) * log((mc + mu_p) / (mc + wM))
                    );

            //std::cout << "etaApert = " << etaApert << std::endl;
            //std::cout << "etaAsoft = " << etaAsoft << std::endl;
            //std::cout << "etaAspec = " << etaAspec << std::endl;

            return power_of<2>(etaApert) - etaAsoft + etaAspec;
        }

        inline double deltaA() const
        {
            const double mb = model->m_b_kin(mu), mb2 = mb * mb;
            const double mc = model->m_c_kin(mu), mc2 = mc * mc, mc3 = mc2 * mc;

            const double deltaA2 = mu2_pi / 4.0 * (1.0 / mc2 + 1.0 / mb2 + 2.0 / (3.0 * mb * mc));
            const double deltaA3 = rho3_D / (4.0 * mc3) + rho3_D / (12.0 * mb) * (1.0 / mc2 + 1.0 / (mc * mb) + 3.0 / mb2);

            return deltaA2 + deltaA3;
        }

        double xiV() const
        {
            // cf. [U2003], eq. (27)

            const double mb = model->m_b_kin(mu), mb2 = mb * mb;
            const double mc = model->m_c_kin(mu), mc2 = mc * mc;

            const double a_s = model->alpha_s(sqrt(mb * mc));
            const double mu = this->mu(), mu2 = mu * mu;
            const double wb = std::sqrt(mu2 + mb2), wb3 = wb * wb * wb;
            const double wc = std::sqrt(mu2 + mc2), wc3 = wc * wc * wc;

            //std::cout << "alpha_s = " << a_s << std::endl;
            //std::cout << "m_b     = " << mb << std::endl;
            //std::cout << "m_c     = " << mc << std::endl;


            const double xiVnlo1 =
                (3.0 * mb2 + 2.0 * mc * mb + 3.0 * mc2) / (2.0 * (mb2 - mc2)) * log((mu + wb) / (mu + wc))
                - 2.0;
            const double xiVnlo2 =
                4.0 / (3.0 * mu2) * (mc * wb - mb * wc) / (mb - mc)
                + 2.0 / 3.0 * (mc / wb - mb / wc) / (mb - mc)
                - 1.0 / 3.0 * (wb / mb - wc / mc) / (mb - mc)
                + 2.0 * mc * mb / (mb2 - mc2) * (1.0 / wb - 1.0 / wc)
                + 1.0 / (6.0 * (mc + mb)) * (wc / mc * (3.0 - mb / mc) + wb / mb * (3.0 - mc / mb))
                + 4.0 * mc * mb / (3.0 * (mc + mb)) * (mb / wb3 + mc / wc3)
                + mu / 6.0 * power_of<2>(1.0 / mc - 1.0 / mb)
                - 2.0 * mu2 / (3.0 * mc * mb) * (mb2 / wc + mc2 / wb) / (mb2 - mc2);

            return 1.0 + 2.0 * a_s / (3.0 * M_PI) * (xiVnlo1 - mu * xiVnlo2);
        }

        inline double deltaV() const
        {
            const double mb = model->m_b_kin(mu);
            const double mc = model->m_c_kin(mu);

            const double deltaV2 = mu2_pi / 4.0 * power_of<2>(1.0 / mc - 1.0 / mb);
            const double deltaV3 = rho3_D / 4.0 * power_of<2>(1.0 / mc - 1.0 / mb) * (1.0 / mc + 1.0 / mb);

            //std::cout << "mu2_pi  = " << mu2_pi << std::endl;
            //std::cout << "rho3_D  = " << rho3_D << std::endl;
            //std::cout << "deltaV2 = " << deltaV2 << std::endl;
            //std::cout << "deltaV3 = " << deltaV3 << std::endl;

            return deltaV2 + deltaV3;
        }

        // exclusive results for the inelastic contributions
        double f_inel() const
        {
            const double q2max = power_of<2>(m_Lambda_b - m_Lambda_c_2625);
            const double r     = power_of<2>((m_Lambda_b + m_Lambda_c_2625) / (m_Lambda_b - m_Lambda_c_2625));
            const double fT    = ff_2595->f_time_v(q2max);
            const double f0    = ff_2595->f_long_v(q2max);
            const double fP    = ff_2595->f_perp_v(q2max);
            const double F12T  = ff_2625->f_time12_v(q2max);
            const double F120  = ff_2625->f_long12_v(q2max);
            const double F12P  = ff_2625->f_perp12_v(q2max);
            const double F32P  = ff_2625->f_perp32_v(q2max);

            // note the normalization N_A = 1.0 in [MvD2015].
            const double f_inel_2595 =  1.0 / 3.0 * (power_of<2>(f0) + r * power_of<2>(fT) + 2.0 * power_of<2>(fP));
            const double f_inel_2625 =  2.0 / 3.0 * (power_of<2>(F120) + r * power_of<2>(F12T) + 2.0 * power_of<2>(F12P) + 6.0 * power_of<2>(F32P));

            return f_inel_2595 + f_inel_2625;
        }

        double g_inel() const
        {
            const double q2max = power_of<2>(m_Lambda_b - m_Lambda_c_2625);
            const double r     = power_of<2>((m_Lambda_b + m_Lambda_c_2625) / (m_Lambda_b - m_Lambda_c_2625));
            const double gT    = ff_2595->f_time_a(q2max);
            const double g0    = ff_2595->f_long_a(q2max);
            const double gP    = ff_2595->f_perp_a(q2max);
            const double G12T  = ff_2625->f_time12_a(q2max);
            const double G120  = ff_2625->f_long12_a(q2max);
            const double G12P  = ff_2625->f_perp12_a(q2max);
            const double G32P  = ff_2625->f_perp32_a(q2max);

            // note the normalization N_A = 3.0 in [MvD2015].
            const double g_inel_2595 = 1.0 / 9.0 * (power_of<2>(g0) + r * power_of<2>(gT) + 2.0 * power_of<2>(gP));
            const double g_inel_2625 = 2.0 / 9.0 * (power_of<2>(G120) + r * power_of<2>(G12T) + 2.0 * power_of<2>(G12P) + 6.0 * power_of<2>(G32P));

            return g_inel_2595 + g_inel_2625;
        }
    };

    const std::vector<OptionSpecification>
    Implementation<ZeroRecoilSumRule<LambdaBToC>>::options
    {
    };

    ZeroRecoilSumRule<LambdaBToC>::ZeroRecoilSumRule(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<ZeroRecoilSumRule<LambdaBToC>>(new Implementation<ZeroRecoilSumRule<LambdaBToC>>(parameters, options, *this))
    {
    }

    ZeroRecoilSumRule<LambdaBToC>::~ZeroRecoilSumRule()
    {
    }

    double
    ZeroRecoilSumRule<LambdaBToC>::axialvector_current() const
    {
        return _imp->xiA() - _imp->deltaA();
    }

    double
    ZeroRecoilSumRule<LambdaBToC>::vector_current() const
    {
        return _imp->xiV() - _imp->deltaV();
    }

    double
    ZeroRecoilSumRule<LambdaBToC>::axialvector_current_inel() const
    {
        return _imp->g_inel();
    }

    double
    ZeroRecoilSumRule<LambdaBToC>::vector_current_inel() const
    {
        return _imp->f_inel();
    }

    const std::set<ReferenceName>
    ZeroRecoilSumRule<LambdaBToC>::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    ZeroRecoilSumRule<LambdaBToC>::begin_options()
    {
        return Implementation<ZeroRecoilSumRule<LambdaBToC>>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    ZeroRecoilSumRule<LambdaBToC>::end_options()
    {
        return Implementation<ZeroRecoilSumRule<LambdaBToC>>::options.cend();
    }
}
