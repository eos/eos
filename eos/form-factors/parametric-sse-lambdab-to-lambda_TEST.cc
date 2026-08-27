/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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

#include <eos/form-factors/parametric-sse-impl-onehalfplus-to-onehalfplus.hh>
#include <eos/maths/power-of.hh>

#include <test/test.hh>

using namespace test;
using namespace eos;

class SSEOneHalfPlusToOneHalfPlusFormFactorsTest : public TestCase
{
    public:
        SSEOneHalfPlusToOneHalfPlusFormFactorsTest() :
            TestCase("sse_onehalfplus_to_onehalfplus_form_factors_test")
        {
        }

        virtual void
        run() const
        {
            /* Lambda_b -> Lambda */
            {
                Parameters                                             p  = Parameters::Defaults();
                std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> ff = FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Lambda_b->Lambda::SSE", p, Options{});
                TEST_CHECK(ff.get() != nullptr);

                // arbitrary, non-zero coefficients: the equations of motion must hold
                // identically, regardless of their values.
                p["Lambda_b->Lambda::alpha^(0,V)_0@SSE"]    = 0.10;
                p["Lambda_b->Lambda::alpha^(0,V)_1@SSE"]    = -0.20;
                p["Lambda_b->Lambda::alpha^(0,V)_2@SSE"]    = 0.30;
                p["Lambda_b->Lambda::alpha^(perp,V)_0@SSE"] = 0.15;
                p["Lambda_b->Lambda::alpha^(perp,V)_1@SSE"] = -0.25;
                p["Lambda_b->Lambda::alpha^(perp,V)_2@SSE"] = 0.35;
                p["Lambda_b->Lambda::alpha^(t,V)_1@SSE"]    = 0.40;
                p["Lambda_b->Lambda::alpha^(t,V)_2@SSE"]    = -0.45;

                p["Lambda_b->Lambda::alpha^(0,A)_0@SSE"]    = 0.12;
                p["Lambda_b->Lambda::alpha^(0,A)_1@SSE"]    = -0.22;
                p["Lambda_b->Lambda::alpha^(0,A)_2@SSE"]    = 0.32;
                p["Lambda_b->Lambda::alpha^(t,A)_1@SSE"]    = 0.42;
                p["Lambda_b->Lambda::alpha^(t,A)_2@SSE"]    = -0.47;
                p["Lambda_b->Lambda::alpha^(perp,A)_1@SSE"] = 0.52;
                p["Lambda_b->Lambda::alpha^(perp,A)_2@SSE"] = -0.14;

                p["Lambda_b->Lambda::alpha^(0,T)_0@SSE"]    = 0.11;
                p["Lambda_b->Lambda::alpha^(0,T)_1@SSE"]    = -0.21;
                p["Lambda_b->Lambda::alpha^(0,T)_2@SSE"]    = 0.31;
                p["Lambda_b->Lambda::alpha^(perp,T)_1@SSE"] = 0.41;
                p["Lambda_b->Lambda::alpha^(perp,T)_2@SSE"] = -0.16;

                p["Lambda_b->Lambda::alpha^(perp,T5)_0@SSE"] = 0.13;
                p["Lambda_b->Lambda::alpha^(perp,T5)_1@SSE"] = -0.23;
                p["Lambda_b->Lambda::alpha^(perp,T5)_2@SSE"] = 0.33;
                p["Lambda_b->Lambda::alpha^(0,T5)_1@SSE"]    = 0.43;
                p["Lambda_b->Lambda::alpha^(0,T5)_2@SSE"]    = -0.18;

                // t_- computed from the masses used by the parametrisation
                const double m_1 = p["mass::Lambda_b@HME"].evaluate();
                const double m_2 = p["mass::Lambda@HME"].evaluate();
                const double tm  = power_of<2>(m_1 - m_2);

                // key regression assertions: the equations of motion must be enforced exactly.
                TEST_CHECK_RELATIVE_ERROR(ff->f_time_v(0.0), ff->f_long_v(0.0), 1.0e-12);
                TEST_CHECK_RELATIVE_ERROR(ff->f_time_a(0.0), ff->f_long_a(0.0), 1.0e-12);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_a(tm), ff->f_long_a(tm), 1.0e-12);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_t(0.0), ff->f_perp_t5(0.0), 1.0e-12);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_t5(tm), ff->f_perp_t5(tm), 1.0e-12);

                // pin the parametrisation to reference values (the z expansion is not
                // shifted, so z(0) != 0)
                static const double eps = 1.0e-6;
                TEST_CHECK_RELATIVE_ERROR(ff->f_time_v(0.0), 0.080724325011045192, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_v(0.0), 0.080724325011045192, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_v(0.0), 0.12556397490158280, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_time_a(0.0), 0.097285779093973457, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_a(0.0), 0.097285779093973457, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_a(0.0), 0.28455837646939042, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_t(0.0), 0.089692254989152731, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_t(0.0), 0.10618053613201273, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_t5(0.0), 0.27320744730468194, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_t5(0.0), 0.10618053613201273, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->f_time_v(10.0), 0.074536170153152739, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_v(10.0), 0.14357601817985308, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_v(10.0), 0.21737846735415511, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_time_a(10.0), 0.099421791251862179, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_a(10.0), 0.16034124175343792, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_a(10.0), 0.33704807880981136, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_t(10.0), 0.15833650801471347, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_t(10.0), 0.10904830756139758, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_t5(10.0), 0.33377893554379762, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_t5(10.0), 0.17406378901236480, eps);
            }

            /* Lambda_c -> neutron */
            {
                Parameters                                             p  = Parameters::Defaults();
                std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> ff = FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Lambda_c->neutron::SSE", p, Options{});
                TEST_CHECK(ff.get() != nullptr);

                p["Lambda_c->neutron::alpha^(0,V)_0@SSE"]    = 0.10;
                p["Lambda_c->neutron::alpha^(0,V)_1@SSE"]    = -0.20;
                p["Lambda_c->neutron::alpha^(0,V)_2@SSE"]    = 0.30;
                p["Lambda_c->neutron::alpha^(perp,V)_0@SSE"] = 0.15;
                p["Lambda_c->neutron::alpha^(perp,V)_1@SSE"] = -0.25;
                p["Lambda_c->neutron::alpha^(perp,V)_2@SSE"] = 0.35;
                p["Lambda_c->neutron::alpha^(t,V)_1@SSE"]    = 0.40;
                p["Lambda_c->neutron::alpha^(t,V)_2@SSE"]    = -0.45;

                p["Lambda_c->neutron::alpha^(0,A)_0@SSE"]    = 0.12;
                p["Lambda_c->neutron::alpha^(0,A)_1@SSE"]    = -0.22;
                p["Lambda_c->neutron::alpha^(0,A)_2@SSE"]    = 0.32;
                p["Lambda_c->neutron::alpha^(t,A)_1@SSE"]    = 0.42;
                p["Lambda_c->neutron::alpha^(t,A)_2@SSE"]    = -0.47;
                p["Lambda_c->neutron::alpha^(perp,A)_1@SSE"] = 0.52;
                p["Lambda_c->neutron::alpha^(perp,A)_2@SSE"] = -0.14;

                p["Lambda_c->neutron::alpha^(0,T)_0@SSE"]    = 0.11;
                p["Lambda_c->neutron::alpha^(0,T)_1@SSE"]    = -0.21;
                p["Lambda_c->neutron::alpha^(0,T)_2@SSE"]    = 0.31;
                p["Lambda_c->neutron::alpha^(perp,T)_1@SSE"] = 0.41;
                p["Lambda_c->neutron::alpha^(perp,T)_2@SSE"] = -0.16;

                p["Lambda_c->neutron::alpha^(perp,T5)_0@SSE"] = 0.13;
                p["Lambda_c->neutron::alpha^(perp,T5)_1@SSE"] = -0.23;
                p["Lambda_c->neutron::alpha^(perp,T5)_2@SSE"] = 0.33;
                p["Lambda_c->neutron::alpha^(0,T5)_1@SSE"]    = 0.43;
                p["Lambda_c->neutron::alpha^(0,T5)_2@SSE"]    = -0.18;

                const double m_1 = p["mass::Lambda_c@HME"].evaluate();
                const double m_2 = p["mass::neutron@HME"].evaluate();
                const double tm  = power_of<2>(m_1 - m_2);

                TEST_CHECK_RELATIVE_ERROR(ff->f_time_v(0.0), ff->f_long_v(0.0), 1.0e-12);
                TEST_CHECK_RELATIVE_ERROR(ff->f_time_a(0.0), ff->f_long_a(0.0), 1.0e-12);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_a(tm), ff->f_long_a(tm), 1.0e-12);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_t(0.0), ff->f_perp_t5(0.0), 1.0e-12);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_t5(tm), ff->f_perp_t5(tm), 1.0e-12);

                static const double eps = 1.0e-6;
                TEST_CHECK_RELATIVE_ERROR(ff->f_time_v(0.0), 0.086657794398335666, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_v(0.0), 0.086657794398335666, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_v(0.0), 0.13318089812831060, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_time_a(0.0), 0.10737396864379092, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_a(0.0), 0.10737396864379091, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_a(0.0), 0.20091151956892409, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_t(0.0), 0.09596241514433064, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_t(0.0), 0.11678190188915193, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_t5(0.0), 0.20020728514670313, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_t5(0.0), 0.11678190188915193, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->f_time_v(1.0), 0.078199650199289997, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_v(1.0), 0.13193436324656799, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_v(1.0), 0.19813768090430819, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_time_a(1.0), 0.11671092472376722, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_a(1.0), 0.14431285647218686, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_a(1.0), 0.20426410789082075, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_t(1.0), 0.14517502677811603, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_t(1.0), 0.11736486650533799, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_t5(1.0), 0.21030409590211466, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_t5(1.0), 0.15635394602160990, eps);
            }

            /* Lambda_c -> proton (rare c -> u) */
            {
                Parameters                                             p  = Parameters::Defaults();
                std::shared_ptr<FormFactors<OneHalfPlusToOneHalfPlus>> ff = FormFactorFactory<OneHalfPlusToOneHalfPlus>::create("Lambda_c->proton::SSE", p, Options{});
                TEST_CHECK(ff.get() != nullptr);

                p["Lambda_c->proton::alpha^(0,V)_0@SSE"]    = 0.10;
                p["Lambda_c->proton::alpha^(0,V)_1@SSE"]    = -0.20;
                p["Lambda_c->proton::alpha^(0,V)_2@SSE"]    = 0.30;
                p["Lambda_c->proton::alpha^(perp,V)_0@SSE"] = 0.15;
                p["Lambda_c->proton::alpha^(perp,V)_1@SSE"] = -0.25;
                p["Lambda_c->proton::alpha^(perp,V)_2@SSE"] = 0.35;
                p["Lambda_c->proton::alpha^(t,V)_1@SSE"]    = 0.40;
                p["Lambda_c->proton::alpha^(t,V)_2@SSE"]    = -0.45;

                p["Lambda_c->proton::alpha^(0,A)_0@SSE"]    = 0.12;
                p["Lambda_c->proton::alpha^(0,A)_1@SSE"]    = -0.22;
                p["Lambda_c->proton::alpha^(0,A)_2@SSE"]    = 0.32;
                p["Lambda_c->proton::alpha^(t,A)_1@SSE"]    = 0.42;
                p["Lambda_c->proton::alpha^(t,A)_2@SSE"]    = -0.47;
                p["Lambda_c->proton::alpha^(perp,A)_1@SSE"] = 0.52;
                p["Lambda_c->proton::alpha^(perp,A)_2@SSE"] = -0.14;

                p["Lambda_c->proton::alpha^(0,T)_0@SSE"]    = 0.11;
                p["Lambda_c->proton::alpha^(0,T)_1@SSE"]    = -0.21;
                p["Lambda_c->proton::alpha^(0,T)_2@SSE"]    = 0.31;
                p["Lambda_c->proton::alpha^(perp,T)_1@SSE"] = 0.41;
                p["Lambda_c->proton::alpha^(perp,T)_2@SSE"] = -0.16;

                p["Lambda_c->proton::alpha^(perp,T5)_0@SSE"] = 0.13;
                p["Lambda_c->proton::alpha^(perp,T5)_1@SSE"] = -0.23;
                p["Lambda_c->proton::alpha^(perp,T5)_2@SSE"] = 0.33;
                p["Lambda_c->proton::alpha^(0,T5)_1@SSE"]    = 0.43;
                p["Lambda_c->proton::alpha^(0,T5)_2@SSE"]    = -0.18;

                const double m_1 = p["mass::Lambda_c@HME"].evaluate();
                const double m_2 = p["mass::proton@HME"].evaluate();
                const double tm  = power_of<2>(m_1 - m_2);

                TEST_CHECK_RELATIVE_ERROR(ff->f_time_v(0.0), ff->f_long_v(0.0), 1.0e-12);
                TEST_CHECK_RELATIVE_ERROR(ff->f_time_a(0.0), ff->f_long_a(0.0), 1.0e-12);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_a(tm), ff->f_long_a(tm), 1.0e-12);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_t(0.0), ff->f_perp_t5(0.0), 1.0e-12);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_t5(tm), ff->f_perp_t5(tm), 1.0e-12);

                static const double eps = 1.0e-6;
                TEST_CHECK_RELATIVE_ERROR(ff->f_time_v(0.0), 0.086627193573833236, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_v(0.0), 0.086627193573833236, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_v(0.0), 0.13314190305031082, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_time_a(0.0), 0.10734565008594721, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_a(0.0), 0.10734565008594722, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_a(0.0), 0.20111668997072799, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_t(0.0), 0.095930135469128749, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_t(0.0), 0.11675220536759258, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_t5(0.0), 0.20038583553509975, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_t5(0.0), 0.11675220536759258, eps);

                TEST_CHECK_RELATIVE_ERROR(ff->f_time_v(1.0), 0.078177157509018394, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_v(1.0), 0.13188212066412514, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_v(1.0), 0.19807232908825209, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_time_a(1.0), 0.11683494785533531, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_a(1.0), 0.14427102658905044, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_a(1.0), 0.20451520969917117, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_t(1.0), 0.14512016234895053, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_t(1.0), 0.11733203468586734, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_long_t5(1.0), 0.21052396646754837, eps);
                TEST_CHECK_RELATIVE_ERROR(ff->f_perp_t5(1.0), 0.15631021203973280, eps);
            }
        }
} sse_onehalfplus_to_onehalfplus_form_factors_test;
