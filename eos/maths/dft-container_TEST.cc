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

#include <test/test.hh>
#include <eos/maths/dft-container.hh>
#include <eos/utils/exception.hh>

using namespace test;
using namespace eos;

class DftContainerDoubleRankOneTest :
    public TestCase
{
    public:
        DftContainerDoubleRankOneTest() :
            TestCase("dft_container_double_rank_one_test")
        {
        }

        virtual void run() const
        {
            // standard initialisation and access for a rank 1 container
            {
                dft::Container<double, 1> container({10});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    container({i}) = static_cast<double>(i);
                }

                const double eps = std::numeric_limits<double>::epsilon();
                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(container({i}), static_cast<double>(i), eps);
                }
            }

            // copy constructor and assignment operator
            {
                dft::Container<double, 1> container({10});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    container({i}) = static_cast<double>(i);
                }

                dft::Container<double, 1> copy(container);
                const double eps = std::numeric_limits<double>::epsilon();
                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(copy({i}), static_cast<double>(i), eps);
                }

                dft::Container<double, 1> assigned = container;
                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(assigned({i}), static_cast<double>(i), eps);
                }
            }

            // move constructor and assignment operator
            {
                dft::Container<double, 1> container({10});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    container({i}) = static_cast<double>(i);
                }

                dft::Container<double, 1> moved(std::move(container));
                const double eps = std::numeric_limits<double>::epsilon();
                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(moved({i}), static_cast<double>(i), eps);
                }

                dft::Container<double, 1> assigned = std::move(moved);
                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(assigned({i}), static_cast<double>(i), eps);
                }
            }

            // assignment operator: self-assignment must not cause issues
            {
                dft::Container<double, 1> container({10});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    container({i}) = static_cast<double>(i);
                }

                container = container;

                const double eps = std::numeric_limits<double>::epsilon();
                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(container({i}), static_cast<double>(i), eps);
                }
            }

            // assignment operator: check result of copy assignment
            {
                dft::Container<double, 1> container({10});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    container({i}) = static_cast<double>(i);
                }

                dft::Container<double, 1> assigned({10});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    assigned({i}) = static_cast<double>(10 - i);
                }

                const double eps = std::numeric_limits<double>::epsilon();
                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(container({i}), static_cast<double>(i), eps);
                    TEST_CHECK_NEARLY_EQUAL(assigned({i}), static_cast<double>(10 - i), eps);
                }

                assigned = container;

                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(assigned({i}), static_cast<double>(i), eps);
                }
            }

            // assignment operator: self-move-assignment must not cause issues
            {
                dft::Container<double, 1> container({10});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    container({i}) = static_cast<double>(i);
                }

                container = std::move(container);

                const double eps = std::numeric_limits<double>::epsilon();
                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(container({i}), static_cast<double>(i), eps);
                }
            }

            // assignment operator: check result of move assignment
            {
                dft::Container<double, 1> container({10});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    container({i}) = static_cast<double>(i);
                }

                dft::Container<double, 1> assigned({10});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    assigned({i}) = static_cast<double>(10 - i);
                }

                const double eps = std::numeric_limits<double>::epsilon();
                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(container({i}), static_cast<double>(i), eps);
                    TEST_CHECK_NEARLY_EQUAL(assigned({i}), static_cast<double>(10 - i), eps);
                }

                assigned = std::move(container);

                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(assigned({i}), static_cast<double>(i), eps);
                }
            }

            // out-of-bounds access: must throw
            {
                dft::Container<double, 1> container({10});
                TEST_CHECK_THROWS(std::out_of_range, container({10}));
            }
        }
} dft_container_double_rank_one_test;

class DftContainerComplexDoubleRankOneTest :
    public TestCase
{
    public:
        DftContainerComplexDoubleRankOneTest() :
            TestCase("dft_container_complex_double_rank_one_test")
        {
        }

        virtual void run() const
        {
            // standard initialisation and access for a rank 1 container
            {
                dft::Container<std::complex<double>, 1> container({10});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    container({i}) = static_cast<std::complex<double>>(static_cast<double>(i));
                }

                const double eps = std::numeric_limits<double>::epsilon();
                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(container({i}), static_cast<std::complex<double>>(static_cast<double>(i)), eps);
                }
            }

            // copy constructor and assignment operator
            {
                dft::Container<std::complex<double>, 1> container({10});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    container({i}) = std::complex<double>(static_cast<double>(i), -static_cast<double>(i + 1));
                }

                dft::Container<std::complex<double>, 1> copy(container);
                const double eps = std::numeric_limits<double>::epsilon();
                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(std::real(copy({i})), +static_cast<double>(i), eps);
                    TEST_CHECK_NEARLY_EQUAL(std::imag(copy({i})), -static_cast<double>(i + 1), eps);
                }

                dft::Container<std::complex<double>, 1> assigned = container;
                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(std::real(assigned({i})), +static_cast<double>(i), eps);
                    TEST_CHECK_NEARLY_EQUAL(std::imag(assigned({i})), -static_cast<double>(i + 1), eps);
                }
            }

            // move constructor and assignment operator
            {
                dft::Container<std::complex<double>, 1> container({10});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    container({i}) = std::complex<double>(static_cast<double>(i), -static_cast<double>(i + 1));
                }

                dft::Container<std::complex<double>, 1> moved(std::move(container));
                const double eps = std::numeric_limits<double>::epsilon();
                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(std::real(moved({i})), +static_cast<double>(i), eps);
                    TEST_CHECK_NEARLY_EQUAL(std::imag(moved({i})), -static_cast<double>(i + 1), eps);
                }

                dft::Container<std::complex<double>, 1> assigned = std::move(moved);
                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(std::real(assigned({i})), +static_cast<double>(i), eps);
                    TEST_CHECK_NEARLY_EQUAL(std::imag(assigned({i})), -static_cast<double>(i + 1), eps);
                }
            }

            // assignment operator: self-assignment must not cause issues
            {
                dft::Container<std::complex<double>, 1> container({10});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    container({i}) = std::complex<double>(static_cast<double>(i), -static_cast<double>(i + 1));
                }

                container = container;

                const double eps = std::numeric_limits<double>::epsilon();
                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(std::real(container({i})), +static_cast<double>(i), eps);
                    TEST_CHECK_NEARLY_EQUAL(std::imag(container({i})), -static_cast<double>(i + 1), eps);
                }
            }

            // assignment operator: check result of copy assignment
            {
                dft::Container<std::complex<double>, 1> container({10});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    container({i}) = std::complex<double>(static_cast<double>(i), -static_cast<double>(i + 1));
                }

                dft::Container<std::complex<double>, 1> assigned({10});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    assigned({i}) = std::complex<double>(static_cast<double>(10 - i), -static_cast<double>(10 - i + 1));
                }

                const double eps = std::numeric_limits<double>::epsilon();
                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(std::real(container({i})), +static_cast<double>(i), eps);
                    TEST_CHECK_NEARLY_EQUAL(std::imag(container({i})), -static_cast<double>(i + 1), eps);
                    TEST_CHECK_NEARLY_EQUAL(std::real(assigned({i})), +static_cast<double>(10 - i), eps);
                    TEST_CHECK_NEARLY_EQUAL(std::imag(assigned({i})), -static_cast<double>(10 - i + 1), eps);
                }

                assigned = container;

                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(std::real(assigned({i})), +static_cast<double>(i), eps);
                    TEST_CHECK_NEARLY_EQUAL(std::imag(assigned({i})), -static_cast<double>(i + 1), eps);
                }
            }

            // assignment operator: self-move-assignment must not cause issues
            {
                dft::Container<std::complex<double>, 1> container({10});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    container({i}) = std::complex<double>(static_cast<double>(i), -static_cast<double>(i + 1));
                }

                container = std::move(container);

                const double eps = std::numeric_limits<double>::epsilon();
                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(std::real(container({i})), +static_cast<double>(i), eps);
                    TEST_CHECK_NEARLY_EQUAL(std::imag(container({i})), -static_cast<double>(i + 1), eps);
                }
            }

            // assignment operator: check result of move assignment
            {
                dft::Container<std::complex<double>, 1> container({10});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    container({i}) = std::complex<double>(static_cast<double>(i), -static_cast<double>(i + 1));
                }

                dft::Container<std::complex<double>, 1> assigned({10});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    assigned({i}) = std::complex<double>(static_cast<double>(10 - i), -static_cast<double>(10 - i + 1));
                }

                const double eps = std::numeric_limits<double>::epsilon();
                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(std::real(container({i})), +static_cast<double>(i), eps);
                    TEST_CHECK_NEARLY_EQUAL(std::imag(container({i})), -static_cast<double>(i + 1), eps);
                    TEST_CHECK_NEARLY_EQUAL(std::real(assigned({i})), +static_cast<double>(10 - i), eps);
                    TEST_CHECK_NEARLY_EQUAL(std::imag(assigned({i})), -static_cast<double>(10 - i + 1), eps);
                }

                assigned = std::move(container);

                for (std::size_t i = 0; i < 10; ++i)
                {
                    TEST_CHECK_NEARLY_EQUAL(std::real(assigned({i})), +static_cast<double>(i), eps);
                    TEST_CHECK_NEARLY_EQUAL(std::imag(assigned({i})), -static_cast<double>(i + 1), eps);
                }
            }

            // out-of-bounds access: must throw
            {
                dft::Container<std::complex<double>, 1> container({10});
                TEST_CHECK_THROWS(std::out_of_range, container({10}));
            }
        }
} dft_container_complex_double_rank_one_test;

class DftContainerDoubleRanksTest : public TestCase
{
    public:
        DftContainerDoubleRanksTest() :
            TestCase("dft_container_double_ranks_test")
        {
        }

        virtual void run() const
        {
            // test rank 2 container layout: FFTW3 uses row-major order, so the last index should be contiguous in memory
            {
                dft::Container<double, 2> container({10, 20});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    for (std::size_t j = 0; j < 20; ++j)
                    {
                        container({i, j}) = static_cast<double>(i * 20 + j);
                    }
                }

                const double eps = std::numeric_limits<double>::epsilon();
                double * data = container.data();
                for (std::size_t i = 0; i < 10; ++i)
                {
                    for (std::size_t j = 0; j < 20; ++j)
                    {
                        TEST_CHECK_NEARLY_EQUAL(data[i * 20 + j], static_cast<double>(i * 20 + j), eps);
                    }
                }
            }

            // test rank 3 container layout: FFTW3 uses row-major order, so the last index should be contiguous in memory
            {
                dft::Container<double, 3> container({10, 20, 30});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    for (std::size_t j = 0; j < 20; ++j)
                    {
                        for (std::size_t k = 0; k < 30; ++k)
                        {
                            container({i, j, k}) = static_cast<double>(i * 600 + j * 30 + k);
                        }
                    }
                }

                const double eps = std::numeric_limits<double>::epsilon();
                double * data = container.data();
                for (std::size_t i = 0; i < 10; ++i)
                {
                    for (std::size_t j = 0; j < 20; ++j)
                    {
                        for (std::size_t k = 0; k < 30; ++k)
                        {
                            TEST_CHECK_NEARLY_EQUAL(data[i * 600 + j * 30 + k], static_cast<double>(i * 600 + j * 30 + k), eps);
                        }
                    }
                }
            }

            // test rank 4 container layout: FFTW3 uses row-major order, so the last index should be contiguous in memory
            {
                dft::Container<double, 4> container({10, 20, 30, 40});
                for (std::size_t i = 0; i < 10; ++i)
                {
                    for (std::size_t j = 0; j < 20; ++j)
                    {
                        for (std::size_t k = 0; k < 30; ++k)
                        {
                            for (std::size_t l = 0; l < 40; ++l)
                            {
                                container({i, j, k, l}) = static_cast<double>(i * 24000 + j * 1200 + k * 40 + l);
                            }
                        }
                    }
                }

                const double eps = std::numeric_limits<double>::epsilon();
                double * data = container.data();
                for (std::size_t i = 0; i < 10; ++i)
                {
                    for (std::size_t j = 0; j < 20; ++j)
                    {
                        for (std::size_t k = 0; k < 30; ++k)
                        {
                            for (std::size_t l = 0; l < 40; ++l)
                            {
                                TEST_CHECK_NEARLY_EQUAL(data[i * 24000 + j * 1200 + k * 40 + l], static_cast<double>(i * 24000 + j * 1200 + k * 40 + l), eps);
                            }
                        }
                    }
                }
            }
        }
} dft_container_double_ranks_test;

class DftContainerComplexDoubleMultiplicationTest :
    public TestCase
{
    public:
        DftContainerComplexDoubleMultiplicationTest() :
            TestCase("dft_container_complex_double_multiplication_test")
        {
        }

        virtual void run() const
        {
            const double eps = std::numeric_limits<double>::epsilon();

            // operator*=: result is correct
            {
                dft::Container<std::complex<double>, 1> lhs({4}), rhs({4});
                lhs({0}) = { 1.0,  2.0}; rhs({0}) = { 3.0,  4.0};
                lhs({1}) = { 5.0,  6.0}; rhs({1}) = { 7.0,  8.0};
                lhs({2}) = { 9.0, 10.0}; rhs({2}) = {11.0, 12.0};
                lhs({3}) = {13.0, 14.0}; rhs({3}) = {15.0, 16.0};

                lhs *= rhs;

                TEST_CHECK_NEARLY_EQUAL(std::real(lhs({0})), 1.0 * 3.0 - 2.0 * 4.0,   eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(lhs({0})), 1.0 * 4.0 + 2.0 * 3.0,   eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(lhs({1})), 5.0 * 7.0 - 6.0 * 8.0,   eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(lhs({1})), 5.0 * 8.0 + 6.0 * 7.0,   eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(lhs({2})), 9.0 * 11.0 - 10.0 * 12.0, eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(lhs({2})), 9.0 * 12.0 + 10.0 * 11.0, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(lhs({3})), 13.0 * 15.0 - 14.0 * 16.0, eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(lhs({3})), 13.0 * 16.0 + 14.0 * 15.0, eps);
            }

            // operator*: result is correct and operands are unchanged
            {
                dft::Container<std::complex<double>, 1> lhs({2}), rhs({2});
                lhs({0}) = {1.0, 2.0}; rhs({0}) = {3.0, 4.0};
                lhs({1}) = {5.0, 6.0}; rhs({1}) = {7.0, 8.0};

                dft::Container<std::complex<double>, 1> result = lhs * rhs;

                TEST_CHECK_NEARLY_EQUAL(std::real(result({0})), 1.0 * 3.0 - 2.0 * 4.0, eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(result({0})), 1.0 * 4.0 + 2.0 * 3.0, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(result({1})), 5.0 * 7.0 - 6.0 * 8.0, eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(result({1})), 5.0 * 8.0 + 6.0 * 7.0, eps);

                // operands must be unchanged
                TEST_CHECK_NEARLY_EQUAL(std::real(lhs({0})), 1.0, eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(lhs({0})), 2.0, eps);
                TEST_CHECK_NEARLY_EQUAL(std::real(rhs({0})), 3.0, eps);
                TEST_CHECK_NEARLY_EQUAL(std::imag(rhs({0})), 4.0, eps);
            }

            // operator*=: dimension mismatch must throw
            {
                dft::Container<std::complex<double>, 1> lhs({4}), rhs({5});
                TEST_CHECK_THROWS(eos::InternalError, lhs *= rhs);
            }
        }
} dft_container_complex_double_multiplication_test;
