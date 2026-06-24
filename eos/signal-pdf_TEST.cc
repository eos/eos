/*
 * Copyright (c) 2023-2026 Danny van Dyk
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

#include <eos/signal-pdf.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/units.hh>

#include <test/test.hh>

#include <cmath>

using namespace test;
using namespace eos;

class SignalPDFTest : public TestCase
{
    public:
        SignalPDFTest() :
            TestCase("signal_pdf_test")
        {
        }

        virtual void
        run() const
        {
            // test that the list of signal pdfs is non-empty
            {
                const auto & signal_pdfs = SignalPDFs();

                TEST_CHECK(std::distance(signal_pdfs.begin(), signal_pdfs.end()) > 0);
            }
        }

} signal_pdf_test;

class SignalPDFRuntimeInsertionTest : public TestCase
{
    public:
        SignalPDFRuntimeInsertionTest() :
            TestCase("signal_pdf_runtime_insertion_test")
        {
        }

        virtual void
        run() const
        {
            // insert a new signal PDF at run time that reuses the observables backing the
            // built-in 'TestLegendre1D::P(z)' PDF: pdf(z) = z (4 - z), int_0^4 dz pdf(z) = 32 / 3
            {
                auto signal_pdfs = SignalPDFs();

                // referencing an unknown numerator or normalization must fail fast
                TEST_CHECK_THROWS(UnknownObservableError,
                                  signal_pdfs.insert("TestLegendre1D::RuntimeP(z)",
                                                     "runtime-inserted test PDF",
                                                     Options{},
                                                     "TestLegendre1D::Unknown(z)",
                                                     { "z" },
                                                     "TestLegendre1D::NormalizationPDF(z)",
                                                     { "z_min", "z_max" }));

                TEST_CHECK_THROWS(UnknownObservableError,
                                  signal_pdfs.insert("TestLegendre1D::RuntimeP(z)",
                                                     "runtime-inserted test PDF",
                                                     Options{},
                                                     "TestLegendre1D::UnnormalizedPDF(z)",
                                                     { "z" },
                                                     "TestLegendre1D::Unknown(z)",
                                                     { "z_min", "z_max" }));

                // a valid insertion referencing the existing test observables
                signal_pdfs.insert("TestLegendre1D::RuntimeP(z)",
                                   "runtime-inserted test PDF",
                                   Options{},
                                   "TestLegendre1D::UnnormalizedPDF(z)",
                                   { "z" },
                                   "TestLegendre1D::NormalizationPDF(z)",
                                   { "z_min", "z_max" });

                // the new PDF must be visible through the container ...
                TEST_CHECK(signal_pdfs["TestLegendre1D::RuntimeP(z)"] != nullptr);

                // ... and through a freshly created container (the registry is shared)
                TEST_CHECK(SignalPDFs()["TestLegendre1D::RuntimeP(z)"] != nullptr);

                Parameters p = Parameters::Defaults();
                Kinematics k{
                    {     "z", 2.0 },
                    { "z_min", 0.0 },
                    { "z_max", 4.0 }
                };
                Options o;

                auto pdf = SignalPDF::make("TestLegendre1D::RuntimeP(z)", p, k, o);
                TEST_CHECK(pdf != nullptr);
                TEST_CHECK_EQUAL(pdf->name(), QualifiedName("TestLegendre1D::RuntimeP(z)"));

                // evaluate() returns the logarithm of the unnormalized PDF; pdf(z = 2) = 4
                TEST_CHECK_NEARLY_EQUAL(pdf->evaluate(), std::log(4.0), 1.0e-10);

                // normalization() returns the logarithm of the normalization integral = 32 / 3
                TEST_CHECK_NEARLY_EQUAL(pdf->normalization(), std::log(32.0 / 3.0), 1.0e-10);
            }
        }

} signal_pdf_runtime_insertion_test;
