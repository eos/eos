/*
 * Copyright (c) 2023 Danny van Dyk
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
#include <eos/utils/options.hh>
#include <eos/utils/units.hh>

#include <test/test.hh>

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
