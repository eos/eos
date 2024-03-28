/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Méril Reboud
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
#include <eos/maths/complex.hh>
#include <eos/nonleptonic-amplitudes/topological-amplitudes.hh>

using namespace test;
using namespace eos;

class TopologicalAmplitudesTest :
    public TestCase
{
    public:
        TopologicalAmplitudesTest() :
            TestCase("topological_amplitudes_test")
        {
        }

        virtual void run() const
        {
            Parameters p = Parameters::Defaults();
            p["nonleptonic::T@Topological"]      = 0.1;
            p["nonleptonic::C@Topological"]      = 0.1;
            p["nonleptonic::A@Topological"]      = 0.1;
            p["nonleptonic::E@Topological"]      = 0.1;
            p["nonleptonic::TES@Topological"]    = 0.1;
            p["nonleptonic::TAS@Topological"]    = 0.1;
            p["nonleptonic::TS@Topological"]     = 0.1;
            p["nonleptonic::TPA@Topological"]    = 0.1;
            p["nonleptonic::TP@Topological"]     = 0.1;
            p["nonleptonic::TSS@Topological"]    = 0.1;
            p["nonleptonic::P@Topological"]      = 0.1;
            p["nonleptonic::PT@Topological"]     = 0.1;
            p["nonleptonic::S@Topological"]      = 0.1;
            p["nonleptonic::PC@Topological"]     = 0.1;
            p["nonleptonic::PTA@Topological"]    = 0.1;
            p["nonleptonic::PA@Topological"]     = 0.1;
            p["nonleptonic::PTE@Topological"]    = 0.1;
            p["nonleptonic::PAS@Topological"]    = 0.1;
            p["nonleptonic::PSS@Topological"]    = 0.1;
            p["nonleptonic::PES@Topological"]    = 0.1;

            Options o
            {
                { "representation", "topological" },
                { "q", "d" },
                { "P1", "pi^+" },
                { "P2", "pi^+" }
            };

            TopologicalRepresentation<PToPP> d(p, o);

            static const double eps = 1.0e-5;

            TEST_CHECK_NEARLY_EQUAL(d.amplitude(), complex<double>(0., 0.), eps);
        }

} topological_amplitudes_test;
