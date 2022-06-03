/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Danny van Dyk
 * Copyright (c) 2022 Philip Lüghausen
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
#include <eos/form-factors/b-lcdas-param.hh>

#include <eos/models/model.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;
using namespace b_lcdas;

class ParamTest :
    public TestCase
{
    public:
        ParamTest() :
            TestCase("b_lcdas_param_test")
        {
        }

        virtual void run() const
        {
            Parameters p = Parameters::Defaults();

            static const double eps = 1e-5;
            {
                Param lcda(p, Options{ { "q", "u" } });
                TEST_CHECK_NEARLY_EQUAL(lcda.phi_plus(1.0), 0.0, eps);
            }
        }
} b_lcdas_param_test;
