/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017, 2019 Danny van Dyk
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

#include <eos/statistics/test-statistic.hh>

namespace eos
{
    namespace test_statistics
    {
        class Empty
        {
            public:
                Empty();

                ~Empty();

                virtual void output(std::ostream & stream) const;
        };

        class ChiSquare
        {
            public:
                double chi2;
                int dof;
                double signed_chi;

                ChiSquare(const double & chi2, const int & dof, const double & signed_chi);

                ChiSquare(const double & chi2, const int & dof);

                ~ChiSquare();

                virtual void output(std::ostream & stream) const;
        };
    }
}
