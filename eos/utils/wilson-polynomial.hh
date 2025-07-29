/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2025 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_WILSON_POLYNOMIAL_HH
#define EOS_GUARD_EOS_UTILS_WILSON_POLYNOMIAL_HH 1

#include <eos/observable.hh>
#include <eos/utils/expression-fwd.hh>
#include <eos/utils/expression.hh>

#include <list>
#include <string>

namespace eos
{
    exp::Expression make_polynomial(const ObservablePtr &, const std::list<std::string> &);

    exp::Expression make_polynomial_ratio(const ObservablePtr & numerator, const ObservablePtr & denominator, const std::list<std::string> & coefficients);
} // namespace eos

#endif
