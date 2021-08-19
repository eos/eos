/*
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

#include <eos/utils/exception.hh>
#include <eos/utils/expression.hh>
#include <eos/utils/stringify.hh>

namespace eos::exp
{
    double BinaryExpression::sum(const double & a, const double & b)        { return a + b; }
    double BinaryExpression::difference(const double & a, const double & b) { return a - b; }
    double BinaryExpression::product(const double & a, const double & b)    { return a * b; }
    double BinaryExpression::ratio(const double & a, const double & b)      { return a / b; }

    BinaryExpression::func
    BinaryExpression::Method(char op)
    {
        switch(op) {
            case '+': return BinaryExpression::sum;
            case '-': return BinaryExpression::difference;
            case '*': return BinaryExpression::product;
            case '/': return BinaryExpression::ratio;
            default:
                InternalError("Unknown binary operator '" + stringify(op) + "' encountered");
                return nullptr;
        }
    }
}
