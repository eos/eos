/*
 * Copyright (c) 2021      Méril Reboud
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

#include <eos/utils/exception.hh>
#include <eos/utils/expression.hh>
#include <eos/utils/stringify.hh>

#include <algorithm>
#include <cmath>

namespace eos::exp
{
    ExpressionError::ExpressionError(const std::string & msg) :
        eos::Exception("Invalid expression statement (" + msg + ")")
    {
    }

    double
    BinaryExpression::sum(const double & a, const double & b)
    {
        return a + b;
    }

    double
    BinaryExpression::difference(const double & a, const double & b)
    {
        return a - b;
    }

    double
    BinaryExpression::product(const double & a, const double & b)
    {
        return a * b;
    }

    double
    BinaryExpression::ratio(const double & a, const double & b)
    {
        return a / b;
    }

    double
    BinaryExpression::power(const double & a, const double & b)
    {
        return pow(a, b);
    }

    BinaryExpression::func
    BinaryExpression::Method(char op)
    {
        switch (op)
        {
            case '+': return BinaryExpression::sum;
            case '-': return BinaryExpression::difference;
            case '*': return BinaryExpression::product;
            case '/': return BinaryExpression::ratio;
            case '^': return BinaryExpression::power;
            default:  InternalError("Unknown binary operator '" + stringify(op) + "' encountered"); return nullptr;
        }
    }

    namespace
    {
        double
        exp_function(std::span<const double> x)
        {
            return std::exp(x[0]);
        }

        double
        sin_function(std::span<const double> x)
        {
            return std::sin(x[0]);
        }

        double
        cos_function(std::span<const double> x)
        {
            return std::cos(x[0]);
        }

        double
        atan_function(std::span<const double> x)
        {
            return std::atan(x[0]);
        }

        double
        log_function(std::span<const double> x)
        {
            return std::log(x[0]);
        }

        // Heaviside step function, with theta(0) = 1
        double
        theta_function(std::span<const double> x)
        {
            return x[0] >= 0.0 ? 1.0 : 0.0;
        }

        // Kernel::Gaussian(u, mu, sigma)
        double
        gaussian_kernel(std::span<const double> x)
        {
            const double z = (x[0] - x[1]) / x[2];

            return std::exp(-0.5 * z * z) / (std::sqrt(2.0 * M_PI) * x[2]);
        }
    } // namespace

    const std::map<std::string, FunctionEntry> &
    functions()
    {
        static const std::map<std::string, FunctionEntry> function_table{
            // elementary functions
            {              "exp",   { 1, &exp_function, false } },
            {              "sin",   { 1, &sin_function, false } },
            {              "cos",   { 1, &cos_function, false } },
            {             "atan",  { 1, &atan_function, false } },
            {              "log",   { 1, &log_function, false } },
            {            "theta", { 1, &theta_function, false } },

            // kernels: unit-area densities in their first argument
            { "Kernel::Gaussian", { 3, &gaussian_kernel, true } },
        };

        return function_table;
    }

    FunctionExpression::FunctionExpression(const std::string & f, std::span<const ExpressionPtr> args) :
        f(nullptr),
        fname(f),
        number_of_arguments(args.size())
    {
        const auto & function_table = functions();

        auto it = function_table.find(f);
        if (function_table.end() == it)
        {
            throw ExpressionError("unknown function name " + f);
        }

        if (it->second.arity != args.size())
        {
            throw ExpressionError("function " + f + " expects " + stringify(it->second.arity) + " arguments, but " + stringify(args.size()) + " were given");
        }

        std::copy(args.begin(), args.end(), this->args.begin());
        this->f = it->second.f;
    }
} // namespace eos::exp
