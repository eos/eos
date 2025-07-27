/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2019-2025 Danny van Dyk
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

#include <eos/statistics/goodness-of-fit.hh>
#include <eos/statistics/log-posterior.hh>
#include <eos/statistics/test-statistic-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <map>

namespace
{
    template <class... T_> struct overloads : T_... { using T_::operator()...; };
    // the following is strictly only necessary for C++17, but clang-14 fails to perform
    // Class Template Argument Deduction (CTAD) for the overloads struct even in C++20 mode.
    // Minimum clang version for CTAD support in C++20 mode is clang-17.
    template <class... T_> overloads(T_...) -> overloads<T_...>;
}

namespace eos
{
    template <>
    struct Implementation<GoodnessOfFit>
    {
        LogPosterior log_posterior;

        double total_chi_square;
        int total_degrees_of_freedom;

        std::map<QualifiedName, test_statistics::ChiSquare> chi_squares;

        Constraint * current_constraint;

        Implementation(const LogPosterior & log_posterior) :
            log_posterior(log_posterior),
            total_chi_square(0.0),
            total_degrees_of_freedom(log_posterior.informative_priors() - log_posterior.varied_parameters().size()),
            current_constraint(nullptr)
        {
            compute_test_statistics();
        }

        ~Implementation() = default;

        void compute_test_statistics()
        {
            auto log_likelihood = log_posterior.log_likelihood();
            auto observable_cache = log_likelihood.observable_cache();
            observable_cache.update();

            /* compute the test statistics for each constraint*/
            for (auto c = log_likelihood.begin(), c_end = log_likelihood.end() ; c != c_end ; ++c)
            {
                current_constraint = &(*c);

                /* compute the test statistics for each log-likelihood block */
                for (auto b = c->begin_blocks(), b_end = c->end_blocks() ; b != b_end ; ++b)
                {
                    auto test_statistic = (*b)->primary_test_statistic();

                    // apply special treatments based on its type
                    const auto visitor = overloads
                    {
                        [this] (const test_statistics::Empty &) { },
                        [this] (const test_statistics::ChiSquare & c)
                        {
                            total_chi_square += c.chi2;
                            total_degrees_of_freedom += c.dof;
                            chi_squares.insert(std::make_pair(current_constraint->name(), c));
                        }
                    };
                    std::visit(visitor, test_statistic);
                }
            }
        }
    };

    template <>
    struct WrappedForwardIteratorTraits<GoodnessOfFit::ChiSquareIteratorTag>
    {
        using UnderlyingIterator = std::map<QualifiedName, test_statistics::ChiSquare>::const_iterator;
    };
    template class WrappedForwardIterator<GoodnessOfFit::ChiSquareIteratorTag, const std::pair<const QualifiedName, test_statistics::ChiSquare>>;

    GoodnessOfFit::GoodnessOfFit(const LogPosterior & log_posterior) :
        PrivateImplementationPattern<GoodnessOfFit>(new Implementation<GoodnessOfFit>(log_posterior))
    {
    }

    GoodnessOfFit::~GoodnessOfFit() = default;

    double
    GoodnessOfFit::total_chi_square() const
    {
        return _imp->total_chi_square;
    }

    int
    GoodnessOfFit::total_degrees_of_freedom() const
    {
        return _imp->total_degrees_of_freedom;
    }

    GoodnessOfFit::ChiSquareIterator
    GoodnessOfFit::begin_chi_square() const
    {
        return _imp->chi_squares.begin();
    }

    GoodnessOfFit::ChiSquareIterator
    GoodnessOfFit::end_chi_square() const
    {
        return _imp->chi_squares.end();
    }
}
