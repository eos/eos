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

#ifndef EOS_GUARD_EOS_UTILS_EXPRESSION_OBSERVABLE_HH
#define EOS_GUARD_EOS_UTILS_EXPRESSION_OBSERVABLE_HH 1

#include <eos/observable.hh>
#include <eos/utils/expression-fwd.hh>
#include <eos/utils/expression-parser.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/qualified-name.hh>
#include <eos/utils/units.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

namespace eos
{
    class ExpressionObservable : public Observable
    {
        private:
            QualifiedName _name;

            Parameters _parameters;

            Kinematics _kinematics;

            Options _options;

            eos::exp::ExpressionPtr _expression;

        public:
            ExpressionObservable(const QualifiedName & name, const Parameters & parameters, const Kinematics & kinematics, const Options & options,
                                 const eos::exp::ExpressionPtr & expression);

            ExpressionObservable(const QualifiedName & name, const ObservableCache & cache, const Kinematics & kinematics, const Options & options,
                                 const eos::exp::ExpressionPtr & expression);

            ~ExpressionObservable() override = default;

            [[nodiscard]] double        evaluate() const override;
            [[nodiscard]] ObservablePtr clone() const override;
            [[nodiscard]] ObservablePtr clone(const Parameters & parameters) const override;

            [[nodiscard]] const QualifiedName &
            name() const override
            {
                return _name;
            }

            Parameters
            parameters() override
            {
                return _parameters;
            }

            Kinematics
            kinematics() override
            {
                return _kinematics;
            }

            Options
            options() override
            {
                return _options;
            }

            [[nodiscard]] const eos::exp::ExpressionPtr &
            expression() const
            {
                return _expression;
            }
    };

    class ExpressionObservableEntry : public ObservableEntry
    {
        private:
            QualifiedName _name;

            std::string _latex;

            Unit _unit;

            const eos::exp::ExpressionPtr _expression;

            std::vector<std::string> _kinematics_names;

            Options _forced_options;

            std::vector<OptionSpecification> _option_specifications;

        public:
            ExpressionObservableEntry(const QualifiedName & name, std::string latex, const Unit & unit, eos::exp::ExpressionPtr expression, const Options & forced_options);

            ~ExpressionObservableEntry() override = default;

            [[nodiscard]] ObservableEntry::KinematicVariableIterator begin_kinematic_variables() const override;
            [[nodiscard]] ObservableEntry::KinematicVariableIterator end_kinematic_variables() const override;

            [[nodiscard]] ObservableEntry::OptionIterator begin_options() const override;
            [[nodiscard]] ObservableEntry::OptionIterator end_options() const override;

            [[nodiscard]] ObservablePtr make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const override;

            [[nodiscard]] const QualifiedName &
            name() const override
            {
                return _name;
            }

            [[nodiscard]] const std::string &
            latex() const override
            {
                return _latex;
            }

            [[nodiscard]] const Unit &
            unit() const override
            {
                return _unit;
            }
    };
} // namespace eos

#endif
