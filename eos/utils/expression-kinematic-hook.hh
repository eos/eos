/*
 * Copyright (c) 2026 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_EXPRESSION_KINEMATIC_HOOK_HH
#define EOS_GUARD_EOS_UTILS_EXPRESSION_KINEMATIC_HOOK_HH 1

#include <eos/utils/expression-fwd.hh>
#include <eos/utils/kinematic.hh>

namespace eos::exp
{
    // Computes a derived kinematic variable from an expression in other kinematic variables, constants, and functions.
    class ExpressionKinematicHook : public KinematicHook
    {
        private:
            ExpressionPtr _expression;

        public:
            // Binds the kinematic variables of the expression to the given kinematics.
            ExpressionKinematicHook(const ExpressionPtr & expression, const Kinematics & kinematics);

            ~ExpressionKinematicHook() = default;

            virtual double evaluate() const;

            virtual std::shared_ptr<const KinematicHook> clone(const Kinematics & kinematics) const;
    };
} // namespace eos::exp

#endif
