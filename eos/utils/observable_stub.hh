/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016-2026 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_OBSERVABLE_STUB_HH
#define EOS_GUARD_EOS_UTILS_OBSERVABLE_STUB_HH 1

#include <eos/observable.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/qualified-name.hh>

namespace eos
{
    class ObservableStub : public Observable, public PrivateImplementationPattern<ObservableStub>
    {
        public:
            ObservableStub(const Parameters & parameters, const QualifiedName & name, const Kinematics & kinematics = Kinematics());

            ~ObservableStub() override;

            [[nodiscard]] const QualifiedName & name() const override;

            [[nodiscard]] double evaluate() const override;

            Kinematics kinematics() override;

            Parameters parameters() override;

            Options options() override;

            [[nodiscard]] ObservablePtr clone() const override;

            [[nodiscard]] ObservablePtr clone(const Parameters & parameters) const override;
    };
} // namespace eos

#endif
