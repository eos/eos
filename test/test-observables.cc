/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2026 Danny van Dyk
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

#include <eos/observable.hh>
#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/test-observable.hh>

#include <test/test-observables.hh>

#include <memory>
#include <string>
#include <vector>

namespace eos::test
{
    void
    register_test_observables()
    {
        static bool registered = false;
        if (registered)
        {
            return;
        }
        registered = true;

        auto test_function = [](const eos::Parameters & p, const std::vector<eos::KinematicVariable> & kv, const eos::Options & o)
        {
            using namespace eos;
            return p["mass::c"] * std::stoi(o.get("multiplier"_ok, "1"_ov).str()) * (kv[1] - kv[0]);
        };

        std::shared_ptr<const eos::TestObservableEntry> obs_entry =
                std::make_shared<const eos::TestObservableEntry>("test::obs1", "", eos::Unit::Undefined(), test_function, std::vector<std::string>{ "q2_min", "q2_max" });
        eos::ObservableEntries::instance()->insert_or_assign("test::obs1", obs_entry);
    }
} // namespace eos::test
