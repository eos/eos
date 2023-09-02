/* vim: set sw=4 sts=4 et tw=150 foldmethod=marker : */

/*
 * Copyright (c) 2023 Danny van Dyk
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

#include <eos/observable-impl.hh>
#include <eos/c-decays/dq-to-l-nu.hh>
#include <eos/utils/concrete-cacheable-observable.hh>
#include <eos/utils/concrete_observable.hh>

namespace eos
{
    // Leptonic D decays
    // {{{
    ObservableGroup
    make_dq_to_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $D_q^+\to \ell^+\nu$ decays)",
            R"(The option "l" selects the charged lepton flavor.)",
            {
                make_observable("D_s^+->l^+nu::BR", R"(\mathcal{B}(D_s^+ \to \ell^+\nu))",
                        Unit::None(),
                        &DqToLeptonNeutrino::branching_ratio,
                        std::make_tuple(),
                        Options{ { "q", "s" } }),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    ObservableSection
    make_c_decays_section()
    {
        auto imp = new Implementation<ObservableSection>(
            "Observables in (semi)leptonic $c$-hadron decays",
            "",
            {
                // D_q^+ -> l^+ nu
                make_dq_to_l_nu_group(),
            }
        );

        return ObservableSection(imp);
    }
}
