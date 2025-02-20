/* vim: set sw=4 sts=4 et tw=150 foldmethod=marker : */

/*
 * Copyright (c) 2024 Danny van Dyk
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
#include <eos/s-decays/k-to-l-nu.hh>
#include <eos/utils/concrete-cacheable-observable.hh>
#include <eos/utils/concrete_observable.hh>

namespace eos
{
    // Leptonic K decays
    // {{{

    // K^- -> l nu
    // {{{
    ObservableGroup
    make_k_to_l_nu_group()
    {
        auto imp = new Implementation<ObservableGroup>(
                R"(Observables in $K \to \ell^-\bar{\nu}_\ell$ decays)",
                R"(The option "l" selects the charged lepton flavor.)",
                { make_observable("K->lnu::BR", R"(\mathcal{B}(K^- \to \ell^-\bar{\nu}_\ell))", Unit::None(), &KToLeptonNeutrino::branching_ratio) });

        return ObservableGroup(imp);
    }

    // }}}
    // }}}

    ObservableSection
    make_s_decays_section()
    {
        auto imp = new Implementation<ObservableSection>("Observables in (semi)leptonic $s$-hadron decays",
                                                         "",
                                                         { // K -> l^- nubar
                                                           make_k_to_l_nu_group() });

        return ObservableSection(imp);
    }
} // namespace eos
