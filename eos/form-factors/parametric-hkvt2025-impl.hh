/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Florian Herren
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_HKVT2025_IMPL_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_HKVT2025_IMPL_HH 1

#include <eos/form-factors/parametric-hkvt2025.hh>
#include <eos/maths/integrate.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/diagnostics.hh>

#include <boost/math/special_functions/legendre.hpp>

#include <numeric>
#include <iostream>
namespace eos
{
    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    HKVT2025FormFactorTraits<Process_, PToPP>::resonance_0m_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), "mass::B_u@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s@BSZ2015" }
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    HKVT2025FormFactorTraits<Process_, PToPP>::resonance_1m_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), "mass::B_u^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c^*@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s^*@BSZ2015" }
    };

    template <typename Process_>
    const std::map<std::tuple<QuarkFlavor, QuarkFlavor>, std::string>
    HKVT2025FormFactorTraits<Process_, PToPP>::resonance_1p_names
    {
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::up), "mass::B_u,1@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::down), "mass::B_d,1@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::strange), "mass::B_s,1@BSZ2015" },
        { std::make_tuple(QuarkFlavor::bottom, QuarkFlavor::charm), "mass::B_c,1@BSZ2015" },
        { std::make_tuple(QuarkFlavor::charm,  QuarkFlavor::strange), "mass::D_s,1@BSZ2015" }
    };

    template <typename Process_>
    const std::map<std::string, unsigned>
    HKVT2025FormFactorTraits<Process_, PToPP>::charge_map
    {
        { "+-", 0 },
        { "00", 1 },
        { "+0", 2 }
    };

    template<typename Process_>
    HKVT2025FormFactors<Process_, PToPP>::HKVT2025FormFactors(const Parameters & p, const Options & o) :
        _a_g_1{{{ std::array{ UsedParameter(p[_exp_par_name("g_1", 0, 0, 0)], *this), UsedParameter(p[_exp_par_name("g_1", 0, 1, 0)], *this), UsedParameter(p[_exp_par_name("g_1", 0, 2, 0)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("g_1", 0, 0, 1)], *this), UsedParameter(p[_exp_par_name("g_1", 0, 1, 1)], *this), UsedParameter(p[_exp_par_name("g_1", 0, 2, 1)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("g_1", 0, 0, 2)], *this), UsedParameter(p[_exp_par_name("g_1", 0, 1, 2)], *this), UsedParameter(p[_exp_par_name("g_1", 0, 2, 2)], *this) } } ,
                { std::array{ UsedParameter(p[_exp_par_name("g_1", 1, 0, 0)], *this), UsedParameter(p[_exp_par_name("g_1", 1, 1, 0)], *this), UsedParameter(p[_exp_par_name("g_1", 1, 2, 0)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("g_1", 1, 0, 1)], *this), UsedParameter(p[_exp_par_name("g_1", 1, 1, 1)], *this), UsedParameter(p[_exp_par_name("g_1", 1, 2, 1)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("g_1", 1, 0, 2)], *this), UsedParameter(p[_exp_par_name("g_1", 1, 1, 2)], *this), UsedParameter(p[_exp_par_name("g_1", 1, 2, 2)], *this) } } ,
                { std::array{ UsedParameter(p[_exp_par_name("g_1", 2, 0, 0)], *this), UsedParameter(p[_exp_par_name("g_1", 2, 1, 0)], *this), UsedParameter(p[_exp_par_name("g_1", 2, 2, 0)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("g_1", 2, 0, 1)], *this), UsedParameter(p[_exp_par_name("g_1", 2, 1, 1)], *this), UsedParameter(p[_exp_par_name("g_1", 2, 2, 1)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("g_1", 2, 0, 2)], *this), UsedParameter(p[_exp_par_name("g_1", 2, 1, 2)], *this), UsedParameter(p[_exp_par_name("g_1", 2, 2, 2)], *this) } }
        }},
        _a_g_2{{{ std::array{ UsedParameter(p[_exp_par_name("g_2", 0, 0, 0)], *this), UsedParameter(p[_exp_par_name("g_2", 0, 1, 0)], *this), UsedParameter(p[_exp_par_name("g_2", 0, 2, 0)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("g_2", 0, 0, 1)], *this), UsedParameter(p[_exp_par_name("g_2", 0, 1, 1)], *this), UsedParameter(p[_exp_par_name("g_2", 0, 2, 1)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("g_2", 0, 0, 2)], *this), UsedParameter(p[_exp_par_name("g_2", 0, 1, 2)], *this), UsedParameter(p[_exp_par_name("g_2", 0, 2, 2)], *this) } } ,
                { std::array{ UsedParameter(p[_exp_par_name("g_2", 1, 0, 0)], *this), UsedParameter(p[_exp_par_name("g_2", 1, 1, 0)], *this), UsedParameter(p[_exp_par_name("g_2", 1, 2, 0)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("g_2", 1, 0, 1)], *this), UsedParameter(p[_exp_par_name("g_2", 1, 1, 1)], *this), UsedParameter(p[_exp_par_name("g_2", 1, 2, 1)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("g_2", 1, 0, 2)], *this), UsedParameter(p[_exp_par_name("g_2", 1, 1, 2)], *this), UsedParameter(p[_exp_par_name("g_2", 1, 2, 2)], *this) } } ,
                { std::array{ UsedParameter(p[_exp_par_name("g_2", 2, 0, 0)], *this), UsedParameter(p[_exp_par_name("g_2", 2, 1, 0)], *this), UsedParameter(p[_exp_par_name("g_2", 2, 2, 0)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("g_2", 2, 0, 1)], *this), UsedParameter(p[_exp_par_name("g_2", 2, 1, 1)], *this), UsedParameter(p[_exp_par_name("g_2", 2, 2, 1)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("g_2", 2, 0, 2)], *this), UsedParameter(p[_exp_par_name("g_2", 2, 1, 2)], *this), UsedParameter(p[_exp_par_name("g_2", 2, 2, 2)], *this) } }
        }},
        _a_f_1{{{ std::array{ UsedParameter(p[_exp_par_name("f_1", 0, 0, 0)], *this), UsedParameter(p[_exp_par_name("f_1", 0, 1, 0)], *this), UsedParameter(p[_exp_par_name("f_1", 0, 2, 0)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("f_1", 0, 0, 1)], *this), UsedParameter(p[_exp_par_name("f_1", 0, 1, 1)], *this), UsedParameter(p[_exp_par_name("f_1", 0, 2, 1)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("f_1", 0, 0, 2)], *this), UsedParameter(p[_exp_par_name("f_1", 0, 1, 2)], *this), UsedParameter(p[_exp_par_name("f_1", 0, 2, 2)], *this) } } ,
                { std::array{ UsedParameter(p[_exp_par_name("f_1", 1, 0, 0)], *this), UsedParameter(p[_exp_par_name("f_1", 1, 1, 0)], *this), UsedParameter(p[_exp_par_name("f_1", 1, 2, 0)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("f_1", 1, 0, 1)], *this), UsedParameter(p[_exp_par_name("f_1", 1, 1, 1)], *this), UsedParameter(p[_exp_par_name("f_1", 1, 2, 1)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("f_1", 1, 0, 2)], *this), UsedParameter(p[_exp_par_name("f_1", 1, 1, 2)], *this), UsedParameter(p[_exp_par_name("f_1", 1, 2, 2)], *this) } } ,
                { std::array{ UsedParameter(p[_exp_par_name("f_1", 2, 0, 0)], *this), UsedParameter(p[_exp_par_name("f_1", 2, 1, 0)], *this), UsedParameter(p[_exp_par_name("f_1", 2, 2, 0)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("f_1", 2, 0, 1)], *this), UsedParameter(p[_exp_par_name("f_1", 2, 1, 1)], *this), UsedParameter(p[_exp_par_name("f_1", 2, 2, 1)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("f_1", 2, 0, 2)], *this), UsedParameter(p[_exp_par_name("f_1", 2, 1, 2)], *this), UsedParameter(p[_exp_par_name("f_1", 2, 2, 2)], *this) } }
        }},
        _a_f_2{{{ std::array{ UsedParameter(p[_exp_par_name("f_2", 0, 0, 0)], *this), UsedParameter(p[_exp_par_name("f_2", 0, 1, 0)], *this), UsedParameter(p[_exp_par_name("f_2", 0, 2, 0)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("f_2", 0, 0, 1)], *this), UsedParameter(p[_exp_par_name("f_2", 0, 1, 1)], *this), UsedParameter(p[_exp_par_name("f_2", 0, 2, 1)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("f_2", 0, 0, 2)], *this), UsedParameter(p[_exp_par_name("f_2", 0, 1, 2)], *this), UsedParameter(p[_exp_par_name("f_2", 0, 2, 2)], *this) } } ,
                { std::array{ UsedParameter(p[_exp_par_name("f_2", 1, 0, 0)], *this), UsedParameter(p[_exp_par_name("f_2", 1, 1, 0)], *this), UsedParameter(p[_exp_par_name("f_2", 1, 2, 0)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("f_2", 1, 0, 1)], *this), UsedParameter(p[_exp_par_name("f_2", 1, 1, 1)], *this), UsedParameter(p[_exp_par_name("f_2", 1, 2, 1)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("f_2", 1, 0, 2)], *this), UsedParameter(p[_exp_par_name("f_2", 1, 1, 2)], *this), UsedParameter(p[_exp_par_name("f_2", 1, 2, 2)], *this) } } ,
                { std::array{ UsedParameter(p[_exp_par_name("f_2", 2, 0, 0)], *this), UsedParameter(p[_exp_par_name("f_2", 2, 1, 0)], *this), UsedParameter(p[_exp_par_name("f_2", 2, 2, 0)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("f_2", 2, 0, 1)], *this), UsedParameter(p[_exp_par_name("f_2", 2, 1, 1)], *this), UsedParameter(p[_exp_par_name("f_2", 2, 2, 1)], *this) } ,
                  std::array{ UsedParameter(p[_exp_par_name("f_2", 2, 0, 2)], *this), UsedParameter(p[_exp_par_name("f_2", 2, 1, 2)], *this), UsedParameter(p[_exp_par_name("f_2", 2, 2, 2)], *this) } }
        }},
        _a_F1_1{{{ std::array{ UsedParameter(p[_exp_par_name("F1_1", 0, 0, 0)], *this), UsedParameter(p[_exp_par_name("F1_1", 0, 1, 0)], *this), UsedParameter(p[_exp_par_name("F1_1", 0, 2, 0)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F1_1", 0, 0, 1)], *this), UsedParameter(p[_exp_par_name("F1_1", 0, 1, 1)], *this), UsedParameter(p[_exp_par_name("F1_1", 0, 2, 1)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F1_1", 0, 0, 2)], *this), UsedParameter(p[_exp_par_name("F1_1", 0, 1, 2)], *this), UsedParameter(p[_exp_par_name("F1_1", 0, 2, 2)], *this) } } ,
                 { std::array{ UsedParameter(p[_exp_par_name("F1_1", 1, 0, 0)], *this), UsedParameter(p[_exp_par_name("F1_1", 1, 1, 0)], *this), UsedParameter(p[_exp_par_name("F1_1", 1, 2, 0)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F1_1", 1, 0, 1)], *this), UsedParameter(p[_exp_par_name("F1_1", 1, 1, 1)], *this), UsedParameter(p[_exp_par_name("F1_1", 1, 2, 1)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F1_1", 1, 0, 2)], *this), UsedParameter(p[_exp_par_name("F1_1", 1, 1, 2)], *this), UsedParameter(p[_exp_par_name("F1_1", 1, 2, 2)], *this) } } ,
                 { std::array{ UsedParameter(p[_exp_par_name("F1_1", 2, 0, 0)], *this), UsedParameter(p[_exp_par_name("F1_1", 2, 1, 0)], *this), UsedParameter(p[_exp_par_name("F1_1", 2, 2, 0)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F1_1", 2, 0, 1)], *this), UsedParameter(p[_exp_par_name("F1_1", 2, 1, 1)], *this), UsedParameter(p[_exp_par_name("F1_1", 2, 2, 1)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F1_1", 2, 0, 2)], *this), UsedParameter(p[_exp_par_name("F1_1", 2, 1, 2)], *this), UsedParameter(p[_exp_par_name("F1_1", 2, 2, 2)], *this) } }
        }},
        _a_F1_2{{{ std::array{ UsedParameter(p[_exp_par_name("F1_2", 0, 0, 0)], *this), UsedParameter(p[_exp_par_name("F1_2", 0, 1, 0)], *this), UsedParameter(p[_exp_par_name("F1_2", 0, 2, 0)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F1_2", 0, 0, 1)], *this), UsedParameter(p[_exp_par_name("F1_2", 0, 1, 1)], *this), UsedParameter(p[_exp_par_name("F1_2", 0, 2, 1)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F1_2", 0, 0, 2)], *this), UsedParameter(p[_exp_par_name("F1_2", 0, 1, 2)], *this), UsedParameter(p[_exp_par_name("F1_2", 0, 2, 2)], *this) } } ,
                 { std::array{ UsedParameter(p[_exp_par_name("F1_2", 1, 0, 0)], *this), UsedParameter(p[_exp_par_name("F1_2", 1, 1, 0)], *this), UsedParameter(p[_exp_par_name("F1_2", 1, 2, 0)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F1_2", 1, 0, 1)], *this), UsedParameter(p[_exp_par_name("F1_2", 1, 1, 1)], *this), UsedParameter(p[_exp_par_name("F1_2", 1, 2, 1)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F1_2", 1, 0, 2)], *this), UsedParameter(p[_exp_par_name("F1_2", 1, 1, 2)], *this), UsedParameter(p[_exp_par_name("F1_2", 1, 2, 2)], *this) } } ,
                 { std::array{ UsedParameter(p[_exp_par_name("F1_2", 2, 0, 0)], *this), UsedParameter(p[_exp_par_name("F1_2", 2, 1, 0)], *this), UsedParameter(p[_exp_par_name("F1_2", 2, 2, 0)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F1_2", 2, 0, 1)], *this), UsedParameter(p[_exp_par_name("F1_2", 2, 1, 1)], *this), UsedParameter(p[_exp_par_name("F1_2", 2, 2, 1)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F1_2", 2, 0, 2)], *this), UsedParameter(p[_exp_par_name("F1_2", 2, 1, 2)], *this), UsedParameter(p[_exp_par_name("F1_2", 2, 2, 2)], *this) } }
        }},
        _a_F2_1{{{ std::array{ UsedParameter(p[_exp_par_name("F2_1", 0, 0, 0)], *this), UsedParameter(p[_exp_par_name("F2_1", 0, 1, 0)], *this), UsedParameter(p[_exp_par_name("F2_1", 0, 2, 0)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F2_1", 0, 0, 1)], *this), UsedParameter(p[_exp_par_name("F2_1", 0, 1, 1)], *this), UsedParameter(p[_exp_par_name("F2_1", 0, 2, 1)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F2_1", 0, 0, 2)], *this), UsedParameter(p[_exp_par_name("F2_1", 0, 1, 2)], *this), UsedParameter(p[_exp_par_name("F2_1", 0, 2, 2)], *this) } } ,
                 { std::array{ UsedParameter(p[_exp_par_name("F2_1", 1, 0, 0)], *this), UsedParameter(p[_exp_par_name("F2_1", 1, 1, 0)], *this), UsedParameter(p[_exp_par_name("F2_1", 1, 2, 0)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F2_1", 1, 0, 1)], *this), UsedParameter(p[_exp_par_name("F2_1", 1, 1, 1)], *this), UsedParameter(p[_exp_par_name("F2_1", 1, 2, 1)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F2_1", 1, 0, 2)], *this), UsedParameter(p[_exp_par_name("F2_1", 1, 1, 2)], *this), UsedParameter(p[_exp_par_name("F2_1", 1, 2, 2)], *this) } } ,
                 { std::array{ UsedParameter(p[_exp_par_name("F2_1", 2, 0, 0)], *this), UsedParameter(p[_exp_par_name("F2_1", 2, 1, 0)], *this), UsedParameter(p[_exp_par_name("F2_1", 2, 2, 0)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F2_1", 2, 0, 1)], *this), UsedParameter(p[_exp_par_name("F2_1", 2, 1, 1)], *this), UsedParameter(p[_exp_par_name("F2_1", 2, 2, 1)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F2_1", 2, 0, 2)], *this), UsedParameter(p[_exp_par_name("F2_1", 2, 1, 2)], *this), UsedParameter(p[_exp_par_name("F2_1", 2, 2, 2)], *this) } }
        }},
        _a_F2_2{{{ std::array{ UsedParameter(p[_exp_par_name("F2_2", 0, 0, 0)], *this), UsedParameter(p[_exp_par_name("F2_2", 0, 1, 0)], *this), UsedParameter(p[_exp_par_name("F2_2", 0, 2, 0)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F2_2", 0, 0, 1)], *this), UsedParameter(p[_exp_par_name("F2_2", 0, 1, 1)], *this), UsedParameter(p[_exp_par_name("F2_2", 0, 2, 1)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F2_2", 0, 0, 2)], *this), UsedParameter(p[_exp_par_name("F2_2", 0, 1, 2)], *this), UsedParameter(p[_exp_par_name("F2_2", 0, 2, 2)], *this) } } ,
                 { std::array{ UsedParameter(p[_exp_par_name("F2_2", 1, 0, 0)], *this), UsedParameter(p[_exp_par_name("F2_2", 1, 1, 0)], *this), UsedParameter(p[_exp_par_name("F2_2", 1, 2, 0)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F2_2", 1, 0, 1)], *this), UsedParameter(p[_exp_par_name("F2_2", 1, 1, 1)], *this), UsedParameter(p[_exp_par_name("F2_2", 1, 2, 1)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F2_2", 1, 0, 2)], *this), UsedParameter(p[_exp_par_name("F2_2", 1, 1, 2)], *this), UsedParameter(p[_exp_par_name("F2_2", 1, 2, 2)], *this) } } ,
                 { std::array{ UsedParameter(p[_exp_par_name("F2_2", 2, 0, 0)], *this), UsedParameter(p[_exp_par_name("F2_2", 2, 1, 0)], *this), UsedParameter(p[_exp_par_name("F2_2", 2, 2, 0)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F2_2", 2, 0, 1)], *this), UsedParameter(p[_exp_par_name("F2_2", 2, 1, 1)], *this), UsedParameter(p[_exp_par_name("F2_2", 2, 2, 1)], *this) } ,
                   std::array{ UsedParameter(p[_exp_par_name("F2_2", 2, 0, 2)], *this), UsedParameter(p[_exp_par_name("F2_2", 2, 1, 2)], *this), UsedParameter(p[_exp_par_name("F2_2", 2, 2, 2)], *this) } }
        }},
        _traits(HKVT2025FormFactorTraits<Process_, PToPP>(p)),
        _mP1(_traits.m_P1),
        _mP2(_traits.m_P2),
        _mP3(_traits.m_P3),
        _numWaves(Process_::num_waves),
        opt_I(o, options, "I"_ok),
        opt_L(o, options, "L"_ok),
        opt_C(o, "C"_ok, { "+-", "+0", "00" }),
        opt_int_points(o, "integration-points"_ok, {"256", "512", "1024", "2048", "4096", "8192", "16384"}, "4096"),
        _int_points(destringify<int>(opt_int_points.value())),
        scattering_amplitudes(ScatteringAmplitudeFactory<PPToPP>::create("pipi->pipi::" + o.get("scattering-amplitudes"_ok, "HKvT2025"), p, o)),
        _charge(_traits.charge_map.at(opt_C.value()))
    {
        _PartialWave_switch[0] = (opt_L.value() && PartialWave::S);
        _PartialWave_switch[1] = (opt_L.value() && PartialWave::P);
        _PartialWave_switch[2] = (opt_L.value() && PartialWave::D);
        _PartialWave_switch[3] = (opt_L.value() && PartialWave::F);

        _Isospin_switch[0] = (opt_I.value() && Isospin::zero) | (opt_I.value() && Isospin::onehalf);
        _Isospin_switch[1] = (opt_I.value() && Isospin::one) | (opt_I.value() && Isospin::threehalves);
    }

    template<typename Process_>
    HKVT2025FormFactors<Process_, PToPP>::~HKVT2025FormFactors() = default;

    template<typename Process_>
    FormFactors<PToPP> *
    HKVT2025FormFactors<Process_, PToPP>::make(const Parameters & parameters, const Options & options)
    {
        return new HKVT2025FormFactors(parameters, options);
    }

    template<typename Process_>
    QualifiedName
    HKVT2025FormFactors<Process_, PToPP>::_exp_par_name(const std::string & ff_name, unsigned wave, unsigned z_order, unsigned y_order) const
    {
        return QualifiedName(stringify(Process_::label) + "::a^" + ff_name + "_" + stringify(wave) + "_" + stringify(z_order) + "_" + stringify(y_order) + "@HKvT2025");
    }

    template<typename Process_>
    double
    HKVT2025FormFactors<Process_, PToPP>::_phi(const int & l, const double & k2, const double & q2, const double & q2p, const double & chi,
                        const unsigned a, const unsigned b, const double N) const
    {
        // Note: The factor 1 / 2 / M_PI at the end originates from the change of q^2 -> z and the 1 / M_PI in front of the dispersive integral
        const double norm = std::sqrt(N / 256 / M_PI / M_PI / M_PI / chi / (2 * l + 1) / 2 / M_PI);

        const double z = _traits.calc_z(q2, q2p, _traits.q20);

        // kinematic_q2p, kinematic_q2m depend on s
        const double kinematic_q2p = power_of<2>(_mP1 + std::sqrt(k2));
        const double kinematic_q2m = power_of<2>(_mP1 - std::sqrt(k2));

        // set Q^2 to 0
        const double q2_term = 1 / (2 * (q2p + std::sqrt(q2p) * std::sqrt(q2p - q2)) - q2);
        const double lambda_q2_term = (kinematic_q2p - q2) * power_of<2>(std::sqrt(q2p - q2) + std::sqrt(q2p - kinematic_q2m));
        const double sqrtjac_q2 = std::sqrt(4 * (1 + z) * (_traits.q20 - q2p) / power_of<3>(z - 1));

        return norm * sqrtjac_q2 * pow(lambda_q2_term, (2 * l + 1 - 2 * a) * 0.25) * pow(q2_term, b * 0.5);
    }

    template<typename Process_>
    inline double
    HKVT2025FormFactors<Process_, PToPP>::_phi_g(const double & q2, const double & k2, const unsigned & l) const
    {
        if (l > 0)
        {
            return _phi(l, k2, q2, _traits.q2p_v, _traits.chi_1m_v, 0, 4, l * (l + 1) / 48.0);
        }
        else
        {
            return 1.0;
        }
    }

    template<typename Process_>
    inline double
    HKVT2025FormFactors<Process_, PToPP>::_phi_f(const double & q2, const double & k2, const unsigned & l) const
    {
        if (l > 0)
        {
            return _phi(l, k2, q2, _traits.q2p_a, _traits.chi_1p_a, 1, 4, l * (l + 1) / 3.0);
        }
        else
        {
            return 1.0;
        }
    }

    template<typename Process_>
    inline double
    HKVT2025FormFactors<Process_, PToPP>::_phi_F1(const double & q2, const double & k2, const unsigned & l) const
    {
        if (l > 0)
        {
            return _phi(l, k2, q2, _traits.q2p_a, _traits.chi_1p_a, 1, 5, 1.0 / 12.0);
        }
        else
        {
            return _phi(0, k2, q2, _traits.q2p_a, _traits.chi_1p_a, -1, 5, 1.0 / 12.0);
        }
    }

    template<typename Process_>
    inline double
    HKVT2025FormFactors<Process_, PToPP>::_phi_F2(const double & q2, const double & k2, const unsigned & l) const
    {
        if (l > 0)
        {
            return _phi(l, k2, q2, _traits.q2p_a, _traits.chi_0m_a, 0, 4, 1.0);
        }
        else
        {
            return _phi(0, k2, q2, _traits.q2p_a, _traits.chi_0m_a, 0, 4, 1.0);
        }
    }

    template <typename Process_>
    complex<double>
    HKVT2025FormFactors<Process_, PToPP>::g_tilde(const double & q2, const double & k2, const unsigned & l, const bool & iso) const
    {
        const double isospin_factor = (iso) ? std::sqrt(Process_::eta_2[l]) : std::sqrt(Process_::eta_1[l]);
        const IsospinRepresentation rep = (iso) ? Process_::rep_2 : Process_::rep_1;
        if ( isospin_factor == 0.0 )
                return 0.0;

        // Partial-wave independent pieces
        const double blaschke      = (power_of<2>(_traits.m_R_1m) < _traits.q2p_v ? _traits.calc_z(q2, _traits.q2p_v, power_of<2>(_traits.m_R_1m)) : 1.0);
        const complex<double> y = iso ? _traits.calc_y(complex<double>(k2), complex<double>(_traits.k2in1), complex<double>(_traits.k20)) : _traits.calc_y(complex<double>(k2), complex<double>(_traits.k2in0), complex<double>(_traits.k20));
        const double z             = _traits.calc_z(q2, _traits.q2p_v, _traits.q20);
        const auto   polynomials   = _traits.orthonormal_polynomials_v(z, k2);
        const auto   y_polynomials = _traits.threshold_improved_polynomials(y, l);

        complex<double> res = 0.0;

        const double phi = _phi_g(q2, k2, l);

        if (iso)
        {
            for (unsigned i = 0 ; i < 3 ; i++)
            {
                complex<double> tmpres = 0.0;
                for (unsigned j = 0 ; j < 3 ; j++)
                {
                    tmpres += double(_a_g_2[l][i][j]) * polynomials[j];
                }
                res += tmpres * y_polynomials[i];
            }
        }
        else
        {
            for (unsigned i = 0 ; i < 3 ; i++)
            {
                complex<double> tmpres = 0.0;
                for (unsigned j = 0 ; j < 3 ; j++)
                {
                    tmpres += double(_a_g_1[l][i][j]) * polynomials[j];
                }
                res += tmpres * y_polynomials[i];
            }
        }

        return scattering_amplitudes->isospin_breaking(k2, l, rep) * scattering_amplitudes->omnes_factor(k2, l, rep) * res / blaschke / phi / isospin_factor;
    }

    template <typename Process_>
    complex<double>
    HKVT2025FormFactors<Process_, PToPP>::v_perp(const double & q2, const double & k2, const unsigned & l, const bool & iso) const
    {
        const double isospin_factor = (iso) ? Process_::lambda_2[l] : Process_::lambda_1[l];
        if ( isospin_factor == 0.0 )
                return 0.0;

        const double kinpref = std::pow(_traits.kappa(k2, q2), l - 1); // Taken into account in matching to helicity FFs: / std::sqrt(k2);

        return isospin_factor * kinpref * g_tilde(q2, k2, l, iso);
    }

    template <typename Process_>
    complex<double>
    HKVT2025FormFactors<Process_, PToPP>::f_tilde(const double & q2, const double & k2, const unsigned & l, const bool & iso) const
    {
        const double isospin_factor = (iso) ? std::sqrt(Process_::eta_2[l]) : std::sqrt(Process_::eta_1[l]);
        const IsospinRepresentation rep = (iso) ? Process_::rep_2 : Process_::rep_1;
        if ( isospin_factor == 0.0 )
                return 0.0;

        // Partial-wave independent pieces
        const double blaschke      = (power_of<2>(_traits.m_R_1p) < _traits.q2p_a ? _traits.calc_z(q2, _traits.q2p_a, power_of<2>(_traits.m_R_1p)) : 1.0);
        const complex<double> y = iso ? _traits.calc_y(complex<double>(k2), complex<double>(_traits.k2in1), complex<double>(_traits.k20)) : _traits.calc_y(complex<double>(k2), complex<double>(_traits.k2in0), complex<double>(_traits.k20));
        const double z             = _traits.calc_z(q2, _traits.q2p_a, _traits.q20);
        const auto   polynomials   = _traits.orthonormal_polynomials_a(z, k2);
        const auto   y_polynomials = _traits.threshold_improved_polynomials(y, l);

        complex<double> res = 0.0;

        const double phi = _phi_f(q2, k2, l);

        if (iso)
        {
            for (unsigned i = 0 ; i < 3 ; i++)
            {
                complex<double> tmpres = 0.0;
                for (unsigned j = 0 ; j < 3 ; j++)
                {
                    tmpres += double(_a_f_2[l][i][j]) * polynomials[j];
                }
                res += tmpres * y_polynomials[i];
            }
        }
        else
        {
            for (unsigned i = 0 ; i < 3 ; i++)
            {
                complex<double> tmpres = 0.0;
                for (unsigned j = 0 ; j < 3 ; j++)
                {
                    tmpres += double(_a_f_1[l][i][j]) * polynomials[j];
                }
                res += tmpres * y_polynomials[i];
            }
        }

        return scattering_amplitudes->isospin_breaking(k2, l, rep) * scattering_amplitudes->omnes_factor(k2, l, rep) * res / blaschke / phi / isospin_factor;
    }

    template <typename Process_>
    complex<double>
    HKVT2025FormFactors<Process_, PToPP>::a_par(const double & q2, const double & k2, const unsigned & l, const bool & iso) const
    {
        const double isospin_factor = (iso) ? Process_::lambda_2[l] : Process_::lambda_1[l];
        if ( isospin_factor == 0.0 )
                return 0.0;

        // Note, the factor sqrt(k2) is removed when matching the helicity amplitudes
        const double kinpref = std::pow(_traits.kappa(k2, q2), l - 1); // / std::sqrt(k2);

        return isospin_factor * kinpref * f_tilde(q2, k2, l, iso);
    }

    template <typename Process_>
    complex<double>
    HKVT2025FormFactors<Process_, PToPP>::F1_tilde(const double & q2, const double & k2, const unsigned & l, const bool & iso) const
    {
        const double isospin_factor = (iso) ? std::sqrt(Process_::eta_2[l]) : std::sqrt(Process_::eta_1[l]);
        const IsospinRepresentation rep = (iso) ? Process_::rep_2 : Process_::rep_1;
        if ( isospin_factor == 0.0 )
                return 0.0;

        // Partial-wave independent pieces
        const double blaschke      = (power_of<2>(_traits.m_R_1p) < _traits.q2p_a ? _traits.calc_z(q2, _traits.q2p_a, power_of<2>(_traits.m_R_1p)) : 1.0);
        const complex<double> y = iso ? _traits.calc_y(complex<double>(k2), complex<double>(_traits.k2in1), complex<double>(_traits.k20)) : _traits.calc_y(complex<double>(k2), complex<double>(_traits.k2in0), complex<double>(_traits.k20));
        const double z             = _traits.calc_z(q2, _traits.q2p_a, _traits.q20);
        const auto   polynomials   = _traits.orthonormal_polynomials_a(z, k2);
        const auto   y_polynomials = _traits.threshold_improved_polynomials(y, l);

        complex<double> res = 0.0;

        const double phi = _phi_F1(q2, k2, l);

        if (iso)
        {
            for (unsigned i = 0 ; i < 3 ; i++)
            {
                complex<double> tmpres = 0.0;
                for (unsigned j = 0 ; j < 3 ; j++)
                {
                    tmpres += double(_a_F1_2[l][i][j]) * polynomials[j];
                }
                res += tmpres * y_polynomials[i];
            }
        }
        else
        {
            for (unsigned i = 0 ; i < 3 ; i++)
            {
                complex<double> tmpres = 0.0;
                for (unsigned j = 0 ; j < 3 ; j++)
                {
                    tmpres += double(_a_F1_1[l][i][j]) * polynomials[j];
                }
                res += tmpres * y_polynomials[i];
            }
        }

        return scattering_amplitudes->isospin_breaking(k2, l, rep) * scattering_amplitudes->omnes_factor(k2, l, rep) * res / blaschke / phi / isospin_factor;
    }

    template <typename Process_>
    complex<double>
    HKVT2025FormFactors<Process_, PToPP>::a_0(const double & q2, const double & k2, const unsigned & l, const bool & iso) const
    {
        const double isospin_factor = (iso) ? Process_::lambda_2[l] : Process_::lambda_1[l];
        if ( isospin_factor == 0.0 )
                return 0.0;

        // Note, we absorb one power of sqrt(lambda_q3) by reducing the power of kappa and one when matching the helicity amplitudes
        const double kinpref = (l == 0) ? 1.0 : std::pow(_traits.kappa(k2, q2), l - 1) * std::sqrt((k2 - power_of<2>(_mP2 + _mP3)) * (k2 - power_of<2>(_mP2 - _mP3))) / k2;
        return isospin_factor * kinpref * F1_tilde(q2, k2, l, iso);
    }

    template <typename Process_>
    complex<double>
    HKVT2025FormFactors<Process_, PToPP>::F2_tilde(const double & q2, const double & k2, const unsigned & l, const bool & iso) const
    {
        const double isospin_factor = (iso) ? std::sqrt(Process_::eta_2[l]) : std::sqrt(Process_::eta_1[l]);
        const IsospinRepresentation rep = (iso) ? Process_::rep_2 : Process_::rep_1;
        if ( isospin_factor == 0.0 )
                return 0.0;

        // Partial-wave independent pieces
        const double blaschke      = (power_of<2>(_traits.m_R_0m) < _traits.q2p_a ? _traits.calc_z(q2, _traits.q2p_a, power_of<2>(_traits.m_R_0m)) : 1.0);
        const complex<double> y = iso ? _traits.calc_y(complex<double>(k2), complex<double>(_traits.k2in1), complex<double>(_traits.k20)) : _traits.calc_y(complex<double>(k2), complex<double>(_traits.k2in0), complex<double>(_traits.k20));
        const double z             = _traits.calc_z(q2, _traits.q2p_a, _traits.q20);
        const auto   polynomials   = _traits.orthonormal_polynomials_a(z, k2);
        const auto   y_polynomials = _traits.threshold_improved_polynomials(y, l);

        complex<double> res = 0.0;

        const double phi = _phi_F2(q2, k2, l);

        if (iso)
        {
            for (unsigned i = 0 ; i < 3 ; i++)
            {
                complex<double> tmpres = 0.0;
                for (unsigned j = 0 ; j < 3 ; j++)
                {
                    tmpres += double(_a_F2_2[l][i][j]) * polynomials[j];
                }
                res += tmpres * y_polynomials[i];
            }
        }
        else
        {
            for (unsigned i = 0 ; i < 3 ; i++)
            {
                complex<double> tmpres = 0.0;
                for (unsigned j = 0 ; j < 3 ; j++)
                {
                    tmpres += double(_a_F2_1[l][i][j]) * polynomials[j];
                }
                res += tmpres * y_polynomials[i];
            }
        }

        return scattering_amplitudes->isospin_breaking(k2, l, rep) * scattering_amplitudes->omnes_factor(k2, l, rep) * res / blaschke / phi / isospin_factor;
    }

    template <typename Process_>
    complex<double>
    HKVT2025FormFactors<Process_, PToPP>::a_t(const double & q2, const double & k2, const unsigned & l, const bool & iso) const
    {
        const double isospin_factor = (iso) ? Process_::lambda_2[l] : Process_::lambda_1[l];
        if ( isospin_factor == 0.0 )
                return 0.0;

        const double kinpref = (l == 0) ? 1.0 : std::pow(_traits.kappa(k2, q2), l);

        return isospin_factor * kinpref * F2_tilde(q2, k2, l, iso);
    }

    template <typename Process_>
    double
    HKVT2025FormFactors<Process_, PToPP>::_unitarity_integrand_0m(const double & k2, const unsigned & l, const bool & iso) const
    {
        const double isospin_factor = (iso) ? std::sqrt(Process_::eta_2[l]) : std::sqrt(Process_::eta_1[l]);
        const IsospinRepresentation rep = (iso) ? Process_::rep_2 : Process_::rep_1;
        if ( isospin_factor == 0.0 )
                return 0.0;

        double res = 0.0;
        double kinpref = std::pow(std::sqrt((k2 - power_of<2>(_mP2 + _mP3))*(k2 - power_of<2>(_mP2 - _mP3))) / k2, 2 * l + 1);
        const complex<double> y = iso ? _traits.calc_y(complex<double>(k2), complex<double>(_traits.k2in1), complex<double>(_traits.k20)) : _traits.calc_y(complex<double>(k2), complex<double>(_traits.k2in0), complex<double>(_traits.k20));
        const std::array<complex<double>, 3> y_polynomials  = _traits.threshold_improved_polynomials(y, l);

        for (unsigned i = 0 ; i < 3 ; i++)
        {
            if (!iso)
            {
                res += std::norm(double(_a_F2_1[l][0][i]) * y_polynomials[0] + double(_a_F2_1[l][1][i]) * y_polynomials[1] + double(_a_F2_1[l][2][i]) * y_polynomials[2]);
            }
            else
            {
                res += std::norm(double(_a_F2_2[l][0][i]) * y_polynomials[0] + double(_a_F2_2[l][1][i]) * y_polynomials[1] + double(_a_F2_2[l][2][i]) * y_polynomials[2]);
            }
        }

        return res * kinpref * std::norm(scattering_amplitudes->isospin_breaking(k2, l, rep) * scattering_amplitudes->omnes_factor(k2, l, rep));
    }

    template <typename Process_>
    double
    HKVT2025FormFactors<Process_, PToPP>::saturation_0m_a() const
    {
        std::function<double(const double &)> integrand = [&] (const double & t)
        {
            double contrib = 0.0;
            for (unsigned l = 0 ; l < 3 ; l++)
            {
                contrib += _unitarity_integrand_0m(1.0 / t, l, false) + _unitarity_integrand_0m(1.0 / t, l, true);
            }
            return contrib / power_of<2>(t);
        };

        double res = integrate1D(integrand, _int_points, 1e-5, 1.0 / power_of<2>(_mP2 + _mP3));

        return res;
    }

    template <typename Process_>
    double
    HKVT2025FormFactors<Process_, PToPP>::_unitarity_integrand_1p(const double & k2, const unsigned & l, const bool & iso) const
    {
        const double isospin_factor = (iso) ? std::sqrt(Process_::eta_2[l]) : std::sqrt(Process_::eta_1[l]);
        const IsospinRepresentation rep = (iso) ? Process_::rep_2 : Process_::rep_1;
        if ( isospin_factor == 0.0 )
                return 0.0;

        double res_f = 0.0;
        double res_F1 = 0.0;
        double kinpref_f = k2 * std::pow(std::sqrt((k2 - power_of<2>(_mP2 + _mP3))*(k2 - power_of<2>(_mP2 - _mP3))) / k2, 2 * l + 1);
        double kinpref_F1 = std::pow(std::sqrt((k2 - power_of<2>(_mP2 + _mP3))*(k2 - power_of<2>(_mP2 - _mP3))) / k2, 2 * l + 1);
        const complex<double> y = iso ? _traits.calc_y(complex<double>(k2), complex<double>(_traits.k2in1), complex<double>(_traits.k20)) : _traits.calc_y(complex<double>(k2), complex<double>(_traits.k2in0), complex<double>(_traits.k20));
        const std::array<complex<double>, 3> y_polynomials  = _traits.threshold_improved_polynomials(y, l);

        for (unsigned i = 0 ; i < 3 ; i++)
        {
            if (!iso)
            {
                res_f += std::norm(double(_a_f_1[l][0][i]) * y_polynomials[0] + double(_a_f_1[l][1][i]) * y_polynomials[1] + double(_a_f_1[l][2][i]) * y_polynomials[2]);
                res_F1 += std::norm(double(_a_F1_1[l][0][i]) * y_polynomials[0] + double(_a_F1_1[l][1][i]) * y_polynomials[1] + double(_a_F1_1[l][2][i]) * y_polynomials[2]);
            }
            else
            {
                res_f += std::norm(double(_a_f_2[l][0][i]) * y_polynomials[0] + double(_a_f_2[l][1][i]) * y_polynomials[1] + double(_a_f_2[l][2][i]) * y_polynomials[2]);
                res_F1 += std::norm(double(_a_F1_2[l][0][i]) * y_polynomials[0] + double(_a_F1_2[l][1][i]) * y_polynomials[1] + double(_a_F1_2[l][2][i]) * y_polynomials[2]);
            }
        }

        return (res_f * kinpref_f + res_F1 * kinpref_F1) * std::norm(scattering_amplitudes->isospin_breaking(k2, l, rep) * scattering_amplitudes->omnes_factor(k2, l, rep));
    }

    template <typename Process_>
    double
    HKVT2025FormFactors<Process_, PToPP>::saturation_1p_a() const
    {
        std::function<double(const double &)> integrand = [&] (const double & t)
        {
            double contrib = 0.0;
            for (unsigned l = 0 ; l < 3 ; l++)
            {
                contrib += _unitarity_integrand_1p(1.0 / t, l, false) + _unitarity_integrand_1p(1.0 / t, l, true);
            }
            return contrib / power_of<2>(t);
        };

        double res = integrate1D(integrand, _int_points, 1e-5, 1.0 / power_of<2>(_mP2 + _mP3));

        return res;
    }

    template <typename Process_>
    double
    HKVT2025FormFactors<Process_, PToPP>::_unitarity_integrand_1m(const double & k2, const unsigned & l, const bool & iso) const
    {
        const double isospin_factor = (iso) ? std::sqrt(Process_::eta_2[l]) : std::sqrt(Process_::eta_1[l]);
        const IsospinRepresentation rep = (iso) ? Process_::rep_2 : Process_::rep_1;
        if ( isospin_factor == 0.0 )
                return 0.0;

        double res = 0.0;
        double kinpref = k2 * std::pow(std::sqrt((k2 - power_of<2>(_mP2 + _mP3))*(k2 - power_of<2>(_mP2 - _mP3))) / k2, 2 * l + 1);
        const complex<double> y = iso ? _traits.calc_y(complex<double>(k2), complex<double>(_traits.k2in1), complex<double>(_traits.k20)) : _traits.calc_y(complex<double>(k2), complex<double>(_traits.k2in0), complex<double>(_traits.k20));
        const std::array<complex<double>, 3> y_polynomials  = _traits.threshold_improved_polynomials(y, l);

        for (unsigned i = 0 ; i < 3 ; i++)
        {
            if (!iso)
            {
                res += std::norm(double(_a_g_1[l][0][i]) * y_polynomials[0] + double(_a_g_1[l][1][i]) * y_polynomials[1] + double(_a_g_1[l][2][i]) * y_polynomials[2]);
            }
            else
            {
                res += std::norm(double(_a_g_2[l][0][i]) * y_polynomials[0] + double(_a_g_2[l][1][i]) * y_polynomials[1] + double(_a_g_2[l][2][i]) * y_polynomials[2]);
            }
        }

        return res * kinpref * std::norm(scattering_amplitudes->isospin_breaking(k2, l, rep) * scattering_amplitudes->omnes_factor(k2, l, rep));
    }

    template <typename Process_>
    double
    HKVT2025FormFactors<Process_, PToPP>::saturation_1m_v() const
    {
        std::function<double(const double &)> integrand = [&] (const double & t)
        {
            double contrib = 0.0;
            for (unsigned l = 0 ; l < 3 ; l++)
            {
                contrib += _unitarity_integrand_1m(1.0 / t, l, false) + _unitarity_integrand_1m(1.0 / t, l, true);
            }
            return contrib / power_of<2>(t);
        };

        double res = integrate1D(integrand, _int_points, 1e-5, 1.0 / power_of<2>(_mP2 + _mP3));

        return res;
    }

    template <typename Process_>
    complex<double>
    HKVT2025FormFactors<Process_, PToPP>::f_perp(const double & q2, const double & k2, const double & z) const
    {
        std::array<complex<double>, 4> partial_waves = f_perp(q2, k2);

        return partial_waves[1] + 3.0 * z * partial_waves[2];
    }

    template <typename Process_>
    complex<double>
    HKVT2025FormFactors<Process_, PToPP>::f_para(const double & q2, const double & k2, const double & z) const
    {
        std::array<complex<double>, 4> partial_waves = f_para(q2, k2);

        return partial_waves[1] + 3.0 * z * partial_waves[2];
    }

    template <typename Process_>
    complex<double>
    HKVT2025FormFactors<Process_, PToPP>::f_long(const double & q2, const double & k2, const double & z) const
    {
        std::array<complex<double>, 4> partial_waves = f_long(q2, k2);

        return partial_waves[0] + z * partial_waves[1] + 0.5 * (3.0 * z * z - 1.0) * partial_waves[2];
    }

    template <typename Process_>
    complex<double>
    HKVT2025FormFactors<Process_, PToPP>::f_time(const double & q2, const double & k2, const double & z) const
    {
        std::array<complex<double>, 4> partial_waves = f_time(q2, k2);

        return partial_waves[0] + z * partial_waves[1] + 0.5 * (3.0 * z * z - 1.0) * partial_waves[2];
    }

    template <typename Process_>
    std::array<complex<double>, 4>
    HKVT2025FormFactors<Process_, PToPP>::f_perp(const double & q2, const double & k2) const
    {
        // Factor std::sqrt(k2) already taken into account
        const double lam = _traits.lamB(k2, q2);
        if (lam <= 0.0) return {0.0, 0.0, 0.0, 0.0};

        const double pref = -std::sqrt(lam) / 4.0;
        std::array<complex<double>, 4> res = {0.0, 0.0, 0.0, 0.0};

        for ( unsigned l = 1 ; l < 3 ; l++ )
        {
            if ( !_PartialWave_switch[l] ) continue;

            if ( _Isospin_switch[0] )
            {
                res[l] += Process_::IsoToPhys[_charge][0] * v_perp(q2, k2, l, false);
            }

            if ( _Isospin_switch[1] )
            {
                res[l] += Process_::IsoToPhys[_charge][1] * v_perp(q2, k2, l, true);
            }

            res[l] *= pref / std::sqrt(2.0 * l + 1.0);
        }

        return res;
    }

    template <typename Process_>
    std::array<complex<double>, 4>
    HKVT2025FormFactors<Process_, PToPP>::f_para(const double & q2, const double & k2) const
    {
        const double pref = 1.0; // Already taken into account: std::sqrt(k2)
        std::array<complex<double>, 4> res = {0.0, 0.0, 0.0, 0.0};

        for ( unsigned l = 1 ; l < 3 ; l++ )
        {
            if ( !_PartialWave_switch[l] ) continue;

            if ( _Isospin_switch[0] )
            {
                res[l] += Process_::IsoToPhys[_charge][0] * a_par(q2, k2, l, false);
            }

            if ( _Isospin_switch[1] )
            {
                res[l] += Process_::IsoToPhys[_charge][1] * a_par(q2, k2, l, true);
            }

            res[l] *= pref / std::sqrt(2.0 * l + 1.0);
        }

        return res;
    }

    template <typename Process_>
    std::array<complex<double>, 4>
    HKVT2025FormFactors<Process_, PToPP>::f_long(const double & q2, const double & k2) const
    {
        const double lam = _traits.lamB(k2, q2);
        const double prefS = (lam > 0.0) ? std::sqrt(lam / q2) / 2.0 : 0.0;
        const double pref = 1.0 / 2.0 / std::sqrt(q2);
        std::array<complex<double>, 4> res = {0.0, 0.0, 0.0, 0.0};

        if ( _PartialWave_switch[0] ) {
            if ( _Isospin_switch[0] )
            {
                res[0] += Process_::IsoToPhys[_charge][0] * a_0(q2, k2, 0, false);
            }

            if ( _Isospin_switch[1]  )
            {
                res[0] += Process_::IsoToPhys[_charge][1] * a_0(q2, k2, 0, true);
            }

            res[0] *= prefS;
        }

        for ( unsigned l = 1 ; l < 3 ; l++ )
        {
            if ( !_PartialWave_switch[l] ) continue;

            if ( _Isospin_switch[0] )
            {
                res[l] += Process_::IsoToPhys[_charge][0] * a_0(q2, k2, l, false);
            }

            if ( _Isospin_switch[1] )
            {
                res[l] += Process_::IsoToPhys[_charge][1] * a_0(q2, k2, l, true);
            }

            res[l] *= pref / std::sqrt(2.0 * l + 1.0);
        }

        return res;
    }

    template <typename Process_>
    std::array<complex<double>, 4>
    HKVT2025FormFactors<Process_, PToPP>::f_time(const double & q2, const double & k2) const
    {
        const double pref = -std::sqrt(1.0 / q2);
        std::array<complex<double>, 4> res = {0.0, 0.0, 0.0, 0.0};

        for ( unsigned l = 0 ; l < 3 ; l++ )
        {
            if ( !_PartialWave_switch[l] ) continue;

            if ( _Isospin_switch[0] )
            {
                res[l] += Process_::IsoToPhys[_charge][0] * a_t(q2, k2, l, false);
            }

            if( _Isospin_switch[1] )
            {
                res[l] += Process_::IsoToPhys[_charge][1] * a_t(q2, k2, l, true);
            }

            res[l] *= pref / std::sqrt(2.0 * l + 1.0);
        }

        return res;
    }

    template <typename Process_>
    Diagnostics
    HKVT2025FormFactors<Process_, PToPP>::diagnostics() const
    {
        Diagnostics results;

        results.add({ _traits.calc_y(4 * 0.135 * 0.135,  _traits.k2in1, _traits.k20), "y(k2 = 4*0.135^2)" });
        results.add({ _traits.calc_y(0.1,                _traits.k2in1, _traits.k20), "y(k2 = 0.1)"       });

        results.add({ _traits.calc_z(0.0,  _traits.q2p_a, _traits.q20), "z_a(q2 =  0)" });
        results.add({ _traits.calc_z(0.0,  _traits.q2p_v, _traits.q20), "z_v(q2 =  0)" });
        results.add({ _traits.calc_z(10.0, _traits.q2p_a, _traits.q20), "z_a(q2 = 10)" });
        results.add({ _traits.calc_z(10.0, _traits.q2p_v, _traits.q20), "z_v(q2 = 10)" });

        {
            const auto & [p0, p1, p2] = _traits.orthonormal_polynomials_v(0.0, 0.1);
            results.add({ p0,              "p_0(z = 0.0, k2 = 0.1)" });
            results.add({ p1,              "p_1(z = 0.0, k2 = 0.1)" });
            results.add({ p2,              "p_2(z = 0.0, k2 = 0.1)" });
        }

        {
            const auto & [p0, p1, p2] = _traits.orthonormal_polynomials_v(_traits.calc_z(10.0, _traits.q2p_v, _traits.q20), 0.1);
            results.add({ p0,              "p_0(z = z(q2 = 10, k2 = 0.1))" });
            results.add({ p1,              "p_1(z = z(q2 = 10, k2 = 0.1))" });
            results.add({ p2,              "p_2(z = z(q2 = 10, k2 = 0.1))" });
        }

        {
            const auto & [p0, p1, p2] = _traits.threshold_improved_polynomials(_traits.calc_y(0.5, _traits.k2in1, _traits.k20), 0);
            results.add({ p0,              "p_0(y = y(k2 = 0.5), 0)" });
            results.add({ p1,              "p_1(y = y(k2 = 0.5), 0)" });
            results.add({ p2,              "p_2(y = y(k2 = 0.5), 0)" });
        }

        {
            const auto & [p0, p1, p2] = _traits.threshold_improved_polynomials(_traits.calc_y(0.5, _traits.k2in1, _traits.k20), 1);
            results.add({ p0,              "p_0(y = y(k2 = 0.5), 1)" });
            results.add({ p1,              "p_1(y = y(k2 = 0.5), 1)" });
            results.add({ p2,              "p_2(y = y(k2 = 0.5), 1)" });
        }

        {
            const auto & [p0, p1, p2] = _traits.threshold_improved_polynomials(_traits.calc_y(0.5, _traits.k2in1, _traits.k20), 2);
            results.add({ p0,              "p_0(y = y(k2 = 0.5), 2)" });
            results.add({ p1,              "p_1(y = y(k2 = 0.5), 2)" });
            results.add({ p2,              "p_2(y = y(k2 = 0.5), 2)" });
        }

        {
            results.add({ _phi_g(-2.0, 0.1, 1),     "phi_g(z = z(q2 = -2.0), y = y(k2 = 0.1), l = 1)" });
            results.add({ _phi_g( 1.0, 0.2, 1),     "phi_g(z = z(q2 =  1.0), y = y(k2 = 0.2), l = 1)" });
            results.add({ _phi_g( 4.0, 0.1, 2),     "phi_g(z = z(q2 =  4.0), y = y(k2 = 0.1), l = 2)" });

            results.add({ _phi_f(-2.0, 0.1, 1),     "phi_f(z = z(q2 = -2.0), y = y(k2 = 0.1), l = 1)" });
            results.add({ _phi_f( 1.0, 0.2, 1),     "phi_f(z = z(q2 =  1.0), y = y(k2 = 0.2), l = 1)" });
            results.add({ _phi_f( 4.0, 0.1, 2),     "phi_f(z = z(q2 =  4.0), y = y(k2 = 0.1), l = 2)" });

            results.add({ _phi_F1(-2.0, 0.1, 1),    "phi_F1(z = z(q2 = -2.0), y = y(k2 = 0.1), l = 1)" });
            results.add({ _phi_F1( 1.0, 0.2, 1),    "phi_F1(z = z(q2 =  1.0), y = y(k2 = 0.2), l = 1)" });
            results.add({ _phi_F1( 4.0, 0.1, 2),    "phi_F1(z = z(q2 =  4.0), y = y(k2 = 0.1), l = 2)" });
            results.add({ _phi_F1( 4.0, 0.1, 0),    "phi_F1(z = z(q2 =  4.0), y = y(k2 = 0.1), l = 0)" });

            results.add({ _phi_F2(-2.0, 0.1, 1),    "phi_F2(z = z(q2 = -2.0), y = y(k2 = 0.1), l = 1)" });
            results.add({ _phi_F2( 1.0, 0.2, 1),    "phi_F2(z = z(q2 =  1.0), y = y(k2 = 0.2), l = 1)" });
            results.add({ _phi_F2( 4.0, 0.1, 2),    "phi_F2(z = z(q2 =  4.0), y = y(k2 = 0.1), l = 2)" });
            results.add({ _phi_F2( 4.0, 0.1, 0),    "phi_F2(z = z(q2 =  4.0), y = y(k2 = 0.1), l = 0)" });
        }

        return results;
    }


    template<typename Process_>
    const std::set<ReferenceName> HKVT2025FormFactors<Process_, PToPP>::references
    {
        "HKvT:2025A"_rn
    };

    template<typename Process_>
    const std::vector<OptionSpecification> HKVT2025FormFactors<Process_, PToPP>::options
    {
        { "I"_ok, { "0|1", "1/2|3/2" }, "0|1" }, // We can handle even or odd cases here, but no mixture
        { "C"_ok, { "+-", "00", "+0" }, "+-" },
        { "L"_ok, { "S|P|D" }, "S|P|D" } // We currently do not support F waves here, as the corresponding y-polynomials are unknown
    };

    template<typename Process_>
    std::vector<OptionSpecification>::const_iterator
    HKVT2025FormFactors<Process_, PToPP>::begin_options()
    {
        return options.cbegin();
    }

    template<typename Process_>
    std::vector<OptionSpecification>::const_iterator
    HKVT2025FormFactors<Process_, PToPP>::end_options()
    {
        return options.cend();
    }
}
#endif
