/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
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

#include <src/constraint.hh>
#include <src/utils/exception.hh>
#include <src/utils/log_likelihood.hh>
#include <src/utils/observable_set.hh>
#include <src/utils/power_of.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/wrapped_forward_iterator-impl.hh>

#include <cmath>
#include <map>
#include <vector>

namespace eos
{
    UnknownConstraintError::UnknownConstraintError(const std::string & name) :
        Exception("Constraint '" + name + "' is unknown")
    {
    }

    struct GaussianConstraintTemplate
    {
        std::string observable;

        Kinematics kinematics;

        Options options;

        double central, sigma_hi_stat, sigma_lo_stat, sigma_hi_sys, sigma_lo_sys;
    };

    namespace templates
    {
        /*
         * LHCb Collaboration
         *
         * Data taken from talk by M. Patel at EPS-HEP 2011
         */
        // BR in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_1_to_6_LHCb_2011
        {
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.9500e-7, +0.3000e-7, -0.3000e-7, +0.1000e-7, -0.1000e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_LHCb_2011
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.0738e-7, +0.1830e-7, -0.1830e-7, +0.0546e-7, -0.0546e-7
        };
        // BR in [16.00, 19.00] (in the preliminary results, LHCb only integrated up to exactly 19.00 GeV^2!)
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_LHCb_2011
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.4400e-7, +0.2400e-7, -0.2400e-7, +0.0800e-7, -0.0800e-7
        };
        // A_FB in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_LHCb_2011
        {
            "B->K^*ll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.10, +0.14, -0.14, +0.05, -0.05
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_LHCb_2011
        {
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.50, +0.09, -0.06, +0.03, -0.03
        };
        // A_FB in [16.00, 19.00] (in the preliminary results, LHCb only integrated up to exactly 19.00 GeV^2!)
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_LHCb_2011
        {
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.10, +0.13, -0.13, +0.06, -0.06
        };
        // F_L in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_1_to_6_LHCb_2011
        {
            "B->K^*ll::F_Lavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.57, +0.11, -0.10, +0.03, -0.03
        };
        // F_L in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_LHCb_2011
        {
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.33, +0.11, -0.08, +0.04, -0.04
        };
        // F_L in [16.00, 19.00] (in the preliminary results, LHCb only integrated up to exactly 19.00 GeV^2!)
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_LHCb_2011
        {
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.28, +0.10, -0.09, +0.04, -0.04
        };
    }

    Constraint
    make_gaussian_constraint(const std::string & name, const Options & options, const GaussianConstraintTemplate & t)
    {
        Parameters parameters(Parameters::Defaults());
        ObservableCache cache(parameters);

        ObservablePtr observable = Observable::make(t.observable, parameters, t.kinematics, t.options + options);
        if (! observable.get())
            throw InternalError("make_gaussian_constraint: " + name + ": '" + t.observable + "' is not a valid observable name");

        double min = 0.0, max = 0.0;
        if ("asymmetric+quadratic" == options.get("uncertainty", "asymmetric+quadratic"))
        {
            min = t.central - std::sqrt(power_of<2>(t.sigma_lo_stat) + power_of<2>(t.sigma_lo_sys));
            max = t.central - std::sqrt(power_of<2>(t.sigma_hi_stat) + power_of<2>(t.sigma_hi_sys));
        }

        LogLikelihoodBlockPtr block = LogLikelihoodBlock::Gaussian(cache, observable, min, t.central, max);

        return Constraint(name, { observable }, { block });
    }

    /* Constraint */
    template class WrappedForwardIterator<Constraint::BlockIteratorTag, LogLikelihoodBlockPtr>;
    template class WrappedForwardIterator<Constraint::ObservableIteratorTag, ObservablePtr>;

    template <>
    struct Implementation<Constraint>
    {
        std::string name;

        ObservableSet observables;

        std::vector<LogLikelihoodBlockPtr> blocks;

        Implementation(const std::string & name,
                const std::initializer_list<ObservablePtr> & observables,
                const std::initializer_list<LogLikelihoodBlockPtr> & blocks) :
            name(name),
            blocks(blocks)
        {
            for (auto o = observables.begin(), o_end = observables.end() ; o != o_end ; ++o)
            {
                this->observables.add(*o);
            }
        }
    };

    Constraint::Constraint(const std::string & name,
                const std::initializer_list<ObservablePtr> & observables,
                const std::initializer_list<LogLikelihoodBlockPtr> & blocks) :
        PrivateImplementationPattern<Constraint>(new Implementation<Constraint>(name, observables, blocks))
    {
    }

    Constraint::~Constraint()
    {
    }

    std::string
    Constraint::name() const
    {
        return _imp->name;
    }

    Constraint::BlockIterator
    Constraint::begin_blocks() const
    {
        return BlockIterator(_imp->blocks.begin());
    }

    Constraint::BlockIterator
    Constraint::end_blocks() const
    {
        return BlockIterator(_imp->blocks.end());
    }

    Constraint::ObservableIterator
    Constraint::begin_observables() const
    {
        return ObservableIterator(_imp->observables.begin());
    }

    Constraint::ObservableIterator
    Constraint::end_observables() const
    {
        return ObservableIterator(_imp->observables.end());
    }

    Constraint
    Constraint::make(const std::string & name, const Options & options)
    {
        typedef std::function<Constraint (const std::string &, const Options & options)> ConstraintFactory;
        static const std::map<std::string, ConstraintFactory> factories
        {
#ifdef MAKEGAUSS
#  undef MAKEGAUSS
#endif
#define MAKEGAUSS(n, t) \
            { \
                n, \
                std::bind(&make_gaussian_constraint, std::placeholders::_1, std::placeholders::_2, t) \
            }
            /* LHCb 2011 */
            MAKEGAUSS("B^0->K^*0mu^+mu^-::BR[1.00,6.00]@LHCb-2011",            templates::Bzero_to_Kstarzero_dimuon_BR_1_to_6_LHCb_2011),
            MAKEGAUSS("B^0->K^*0mu^+mu^-::BR[14.18,16.00]@LHCb-2011",          templates::Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_LHCb_2011),
            MAKEGAUSS("B^0->K^*0mu^+mu^-::BR[16.00,19.21]@LHCb-2011",          templates::Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_LHCb_2011),
            MAKEGAUSS("B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@LHCb-2011",          templates::Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_LHCb_2011),
            MAKEGAUSS("B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@LHCb-2011",        templates::Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_LHCb_2011),
            MAKEGAUSS("B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@LHCb-2011",        templates::Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_LHCb_2011),
            MAKEGAUSS("B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@LHCb-2011",           templates::Bzero_to_Kstarzero_dimuon_F_L_1_to_6_LHCb_2011),
            MAKEGAUSS("B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@LHCb-2011",         templates::Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_LHCb_2011),
            MAKEGAUSS("B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@LHCb-2011",         templates::Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_LHCb_2011),
#undef MAKEGAUSS
        };

        auto f = factories.find(name);
        if (f == factories.end())
            throw UnknownConstraintError(name);

        return f->second(f->first, options);
    }
}
