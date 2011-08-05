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
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/wrapped_forward_iterator-impl.hh>

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

        double min, central, max;
    };

    namespace templates
    {
        /* [arxiv:1107.3753] */
        static const GaussianConstraintTemplate Bu_to_Kplus_dimuon_br_1_to_6_cdf_2011
        {
            "B->Kll::BR@LargeRecoil",
            Kinematics
            {
                { "s_min", 1.0 },
                { "s_max", 6.0 }
            },
            Options
            {
                { "q", "u"  },
                { "l", "mu" }
            },
            1.1907e-7, 1.4100e-7, 1.6293e-7
        };
        static const GaussianConstraintTemplate Bd_to_Kstarzero_dimuon_br_1_to_6_cdf_2011
        {
            "B->K^*ll::BR@LargeRecoil",
            Kinematics
            {
                { "s_min", 1.0 },
                { "s_max", 6.0 }
            },
            Options
            {
                { "q", "d"  },
                { "l", "mu" }
            },
            1.0023e-7, 1.4200e-7, 1.8377e-7
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

        LogLikelihoodBlockPtr block = LogLikelihoodBlock::Gaussian(cache, observable, t.min, t.central, t.max);

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
            std::make_pair("B_d->K^*0mu^+mu^-::BR[1.0,6.0]@CDF-2011", std::bind(&make_gaussian_constraint, std::placeholders::_1,
                        std::placeholders::_2, templates::Bd_to_Kstarzero_dimuon_br_1_to_6_cdf_2011)),
            std::make_pair("B_u->K^+mu^+mu^-::BR[1.0,6.0]@CDF-2011", std::bind(&make_gaussian_constraint, std::placeholders::_1,
                        std::placeholders::_2, templates::Bu_to_Kplus_dimuon_br_1_to_6_cdf_2011)),
        };

        auto f = factories.find(name);
        if (f == factories.end())
            throw UnknownConstraintError(name);

        return f->second(f->first, options);
    }
}
