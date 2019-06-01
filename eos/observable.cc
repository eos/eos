/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017 Danny van Dyk
 * Copyright (c) 2011 Christian Wacker
 * Copyright (c) 2018, 2019 Ahmet Kokulu
 * Copyright (c) 2018, 2019 Nico Gubernari
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
#include <eos/b-decays/observables.hh>
#include <eos/rare-b-decays/observables.hh>
#include <eos/form-factors/observables.hh>
//#include <eos/utils/concrete_observable.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/observable_stub.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>
#include <eos/observable-impl.hh>

#include <algorithm>
#include <map>

namespace eos
{
    std::vector<ObservableSection>
    make_observable_sections()
    {
        return std::vector<ObservableSection>({
            make_b_decays_section(),
            make_rare_b_decays_section(),
            make_form_factors_section()
        });
    }

    std::map<QualifiedName, ObservableEntryPtr>
    make_observable_entries()
    {
        std::map<QualifiedName, ObservableEntryPtr> observable_entries;

        for (auto && section : make_observable_sections())
        {
            for (auto && group : section)
            {
                observable_entries.insert(group.begin(), group.end());
            }
        }
        return observable_entries;
    }

    ObservablePtr
    Observable::make(const QualifiedName & name, const Parameters & parameters, const Kinematics & kinematics, const Options & _options)
    {
        static const std::map<QualifiedName, ObservableEntryPtr> observable_entries = std::move(make_observable_entries());

        // check if 'name' matches a simple observable
        {
            auto i = observable_entries.find(name);
            if (observable_entries.end() != i)
                return i->second->make(parameters, kinematics, name.options() + _options);
        }

        // check if 'name' matches a parameter
        if (name.options().empty())
        {
            auto i = std::find_if(parameters.begin(), parameters.end(), [&] (const Parameter & p) { return p.name() == name.str(); });
            if (parameters.end() != i)
            {
                return ObservablePtr(new ObservableStub(parameters, name));
            }
        }

        return ObservablePtr();
    }

    /* ObservableEntry */

    ObservableEntry::ObservableEntry()
    {
    }

    ObservableEntry::~ObservableEntry()
    {
    }

    std::ostream &
    ObservableEntry::insert(std::ostream & os) const
    {
        os << "<empty Observable description>" << std::endl;
        return os;
    }

    /* ObservableGroup */

    template <>
    struct WrappedForwardIteratorTraits<ObservableGroup::ObservableIteratorTag>
    {
        typedef std::map<QualifiedName, ObservableEntryPtr>::const_iterator UnderlyingIterator;
    };
    template class WrappedForwardIterator<ObservableGroup::ObservableIteratorTag, const std::pair<const QualifiedName, ObservableEntryPtr>>;

    ObservableGroup::ObservableGroup(Implementation<ObservableGroup> * imp) :
        PrivateImplementationPattern<ObservableGroup>(imp)
    {
    }

    ObservableGroup::~ObservableGroup() = default;

    ObservableGroup::ObservableIterator
    ObservableGroup::begin() const
    {
        return _imp->entries.begin();
    }

    ObservableGroup::ObservableIterator
    ObservableGroup::end() const
    {
        return _imp->entries.end();
    }

    const std::string &
    ObservableGroup::name() const
    {
        return _imp->name;
    }

    const std::string &
    ObservableGroup::description() const
    {
        return _imp->description;
    }

    /* ObservableSection */

    template <>
    struct WrappedForwardIteratorTraits<ObservableSection::GroupIteratorTag>
    {
        typedef std::vector<ObservableGroup>::const_iterator UnderlyingIterator;
    };
    template class WrappedForwardIterator<ObservableSection::GroupIteratorTag, const ObservableGroup &>;

    ObservableSection::ObservableSection(Implementation<ObservableSection> * imp) :
        PrivateImplementationPattern<ObservableSection>(imp)
    {
    }

    ObservableSection::~ObservableSection() = default;

    ObservableSection::GroupIterator
    ObservableSection::begin() const
    {
        return _imp->groups.begin();
    }

    ObservableSection::GroupIterator
    ObservableSection::end() const
    {
        return _imp->groups.end();
    }

    const std::string &
    ObservableSection::name() const
    {
        return _imp->name;
    }

    const std::string &
    ObservableSection::description() const
    {
        return _imp->description;
    }

    /* Observables */

    template <>
    struct WrappedForwardIteratorTraits<Observables::ObservableIteratorTag>
    {
        typedef std::map<QualifiedName, ObservableEntryPtr>::const_iterator UnderlyingIterator;
    };
    template class WrappedForwardIterator<Observables::ObservableIteratorTag, const std::pair<const QualifiedName, ObservableEntryPtr>>;

    template <>
    struct WrappedForwardIteratorTraits<Observables::SectionIteratorTag>
    {
        typedef std::vector<ObservableSection>::const_iterator UnderlyingIterator;
    };
    template class WrappedForwardIterator<Observables::SectionIteratorTag, const ObservableSection &>;

    template<>
    struct Implementation<Observables>
    {
        std::vector<ObservableSection> observable_sections;

        std::map<QualifiedName, ObservableEntryPtr> observable_entries;

        Implementation() :
            observable_sections(make_observable_sections())
        {
            for (auto && section : observable_sections)
            {
                for (auto && group : section)
                {
                    observable_entries.insert(group.begin(), group.end());
                }
            }
        }
    };

    Observables::Observables() :
        PrivateImplementationPattern<Observables>(new Implementation<Observables>())
    {
    }

    Observables::~Observables()
    {
    }

    ObservableEntryPtr
    Observables::operator[] (const QualifiedName & qn) const
    {
        auto i = _imp->observable_entries.find(qn);
        if (i != _imp->observable_entries.end())
            return i->second;

        return ObservableEntryPtr(nullptr);
    }

    Observables::ObservableIterator
    Observables::begin() const
    {
        return ObservableIterator(_imp->observable_entries.begin());
    }

    Observables::ObservableIterator
    Observables::end() const
    {
        return ObservableIterator(_imp->observable_entries.end());
    }

    Observables::SectionIterator
    Observables::begin_sections() const
    {
        return SectionIterator(_imp->observable_sections.begin());
    }

    Observables::SectionIterator
    Observables::end_sections() const
    {
        return SectionIterator(_imp->observable_sections.end());
    }
}
