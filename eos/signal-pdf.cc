/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015-2026 Danny van Dyk
 * Copyright (c) 2019      Ahmet Kokulu
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

#include <eos/b-decays/signal-pdfs.hh>
#include <eos/maths/power-of.hh>
#include <eos/rare-b-decays/signal-pdfs.hh>
#include <eos/signal-pdf-impl.hh>
#include <eos/signal-pdf.hh>
#include <eos/utils/concrete-signal-pdf.hh>
#include <eos/utils/density.hh>
#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/log.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <algorithm>
#include <functional>
#include <map>
#include <ostream>

namespace eos
{
    namespace impl
    {
        std::map<QualifiedName, SignalPDFEntryPtr> signal_pdf_entries;
    }

    SignalPDFEntries::SignalPDFEntries() :
        _entries(&impl::signal_pdf_entries)
    {
        std::vector<std::function<SignalPDFSection()>> section_makers = {
            make_b_decays_pdf_section,
            make_rare_b_decays_pdf_section,
        };

        for (const auto & section_maker : section_makers)
        {
            for (const auto & group : section_maker())
            {
                _entries->insert(group.begin(), group.end());
            }
        }

        // add test entries to the list of available signal PDFs, but avoid adding it via a group/section
        // 1D Legendre PDF
        {
            using std::literals::string_literals::operator""s;

            auto name_and_entry_pair = make_signal_pdf("TestLegendre1D::P(z)",
                                                       "PDF for testing purpose only: 1D PDF based on Legendre polynomials as a function of the variable $z$.",
                                                       Options{},
                                                       "TestLegendre1D::UnnormalizedPDF(z)",
                                                       std::make_tuple("z"s),
                                                       "TestLegendre1D::NormalizationPDF(z)",
                                                       std::make_tuple("z_min"s, "z_max"s));

            _entries->insert(name_and_entry_pair);
        }
    }

    SignalPDFEntries::~SignalPDFEntries() = default;

    void
    SignalPDFEntries::insert_or_assign(const QualifiedName & key, const std::shared_ptr<const SignalPDFEntry> & value)
    {
        auto result = _entries->insert_or_assign(key, value);

        if (! result.second)
        {
            Log::instance()->message("[SignalPDFEntries.insert_or_assign]", ll_warning) << "Entry for signal PDF " << key.str() << " has been replaced.";
        }
    }

    SignalPDFNameError::SignalPDFNameError(const std::string & name) :
        Exception("SignalPDF name '" + name + "' is malformed")
    {
    }

    SignalPDFEntry::SignalPDFEntry() {}

    SignalPDFEntry::~SignalPDFEntry() = default;

    std::ostream &
    SignalPDFEntry::insert(std::ostream & os) const
    {
        os << "<empty SignalPDF description>" << std::endl;
        return os;
    }

    SignalPDFPtr
    SignalPDF::make(const QualifiedName & name, const Parameters & parameters, const Kinematics & kinematics, const Options & _options)
    {
        const std::map<QualifiedName, SignalPDFEntryPtr> & signal_pdf_entries = SignalPDFEntries::instance()->entries();

        // check if 'name' matches any of the implemented Signal PDFs
        {
            auto i = signal_pdf_entries.find(name);
            if (signal_pdf_entries.end() != i)
            {
                return i->second->make(parameters, kinematics, name.options() + _options);
            }
        }

        return SignalPDFPtr();
    }

    /* SignalPDFSections */

    class SignalPDFSections : public InstantiationPolicy<SignalPDFSections, Singleton>
    {
        private:
            std::vector<SignalPDFSection> _sections;

            SignalPDFSections()
            {
                // ensure that the observable entries have been generated already
                auto entries = std::distance(SignalPDFEntries::instance()->entries().begin(), SignalPDFEntries::instance()->entries().end());
                Log::instance()->message("SignalPDFSections::SignalPDFSections()", ll_debug) << "Total number of registered signal PDFs: " << entries;

                _sections = std::vector<SignalPDFSection>({
                    make_b_decays_pdf_section(),
                    make_rare_b_decays_pdf_section(),
                });
            }

            ~SignalPDFSections() = default;

        public:
            friend class InstantiationPolicy<SignalPDFSections, Singleton>;

            const std::vector<SignalPDFSection> &
            sections() const
            {
                return _sections;
            }
    };

    /* SignalPDFGroup */

    template <> struct WrappedForwardIteratorTraits<SignalPDFGroup::SignalPDFIteratorTag>
    {
            using UnderlyingIterator = std::map<QualifiedName, SignalPDFEntryPtr>::const_iterator;
    };
    template class WrappedForwardIterator<SignalPDFGroup::SignalPDFIteratorTag, const std::pair<const QualifiedName, SignalPDFEntryPtr>>;

    SignalPDFGroup::SignalPDFGroup(Implementation<SignalPDFGroup> * imp) :
        PrivateImplementationPattern<SignalPDFGroup>(imp)
    {
    }

    SignalPDFGroup::~SignalPDFGroup() = default;

    SignalPDFGroup::SignalPDFIterator
    SignalPDFGroup::begin() const
    {
        return _imp->entries.begin();
    }

    SignalPDFGroup::SignalPDFIterator
    SignalPDFGroup::end() const
    {
        return _imp->entries.end();
    }

    const std::string &
    SignalPDFGroup::name() const
    {
        return _imp->name;
    }

    const std::string &
    SignalPDFGroup::description() const
    {
        return _imp->description;
    }

    /* SignalPDFSection */

    template <> struct WrappedForwardIteratorTraits<SignalPDFSection::GroupIteratorTag>
    {
            using UnderlyingIterator = std::vector<SignalPDFGroup>::const_iterator;
    };
    template class WrappedForwardIterator<SignalPDFSection::GroupIteratorTag, const SignalPDFGroup &>;

    SignalPDFSection::SignalPDFSection(Implementation<SignalPDFSection> * imp) :
        PrivateImplementationPattern<SignalPDFSection>(imp)
    {
    }

    SignalPDFSection::~SignalPDFSection() = default;

    SignalPDFSection::GroupIterator
    SignalPDFSection::begin() const
    {
        return _imp->groups.begin();
    }

    SignalPDFSection::GroupIterator
    SignalPDFSection::end() const
    {
        return _imp->groups.end();
    }

    const std::string &
    SignalPDFSection::name() const
    {
        return _imp->name;
    }

    const std::string &
    SignalPDFSection::description() const
    {
        return _imp->description;
    }

    /* SignalPDFs */

    template <> struct WrappedForwardIteratorTraits<SignalPDFs::SignalPDFIteratorTag>
    {
            using UnderlyingIterator = std::map<QualifiedName, SignalPDFEntryPtr>::const_iterator;
    };
    template class WrappedForwardIterator<SignalPDFs::SignalPDFIteratorTag, const std::pair<const QualifiedName, SignalPDFEntryPtr>>;

    template <> struct WrappedForwardIteratorTraits<SignalPDFs::SectionIteratorTag>
    {
            using UnderlyingIterator = std::vector<SignalPDFSection>::const_iterator;
    };
    template class WrappedForwardIterator<SignalPDFs::SectionIteratorTag, const SignalPDFSection &>;

    template <> struct Implementation<SignalPDFs>
    {
            std::vector<SignalPDFSection> signal_pdf_sections;

            Implementation() :
                signal_pdf_sections(SignalPDFSections::instance()->sections())
            {
            }
    };

    SignalPDFs::SignalPDFs() :
        PrivateImplementationPattern<SignalPDFs>(new Implementation<SignalPDFs>())
    {
    }

    SignalPDFs::~SignalPDFs() {}

    SignalPDFEntryPtr
    SignalPDFs::operator[] (const QualifiedName & qn) const
    {
        const auto & signal_pdf_entries = SignalPDFEntries::instance()->entries();

        auto i = signal_pdf_entries.find(qn);
        if (i != signal_pdf_entries.end())
        {
            return i->second;
        }

        return SignalPDFEntryPtr(nullptr);
    }

    SignalPDFs::SignalPDFIterator
    SignalPDFs::begin() const
    {
        const auto & signal_pdf_entries = SignalPDFEntries::instance()->entries();

        return SignalPDFIterator(signal_pdf_entries.begin());
    }

    SignalPDFs::SignalPDFIterator
    SignalPDFs::end() const
    {
        const auto & signal_pdf_entries = SignalPDFEntries::instance()->entries();

        return SignalPDFIterator(signal_pdf_entries.end());
    }

    SignalPDFs::SectionIterator
    SignalPDFs::begin_sections() const
    {
        return SectionIterator(_imp->signal_pdf_sections.begin());
    }

    SignalPDFs::SectionIterator
    SignalPDFs::end_sections() const
    {
        return SectionIterator(_imp->signal_pdf_sections.end());
    }

    void
    SignalPDFs::insert(const QualifiedName & name, const std::string & description, const Options & options, const QualifiedName & numerator,
                       const std::vector<std::string> & numerator_kinematic_names, const QualifiedName & normalization,
                       const std::vector<std::string> & normalization_kinematic_names) const
    {
        const auto & observable_entries = ObservableEntries::instance()->entries();

        // the numerator and normalization must reference known observables; fail fast otherwise
        if (observable_entries.end() == observable_entries.find(numerator))
        {
            throw UnknownObservableError("Cannot create SignalPDF '" + name.str() + "': its numerator '" + numerator.str() + "' is not a known observable");
        }

        if (observable_entries.end() == observable_entries.find(normalization))
        {
            throw UnknownObservableError("Cannot create SignalPDF '" + name.str() + "': its normalization '" + normalization.str() + "' is not a known observable");
        }

        // each sampling variable 'v' (a numerator kinematic variable) should have matching bounds 'v_min' and 'v_max'
        // among the normalization kinematic variables; the Python sampling layer relies on this convention, so warn if it is not met
        for (const auto & variable : numerator_kinematic_names)
        {
            for (const auto & bound : { variable + "_min", variable + "_max" })
            {
                if (normalization_kinematic_names.end() == std::find(normalization_kinematic_names.begin(), normalization_kinematic_names.end(), bound))
                {
                    Log::instance()->message("[SignalPDFs.insert]", ll_warning)
                            << "SignalPDF '" << name.str() << "': the normalization is missing the bound '" << bound << "' for the sampling variable '" << variable
                            << "'; sampling from this PDF may not work as expected";
                }
            }
        }

        SignalPDFEntry * entry = new ConcreteSignalPDFEntry(name, description, options, numerator, normalization, numerator_kinematic_names, normalization_kinematic_names);

        SignalPDFEntries::instance()->insert_or_assign(name, std::shared_ptr<const SignalPDFEntry>(entry));
    }
} // namespace eos
