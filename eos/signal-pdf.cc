/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015-2023 Danny van Dyk
 * Copyright (c) 2019 Ahmet Kokulu
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
#include <eos/utils/density.hh>
#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/log.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <functional>
#include <map>
#include <ostream>

namespace eos
{
    namespace test
    {
        // PDF = (1/2 L_0 + 1/3 L_1 + 1/4 L_2) / 2
        class Legendre1DPDF
        {
            public:
                Legendre1DPDF(const Parameters &, const Options &) {}

                double
                pdf(const double & z) const
                {
                    return (9.0 + 8.0 * z + 9.0 * z * z);
                }

                double
                norm(const double & z_min, const double & z_max) const
                {
                    return (9.0 * (z_max - z_min) + 4.0 * (power_of<2>(z_max) - power_of<2>(z_min)) + 3.0 * (power_of<3>(z_max) - power_of<3>(z_min)));
                }

                static const std::string description;
        };

        const std::string Legendre1DPDF::description = "1D PDF up to 2nd order in z; used for unit tests only.";
    } // namespace test
} // namespace eos

#include <eos/utils/concrete-signal-pdf.hh>

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
            auto name_and_entry_pair = make_signal_pdf("Test::Legendre1D",
                                                       Options{},
                                                       &test::Legendre1DPDF::pdf,
                                                       std::make_tuple(KinematicRange{ "z", -1.0, +1.0, "" }),
                                                       &test::Legendre1DPDF::norm,
                                                       std::make_tuple("z_min", "z_max"));

            _entries->insert(name_and_entry_pair);
        }
    }

    SignalPDFEntries::~SignalPDFEntries() = default;

    SignalPDFNameError::SignalPDFNameError(const std::string & name) :
        Exception("SignalPDF name '" + name + "' is malformed")
    {
    }

    SignalPDFEntry::SignalPDFEntry() {}

    SignalPDFEntry::~SignalPDFEntry() = default;

    template class WrappedForwardIterator<SignalPDFEntry::KinematicRangeIteratorTag, const KinematicRange>;

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

            std::map<QualifiedName, SignalPDFEntryPtr> signal_pdf_entries;

            Implementation() :
                signal_pdf_sections(SignalPDFSections::instance()->sections()),
                signal_pdf_entries(SignalPDFEntries::instance()->entries())
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
        auto i = _imp->signal_pdf_entries.find(qn);
        if (i != _imp->signal_pdf_entries.end())
        {
            return i->second;
        }

        return SignalPDFEntryPtr(nullptr);
    }

    SignalPDFs::SignalPDFIterator
    SignalPDFs::begin() const
    {
        return SignalPDFIterator(_imp->signal_pdf_entries.begin());
    }

    SignalPDFs::SignalPDFIterator
    SignalPDFs::end() const
    {
        return SignalPDFIterator(_imp->signal_pdf_entries.end());
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
} // namespace eos
