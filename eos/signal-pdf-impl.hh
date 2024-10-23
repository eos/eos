/* vim: set sw=4 sts=4 et tw=150 foldmethod=syntax : */

/*
 * Copyright (c) 2019-2023 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_SIGNAL_PDF_IMPL_HH
#define EOS_GUARD_EOS_SIGNAL_PDF_IMPL_HH 1

#include <eos/signal-pdf.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/stringify.hh>

#include <array>
#include <map>

namespace eos
{
    namespace impl
    {
        extern std::map<QualifiedName, std::shared_ptr<const SignalPDFEntry>> signal_pdf_entries;
    }

    template <> struct Implementation<SignalPDFGroup>
    {
            std::string name;

            std::string description;

            std::map<QualifiedName, SignalPDFEntryPtr> entries;

            Implementation(const std::string & name, const std::string & description, std::initializer_list<std::pair<const QualifiedName, SignalPDFEntryPtr>> && entries) :
                name(name),
                description(description),
                entries(entries)
            {
            }
    };

    template <> struct Implementation<SignalPDFSection>
    {
            std::string name;

            std::string description;

            std::vector<SignalPDFGroup> groups;

            Implementation(const std::string & name, const std::string & description, std::initializer_list<SignalPDFGroup> && groups) :
                name(name),
                description(description),
                groups(groups)
            {
            }
    };

    template <typename Decay_, typename... PDFArgs_, typename... PDFKinematicRanges_, typename... NormArgs_, typename... NormKinematicNames_>
    std::pair<QualifiedName, std::shared_ptr<SignalPDFEntry>>
    make_signal_pdf(const char * name, const Options & default_options, double (Decay_::*pdf)(const PDFArgs_ &...) const,
                    const std::tuple<PDFKinematicRanges_...> & pdf_kinematic_ranges, double (Decay_::*norm)(const NormArgs_ &...) const,
                    const std::tuple<NormKinematicNames_...> & norm_kinematic_names)
    {
        static_assert(sizeof...(PDFArgs_) == sizeof...(PDFKinematicRanges_), "Need as many function arguments for the PDF as kinematics ranges!");
        static_assert(sizeof...(NormArgs_) == sizeof...(NormKinematicNames_), "Need as many function arguments for the normalization as kinematics names!");

        QualifiedName qn(name);

        std::function<double(const Decay_ *, const PDFArgs_ &...)>  pdf_function  = std::mem_fn(pdf);
        std::function<double(const Decay_ *, const NormArgs_ &...)> norm_function = std::mem_fn(norm);

        auto entry_ptr =
                std::shared_ptr<SignalPDFEntry>(make_concrete_signal_pdf_entry(qn, default_options, pdf_function, pdf_kinematic_ranges, norm_function, norm_kinematic_names));

        return std::make_pair(qn, entry_ptr);
    }

    template <typename Decay_, typename... PDFArgs_, typename... PDFKinematicRanges_, typename... NormArgs_, typename... NormKinematicNames_>
    std::pair<QualifiedName, std::shared_ptr<SignalPDFEntry>>
    make_signal_pdf(const char * name, const Options & default_options, double (Decay_::*pdf)(const PDFArgs_ &...) const,
                    const std::tuple<PDFKinematicRanges_...> & pdf_kinematic_ranges, const std::function<double(const Decay_ *, const NormArgs_ &...)> & norm_function,
                    const std::tuple<NormKinematicNames_...> & norm_kinematic_names)
    {
        static_assert(sizeof...(PDFArgs_) == sizeof...(PDFKinematicRanges_), "Need as many function arguments for the PDF as kinematics ranges!");
        static_assert(sizeof...(NormArgs_) == sizeof...(NormKinematicNames_), "Need as many function arguments for the normalization as kinematics names!");

        QualifiedName qn(name);

        std::function<double(const Decay_ *, const PDFArgs_ &...)> pdf_function = std::mem_fn(pdf);

        auto entry_ptr =
                std::shared_ptr<SignalPDFEntry>(make_concrete_signal_pdf_entry(qn, default_options, pdf_function, pdf_kinematic_ranges, norm_function, norm_kinematic_names));

        return std::make_pair(qn, entry_ptr);
    }

    template <> struct WrappedForwardIteratorTraits<SignalPDFEntry::KinematicRangeIteratorTag>
    {
            using UnderlyingIterator = const KinematicRange *;
    };
} // namespace eos

#endif
