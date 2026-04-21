/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015-2026 Danny van Dyk
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

#include <eos/utils/concrete-signal-pdf.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

namespace eos
{
    ConcreteSignalPDF::ConcreteSignalPDF(const QualifiedName & name, const Parameters & parameters, const Kinematics & kinematics, const Options & options,
                                         const QualifiedName & unnormalized_pdf, const QualifiedName & normalization) :
        _name(name),
        _parameters(parameters),
        _kinematics(kinematics),
        _options(options),
        _unnormalized_pdf(Observable::make(unnormalized_pdf, parameters, kinematics, options)),
        _normalization(Observable::make(normalization, parameters, kinematics, options))
    {
        if (_unnormalized_pdf == nullptr)
        {
            throw InternalError("ConcreteSignalPDF: failed to construct unnormalized pdf from " + unnormalized_pdf.str());
        }

        if (_normalization == nullptr)
        {
            throw InternalError("ConcreteSignalPDF: failed to construct normalization from " + normalization.str());
        }
    }

    const QualifiedName &
    ConcreteSignalPDF::name() const
    {
        return _name;
    }

    double
    ConcreteSignalPDF::evaluate() const
    {
        if (auto result = _unnormalized_pdf->evaluate(); result > 0.0) [[likely]]
        {
            return std::log(result);
        }

        return -std::numeric_limits<double>::infinity();
    }

    double
    ConcreteSignalPDF::normalization() const
    {
        if (auto result = _normalization->evaluate(); result > 0.0) [[likely]]
        {
            return std::log(result);
        }

        return -std::numeric_limits<double>::infinity();
    }

    Parameters
    ConcreteSignalPDF::parameters()
    {
        return _parameters;
    }

    Kinematics
    ConcreteSignalPDF::kinematics()
    {
        return _kinematics;
    }

    Options
    ConcreteSignalPDF::options()
    {
        return _options;
    }

    DensityPtr
    ConcreteSignalPDF::clone() const
    {
        return DensityPtr(new ConcreteSignalPDF(_name, _parameters.clone(), _kinematics.clone(), _options, _unnormalized_pdf->name(), _normalization->name()));
    }

    DensityPtr
    ConcreteSignalPDF::clone(const Parameters & parameters) const
    {
        return DensityPtr(new ConcreteSignalPDF(_name, parameters, _kinematics.clone(), _options, _unnormalized_pdf->name(), _normalization->name()));
    }

    Density::Iterator
    ConcreteSignalPDF::begin() const
    {
        return Density::Iterator(_descriptions.cbegin());
    }

    Density::Iterator
    ConcreteSignalPDF::end() const
    {
        return Density::Iterator(_descriptions.cend());
    }

    ConcreteSignalPDFEntry::ConcreteSignalPDFEntry(const QualifiedName & name, const std::string & description, const Options & default_options, const QualifiedName & numerator,
                                                   const QualifiedName & normalization, const std::vector<std::string> & numerator_kinematic_names,
                                                   const std::vector<std::string> & normalization_kinematic_names) :
        _name(name),
        _description(description),
        _default_options(default_options),
        _numerator(numerator),
        _normalization(normalization),
        _numerator_kinematic_names(numerator_kinematic_names),
        _normalization_kinematic_names(normalization_kinematic_names)
    {
    }

    ConcreteSignalPDFEntry::~ConcreteSignalPDFEntry() = default;

    const QualifiedName &
    ConcreteSignalPDFEntry::name() const
    {
        return _name;
    }

    const std::string &
    ConcreteSignalPDFEntry::description() const
    {
        return _description;
    }

    SignalPDFEntry::NumeratorKinematicVariableIterator
    ConcreteSignalPDFEntry::begin_numerator_kinematic_variables() const
    {
        return SignalPDFEntry::NumeratorKinematicVariableIterator(_numerator_kinematic_names.cbegin());
    }

    SignalPDFEntry::NumeratorKinematicVariableIterator
    ConcreteSignalPDFEntry::end_numerator_kinematic_variables() const
    {
        return SignalPDFEntry::NumeratorKinematicVariableIterator(_numerator_kinematic_names.cend());
    }

    SignalPDFEntry::DenominatorKinematicVariableIterator
    ConcreteSignalPDFEntry::begin_denominator_kinematic_variables() const
    {
        return SignalPDFEntry::DenominatorKinematicVariableIterator(_normalization_kinematic_names.cbegin());
    }

    SignalPDFEntry::DenominatorKinematicVariableIterator
    ConcreteSignalPDFEntry::end_denominator_kinematic_variables() const
    {
        return SignalPDFEntry::DenominatorKinematicVariableIterator(_normalization_kinematic_names.cend());
    }

    SignalPDFPtr
    ConcreteSignalPDFEntry::make(const Parameters & parameters, const Kinematics & kinematics, const Options & options) const
    {
        return SignalPDFPtr(new ConcreteSignalPDF(_name, parameters, kinematics, _default_options + options, _numerator, _normalization));
    }

    std::ostream &
    ConcreteSignalPDFEntry::insert(std::ostream & os) const
    {
        return os;
    }

    template <> struct WrappedForwardIteratorTraits<SignalPDFEntry::NumeratorKinematicVariableIteratorTag>
    {
            using UnderlyingIterator = std::vector<std::string>::const_iterator;
    };
    template class WrappedForwardIterator<SignalPDFEntry::NumeratorKinematicVariableIteratorTag, const std::string &>;

    template <> struct WrappedForwardIteratorTraits<SignalPDFEntry::DenominatorKinematicVariableIteratorTag>
    {
            using UnderlyingIterator = std::vector<std::string>::const_iterator;
    };
    template class WrappedForwardIterator<SignalPDFEntry::DenominatorKinematicVariableIteratorTag, const std::string &>;
} // namespace eos
