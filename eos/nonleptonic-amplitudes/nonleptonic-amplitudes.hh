/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Méril Reboud
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

#ifndef EOS_GUARD_EOS_NONLEPTONIC_AMPLITUDES_NONLEPTONIC_AMPLITUDES_HH
#define EOS_GUARD_EOS_NONLEPTONIC_AMPLITUDES_NONLEPTONIC_AMPLITUDES_HH 1

// TODO remove this import statement when the PToPP structure moves to another header
#include <eos/form-factors/mesonic.hh>
#include <eos/nonleptonic-amplitudes/nonleptonic-amplitudes-fwd.hh>
#include <eos/maths/complex.hh>
#include <eos/models/model.hh>
#include <eos/utils/parameters.hh>

#include <array>
#include <map>

namespace eos
{
    class NoSuchNonleptonicAmplitudeError :
        public Exception
    {
        public:
            NoSuchNonleptonicAmplitudeError(const std::string & process, const std::string & tag);
    };


    /* P -> PP transitions */

    // struct PToPP { };

    template <>
    class NonleptonicAmplitudes<PToPP> :
        public virtual ParameterUser
    {
        public:
            virtual ~NonleptonicAmplitudes();

            virtual complex<double> amplitude() const = 0;
            double re_amplitude() const { return real(amplitude()); }
            double im_amplitude() const { return imag(amplitude()); }
            double abs_amplitude() const { return abs(amplitude()); }
            double arg_amplitude() const { return arg(amplitude()); }
    };

    template <>
    class NonleptonicAmplitudeFactory<PToPP>
    {
        public:
            using KeyType = QualifiedName;
            using ValueType = std::function<NonleptonicAmplitudes<PToPP> * (const Parameters &, const Options &)>;

            static const std::map<KeyType, ValueType> amplitudes;

            static std::shared_ptr<NonleptonicAmplitudes<PToPP>> create(const QualifiedName & name, const Parameters & parameters, const Options & options = Options{ });
            static OptionSpecification option_specification(const qnp::Prefix & process);
            static OptionSpecification option_specification();
    };
}
#endif
