/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Florian Herren
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

#ifndef EOS_GUARD_SRC_SCATTERING_SINGLE_CHANNEL_HH
#define EOS_GUARD_SRC_SCATTERING_SINGLE_CHANNEL_HH 1

#include <eos/scattering/scattering-amplitudes-fwd.hh>
#include <eos/maths/complex.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>
#include <eos/utils/qualified-name.hh>
#include <eos/utils/quantum-numbers.hh>
#include <eos/utils/transitions.hh>

#include <map>
#include <memory>

namespace eos
{
    template <>
    class ScatteringAmplitudes<PPToPP> :
        public virtual ParameterUser
    {
        public:
            virtual ~ScatteringAmplitudes();

            virtual complex<double> scattering_amplitude(const double & s, const unsigned & l, const IsospinRepresentation & i) const = 0;
            virtual complex<double> omnes_factor(const double & s, const unsigned & l, const IsospinRepresentation & i) const = 0;
            virtual complex<double> omnes_outer_function(const double & s, const double & sp, const double & s0, const unsigned & npoints, const unsigned & l, const IsospinRepresentation & i) const = 0;

    };

    template <>
    class ScatteringAmplitudeFactory<PPToPP>
    {
        public:
            using KeyType = QualifiedName;
            using ValueType = std::function<ScatteringAmplitudes<PPToPP> * (const Parameters &, const Options &)>;

            static const std::map<KeyType, ValueType> scattering_amplitudes;

            static std::shared_ptr<ScatteringAmplitudes<PPToPP>> create(const QualifiedName & label, const Parameters & parameters, const Options & options = Options{ });
            static OptionSpecification option_specification(const qnp::Prefix & process);
            static OptionSpecification option_specification();
    };

}

#endif
