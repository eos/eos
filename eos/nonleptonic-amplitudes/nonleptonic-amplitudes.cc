/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Méril Reboud
 * Copyright (c) 2025 Danny van Dyk
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

#include <eos/maths/power-of.hh>
#include <eos/nonleptonic-amplitudes/nonleptonic-amplitudes.hh>
#include <eos/nonleptonic-amplitudes/qcdf-amplitudes.hh>
#include <eos/nonleptonic-amplitudes/su3f-amplitudes.hh>
#include <eos/nonleptonic-amplitudes/topological-amplitudes.hh>
#include <eos/nonleptonic-amplitudes/qcdf_coefs.hh>
#include <eos/utils/options-impl.hh>

#include <map>

namespace eos
{
    NoSuchNonleptonicAmplitudeError::NoSuchNonleptonicAmplitudeError(const std::string & process, const std::string & tag) :
        Exception("No nonleptonic amplitude found for process '" + process + "' and tag '" + tag + "'!")
    {
    }

    namespace su3f
    {
        void
        transpose(rank2 & M)
        {
            for (size_t i = 0; i < 3; ++i)
            {
                for (size_t j = i + 1; j < 3; ++j)
                {
                    std::swap(M[i][j], M[j][i]);
                }
            }
        }

        // Note that this matrix is transposed w.r.t [HTX:2021A] to follow the convention M^i_j = M[i][j]
        // clang-format off
        const std::map<LightMeson, std::function<void (const double &, rank2 &)>>
        psd_octet
        {
            { LightMeson::pi0,      [](const double &, rank2 & res) { res = {{{1.0 / sqrt(2.0), 0.0, 0.0}, {0.0, -1.0 / sqrt(2.0), 0.0}, {0.0, 0.0, 0.0}}}; } },
            { LightMeson::piplus,   [](const double &, rank2 & res) { res = {{{0.0,             1.0, 0.0}, {0.0,  0.0,             0.0}, {0.0, 0.0, 0.0}}}; } },
            { LightMeson::piminus,  [](const double &, rank2 & res) { res = {{{0.0,             0.0, 0.0}, {1.0,  0.0,             0.0}, {0.0, 0.0, 0.0}}}; } },
            { LightMeson::K0,       [](const double &, rank2 & res) { res = {{{0.0,             0.0, 0.0}, {0.0,  0.0,             1.0}, {0.0, 0.0, 0.0}}}; } },
            { LightMeson::K0bar,    [](const double &, rank2 & res) { res = {{{0.0,             0.0, 0.0}, {0.0,  0.0,             0.0}, {0.0, 1.0, 0.0}}}; } },
            { LightMeson::KS,       [](const double &, rank2 & res) { res = {{{0.0,             0.0, 0.0}, {0.0,  0.0,             1.0 / sqrt(2.0)}, {0.0, -1.0 / sqrt(2.0), 0.0}}}; } },
            { LightMeson::Kplus,    [](const double &, rank2 & res) { res = {{{0.0,             0.0, 1.0}, {0.0,  0.0,             0.0}, {0.0, 0.0, 0.0}}}; } },
            { LightMeson::Kminus,   [](const double &, rank2 & res) { res = {{{0.0,             0.0, 0.0}, {0.0,  0.0,             0.0}, {1.0, 0.0, 0.0}}}; } },
            { LightMeson::eta,      [](const double & theta_18, rank2 & res)
                                    {
                                        const double c18 = std::cos(theta_18),
                                                     s18 = std::sin(theta_18);
                                        res = {{
                                                {c18 / sqrt(6.0) - s18 / sqrt(3.0), 0.0,                            0.0},
                                                {0.0,                           c18 / sqrt(6.0) - s18 / sqrt(3.0) ,  0.0},
                                                {0.0,                           0.0,                           -2.0 * c18 / sqrt(6.0) - s18 / sqrt(3.0)}
                                        }};
                                    } },
            { LightMeson::etap,     [](const double & theta_18, rank2 & res)
                                    {
                                        const double c18 = std::cos(theta_18),
                                                     s18 = std::sin(theta_18);
                                        res = {{
                                                {s18 / sqrt(6.0) + c18 / sqrt(3.0), 0.0,                            0.0},
                                                {0.0,                           s18 / sqrt(6.0) + c18 / sqrt(3.0),  0.0},
                                                {0.0,                           0.0,                           -2.0 * s18 / sqrt(6.0) + c18 / sqrt(3.0)}
                                        }};
                                    } },
            { LightMeson::etaq,      [](const double &, rank2 & res) { res = {{{1.0 / sqrt(2.0), 0.0, 0.0}, {0.0, 1.0 / sqrt(2.0), 0.0}, {0.0, 0.0, 0.0}}}; } },
            { LightMeson::etas,      [](const double &, rank2 & res) { res = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}}}; } }
        };
        // clang-format on
    } // namespace su3f

    NonleptonicAmplitudes<PToPP>::~NonleptonicAmplitudes(){};

    const std::map<NonleptonicAmplitudeFactory<PToPP>::KeyType, NonleptonicAmplitudeFactory<PToPP>::ValueType> NonleptonicAmplitudeFactory<PToPP>::amplitudes{
        { "B->PP::topological", &TopologicalRepresentation<PToPP>::make },
        {        "B->PP::SU3F",        &SU3FRepresentation<PToPP>::make },
        {        "B->PP::QCDF",        &QCDFRepresentation<PToPP>::make },
        {   "B->PP::QCDFcoefs",          &QCDFCoefficients<PToPP>::make }
    };

    std::shared_ptr<NonleptonicAmplitudes<PToPP>>
    NonleptonicAmplitudeFactory<PToPP>::create(const QualifiedName & name, const Parameters & parameters, const Options & options)
    {
        Context ctx("When creating a P->PP nonleptonic amplitude");

        std::shared_ptr<NonleptonicAmplitudes<PToPP>> result;

        auto i = NonleptonicAmplitudeFactory<PToPP>::amplitudes.find(name);
        if (NonleptonicAmplitudeFactory<PToPP>::amplitudes.end() != i)
        {
            result.reset(i->second(parameters, name.options() + options));
            return result;
        }

        throw NoSuchNonleptonicAmplitudeError(name.prefix_part().str(), name.name_part().str());
        return result;
    }

    OptionSpecification
    NonleptonicAmplitudeFactory<PToPP>::option_specification(const qnp::Prefix & process)
    {
        OptionSpecification result{ "representation"_ok, {}, "" };
        for (const auto & ff : NonleptonicAmplitudeFactory<PToPP>::amplitudes)
        {
            if (process == std::get<0>(ff).prefix_part())
            {
                result.allowed_values.push_back(std::get<0>(ff).name_part().str());
            }
        }

        return result;
    }

    OptionSpecification
    NonleptonicAmplitudeFactory<PToPP>::option_specification()
    {
        std::set<std::string> allowed_values;
        for (const auto & ff : NonleptonicAmplitudeFactory<PToPP>::amplitudes)
        {
            allowed_values.insert(std::get<0>(ff).name_part().str());
        }

        OptionSpecification result{
            "representation"_ok,
            { allowed_values.cbegin(), allowed_values.cend() },
            ""
        };
        return result;
    }

} // namespace eos
