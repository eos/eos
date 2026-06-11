/* vim: set sw=4 sts=4 et foldmethod=marker : */

/*
 * Copyright (c) 2025      Danny van Dyk
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

#include <eos/utils/qualified-name-parts.hh>

#include <boost/python.hpp>

#include <string>
#include <variant>
#include <vector>

#ifndef EOS_PYTHON__EOS_CONVERTERS_HH
#  define EOS_PYTHON__EOS_CONVERTERS_HH 1

namespace impl
{
    struct VariantOptionAllowedValuesConverter
    {
            static PyObject *
            convert(const std::variant<eos::qnp::OptionValue, std::vector<eos::qnp::OptionValue>> & value)
            {
                boost::python::list result;

                if (std::holds_alternative<eos::qnp::OptionValue>(value))
                {
                    result.append(std::get<eos::qnp::OptionValue>(value).str());
                }
                else
                {
                    const auto & vector = std::get<std::vector<eos::qnp::OptionValue>>(value);

                    for (const auto & element : vector)
                    {
                        result.append(element.str());
                    }
                }

                return boost::python::incref(boost::python::object(result).ptr());
            }

            static const PyTypeObject *
            get_pytype()
            {
                return &PyList_Type;
            }
    };

} // namespace impl

#endif // EOS_PYTHON__EOS_CONVERTERS_HH
