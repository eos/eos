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

#include <boost/python.hpp>

#include <string>
#include <variant>
#include <vector>

namespace impl
{
    struct VariantOptionAllowedValuesConverter
    {
            static PyObject *
            convert(const std::variant<std::string, std::vector<std::string>> & value)
            {
                boost::python::list result;

                if (std::holds_alternative<std::string>(value))
                {
                    result.append(std::get<std::string>(value));
                }
                else
                {
                    const auto & vector = std::get<std::vector<std::string>>(value);

                    for (const auto & element : vector)
                    {
                        result.append(element);
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
