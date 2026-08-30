/* vim: set sw=4 sts=4 et foldmethod=marker : */

/*
 * Copyright (c) 2025-2026 Danny van Dyk
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
#include <eos/utils/quantum-numbers.hh>

#include <boost/python.hpp>

#include <map>
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

    struct LeptonFlavorFromPythonStringConverter
    {
            LeptonFlavorFromPythonStringConverter() { boost::python::converter::registry::push_back(&convertible, &construct, boost::python::type_id<eos::LeptonFlavor>()); }

            static const std::map<std::string, eos::LeptonFlavor> &
            values()
            {
                static const std::map<std::string, eos::LeptonFlavor> result{
                    {   "e", eos::LeptonFlavor::electron },
                    {  "mu",     eos::LeptonFlavor::muon },
                    { "tau",    eos::LeptonFlavor::tauon }
                };

                return result;
            }

            static void *
            convertible(PyObject * object)
            {
                if (! PyUnicode_Check(object))
                {
                    return nullptr;
                }

                const char * const string = PyUnicode_AsUTF8(object);

                if ((nullptr == string) || (values().end() == values().find(string)))
                {
                    return nullptr;
                }

                return object;
            }

            static void
            construct(PyObject * object, boost::python::converter::rvalue_from_python_stage1_data * data)
            {
                void * const storage = reinterpret_cast<boost::python::converter::rvalue_from_python_storage<eos::LeptonFlavor> *>(data)->storage.bytes;

                new (storage) eos::LeptonFlavor(values().at(PyUnicode_AsUTF8(object)));

                data->convertible = storage;
            }
    };

} // namespace impl

#endif // EOS_PYTHON__EOS_CONVERTERS_HH
