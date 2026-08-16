/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_YAML_SCHEMA_HH
#define EOS_GUARD_EOS_UTILS_YAML_SCHEMA_HH 1

#include <functional>
#include <string>
#include <vector>

namespace YAML
{
    class Node;
} // namespace YAML

namespace eos
{
    namespace yaml_schema
    {
        /*!
         * The expected YAML node type of a schema field.
         */
        enum class Kind
        {
            Scalar,
            Map,
            Sequence
        };

        /*!
         * Describes one field of an on-disk YAML mapping: its key, expected node
         * type, whether it must be present, and a one-line human-readable
         * description. The same table drives both validation (see validate())
         * and, via a caller's own introspection, generated documentation.
         */
        struct Field
        {
                std::string key;
                Kind        kind;
                bool        required;
                std::string description;
        };

        /*!
         * Checks a YAML mapping node against a schema, invoking @p raise with a
         * human-readable message on the first violation encountered:
         *   - a required field is missing, or
         *   - a present field (required or optional) has the wrong node type.
         *
         * Fields absent from @p n that are not required are silently skipped.
         * Keys present in @p n that are not listed in @p fields are ignored;
         * this is not a "reject unknown keys" validator.
         */
        void validate(const YAML::Node & n, const std::vector<Field> & fields, const std::function<void(const std::string &)> & raise);
    } // namespace yaml_schema
} // namespace eos

#endif
