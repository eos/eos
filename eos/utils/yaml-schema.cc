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

#include <eos/utils/yaml-schema.hh>

#include <yaml-cpp/yaml.h>

namespace eos
{
    namespace yaml_schema
    {
        namespace
        {
            const char *
            kind_word(Kind kind)
            {
                switch (kind)
                {
                    case Kind::Scalar:   return "scalar value";
                    case Kind::Map:      return "map";
                    case Kind::Sequence: return "sequence";
                }

                return "value";
            }

            bool
            matches(Kind kind, const YAML::Node & n)
            {
                switch (kind)
                {
                    case Kind::Scalar:   return YAML::NodeType::Scalar == n.Type();
                    case Kind::Map:      return YAML::NodeType::Map == n.Type();
                    case Kind::Sequence: return YAML::NodeType::Sequence == n.Type();
                }

                return false;
            }
        } // namespace

        void
        validate(const YAML::Node & n, const std::vector<Field> & fields, const std::function<void(const std::string &)> & raise)
        {
            for (const auto & f : fields)
            {
                if (! n[f.key].IsDefined())
                {
                    if (f.required)
                    {
                        raise("required key '" + f.key + "' not specified");
                        return;
                    }

                    continue;
                }

                if (! matches(f.kind, n[f.key]))
                {
                    raise("key '" + f.key + "' not mapped to a " + std::string(kind_word(f.kind)));
                    return;
                }
            }
        }
    } // namespace yaml_schema
} // namespace eos
