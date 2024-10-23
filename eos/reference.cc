/* vim: set sw=4 sts=4 et foldmethod=marker foldmarker={{{,}}} : */

/*
 * Copyright (c) 2019 Danny van Dyk
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

#include <eos/reference.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#include <iostream>
#include <map>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace fs = boost::filesystem;

namespace eos
{
    template <> struct Implementation<Reference>
    {
            ReferenceName name;
            std::string   authors;
            std::string   title;
            std::string   eprint_archive;
            std::string   eprint_id;
            std::string   inspire_id;
    };

    Reference::Reference(Implementation<Reference> * imp) :
        PrivateImplementationPattern<Reference>(imp)
    {
    }

    Reference::~Reference() = default;

    const ReferenceName &
    Reference::name() const
    {
        return _imp->name;
    }

    const std::string &
    Reference::authors() const
    {
        return _imp->authors;
    }

    const std::string &
    Reference::title() const
    {
        return _imp->title;
    }

    const std::string &
    Reference::eprint_archive() const
    {
        return _imp->eprint_archive;
    }

    const std::string &
    Reference::eprint_id() const
    {
        return _imp->eprint_id;
    }

    const std::string &
    Reference::inspire_id() const
    {
        return _imp->inspire_id;
    }

    template <> struct Implementation<References>
    {
            std::map<ReferenceName, std::shared_ptr<const Reference>> reference_map;

            Implementation() { load(); }

            Implementation(const Implementation &) = delete;

            void
            load()
            {
                fs::path base;
                if (std::getenv("EOS_TESTS_REFERENCES"))
                {
                    std::string envvar = std::string(std::getenv("EOS_TESTS_REFERENCES"));
                    base               = fs::system_complete(envvar);
                }
                else if (std::getenv("EOS_HOME"))
                {
                    std::string envvar = std::string(std::getenv("EOS_HOME"));
                    base               = fs::system_complete(envvar);
                }
                else
                {
                    base = fs::system_complete(EOS_DATADIR "/eos/");
                }

                if (! fs::exists(base))
                {
                    throw InternalError("Could not find the directory containing the references file");
                }

                if (! fs::is_directory(base))
                {
                    throw InternalError("Expect '" + base.string() + " to be a directory");
                }

                auto file_path = base / "references.yaml";

                if (! fs::is_regular_file(status(file_path)))
                {
                    throw InternalError("Expect '" + file_path.string() + " to be a regular file");
                }

                std::string file = file_path.string();
                try
                {
                    YAML::Node root_node = YAML::LoadFile(file);

                    // parse the references
                    for (auto && r : root_node)
                    {
                        ReferenceName name(r.first.Scalar());

                        // authors
                        auto authors_node = r.second["authors"];
                        if (! authors_node)
                        {
                            throw ReferencesInputFileNodeError(file, name.str(), "has no entry named 'authors'");
                        }
                        if (YAML::NodeType::Scalar != authors_node.Type())
                        {
                            throw ReferencesInputFileNodeError(file, name.str() + ".authors", "is not a scalar");
                        }
                        auto authors = authors_node.as<std::string>();

                        // title
                        auto title_node = r.second["title"];
                        if (! title_node)
                        {
                            throw ReferencesInputFileNodeError(file, name.str(), "has no entry named 'title'");
                        }
                        if (YAML::NodeType::Scalar != title_node.Type())
                        {
                            throw ReferencesInputFileNodeError(file, name.str() + ".title", "is not a scalar");
                        }
                        auto title = title_node.as<std::string>();

                        // eprint
                        auto        eprint_node = r.second["eprint"];
                        std::string eprint_archive, eprint_id;
                        if (eprint_node)
                        {
                            if (YAML::NodeType::Map != eprint_node.Type())
                            {
                                throw ReferencesInputFileNodeError(file, name.str() + ".eprint", "is not a map");
                            }

                            auto eprint_archive_node = eprint_node["archive"];
                            if (! eprint_archive_node)
                            {
                                throw ReferencesInputFileNodeError(file, name.str() + ".eprint", "has no entry named 'archive'");
                            }
                            if (YAML::NodeType::Scalar != eprint_archive_node.Type())
                            {
                                throw ReferencesInputFileNodeError(file, name.str() + ".eprint.archive", "is not a scalar");
                            }

                            auto eprint_id_node = eprint_node["id"];
                            if (! eprint_id_node)
                            {
                                throw ReferencesInputFileNodeError(file, name.str() + ".eprint", "has no entry named 'id'");
                            }
                            if (YAML::NodeType::Scalar != eprint_id_node.Type())
                            {
                                throw ReferencesInputFileNodeError(file, name.str() + ".eprint.id", "is not a scalar");
                            }

                            eprint_archive = eprint_archive_node.as<std::string>();
                            eprint_id      = eprint_id_node.as<std::string>();
                        }

                        // inspire id
                        auto        inspire_id_node = r.second["inspire-id"];
                        std::string inspire_id;
                        if (inspire_id_node)
                        {
                            if (YAML::NodeType::Scalar != inspire_id_node.Type())
                            {
                                throw ReferencesInputFileNodeError(file, name.str() + ".inspire-id", "is not a scalar");
                            }

                            inspire_id = inspire_id_node.as<std::string>();
                        }

                        // check for duplicates
                        if (reference_map.end() != reference_map.find(name))
                        {
                            throw ReferencesInputDuplicateError(file, name.str());
                        }

                        reference_map[name] =
                                std::shared_ptr<const Reference>(new Reference(new Implementation<Reference>{ name, authors, title, eprint_archive, eprint_id, inspire_id }));
                    }
                }
                catch (ReferenceNameSyntaxError & e)
                {
                    throw ReferencesInputFileParseError(file, e.what());
                }
                catch (std::exception & e)
                {
                    throw ReferencesInputFileParseError(file, e.what());
                }
            }
    };

    template <> struct WrappedForwardIteratorTraits<References::ReferenceIteratorTag>
    {
            using UnderlyingIterator = std::map<ReferenceName, std::shared_ptr<const Reference>>::const_iterator;
    };
    template class WrappedForwardIterator<References::ReferenceIteratorTag, const std::pair<const ReferenceName, std::shared_ptr<const Reference>>>;

    References::References() :
        PrivateImplementationPattern<References>(new Implementation<References>())
    {
    }

    References::~References() {}

    References::ReferenceIterator
    References::begin() const
    {
        return ReferenceIterator(_imp->reference_map.cbegin());
    }

    References::ReferenceIterator
    References::end() const
    {
        return ReferenceIterator(_imp->reference_map.cend());
    }

    std::shared_ptr<const Reference>
    References::operator[] (const ReferenceName & name) const
    {
        auto i = _imp->reference_map.find(name);

        if (_imp->reference_map.end() == i)
        {
            return {};
        }

        return i->second;
    }

    template <> struct WrappedForwardIteratorTraits<ReferenceUser::ReferenceIteratorTag>
    {
            using UnderlyingIterator = std::set<ReferenceName>::const_iterator;
    };
    template class WrappedForwardIterator<ReferenceUser::ReferenceIteratorTag, const ReferenceName>;

    ReferenceUser::ReferenceIterator
    ReferenceUser::begin_references() const
    {
        return ReferenceIterator(_references.begin());
    }

    ReferenceUser::ReferenceIterator
    ReferenceUser::end_references() const
    {
        return ReferenceIterator(_references.end());
    }

    UnknownReferenceError::UnknownReferenceError(const ReferenceName & name) :
        Exception("Reference '" + name.str() + "' is unknown")
    {
    }

    ReferencesInputFileParseError::ReferencesInputFileParseError(const std::string & file, const std::string & msg) throw() :
        Exception("Malformed references file '" + file + "': " + msg)
    {
    }

    ReferencesInputFileNodeError::ReferencesInputFileNodeError(const std::string & file, const std::string & node, const std::string & msg) throw() :
        Exception("Malformed references file '" + file + "': Node '" + node + "' " + msg)
    {
    }

    ReferencesInputDuplicateError::ReferencesInputDuplicateError(const std::string & file, const std::string & node) throw() :
        Exception("Malformed references file '" + file + "': Duplicate entry for reference '" + node + "'")
    {
    }
} // namespace eos
