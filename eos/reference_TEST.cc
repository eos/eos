/* vim: set sw=4 sts=4 et foldmethod=marker foldmarker={{{,}}} : */

/*
 * Copyright (c) 2019-2026 Danny van Dyk
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

#include <test/test.hh>

#include <map>
#include <string>
#include <vector>

using namespace test;
using namespace eos;

class ReferencesTest : public TestCase
{
    public:
        ReferencesTest() :
            TestCase("references_test")
        {
        }

        virtual void
        run() const
        {
            /* Test iterating over references */
            {
                std::cout << "# References :" << std::endl;

                auto references = References();

                for (const auto & r : references)
                {
                    std::cout << "#  " << r.first.str() << ": ";
                    std::cout << std::endl;
                }
                std::cout << std::endl;
            }

            /* Test retrieving Reference by name */
            {
                auto references = References();

                static const std::vector<ReferenceName> names{ "ATLAS:2013A", "BHvD:2012A" };

                for (auto & n : names)
                {
                    std::shared_ptr<const Reference> r;

                    TEST_CHECK_NO_THROW(r = references[n]);

                    TEST_CHECK(r.get() != nullptr);
                }
            }

            /* Test that names, INSPIRE ids, and eprints map one-to-one */
            {
                auto references = References();

                std::map<std::string, ReferenceName> name_of_inspire_id;
                std::map<std::string, ReferenceName> name_of_eprint;
                std::map<std::string, std::string>   eprint_of_inspire_id;
                std::map<std::string, std::string>   inspire_id_of_eprint;

                for (const auto & r : references)
                {
                    const ReferenceName & name       = r.first;
                    const std::string &   inspire_id = r.second->inspire_id();
                    const std::string &   eprint     = r.second->eprint_id();

                    if (! inspire_id.empty())
                    {
                        auto [i, inserted] = name_of_inspire_id.try_emplace(inspire_id, name);
                        TEST_CHECK_MSG(inserted, "references '" + name.str() + "' and '" + i->second.str() + "' share the INSPIRE id '" + inspire_id + "'");
                    }

                    if (! eprint.empty())
                    {
                        auto [i, inserted] = name_of_eprint.try_emplace(eprint, name);
                        TEST_CHECK_MSG(inserted, "references '" + name.str() + "' and '" + i->second.str() + "' share the eprint '" + eprint + "'");
                    }

                    if (inspire_id.empty() || eprint.empty())
                    {
                        continue;
                    }

                    auto [i, i_inserted] = eprint_of_inspire_id.try_emplace(inspire_id, eprint);
                    TEST_CHECK_MSG(i_inserted || (i->second == eprint), "INSPIRE id '" + inspire_id + "' maps to the eprints '" + i->second + "' and '" + eprint + "'");

                    auto [j, j_inserted] = inspire_id_of_eprint.try_emplace(eprint, inspire_id);
                    TEST_CHECK_MSG(j_inserted || (j->second == inspire_id), "eprint '" + eprint + "' maps to the INSPIRE ids '" + j->second + "' and '" + inspire_id + "'");
                }
            }

            /* Test References::has */
            {
                auto references = References();

                TEST_CHECK_EQUAL(references.has("ATLAS:2013A"), true);
                TEST_CHECK_EQUAL(references.has("Nobody:2000A"), false);

                // has() and operator[] must agree
                TEST_CHECK(nullptr != references["ATLAS:2013A"].get());
                TEST_CHECK(nullptr == references["Nobody:2000A"].get());
            }
        }
} references_test;
