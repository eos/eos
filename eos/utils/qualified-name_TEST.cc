/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016 Danny van Dyk
 * Copyright (c) 2016 Rafael Silva Coutinho
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

#include <eos/utils/qualified-name.hh>

#include <test/test.hh>

#include <cmath>
#include <iostream>

using namespace test;
using namespace eos;

class PrefixTest : public TestCase
{
    public:
        PrefixTest() :
            TestCase("prefix_test")
        {
        }

        virtual void
        run() const
        {
            TEST_CHECK_NO_THROW(auto p = qnp::Prefix("B->K^*ll"));
            TEST_CHECK_NO_THROW(auto p = qnp::Prefix("B->B"));

            TEST_CHECK_THROWS(QualifiedNameSyntaxError, auto p = qnp::Prefix("B->K^*ll::"));
            TEST_CHECK_THROWS(QualifiedNameSyntaxError, auto p = qnp::Prefix("B->>B"));
            TEST_CHECK_THROWS(QualifiedNameSyntaxError, auto p = qnp::Prefix(""));
        }
} prefix_test;

class NameTest : public TestCase
{
    public:
        NameTest() :
            TestCase("name_test")
        {
        }

        virtual void
        run() const
        {
            TEST_CHECK_NO_THROW(auto n = qnp::Name("A_FB(s)"));
            TEST_CHECK_NO_THROW(auto n = qnp::Name("dBR/ds"));
            TEST_CHECK_NO_THROW(auto n = qnp::Name("BR(q2,k2,z)"));
            TEST_CHECK_NO_THROW(auto n = qnp::Name("F_perp[sys=0.30]"));

            TEST_CHECK_THROWS(QualifiedNameSyntaxError, auto n = qnp::Name("::BR(s)"));
            TEST_CHECK_THROWS(QualifiedNameSyntaxError, auto n = qnp::Name(""));
        }
} name_test;

class SuffixTest : public TestCase
{
    public:
        SuffixTest() :
            TestCase("suffix_test")
        {
        }

        virtual void
        run() const
        {
            TEST_CHECK_NO_THROW(auto p = qnp::Suffix("LargeRecoil"));
            TEST_CHECK_NO_THROW(auto p = qnp::Suffix("LHCb:2014A"));

            // empty suffixes are permitted
            TEST_CHECK_NO_THROW(auto p = qnp::Suffix());
            TEST_CHECK_NO_THROW(auto p = qnp::Suffix(""));

            TEST_CHECK_THROWS(QualifiedNameSyntaxError, auto p = qnp::Suffix("@LargeRecoil"));
        }
} suffix_test;

class OptionKeyTest : public TestCase
{
    public:
        OptionKeyTest() :
            TestCase("option_key_test")
        {
        }

        virtual void
        run() const
        {
            TEST_CHECK_NO_THROW(auto k = qnp::OptionKey("form-factors"));
            TEST_CHECK_NO_THROW(auto k = qnp::OptionKey("correlator"));
            TEST_CHECK_NO_THROW(auto k = qnp::OptionKey("q"));
            TEST_CHECK_NO_THROW(auto p = qnp::OptionKey("0-recoil"));
            TEST_CHECK_NO_THROW(auto p = qnp::OptionKey("KEY"));

            TEST_CHECK_THROWS(QualifiedNameSyntaxError, auto p = qnp::OptionKey("key1+key2"));
        }
} option_key_test;

class OptionValueTest : public TestCase
{
    public:
        OptionValueTest() :
            TestCase("option_value_test")
        {
        }

        virtual void
        run() const
        {
            TEST_CHECK_NO_THROW(auto p = qnp::OptionValue("KMPW2010"));
            TEST_CHECK_NO_THROW(auto p = qnp::OptionValue("BCvD2016-model1"));
            TEST_CHECK_NO_THROW(auto p = qnp::OptionValue("d"));
            TEST_CHECK_NO_THROW(auto p = qnp::OptionValue("tw2+tw3"));
            TEST_CHECK_NO_THROW(auto p = qnp::OptionValue("6.346"));

            TEST_CHECK_THROWS(QualifiedNameSyntaxError, auto p = qnp::OptionKey("D^*"));
            TEST_CHECK_THROWS(QualifiedNameSyntaxError, auto p = qnp::OptionKey(""));
        }
} option_value_test;

class QualifiedNameTest : public TestCase
{
    public:
        QualifiedNameTest() :
            TestCase("qualified_name_test")
        {
        }

        virtual void
        run() const
        {
            TEST_CHECK_NO_THROW(auto qn = QualifiedName("B->K^*ll::A_FB(s)@LargeRecoil;form-factors=KMPW2010"));
            TEST_CHECK_NO_THROW(auto qn = QualifiedName("B->K^*ll::A_FB(s)@LargeRecoil;form-factors=BSZ2015"));
            TEST_CHECK_NO_THROW(auto qn = QualifiedName("B^0->K^*0mu^+mu^-::A_FB@LHCb:2015-A;form-factors=BSZ2015"));
            TEST_CHECK_NO_THROW(auto qn = QualifiedName("mass::b(MSbar);opt=har"));
            TEST_CHECK_NO_THROW(auto qn = QualifiedName("D->K::f_++f_0@ETM:2017B;form-factors=BFW2010,rescale-factor=6.346"));
            TEST_CHECK_NO_THROW(auto qn = QualifiedName(qnp::Prefix("mass"), qnp::Name("b(MSbar)"), qnp::Suffix("non-empty")));

            TEST_CHECK_EQUAL_STR("B->K^*ll::A_FB(s)@LargeRecoil;form-factors=KMPW2010", QualifiedName("B->K^*ll::A_FB(s)@LargeRecoil;form-factors=KMPW2010").full());
            TEST_CHECK_EQUAL_STR("B->K^*ll::A_FB(s)@LargeRecoil", QualifiedName("B->K^*ll::A_FB(s)@LargeRecoil;form-factors=KMPW2010").str());

            TEST_CHECK_EQUAL_STR("B->K^*ll::A_FB(s)@LargeRecoil;form-factors=KMPW2010,model=WET",
                                 QualifiedName("B->K^*ll::A_FB(s)@LargeRecoil;form-factors=KMPW2010,model=WET").full());
            TEST_CHECK_EQUAL_STR("B->K^*ll::A_FB(s)@LargeRecoil", QualifiedName("B->K^*ll::A_FB(s)@LargeRecoil;form-factors=KMPW2010,model=WET").str());

            TEST_CHECK_EQUAL_STR("mass::b(MSbar)@non-empty", QualifiedName(qnp::Prefix("mass"), qnp::Name("b(MSbar)"), qnp::Suffix("non-empty")).str());

            TEST_CHECK_EQUAL_STR("mass::b(MSbar)", QualifiedName(qnp::Prefix("mass"), qnp::Name("b(MSbar)")).str());

            TEST_CHECK_THROWS(QualifiedNameSyntaxError, auto qn = QualifiedName(""));
        }
} qualified_name_test;
