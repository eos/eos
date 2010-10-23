/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <test/test.hh>
#include <src/utils/stringify.hh>

using namespace test;
using namespace eos;

class StringifyTest :
    public TestCase
{
    public:
        StringifyTest() :
            TestCase("b_to_kstar_dilepton_test")
        {
        }

        virtual void run() const
        {
            TEST_CHECK("foo" == stringify("foo"));
            TEST_CHECK("foobar" == stringify(std::string("foobar")));
            TEST_CHECK("-1" == stringify(-1));
        }
} stringify_test;

