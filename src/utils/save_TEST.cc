/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <test/test.hh>
#include <src/utils/save.hh>

using namespace test;
using namespace eos;

class SaveTest :
    public TestCase
{
    public:
        SaveTest() :
            TestCase("save_test")
        {
        }

        virtual void run() const
        {
            /* bool */
            {
                bool x = false;

                TEST_CHECK_EQUAL(x, false);
                {
                    Save<bool> s(x, false);
                    TEST_CHECK_EQUAL(x, false);
                }
                TEST_CHECK_EQUAL(x, false);
                {
                    Save<bool> s(x, true);
                    TEST_CHECK_EQUAL(x, true);
                }
                TEST_CHECK_EQUAL(x, false);
            }
        }
} save_test;
