/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <test/test.hh>
#include <src/rare-b-decays/charm-loops.hh>

#include <cmath>
#include <iostream>

using namespace test;
using namespace wf;

class TwoLoopTest :
    public TestCase
{
    public:
        TwoLoopTest() :
            TestCase("two_loop_test")
        {
        }

        virtual void run() const
        {
            std::cout << "# Compare to Figures in [S2004], consider a global sign!" << std::endl;
            std::cout << "# s F_1^7(R,I) F_2^7(R,I) F_1^9(R,I) F_2^9(R,I)" << std::endl;
            for (auto i(0) ; i < 40 ; ++i)
            {
                double s = 0.0 + i * (19.21) / 40.0;
                double mu = 4.6, m_b = 4.6;

                complex<double> f17 = CharmLoops::F17(mu, s, m_b), f27 = CharmLoops::F27(mu, s, m_b),
                    f19 = CharmLoops::F19(mu, s, m_b), f29 = CharmLoops::F29(mu, s, m_b);

                std::cout
                    << s << '\t'
                    << real(f17) << ' ' << imag(f17) << '\t'
                    << real(f19) << ' ' << imag(f19) << '\t'
                    << real(f27) << ' ' << imag(f27) << '\t'
                    << real(f29) << ' ' << imag(f29) << '\t'
                    << std::endl;
            }

            std::cout << std::endl;
            std::cout << "# Comparison with Christoph Bobeth at q^2 =14 GeV^2, mu = mb = 4.45 GeV" << std::endl;
            std::cout << "A = " << CharmLoops::A(4.45, 14.0, 4.45) << std::endl;
        }
} two_loop_test;
