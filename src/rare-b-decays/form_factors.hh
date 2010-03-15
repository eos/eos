/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_RARE_B_DECAYS_FORM_FACTORS_HH
#define WFITTER_GUARD_SRC_RARE_B_DECAYS_FORM_FACTORS_HH 1

#include <src/utils/parameters.hh>

#include <string>
#include <tr1/memory>

namespace wf
{
    template <typename Transition_>
    class FormFactors;

    template <typename Transition_>
    class FormFactorFactory;

    struct BToKstar { };

    template <>
    class FormFactors<BToKstar>
    {
        public:
            virtual ~FormFactors();

            virtual double v(const double & s_hat) = 0;

            virtual double a_0(const double & s_hat) = 0;
            virtual double a_1(const double & s_hat) = 0;
            virtual double a_2(const double & s_hat) = 0;
    };

    template <>
    class FormFactorFactory<BToKstar>
    {
        public:
            static std::tr1::shared_ptr<FormFactors<BToKstar>>
                create(const std::string & label, const Parameters & parameters);
    };
}


#endif
