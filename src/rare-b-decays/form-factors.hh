/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_FORM_FACTORS_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_FORM_FACTORS_HH 1

#include <src/utils/parameters.hh>

#include <memory>
#include <string>

namespace eos
{
    template <typename Transition_>
    class FormFactors;

    template <typename Transition_>
    class FormFactorFactory;

    /* Tags */

    /*
     * P -> V transitions
     *
     * P: (Heavy) pseudoscalar meson.
     * V: Light vector meson.
     */
    struct PToV { };

    struct PToP { };

    template <>
    class FormFactors<PToV>
    {
        public:
            virtual ~FormFactors();

            virtual double v(const double & s) const = 0;

            virtual double a_0(const double & s) const = 0;
            virtual double a_1(const double & s) const = 0;
            virtual double a_2(const double & s) const = 0;

            // TODO: dipole form factors
    };

    template <>
    class FormFactorFactory<PToV>
    {
        public:
            static std::shared_ptr<FormFactors<PToV>> create(const std::string & label, const Parameters & parameters);
    };

    template <>
    class FormFactors<PToP>
    {
        public:
            virtual ~FormFactors();

            virtual double f_p(const double & s) const = 0;
            virtual double f_0(const double & s) const = 0;
            virtual double f_t(const double & s) const = 0;
    };

    template <>
    class FormFactorFactory<PToP>
    {
        public:
            static std::shared_ptr<FormFactors<PToP>> create(const std::string & label, const Parameters & parameters);
    };
}


#endif
