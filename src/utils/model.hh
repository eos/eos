/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_UTILS_MODEL_HH
#define EOS_GUARD_SRC_UTILS_MODEL_HH 1

#include <src/utils/private_implementation_pattern.hh>
#include <src/utils/parameters.hh>

#include <complex>

namespace eos
{
    using std::complex;

    class Model
    {
        public:
            virtual ~Model() = 0;

            static std::shared_ptr<Model> make(const std::string & name, const Parameters & parameters);

            /* QCD */
            virtual double alpha_s(const double & mu) const = 0;
            virtual double m_b_msbar(const double & mu) const = 0;
            virtual double m_b_pole() const = 0;
            virtual double m_b_ps(const double & mu_f) const = 0;
            virtual double m_c_msbar(const double & mu) const = 0;
            virtual double m_c_pole() const = 0;

            /* CKM matrix elements */
            virtual complex<double> ckm_cb() const = 0;
            virtual complex<double> ckm_us() const = 0;
            virtual complex<double> ckm_ub() const = 0;
            virtual complex<double> ckm_ts() const = 0;
            virtual complex<double> ckm_tb() const = 0;
    };

    struct NoSuchModelError :
        public Exception
    {
        NoSuchModelError(const std::string & name);
    };
}

#endif
