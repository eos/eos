/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_UTILS_STANDARD_MODEL_HH
#define EOS_GUARD_SRC_UTILS_STANDARD_MODEL_HH 1

#include <src/utils/model.hh>
#include <src/utils/private_implementation_pattern.hh>

namespace eos
{
    class StandardModel :
        public Model,
        public PrivateImplementationPattern<StandardModel>
    {
        public:
            StandardModel(const Parameters &);
            virtual ~StandardModel();

            static std::shared_ptr<Model> make(const Parameters &);

            /* QCD */
            virtual double alpha_s(const double &) const;
            virtual double m_b_msbar(const double & mu) const;
            virtual double m_b_pole() const;
            virtual double m_b_ps(const double & mu_f) const;
            virtual double m_c_msbar(const double & mu) const;
            virtual double m_c_pole() const;

            /* CKM matrix elements */
            virtual complex<double> ckm_cb() const;
            virtual complex<double> ckm_us() const;
            virtual complex<double> ckm_ub() const;
            virtual complex<double> ckm_ts() const;
            virtual complex<double> ckm_tb() const;
    };
}

#endif
