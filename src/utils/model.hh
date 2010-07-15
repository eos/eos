/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_UTILS_MODEL_HH
#define WFITTER_GUARD_SRC_UTILS_MODEL_HH 1

#include <src/utils/private_implementation_pattern.hh>
#include <src/utils/parameters.hh>

#include <complex>

namespace wf
{
    using std::complex;

    class Model
    {
        public:
            virtual ~Model() = 0;

            virtual complex<double> ckm_cb() const = 0;
            virtual complex<double> ckm_us() const = 0;
            virtual complex<double> ckm_ub() const = 0;
            virtual complex<double> ckm_ts() const = 0;
            virtual complex<double> ckm_tb() const = 0;
    };

    class StandardModel :
        public Model,
        public PrivateImplementationPattern<StandardModel>
    {
        public:
            StandardModel(const Parameters &);
            ~StandardModel();

            virtual complex<double> ckm_cb() const;
            virtual complex<double> ckm_us() const;
            virtual complex<double> ckm_ub() const;
            virtual complex<double> ckm_ts() const;
            virtual complex<double> ckm_tb() const;
    };
}

#endif
