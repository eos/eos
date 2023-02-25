/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010-2023 Danny van Dyk
 * Copyright (c) 2018 Ahmet Kokulu
 * Copyright (c) 2018 Christoph Bobeth
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

#ifndef EOS_GUARD_EOS_MODELS_MODEL_HH
#define EOS_GUARD_EOS_MODELS_MODEL_HH 1

#include <eos/models/wilson-coefficients.hh>
#include <eos/utils/options.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/quantum-numbers.hh>

#include <complex>
#include <map>

namespace eos
{
    using std::complex;

    namespace components
    {
        /*!
         * Tags for model components.
         */
        ///@{
        struct CKM;
        struct QCD;
        namespace WET
        {
            struct SBSB;
            struct CBLNu;
            struct UBLNu;
            struct SBNuNu;
            struct DBCU;
            struct SBCU;
        }
        struct DeltaBS1;
        ///@}
    }

    /*!
     * Base classes for individual model components.
     */
    template <typename Tag_> class ModelComponent;

    /*!
     * Base class for the CKM component of models.
     */
    template <> class ModelComponent<components::CKM>
    {
        public:
            /* CKM matrix elements */
            virtual complex<double> ckm_cd() const = 0;
            virtual complex<double> ckm_cs() const = 0;
            virtual complex<double> ckm_cb() const = 0;
            virtual complex<double> ckm_ud() const = 0;
            virtual complex<double> ckm_us() const = 0;
            virtual complex<double> ckm_ub() const = 0;
            virtual complex<double> ckm_td() const = 0;
            virtual complex<double> ckm_ts() const = 0;
            virtual complex<double> ckm_tb() const = 0;
    };

    /*!
     * Base class for the QCD component of models.
     */
    template <> class ModelComponent<components::QCD>
    {
        public:
            /* QCD */
            virtual double alpha_s(const double & mu) const = 0;
            virtual double m_t_msbar(const double & mu) const = 0;
            virtual double m_t_pole() const = 0;
            virtual double m_b_kin(const double & mu_kin) const = 0;
            virtual double m_b_msbar(const double & mu) const = 0;
            virtual double m_b_pole(unsigned int loop_order = 3) const = 0;
            virtual double m_b_ps(const double & mu_f) const = 0;
            virtual double m_c_kin(const double & mu_kin) const = 0;
            virtual double m_c_msbar(const double & mu) const = 0;
            virtual double m_c_pole() const = 0;
            virtual double m_s_msbar(const double & mu) const = 0;
            virtual double m_ud_msbar(const double & mu) const = 0;
            virtual double m_u_msbar(const double & mu) const = 0;
            virtual double m_d_msbar(const double & mu) const = 0;
    };

    /*!
     * Base class for the Delta B = 1 = -Delta S FCNC component of models.
     */
    template <> class ModelComponent<components::DeltaBS1>
    {
        public:
            /* b->s Wilson coefficients */
            virtual WilsonCoefficients<BToS> wilson_coefficients_b_to_s(const double & mu, const std::string & lepton_flavor, const bool & cp_conjugate = false) const = 0;
    };

    /*!
     * Base class for the Delta B = 2 = -Delta S FCNC component of models.
     */
    template <> class ModelComponent<components::WET::SBSB>
    {
        public:
            /* sbar b sbar b Wilson coefficients */
            virtual WilsonCoefficients<wc::SBSB> wet_sbsb() const = 0;
    };

    /*!
     * Base class for the Delta B = 1 = Delta U CC component of models.
     */
    template <> class ModelComponent<components::WET::UBLNu>
    {
        public:
            /* b->u Wilson coefficients */
            virtual WilsonCoefficients<ChargedCurrent> wet_ublnu(LeptonFlavor lepton_flavor, const bool & cp_conjugate = false) const = 0;
    };

    /*!
     * Base class for the Delta B = 1 = Delta C CC component of models.
     */
    template <> class ModelComponent<components::WET::CBLNu>
    {
        public:
            /* b->c Wilson coefficients */
            virtual WilsonCoefficients<ChargedCurrent> wet_cblnu(LeptonFlavor lepton_flavor, const bool & cp_conjugate = false) const = 0;
    };

    /*!
     * Base class for the sbnunu component of models.
     */
    template <> class ModelComponent<components::WET::SBNuNu>
    {
        public:
            /* sbnunu Wilson coefficients */
            virtual WilsonCoefficients<wc::SBNuNu> wet_sbnunu(const bool & cp_conjguate = false) const = 0;
    };

    /*!
     * Base class for the sbcu component of models.
     */
    template <> class ModelComponent<components::WET::SBCU>
    {
        public:
            /* sbcu Wilson coefficients */
            virtual WilsonCoefficients<wc::SBCU> wet_sbcu(const bool & cp_conjguate = false) const = 0;
    };

    /*!
     * Base class for the dbcu component of models.
     */
    template <> class ModelComponent<components::WET::DBCU>
    {
        public:
            /* dbcu Wilson coefficients */
            virtual WilsonCoefficients<wc::DBCU> wet_dbcu(const bool & cp_conjguate = false) const = 0;
    };

    /*!
     * Base class for all models.
     */
    class Model :
        public ParameterUser,
        public virtual ModelComponent<components::CKM>,
        public virtual ModelComponent<components::QCD>,
        public virtual ModelComponent<components::WET::SBSB>,
        public virtual ModelComponent<components::DeltaBS1>,
        public virtual ModelComponent<components::WET::UBLNu>,
        public virtual ModelComponent<components::WET::CBLNu>,
        public virtual ModelComponent<components::WET::SBNuNu>,
        public virtual ModelComponent<components::WET::SBCU>,
        public virtual ModelComponent<components::WET::DBCU>
    {
        public:
            virtual ~Model() = 0;

            using KeyType   = std::string;
            using ValueType = std::function<std::shared_ptr<Model> (const Parameters &, const Options &)>;
            static const std::map<KeyType, ValueType> models;

            static std::shared_ptr<Model> make(const std::string & name, const Parameters & parameters, const Options & options);
            static OptionSpecification option_specification();
    };

    struct NoSuchModelError :
        public Exception
    {
        NoSuchModelError(const std::string & name);
    };
}

#endif
