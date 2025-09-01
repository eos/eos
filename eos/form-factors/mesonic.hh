/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2022 Stephan Kuerten
 * Copyright (c) 2010-2024 Danny van Dyk
 * Copyright (c) 2015 Christoph Bobeth
 * Copyright (c) 2022 Philip Lüghausen
 * Copyright (c) 2010 Christian Wacker
 * Copyright (c) 2024 Matthew J. Kirk
 * Copyright (c) 2025 Florian Herren
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

#ifndef EOS_GUARD_SRC_FORM_FACTORS_MESONIC_HH
#define EOS_GUARD_SRC_FORM_FACTORS_MESONIC_HH 1

#include <eos/form-factors/form-factors-fwd.hh>
#include <eos/maths/complex.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/options.hh>
#include <eos/utils/qualified-name.hh>
#include <eos/utils/transitions.hh>

#include <array>
#include <map>
#include <memory>
#include <string>

namespace eos
{
    template <>
    class FormFactors<PToV> :
        public virtual ParameterUser
    {
        public:
            virtual ~FormFactors();

            virtual double v(const double & q2) const = 0;

            virtual double a_0(const double & q2) const = 0;
            virtual double a_1(const double & q2) const = 0;
            virtual double a_2(const double & q2) const = 0;
            virtual double a_12(const double & q2) const = 0;

            virtual double t_1(const double & q2) const = 0;
            virtual double t_2(const double & q2) const = 0;
            virtual double t_3(const double & q2) const = 0;
            virtual double t_23(const double & q2) const = 0;

            virtual double f_perp(const double & q2) const = 0;
            virtual double f_para(const double & q2) const = 0;
            virtual double f_long(const double & q2) const = 0;

            virtual double f_perp_T(const double & q2) const = 0;
            virtual double f_para_T(const double & q2) const = 0;
            virtual double f_long_T(const double & q2) const = 0;

            // for access in the complex q2 plane
            virtual complex<double> v(const complex<double> & q2) const;

            virtual complex<double> a_0(const complex<double> & q2) const;
            virtual complex<double> a_1(const complex<double> & q2) const;
            virtual complex<double> a_12(const complex<double> & q2) const;
            virtual complex<double> a_2(const complex<double> & q2) const;

            virtual complex<double> t_1(const complex<double> & q2) const;
            virtual complex<double> t_2(const complex<double> & q2) const;
            virtual complex<double> t_23(const complex<double> & q2) const;
    };

    template <>
    class FormFactorFactory<PToV>
    {
        public:
            using KeyType = QualifiedName;
            using ValueType = std::function<FormFactors<PToV> * (const Parameters &, const Options &)>;

            static const std::map<KeyType, ValueType> form_factors;

            static std::shared_ptr<FormFactors<PToV>> create(const QualifiedName & name, const Parameters & parameters, const Options & options = Options{ });
            static OptionSpecification option_specification(const qnp::Prefix & process);
            static OptionSpecification option_specification();
    };

    template <>
    class FormFactors<PToGamma> :
        public virtual ParameterUser
    {
        public:
            virtual ~FormFactors();

            // axial form factor
            virtual double F_A(const double & Egamma) const = 0;

            // vector form factor
            virtual double F_V(const double & Egamma) const = 0;
    };

    template <>
    class FormFactorFactory<PToGamma>
    {
        public:
            using KeyType = QualifiedName;
            using ValueType = std::function<FormFactors<PToGamma> * (const Parameters &, const Options &)>;

            static const std::map<KeyType, ValueType> form_factors;

            static std::shared_ptr<FormFactors<PToGamma>> create(const QualifiedName & name, const Parameters & parameters, const Options & options = Options{ });
            static OptionSpecification option_specification(const qnp::Prefix & process);
            static OptionSpecification option_specification();
    };

    template <>
    class FormFactors<PToGammaOffShell> :
        public virtual ParameterUser
    {
        public:
            virtual ~FormFactors();

            // axial current, superposition of transverse and longitudinal polarizations of both currents
            virtual complex<double> F_1(const double & q2, const double & k2) const = 0;

            // axial current, superposition of transverse and longitudinal polarizations of both currents
            virtual complex<double> F_2(const double & q2, const double & k2) const = 0;

            // axial current, pseudoscalar form factor, time-like polarization of the weak current and longitudinal polarization of the em current
            virtual complex<double> F_3(const double & q2, const double & k2) const = 0;

            // vector current, transverse polarization of the weak current and the off-shell photon
            virtual complex<double> F_4(const double & q2, const double & k2) const = 0;

            inline double arg_F_1(const double & q2, const double & k2) const { return std::arg(F_1(q2, k2)); }
            inline double arg_F_2(const double & q2, const double & k2) const { return std::arg(F_2(q2, k2)); }
            inline double arg_F_3(const double & q2, const double & k2) const { return std::arg(F_3(q2, k2)); }
            inline double arg_F_4(const double & q2, const double & k2) const { return std::arg(F_4(q2, k2)); }
            inline double abs_F_1(const double & q2, const double & k2) const { return std::abs(F_1(q2, k2)); }
            inline double abs_F_2(const double & q2, const double & k2) const { return std::abs(F_2(q2, k2)); }
            inline double abs_F_3(const double & q2, const double & k2) const { return std::abs(F_3(q2, k2)); }
            inline double abs_F_4(const double & q2, const double & k2) const { return std::abs(F_4(q2, k2)); }
    };

    template <>
    class FormFactorFactory<PToGammaOffShell>
    {
        public:
            using KeyType = QualifiedName;
            using ValueType = std::function<FormFactors<PToGammaOffShell> * (const Parameters &, const Options &)>;

            static const std::map<KeyType, ValueType> form_factors;

            static std::shared_ptr<FormFactors<PToGammaOffShell>> create(const QualifiedName & name, const Parameters & parameters, const Options & options = Options{ });
            static OptionSpecification option_specification(const qnp::Prefix & process);
            static OptionSpecification option_specification();
    };

    template <>
    class FormFactors<PToP> :
        public virtual ParameterUser
    {
        public:
            virtual ~FormFactors();

            virtual double f_p(const double & s) const = 0;
            virtual double f_0(const double & s) const = 0;
            virtual double f_t(const double & s) const = 0;
            virtual double f_m(const double & s) const;

            // Conventions of GvDV:2020A eq. (A.5)
            virtual double f_plus_T(const double & s) const = 0;

            virtual double f_p_d1(const double & s) const;
            virtual double f_p_d2(const double & s) const;

            // for access in the complex q2 plane
            virtual complex<double> f_p(const complex<double> & q2) const;
            virtual complex<double> f_0(const complex<double> & q2) const;
            virtual complex<double> f_t(const complex<double> & q2) const;

    };

    template <>
    class FormFactorFactory<PToP>
    {
        public:
            using KeyType = QualifiedName;
            using ValueType = std::function<FormFactors<PToP> * (const Parameters &, const Options &)>;

            static const std::map<KeyType, ValueType> form_factors;

            static std::shared_ptr<FormFactors<PToP>> create(const QualifiedName & label, const Parameters & parameters, const Options & options = Options{ });
            static OptionSpecification option_specification(const qnp::Prefix & process);
            static OptionSpecification option_specification();
    };

    template <>
    class FormFactors<PToPP> :
        public virtual ParameterUser
    {
        public:
            virtual ~FormFactors();

            // Partial waves
            virtual std::array<complex<double>, 4> f_perp(const double & q2, const double & k2) const = 0;
            virtual std::array<complex<double>, 4> f_para(const double & q2, const double & k2) const = 0;
            virtual std::array<complex<double>, 4> f_long(const double & q2, const double & k2) const = 0;
            virtual std::array<complex<double>, 4> f_time(const double & q2, const double & k2) const = 0;

            // form factors
            virtual complex<double> f_perp(const double & q2, const double & k2, const double & z) const = 0;
            virtual complex<double> f_para(const double & q2, const double & k2, const double & z) const = 0;
            virtual complex<double> f_long(const double & q2, const double & k2, const double & z) const = 0;
            virtual complex<double> f_time(const double & q2, const double & k2, const double & z) const = 0;

            double im_f_perp(const double & q2, const double & k2, const double & z) const { return std::imag(f_perp(q2, k2, z)); }
            double im_f_para(const double & q2, const double & k2, const double & z) const { return std::imag(f_para(q2, k2, z)); }
            double im_f_long(const double & q2, const double & k2, const double & z) const { return std::imag(f_long(q2, k2, z)); }
            double im_f_time(const double & q2, const double & k2, const double & z) const { return std::imag(f_time(q2, k2, z)); }
    };

    template <>
    class FormFactorFactory<PToPP>
    {
        public:
            using KeyType = QualifiedName;
            using ValueType = std::function<FormFactors<PToPP> * (const Parameters &, const Options &)>;

            static const std::map<KeyType, ValueType> form_factors;

            static std::shared_ptr<FormFactors<PToPP>> create(const QualifiedName & name, const Parameters & parameters, const Options & options = Options{ });
            static OptionSpecification option_specification(const qnp::Prefix & process);
            static OptionSpecification option_specification();
    };

    template <>
    class FormFactors<VToP> :
        public virtual ParameterUser
    {
        public:
            virtual ~FormFactors();

            virtual double h_vbar(const double & s) const = 0;

            virtual double h_abar_1(const double & s) const = 0;
            virtual double h_abar_2(const double & s) const = 0;
            virtual double h_abar_3(const double & s) const = 0;
    };

    template <>
    class FormFactorFactory<VToP>
    {
        public:
            using KeyType = QualifiedName;
            using ValueType = std::function<FormFactors<VToP> * (const Parameters &, const Options &)>;

            static const std::map<KeyType, ValueType> form_factors;

            static std::shared_ptr<FormFactors<VToP>> create(const QualifiedName & label, const Parameters & parameters, const Options & options = Options{ });
            static OptionSpecification option_specification(const qnp::Prefix & process);
    };

    template <>
    class FormFactors<VToV> :
        public virtual ParameterUser
    {
        public:
            virtual ~FormFactors();

            // vector current
            virtual double h_1(const double & s) const = 0;
            virtual double h_2(const double & s) const = 0;
            virtual double h_3(const double & s) const = 0;
            virtual double h_4(const double & s) const = 0;
            virtual double h_5(const double & s) const = 0;
            virtual double h_6(const double & s) const = 0;

            // axial current
            virtual double h_7(const double & s) const = 0;
            virtual double h_8(const double & s) const = 0;
            virtual double h_9(const double & s) const = 0;
            virtual double h_10(const double & s) const = 0;
    };

    template <>
    class FormFactorFactory<VToV>
    {
        public:
            using KeyType = QualifiedName;
            using ValueType = std::function<FormFactors<VToV> * (const Parameters &, const Options &)>;

            static const std::map<KeyType, ValueType> form_factors;

            static std::shared_ptr<FormFactors<VToV>> create(const QualifiedName & label, const Parameters & parameters, const Options & options = Options{ });
            static OptionSpecification option_specification(const qnp::Prefix & process);
    };


    /*
     * Vacuum -> P P transitions
     */
    template <>
    class FormFactors<VacuumToPP> :
        public virtual ParameterUser
    {
        public:
            virtual ~FormFactors();

            // vector form factor
            virtual complex<double> f_p(const double & q2) const = 0;
            virtual double abs2_f_p(const double & q2) const;
            virtual double arg_f_p(const double & q2) const;

            virtual complex<double> f_p(const complex<double> & q2) const = 0;
            virtual double re_f_p(const double & re_q2, const double & im_q2) const;
            virtual double im_f_p(const double & re_q2, const double & im_q2) const;

            // scalar form factor
            virtual complex<double> f_0(const double & q2) const = 0;
            virtual double abs2_f_0(const double & q2) const;
            virtual double arg_f_0(const double & q2) const;

            virtual complex<double> f_0(const complex<double> & q2) const = 0;
            virtual double re_f_0(const double & re_q2, const double & im_q2) const;
            virtual double im_f_0(const double & re_q2, const double & im_q2) const;

            // tensor form factor
            virtual complex<double> f_t(const double & q2) const = 0;
            virtual double abs2_f_t(const double & q2) const;
            virtual double arg_f_t(const double & q2) const;

            virtual complex<double> f_t(const complex<double> & q2) const = 0;
            virtual double re_f_t(const double & re_q2, const double & im_q2) const;
            virtual double im_f_t(const double & re_q2, const double & im_q2) const;
    };

    template <>
    class FormFactorFactory<VacuumToPP>
    {
        public:
            using KeyType = QualifiedName;
            using ValueType = std::function<FormFactors<VacuumToPP> * (const Parameters &, const Options &)>;

            static const std::map<KeyType, ValueType> form_factors;

            static std::shared_ptr<FormFactors<VacuumToPP>> create(const QualifiedName & label, const Parameters & parameters, const Options & options = Options{ });
            static OptionSpecification option_specification(const qnp::Prefix & process);
    };
}


#endif
