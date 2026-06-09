/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013-2016, 2018 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BCL2008_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BCL2008_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>

namespace eos
{
    /* Form Factors according to [BCL:2008A] */
    template <typename Process_, unsigned K_> class BCL2008FormFactors;

    /*
     * P -> P Form Factors in the simplified series expansion according to
     * [BCL:2008A].
     */
    template <typename Process_, unsigned K_, bool with_tensor_> class BCL2008FormFactorBase;

    template <typename Process_> class BCL2008FormFactorBase<Process_, 3u, false> :
        public FormFactors<typename Process_::Transition>
    {
        private:
            /*
             * Fit parametrisation for P -> P according to [BCL:2008A], eq. (11)
             * with K = 3. Note that we factor out the form factor at q^2 = 0
             * by setting t_0 = 0.0, thus b_k -> b_k / b_0. Note that the last
             * coefficient b_K is fixed by eq. (14).
             */
            UsedParameter _f_plus_0, _b_plus_1, _b_plus_2;
            UsedParameter            _b_zero_1, _b_zero_2, _b_zero_3;

        protected:
            double _z(const double & s) const;

        public:
            BCL2008FormFactorBase(const Parameters & p, const Options &);

            virtual double f_p(const double & s) const;

            virtual double f_0(const double & s) const;

            virtual double f_t(const double &) const;

            virtual double f_plus_T(const double &) const;
    };

    template <typename Process_> class BCL2008FormFactorBase<Process_, 4u, false> :
        public FormFactors<typename Process_::Transition>
    {
        private:
            /*
             * Fit parametrisation for P -> P according to [BCL:2008A], eq. (11)
             * with K = 4. Note that we factor out the form factor at q^2 = 0
             * by setting t_0 = 0.0, thus b_k -> b_k / b_0. Note that the last
             * coefficient b_K is fixed by eq. (14).
             */
            UsedParameter _f_plus_0, _b_plus_1, _b_plus_2, _b_plus_3;
            UsedParameter            _b_zero_1, _b_zero_2, _b_zero_3, _b_zero_4;

        protected:
            double _z(const double & s) const;

        public:
            BCL2008FormFactorBase(const Parameters & p, const Options &);

            virtual double f_p(const double & s) const;

            virtual double f_0(const double & s) const;

            virtual double f_t(const double &) const;

            virtual double f_plus_T(const double &) const;
    };

    template <typename Process_> class BCL2008FormFactorBase<Process_, 5u, false> :
        public FormFactors<typename Process_::Transition>
    {
        private:
            /*
             * Fit parametrisation for P -> P according to [BCL:2008A], eq. (11)
             * with K = 5. Note that we factor out the form factor at q^2 = 0
             * by setting t_0 = 0.0, thus b_k -> b_k / b_0. Note that the last
             * coefficient b_K is fixed by eq. (14).
             */
            UsedParameter _f_plus_0, _b_plus_1, _b_plus_2, _b_plus_3, _b_plus_4;
            UsedParameter            _b_zero_1, _b_zero_2, _b_zero_3, _b_zero_4, _b_zero_5;

        protected:
            double _z(const double & s) const;

        public:
            BCL2008FormFactorBase(const Parameters & p, const Options &);

            virtual double f_p(const double & s) const;

            virtual double f_0(const double & s) const;

            virtual double f_t(const double &) const;

            virtual double f_plus_T(const double &) const;
    };

    template <typename Process_> class BCL2008FormFactorBase<Process_, 3u, true> :
        public BCL2008FormFactorBase<Process_, 3u, false>
    {
        private:
            /*
             * Fit parametrisation for P -> P according to [BCL:2008A], eq. (11)
             * with K = 3. Note that we factor out the form factor at q^2 = 0
             * by setting t_0 = 0.0, thus b_k -> b_k / b_0. Note that the last
             * coefficient b_K is fixed by eq. (14).
             *
             * Tensor form factors only. Vector and scalar form factor in BCL2008FormFactorBase<>.
             */
            UsedParameter _f_t_0,    _b_t_1,    _b_t_2;

        public:
            BCL2008FormFactorBase(const Parameters & p, const Options & o);

            virtual double f_t(const double & s) const;
    };

    template <typename Process_> class BCL2008FormFactorBase<Process_, 4u, true> :
        public BCL2008FormFactorBase<Process_, 4u, false>
    {
        private:
            /*
             * Fit parametrisation for P -> P according to [BCL:2008A], eq. (11)
             * with K = 4. Note that we factor out the form factor at q^2 = 0
             * by setting t_0 = 0.0, thus b_k -> b_k / b_0. Note that the last
             * coefficient b_K is fixed by eq. (14).
             *
             * Tensor form factors only. Vector and scalar form factor in BCL2008FormFactorBase<>.
             */
            UsedParameter _f_t_0,    _b_t_1,    _b_t_2,    _b_t_3;

        public:
            BCL2008FormFactorBase(const Parameters & p, const Options & o);

            virtual double f_t(const double & s) const;
    };

    template <typename Process_> class BCL2008FormFactorBase<Process_, 5u, true> :
        public BCL2008FormFactorBase<Process_, 5u, false>
    {
        private:
            /*
             * Fit parametrisation for P -> P according to [BCL:2008A], eq. (11)
             * with K = 5. Note that we factor out the form factor at q^2 = 0
             * by setting t_0 = 0.0, thus b_k -> b_k / b_0. Note that the last
             * coefficient b_K is fixed by eq. (14).
             *
             * Tensor form factors only. Vector and scalar form factor in BCL2008FormFactorBase<>.
             */
            UsedParameter _f_t_0,    _b_t_1,    _b_t_2,    _b_t_3,    _b_t_4;

        public:
            BCL2008FormFactorBase(const Parameters & p, const Options & o);

            virtual double f_t(const double & s) const;
    };


    template <typename Process_, unsigned K_> class BCL2008FormFactors :
        public BCL2008FormFactorBase<Process_, K_, Process_::uses_tensor_form_factors>
    {
        public:
            BCL2008FormFactors(const Parameters & p, const Options & o);

            static FormFactors<PToP> * make(const Parameters & parameters, const Options & options);
    };

    extern template class BCL2008FormFactors<BToPi, 3u>;
    extern template class BCL2008FormFactors<BToPi, 4u>;
    extern template class BCL2008FormFactors<BToPi, 5u>;

    extern template class BCL2008FormFactors<BToK, 3u>;

    extern template class BCL2008FormFactors<BToD, 3u>;
}

#endif
