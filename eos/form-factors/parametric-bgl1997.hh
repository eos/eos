/* vim: set sw=4 sts=4 et tw=120 foldmethod=syntax : */

/*
 * Copyright (c) 2020-2025 Danny van Dyk
 * Copyright (c) 2020      Nico Gubernari
 * Copyright (c) 2020      Christoph Bobeth
 * Copyright (c) 2025      Maximilian Hoverath
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BGL1997_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BGL1997_HH 1

#include <eos/form-factors/form-factors-fwd.hh>
#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>
#include <eos/models/model.hh>
#include <eos/utils/reference-name.hh>
#include <eos/maths/power-of.hh>
#include <eos/maths/szego-polynomial.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/options.hh>

#include <array>

namespace eos
{
    template <typename Process_, typename Transition_> class BGL1997FormFactorTraits;

    template <typename Process_, typename Transition_> class BGL1997FormFactors;


    template <typename Process_>
    class BGL1997FormFactorTraits<Process_, PToV> :
        public virtual ParameterUser
    {
        public:
            UsedParameter m_B, m_V;
            std::array<UsedParameter, 4> masses_1m;
            std::array<UsedParameter, 4> masses_1p;
            std::array<UsedParameter, 3> masses_0m;
            std::array<UsedParameter, 2> masses_0p;
            UsedParameter chi_1m, chi_0p;
            UsedParameter chi_1p, chi_0m;
            UsedParameter chi_T_1m, chi_T_1p;
            UsedParameter t_0;

            IntegerOption n_bound_states_1m;
            IntegerOption n_bound_states_1p;
            IntegerOption n_bound_states_0m;
            IntegerOption n_bound_states_0p;

            BGL1997FormFactorTraits(const Parameters & p, const Options & o, const std::vector<OptionSpecification> & options) :
                m_B(UsedParameter(p[std::string(Process_::name_B) + "@BSZ2015"], *this)),
                m_V(UsedParameter(p[std::string(Process_::name_V) + "@BSZ2015"], *this)),
                masses_1m{{ UsedParameter(p["mass::B_c^*@BSZ2015"], *this),
                            UsedParameter(p["mass::B_c^*[1]@BSZ2015"], *this),
                            UsedParameter(p["mass::B_c^*[2]@BSZ2015"], *this),
                            UsedParameter(p["mass::B_c^*[3]@BSZ2015"], *this)
                }},
                masses_1p{{ UsedParameter(p["mass::B_c,1@BSZ2015"], *this),
                            UsedParameter(p["mass::B_c,1[1]@BSZ2015"], *this),
                            UsedParameter(p["mass::B_c,1[2]@BSZ2015"], *this),
                            UsedParameter(p["mass::B_c,1[3]@BSZ2015"], *this)
                }},
                masses_0m{{ UsedParameter(p["mass::B_c@BSZ2015"], *this),
                            UsedParameter(p["mass::B_c[1]@BSZ2015"], *this),
                            UsedParameter(p["mass::B_c[2]@BSZ2015"], *this)
                }},
                masses_0p{{ UsedParameter(p["mass::B_c,0@BSZ2015"], *this),
                            UsedParameter(p["mass::B_c,0[1]@BSZ2015"], *this)
                }},
                chi_1m(UsedParameter(p["b->c::chiOPE[1^-_V]"], *this)),
                chi_0p(UsedParameter(p["b->c::chiOPE[0^+_V]"], *this)),
                chi_1p(UsedParameter(p["b->c::chiOPE[1^+_A]"], *this)),
                chi_0m(UsedParameter(p["b->c::chiOPE[0^-_A]"], *this)),
                chi_T_1m(UsedParameter(p["b->c::chiOPE[1^-_T]"], *this)),
                chi_T_1p(UsedParameter(p["b->c::chiOPE[1^+_T5]"], *this)),
                t_0(UsedParameter(p["B->D^*::t_0@BGL1997"], *this)),
                n_bound_states_1m(o, options, "n-bound-states-1m"_ok),
                n_bound_states_1p(o, options, "n-bound-states-1p"_ok),
                n_bound_states_0m(o, options, "n-bound-states-0m"_ok),
                n_bound_states_0p(o, options, "n-bound-states-0p"_ok)
            {
            }

            double tp() const
            {
                return power_of<2>(m_B + m_V);
            }

            double tm() const
            {
                return power_of<2>(m_B - m_V);
            }

            complex<double> _z(const complex<double> & s, const complex<double> & s_0, const complex<double> & s_p) const
            {
                return (std::sqrt(s_p - s) - std::sqrt(s_p - s_0)) / (std::sqrt(s_p - s) + std::sqrt(s_p - s_0));
            }

            double _z(const double & s, const double & s_0, const double & s_p) const
            {
                if (s > s_p)
                    throw InternalError("The real conformal mapping is used above threshold: " + stringify(s) + " > " + stringify(s_p));

                return real(_z(complex<double>(s, 0.0), complex<double>(s_0, 0.0), complex<double>(s_p, 0.0)));
            }

            double blaschke_1m(const double & s) const
            {
                // bound states for 1^-
                double blaschke = 1.0;
                for (int i = 0; i < n_bound_states_1m.value(); ++i)
                {
                    if (masses_1m[i]() * masses_1m[i]() <= tp())
                    {
                        blaschke *= _z(s, masses_1m[i]() * masses_1m[i](), tp());
                    }
                }
                return blaschke;
            }

            double blaschke_1p(const double & s) const
            {
                // bound states for 1^+
                double blaschke = 1.0;
                for (int i = 0; i < n_bound_states_1p.value(); ++i)
                {
                    if (masses_1p[i]() * masses_1p[i]() <= tp())
                    {
                        blaschke *= _z(s, masses_1p[i]() * masses_1p[i](), tp());
                    }
                }
                return blaschke;
            }

            double blaschke_0m(const double & s) const
            {
                // bound states for 0^-
                double blaschke = 1.0;
                for (int i = 0; i < n_bound_states_0m.value(); ++i)
                {
                    if (masses_0m[i]() * masses_0m[i]() <= tp())
                    {
                        blaschke *= _z(s, masses_0m[i]() * masses_0m[i](), tp());
                    }
                }
                return blaschke;
            }

            double blaschke_0p(const double & s) const
            {
                // bound states for 0^+
                double blaschke = 1.0;
                for (int i = 0; i < n_bound_states_0p.value(); ++i)
                {
                    if (masses_0p[i]() * masses_0p[i]() <= tp())
                    {
                        blaschke *= _z(s, masses_0p[i]() * masses_0p[i](), tp());
                    }
                }
                return blaschke;
            }
    };

    template <typename Process_>
    class BGL1997FormFactorTraits<Process_, PToP> :
        public virtual ParameterUser
    {
        public:
            UsedParameter m_B, m_P;
            std::array<UsedParameter, 4> masses_1m;
            std::array<UsedParameter, 2> masses_0p;
            UsedParameter chi_1m, chi_0p;
            UsedParameter chi_T_1m;
            UsedParameter t_0;

            IntegerOption n_bound_states_1m;
            IntegerOption n_bound_states_0p;

            BGL1997FormFactorTraits(const Parameters & p, const Options & o, const std::vector<OptionSpecification> & options) :
                m_B(UsedParameter(p[std::string(Process_::name_B) + "@BSZ2015"], *this)),
                m_P(UsedParameter(p[std::string(Process_::name_P) + "@BSZ2015"], *this)),
                masses_1m{{ UsedParameter(p["mass::B_c^*@BSZ2015"], *this),
                            UsedParameter(p["mass::B_c^*[1]@BSZ2015"], *this),
                            UsedParameter(p["mass::B_c^*[2]@BSZ2015"], *this),
                            UsedParameter(p["mass::B_c^*[3]@BSZ2015"], *this)
                }},
                masses_0p{{ UsedParameter(p["mass::B_c,0@BSZ2015"], *this),
                            UsedParameter(p["mass::B_c,0[1]@BSZ2015"], *this)
                }},
                chi_1m(UsedParameter(p["b->c::chiOPE[1^-_V]"], *this)),
                chi_0p(UsedParameter(p["b->c::chiOPE[0^+_V]"], *this)),
                chi_T_1m(UsedParameter(p["b->c::chiOPE[1^-_T]"], *this)),
                t_0(UsedParameter(p["B->D::t_0@BGL1997"], *this)),
                n_bound_states_1m(o, options, "n-bound-states-1m"_ok),
                n_bound_states_0p(o, options, "n-bound-states-0p"_ok)
            {
            }

            double tp() const
            {
                return power_of<2>(m_B + m_P);
            }

            double tm() const
            {
                return power_of<2>(m_B - m_P);
            }

            complex<double> _z(const complex<double> & s, const complex<double> & s_0, const complex<double> & s_p) const
            {
                return (std::sqrt(s_p - s) - std::sqrt(s_p - s_0)) / (std::sqrt(s_p - s) + std::sqrt(s_p - s_0));
            }

            double _z(const double & s, const double & s_0, const double & s_p) const
            {
                if (s > s_p)
                    throw InternalError("The real conformal mapping is used above threshold: " + stringify(s) + " > " + stringify(s_p));

                return real(_z(complex<double>(s, 0.0), complex<double>(s_0, 0.0), complex<double>(s_p, 0.0)));
            }

            double blaschke_1m(const double & s) const
            {
                // bound states for 1^-
                double blaschke = 1.0;
                for (int i = 0; i < n_bound_states_1m.value(); ++i)
                {
                    if (masses_1m[i]() * masses_1m[i]() <= tp())
                    {
                        blaschke *= _z(s, masses_1m[i]() * masses_1m[i](), tp());
                    }
                }
                return blaschke;
            }

            double blaschke_0p(const double & s) const
            {
                // bound states for 0^+
                double blaschke = 1.0;
                for (int i = 0; i < n_bound_states_0p.value(); ++i)
                {
                    if (masses_0p[i]() * masses_0p[i]() <= tp())
                    {
                        blaschke *= _z(s, masses_0p[i]() * masses_0p[i](), tp());
                    }
                }
                return blaschke;
            }
    };

    template <typename Process_> class BGL1997FormFactors<Process_, PToV> :
        public FormFactors<PToV>
    {
        private:
            std::array<UsedParameter, 4> _a_g, _a_f;
            std::array<UsedParameter, 3> _a_F1, _a_F2;
            std::array<UsedParameter, 4> _a_T1;
            std::array<UsedParameter, 3> _a_T2, _a_T23;

            const BGL1997FormFactorTraits<Process_, PToV> _traits;

            const UsedParameter & _mB;
            const UsedParameter & _mV;
            const UsedParameter & t_0;

            static std::string _par_name(const std::string & ff_name);

        public:
            BGL1997FormFactors(const Parameters &, const Options &);
            ~BGL1997FormFactors();

            static FormFactors<PToV> * make(const Parameters & parameters, const Options & options);

            double _phi(const double & s, const double & s_0, const double & K, const unsigned & a, const unsigned & b, const unsigned & c, const                      double & chi) const;

            double g(const double & s) const;
            double f(const double & s) const;
            double F1(const double & s) const;
            double F2(const double & s) const;

            double a_F1_0() const;
            double a_F2_0() const;
            double a_T2_0() const;
            double a_T23_0() const;

            virtual double v(const double & s) const;
            virtual double a_0(const double & s) const;
            virtual double a_1(const double & s) const;
            virtual double a_2(const double & s) const;
            virtual double a_12(const double & s) const;
            virtual double t_1(const double & s) const;
            virtual double t_2(const double & s) const;
            virtual double t_3(const double & s) const;
            virtual double t_23(const double & s) const;

            virtual double f_perp(const double & s) const;
            virtual double f_para(const double & s) const;
            virtual double f_long(const double & s) const;

            virtual double f_perp_T(const double & s) const;
            virtual double f_para_T(const double & s) const;
            virtual double f_long_T(const double & s) const;

            /*!
             * References used in the computation of our (pseudo)observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our (pseudo)observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
            static const std::vector<OptionSpecification> _options;
    };
    extern template class BGL1997FormFactors<BToDstar, PToV>;

    template <typename Process_> class BGL1997FormFactors<Process_, PToP> :
        public FormFactors<PToP>
    {
        private:
            std::array<UsedParameter, 4> _a_f_p, _a_f_0, _a_f_t;

            const BGL1997FormFactorTraits<Process_, PToP> _traits;

            const UsedParameter & _mB;
            const UsedParameter & _mP;
            const UsedParameter & t_0;

            static std::string _par_name(const std::string & ff_name);

        public:
            BGL1997FormFactors(const Parameters &, const Options &);
            ~BGL1997FormFactors();

            static FormFactors<PToP> * make(const Parameters & parameters, const Options & options);

            double _phi(const double & s, const double & s_0, const double & K, const unsigned & a, const unsigned & b, const unsigned & c, const                       double & chi) const;

            virtual double f_p(const double & s) const;
            virtual double f_0(const double & s) const;
            virtual double f_t(const double & s) const;

            virtual double f_plus_T(const double & s) const;

            /*!
             * References used in the computation of our (pseudo)observables.
             */
            static const std::set<ReferenceName> references;

            /*!
             * Options used in the computation of our (pseudo)observables.
             */
            static std::vector<OptionSpecification>::const_iterator begin_options();
            static std::vector<OptionSpecification>::const_iterator end_options();
            static const std::vector<OptionSpecification> _options;
    };
    extern template class BGL1997FormFactors<BToD, PToP>;
}

#endif
