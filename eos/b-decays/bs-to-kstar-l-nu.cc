/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013-2016,2021 Danny van Dyk
 * Copyright (c) 2013 Bastian Müller
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

#include <eos/b-decays/b-to-l-nu.hh>
#include <eos/b-decays/b-to-psd-l-nu.hh>
#include <eos/b-decays/b-to-v-l-nu.hh>
#include <eos/b-decays/bs-to-kstar-l-nu.hh>
#include <eos/form-factors/form-factors.hh>
#include <eos/models/model.hh>
#include <eos/utils/destringify.hh>
#include <eos/maths/integrate-impl.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <cmath>
#include <functional>

namespace eos
{
    using namespace eos::btovlnu;
    using std::norm;

    /*
     * Decay: B_s -> K^* l nubar, cf. [FMvD2015]
     */
    template <>
    struct Implementation<BsToKstarLeptonNeutrino>
    {
        std::shared_ptr<Model> model;

        UsedParameter hbar;

        UsedParameter m_Bs;

        UsedParameter m_Kstar;

        UsedParameter m_b_MSbar;

        LeptonFlavorOption opt_l;

        UsedParameter m_l;

        UsedParameter mu;

        UsedParameter g_fermi;

        UsedParameter tau;

        std::shared_ptr<FormFactors<PToV>> form_factors;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            hbar(p["QM::hbar"], u),
            m_Bs(p["mass::B_s"], u),
            m_Kstar(p["mass::K_u^*"], u),
            m_b_MSbar(p["mass::b(MSbar)"], u),
            opt_l(o, options, "l"),
            m_l(p["mass::" + opt_l.str()], u),
            mu(p["ub" + opt_l.str() + "nu" + opt_l.str() + "::mu"], u),
            g_fermi(p["WET::G_Fermi"], u),
            tau(p["life_time::B_s"], u),
            form_factors(FormFactorFactory<PToV>::create("B_s->K^*::" + o.get("form-factors", "FMvD2015"), p, o))
        {
            Context ctx("When constructing B_s->K*lnu observable");

            u.uses(*form_factors);
            u.uses(*model);
        }

        inline double beta_l(const double & s) const
        {
            return (s - m_l * m_l) / s;
        }

        double norm(const double & s) const
        {
            return std::sqrt(power_of<2>(g_fermi()) / 3.0 / 1024 / power_of<3>(M_PI * m_Bs())
                    * std::norm(model->ckm_ub()) * s * beta_l(s)
                    * std::sqrt(lambda(m_Bs * m_Bs, m_Kstar * m_Kstar, s))); // cf. [FMvD2015], Eq. (10), p. 3
        }

        Amplitudes amplitudes(const double & s) const
        {
            static const double sqrt2 = sqrt(2.0);

            WilsonCoefficients<ChargedCurrent> wc = model->wet_ublnu(opt_l.value(), false);
            double m_Bs2 = m_Bs * m_Bs;
            double sqrts = sqrt(s), lam = lambda(m_Bs * m_Bs, m_Kstar * m_Kstar, s), sqrtlam = sqrt(lam);
            double N = this->norm(s);

            // form factors as in [FMvD2015], p. 11, eqs. (A5) and (A6)
            //double Ftime = sqrtlam / m_Bs2 * form_factors->a_0(s);
            double Flong = 8.0 * m_Kstar / m_Bs * form_factors->a_12(s);
            double Fperp = sqrt2 * sqrtlam / m_Bs / (m_Bs + m_Kstar) * form_factors->v(s);
            double Fpara = sqrt2 * (m_Bs + m_Kstar) / m_Bs * form_factors->a_1(s);

            // compute amplitudes, [FMvD2015], p. 13, Eqs. (B9) and (B10)
            Amplitudes result;
            result.a_long_left = -4.0 * N * m_Bs2 / sqrts * (wc.cvl() - wc.cvr()) * Flong;
            result.a_perp_left = +4.0 * N * m_Bs * (wc.cvl() + wc.cvr()) * Fperp;
            result.a_para_left = -4.0 * N * m_Bs * (wc.cvl() - wc.cvr()) * Fpara;
            result.a_time_left = 0.0;//-4.0 * N * m_Bs2 / m_b * (wc.csl() - wc.csr()) * Ftime;
            result.a_paraperp  = 0.0;//+8.0 * N * m_Bs * wc.ct() * FTlong;
            result.a_longpara  = 0.0;//+4.0 * sqrt(2.0) * N * m_Bs * wc.ct() * FTpara;
            result.a_timeperp  = 0.0;//+4.0 * sqrt(2.0) * N * m_Bs * wc.ct() * FTperp;

            return result;
        }

        std::array<double, 12> differential_angular_coefficients_array(const double & s) const
        {
            return angular_coefficients_array(amplitudes(s));
        }

        AngularCoefficients differential_angular_coefficients(const double & s) const
        {
            return AngularCoefficients(angular_coefficients_array(amplitudes(s)));
        }

        AngularCoefficients integrated_angular_coefficients(const double & s_min, const double & s_max) const
        {
            std::function<std::array<double, 12> (const double &)> integrand =
                    std::bind(&Implementation<BsToKstarLeptonNeutrino>::differential_angular_coefficients_array, this, std::placeholders::_1);
            std::array<double, 12> integrated_angular_coefficients_array = integrate1D(integrand, 64, s_min, s_max);

            return AngularCoefficients(integrated_angular_coefficients_array);
        }

        double Ftime(const double & s)
        {
            double lam = lambda(m_Bs * m_Bs, m_Kstar * m_Kstar, s), sqrtlam = sqrt(lam);
            return -sqrtlam / m_Bs * form_factors->a_0(s);
        }

        double Flong(const double & s)
        {
            return 8.0 * m_Kstar * form_factors->a_12(s);
        }

        double Fperp(const double & s)
        {
            double lam = lambda(m_Bs * m_Bs, m_Kstar * m_Kstar, s), sqrtlam = sqrt(lam);
            return sqrt(2.0) * sqrtlam / (m_Bs + m_Kstar) * form_factors->v(s);
        }

        double Fpara(const double & s)
        {
            return sqrt(2.0) * (m_Bs + m_Kstar) * form_factors->a_1(s);
        }
    };

    const std::vector<OptionSpecification>
    Implementation<BsToKstarLeptonNeutrino>::options
    {
        Model::option_specification(),
        FormFactorFactory<PToV>::option_specification(),
        { "l", { "e", "mu", "tau" }, "mu" }
    };

    BsToKstarLeptonNeutrino::BsToKstarLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BsToKstarLeptonNeutrino>(new Implementation<BsToKstarLeptonNeutrino>(parameters, options, *this))
    {

    }

    BsToKstarLeptonNeutrino::~BsToKstarLeptonNeutrino()
    {
    }

    double
    BsToKstarLeptonNeutrino::Ftime(const double & s) const
    {
        return _imp->Ftime(s);
    }

    double
    BsToKstarLeptonNeutrino::Flong(const double & s) const
    {
        return _imp->Flong(s);
    }

    double
    BsToKstarLeptonNeutrino::Fpara(const double & s) const
    {
        return _imp->Fpara(s);
    }

    double
    BsToKstarLeptonNeutrino::Fperp(const double & s) const
    {
        return _imp->Fperp(s);
    }

    double
    BsToKstarLeptonNeutrino::differential_branching_ratio(const double & s) const
    {
        return differential_decay_width(s) * _imp->tau() / _imp->hbar();
    }

    double
    BsToKstarLeptonNeutrino::differential_decay_width(const double & s) const
    {
        return decay_width(_imp->differential_angular_coefficients(s));
    }

    double
    BsToKstarLeptonNeutrino::differential_forward_backward_asymmetry(const double & s) const
    {
        // cf. [BHvD2010], p. 6, eq. (2.8)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j6s / decay_width(a_c);
    }

    double
    BsToKstarLeptonNeutrino::differential_transverse_asymmetry_2(const double & s) const
    {
        // cf. [BHvD2010], p. 6, eq. (2.10)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 0.5 * a_c.j3 / a_c.j2s;
    }

    double
    BsToKstarLeptonNeutrino::differential_transverse_asymmetry_3(const double & s) const
    {
        // cf. [BHvD2010], p. 6, eq. (2.11)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return sqrt((4.0 * power_of<2>(a_c.j4) + power_of<2>(_imp->beta_l(s) * a_c.j7)) / (-2.0 * a_c.j2c * (2.0 * a_c.j2s + a_c.j3)));
    }

    double
    BsToKstarLeptonNeutrino::differential_transverse_asymmetry_4(const double & s) const
    {
        // cf. [BHvD2010], p. 6, eq. (2.12)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return sqrt((power_of<2>(_imp->beta_l(s) * a_c.j5) + 4.0 * power_of<2>(a_c.j8)) / (4.0 * power_of<2>(a_c.j4) + power_of<2>(_imp->beta_l(s) * a_c.j7)));
    }

    double
    BsToKstarLeptonNeutrino::differential_transverse_asymmetry_5(const double & s) const
    {
        // cf. [BS2011], eq. (34), p. 9 for the massless case
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return std::sqrt(16.0 * power_of<2>(a_c.j2s) - power_of<2>(a_c.j6s) - 4.0 * (power_of<2>(a_c.j3) + power_of<2>(a_c.j9)))
            / 8.0 / a_c.j2s;
    }

    double
    BsToKstarLeptonNeutrino::differential_transverse_asymmetry_re(const double & s) const
    {
        // cf. [BS2011], eq. (38), p. 10
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 0.25 * _imp->beta_l(s) * a_c.j6s / a_c.j2s;
    }

    double
    BsToKstarLeptonNeutrino::differential_transverse_asymmetry_im(const double & s) const
    {
        // cf. [BS2011], eq. (30), p. 8
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 0.5 * a_c.j9 / a_c.j2s;
    }

    double
    BsToKstarLeptonNeutrino::differential_longitudinal_polarisation(const double & s) const
    {
        // cf. [BHvD2012], p. 5, eq. (3.15)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return (a_c.j1c - a_c.j2c / 3.0) / decay_width(a_c);
    }

    double
    BsToKstarLeptonNeutrino::differential_transversal_polarisation(const double & s) const
    {
        // cf. [BHvD2012], p. 5, eq. (3.14)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 2.0 * (a_c.j1s - a_c.j2s / 3.0) / decay_width(a_c);
    }

    double
    BsToKstarLeptonNeutrino::differential_h_1(const double & s) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.13)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return sqrt(2.0) * a_c.j4 / sqrt(-a_c.j2c * (2.0 * a_c.j2s - a_c.j3));
    }

    double
    BsToKstarLeptonNeutrino::differential_h_2(const double & s) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.14)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return _imp->beta_l(s) * a_c.j5 / sqrt(-2.0 * a_c.j2c * (2.0 * a_c.j2s + a_c.j3));
    }

    double
    BsToKstarLeptonNeutrino::differential_h_3(const double & s) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.15)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return _imp->beta_l(s) * a_c.j6s / (2.0 * sqrt(power_of<2>(2.0 * a_c.j2s) - power_of<2>(a_c.j3)));
    }

    double
    BsToKstarLeptonNeutrino::differential_h_4(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return sqrt(2.0) * a_c.j8 / sqrt(-a_c.j2c * (2.0 * a_c.j2s + a_c.j3));
    }

    double
    BsToKstarLeptonNeutrino::differential_h_5(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return -a_c.j9 / sqrt(power_of<2>(2.0 * a_c.j2s) + power_of<2>(a_c.j3));
    }

    double
    BsToKstarLeptonNeutrino::integrated_decay_width(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return decay_width(a_c);
    }

    double
    BsToKstarLeptonNeutrino::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return integrated_decay_width(s_min, s_max) * _imp->tau() / _imp->hbar();
    }

    double
    BsToKstarLeptonNeutrino::integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.8), p. 6
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);

        return a_c.j6s / decay_width(a_c);
    }

    double
    BsToKstarLeptonNeutrino::integrated_longitudinal_polarisation(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2012], p. 5, eq. (3.15)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return (a_c.j1c - a_c.j2c / 3.0) / decay_width(a_c);
    }

    double
    BsToKstarLeptonNeutrino::integrated_transversal_polarisation(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2012], p. 5, eq. (3.14)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 2.0 * (a_c.j1s - a_c.j2s / 3.0) / decay_width(a_c);
    }

    double
    BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_2(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.10), p. 6
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 0.5 * a_c.j3 / a_c.j2s;
    }

    double
    BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_3(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.11), p. 6
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);

        return sqrt((4.0 * power_of<2>(a_c.j4) + power_of<2>(a_c.j7)) / (-2.0 * a_c.j2c * (2.0 * a_c.j2s + a_c.j3)));
    }

    double
    BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_4(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.12), p. 6
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);

        return sqrt((power_of<2>(a_c.j5) + 4.0 * power_of<2>(a_c.j8)) / (4.0 * power_of<2>(a_c.j4) + power_of<2>(a_c.j7)));
    }

    double
    BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_5(const double & s_min, const double & s_max) const
    {
        // cf. [BS2011], eq. (34), p. 9 for the massless case
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return std::sqrt(16.0 * power_of<2>(a_c.j2s) - power_of<2>(a_c.j6s) - 4.0 * (power_of<2>(a_c.j3) + power_of<2>(a_c.j9)))
            / 8.0 / a_c.j2s;
    }

    double
    BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_re(const double & s_min, const double & s_max) const
    {
        // cf. [BS2011], eq. (38), p. 10
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 0.25 * a_c.j6s / a_c.j2s;
    }

    double
    BsToKstarLeptonNeutrino::integrated_transverse_asymmetry_im(const double & s_min, const double & s_max) const
    {
        // cf. [BS2011], eq. (30), p. 8
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 0.5 * a_c.j9 / a_c.j2s;
    }

    double
    BsToKstarLeptonNeutrino::integrated_h_1(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.13)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return sqrt(2.0) * a_c.j4 / sqrt(-a_c.j2c * (2.0 * a_c.j2s - a_c.j3));
    }

    double
    BsToKstarLeptonNeutrino::integrated_h_2(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.14)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return  a_c.j5 / sqrt(-2.0 * a_c.j2c * (2.0 * a_c.j2s + a_c.j3));
    }

    double
    BsToKstarLeptonNeutrino::integrated_h_3(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.15)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j6s / (2.0 * sqrt(power_of<2>(2.0 * a_c.j2s) - power_of<2>(a_c.j3)));
    }

    double
    BsToKstarLeptonNeutrino::integrated_h_4(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return sqrt(2.0) * a_c.j8 / sqrt(-a_c.j2c * (2.0 * a_c.j2s + a_c.j3));
    }

    double
    BsToKstarLeptonNeutrino::integrated_h_5(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return -a_c.j9 / sqrt(power_of<2>(2.0 * a_c.j2s) + power_of<2>(a_c.j3));
    }

    double
    BsToKstarLeptonNeutrino::integrated_s_1s(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j1s / decay_width(a_c);
    }

    double
    BsToKstarLeptonNeutrino::integrated_s_1c(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j1c / decay_width(a_c);
    }

    double
    BsToKstarLeptonNeutrino::integrated_s_2s(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j2s / decay_width(a_c);
    }

    double
    BsToKstarLeptonNeutrino::integrated_s_2c(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j2c / decay_width(a_c);
    }

    double
    BsToKstarLeptonNeutrino::integrated_s_3(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j3 / decay_width(a_c);
    }

    double
    BsToKstarLeptonNeutrino::integrated_s_4(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j4 / decay_width(a_c);
    }

    double
    BsToKstarLeptonNeutrino::integrated_s_5(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j5 / decay_width(a_c);
    }

    double
    BsToKstarLeptonNeutrino::integrated_s_6s(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j6s / decay_width(a_c);
    }

    double
    BsToKstarLeptonNeutrino::four_differential_decay_width(const double & s, const double & c_theta_l, const double & c_theta_k, const double & phi) const
    {
        // compute d^4 Gamma, cf. [BHvD2010], p. 5, Eq. (2.6)
        // Cosine squared of the angles
        double c_theta_k_2 = c_theta_k * c_theta_k;
        double c_theta_l_2 = c_theta_l * c_theta_l;
        double c_phi = cos(phi);
        // Sine squared of the angles
        double s_theta_k_2 = 1.0 - c_theta_k_2;
        double s_theta_l_2 = 1.0 - c_theta_l_2;
        // Sine of the angles
        double s_theta_k = sqrt(s_theta_k_2);
        double s_theta_l = sqrt(s_theta_l_2);
        double s_phi = sin(phi);
        // Cosine of twice the angle
        //double c_2_theta_k = 2.0 * c_theta_k_2 - 1.0;
        double c_2_theta_l = 2.0 * c_theta_l_2 - 1.0;
        double c_2_phi = cos(2.0 * phi);
        // Sine of twice the angle
        double s_2_theta_k = 2.0 * s_theta_k * c_theta_k;
        double s_2_theta_l = 2.0 * s_theta_l * c_theta_l;
        double s_2_phi = sin(2.0 * phi);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        double Gamma = decay_width(_imp->integrated_angular_coefficients(0.02, 19.71));

        double result = 3.0 / 8.0 / M_PI * (
                 a_c.j1s + (a_c.j1c - a_c.j1s) * c_theta_k_2
                +  (a_c.j2s + (a_c.j2c - a_c.j2s) * c_theta_k_2) * c_2_theta_l
                +  a_c.j3 * s_theta_k_2 * s_theta_l_2 * c_2_phi
                +  a_c.j4 * s_2_theta_k * s_2_theta_l * c_phi
                +  a_c.j5 * s_2_theta_k * s_theta_l * c_phi
                +  (a_c.j6s * s_theta_k_2 + a_c.j6c * c_theta_k_2) * c_theta_l
                +  a_c.j7 * s_2_theta_k * s_theta_l * s_phi
                +  a_c.j8 * s_2_theta_k * s_2_theta_l * s_phi
                +  a_c.j9 * s_theta_k_2 * s_theta_l_2 * s_2_phi
                ) / Gamma;

        return result;
    }

    std::vector<OptionSpecification>::const_iterator
    BsToKstarLeptonNeutrino::begin_options()
    {
        return Implementation<BsToKstarLeptonNeutrino>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BsToKstarLeptonNeutrino::end_options()
    {
        return Implementation<BsToKstarLeptonNeutrino>::options.cend();
    }

    template <>
    struct Implementation<BsToKstarLeptonNeutrinoRatios>
    {
        UsedParameter hbar;

        UsedParameter tau;

        BsToKstarLeptonNeutrino bstokstarlnu;

        BToLeptonNeutrino btolnu;

        BToPseudoscalarLeptonNeutrino btopilnu;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            hbar(p["QM::hbar"], u),
            tau(p["life_time::B_d"], u),
            bstokstarlnu(p, o + Options{ std::make_pair("l", "mu")}),
            btolnu(p, o + Options{ std::make_pair("l", "tau") }),
            btopilnu(p, o + Options{ { "l", "mu" }, { "U", "u" }, { "q", "d" } })
        {
            u.uses(bstokstarlnu);
            u.uses(btolnu);
            u.uses(btopilnu);
        }
    };

    const std::vector<OptionSpecification>
    Implementation<BsToKstarLeptonNeutrinoRatios>::options
    {
    };

    BsToKstarLeptonNeutrinoRatios::BsToKstarLeptonNeutrinoRatios(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BsToKstarLeptonNeutrinoRatios>(new Implementation<BsToKstarLeptonNeutrinoRatios>(parameters, options, *this))
    {
    }

    BsToKstarLeptonNeutrinoRatios::~BsToKstarLeptonNeutrinoRatios()
    {
    }

    double
    BsToKstarLeptonNeutrinoRatios::ratio_long() const
    {
        AngularCoefficients a_c = _imp->bstokstarlnu._imp->integrated_angular_coefficients(0.02, 19.71);

        return (a_c.j1c - a_c.j2c / 3.0) / _imp->btolnu.decay_width();
    }

    double
    BsToKstarLeptonNeutrinoRatios::ratio_para() const
    {
        AngularCoefficients a_c = _imp->bstokstarlnu._imp->integrated_angular_coefficients(0.02, 19.71);

        return 4.0 / 9.0 * (2.0 * a_c.j1s - 3.0 * a_c.j3) / _imp->btolnu.decay_width();
    }

    double
    BsToKstarLeptonNeutrinoRatios::ratio_perp() const
    {
        AngularCoefficients a_c = _imp->bstokstarlnu._imp->integrated_angular_coefficients(0.02, 19.71);

        return 4.0 / 9.0 * (2.0 * a_c.j1s + 3.0 * a_c.j3) * _imp->tau() / _imp->hbar() / _imp->btopilnu.integrated_branching_ratio(0.02, 12.0);
    }

    const std::string
    BsToKstarLeptonNeutrino::description = "\
The decay B_s->K^* l nubar, where l=e,mu is a light lepton.";

    const std::string
    BsToKstarLeptonNeutrino::kinematics_description_s = "\
The invariant mass of the l-nubar pair in GeV^2.";

    const std::string
    BsToKstarLeptonNeutrino::kinematics_description_c_theta_l = "\
The cosine of the charged lepton l's helicity angle theta_l in the l-nubar_l rest frame.";

    const std::string
    BsToKstarLeptonNeutrino::kinematics_description_c_theta_k = "\
The cosine of the kaon's helicity angle theta_k in the K-pi rest frame of the decaying K^*.";

    const std::string
    BsToKstarLeptonNeutrino::kinematics_description_phi = "\
The azimuthal angle between the l-nubar plane and the K-pi plane.";

    const std::set<ReferenceName>
    BsToKstarLeptonNeutrino::references
    {
    };

    const std::set<ReferenceName>
    BsToKstarLeptonNeutrinoRatios::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    BsToKstarLeptonNeutrinoRatios::begin_options()
    {
        return Implementation<BsToKstarLeptonNeutrinoRatios>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BsToKstarLeptonNeutrinoRatios::end_options()
    {
        return Implementation<BsToKstarLeptonNeutrinoRatios>::options.cend();
    }
}
