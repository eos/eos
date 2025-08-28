/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2021 Méril Reboud
 * Copyright (c) 2025 Danny van Dyk
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

#include <eos/maths/integrate.hh>
#include <eos/maths/integrate-impl.hh>
#include <eos/rare-b-decays/bs-to-phi-ll-base.hh>
#include <eos/rare-b-decays/bs-to-phi-ll-bfs2004.hh>
#include <eos/rare-b-decays/bs-to-phi-ll-gvdv2020.hh>
#include <eos/rare-b-decays/bs-to-phi-ll-naive.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

namespace eos
{
    struct BsToPhiDilepton::AngularCoefficients
    {
        double j1s, j1c;
        double j2s, j2c;
        double j3;
        double j4;
        double j5;
        double j6s, j6c;
        double j7;
        double j8;
        double j9;

        AngularCoefficients()
        {
        }

        AngularCoefficients(const std::array<double, 12> & a) :
            j1s(a[0]),
            j1c(a[1]),
            j2s(a[2]),
            j2c(a[3]),
            j3(a[4]),
            j4(a[5]),
            j5(a[6]),
            j6s(a[7]),
            j6c(a[8]),
            j7(a[9]),
            j8(a[10]),
            j9(a[11])
        {
        }
    };

    /*!
     * Implementation for the decay @f$\bar{B_s} \to \phi \ell^+ \ell^-@f$.
     */
    template <>
    struct Implementation<BsToPhiDilepton>
    {
        std::shared_ptr<BsToPhiDilepton::AmplitudeGenerator> amplitude_generator;

        std::shared_ptr<Model> model;

        LeptonFlavorOption opt_l;

        UsedParameter hbar;
        UsedParameter m_l;
        UsedParameter tau;
        UsedParameter mu;
        UsedParameter phiBs;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model"_ok, "WET"), p, o)),
            opt_l(o, options, "l"_ok),
            hbar(p["QM::hbar"], u),
            m_l(p["mass::" + opt_l.str()], u),
            tau(p["life_time::B_s"], u),
            mu(p["sb" + opt_l.str() + opt_l.str() + "::mu"], u),
            phiBs(p["B_s::q_over_p_phase"], u)
        {
            Context ctx("When constructing Bs->Phill observables");

            std::string tag = o.get("tag"_ok, "");

            if ("BFS2004" == tag)
            {
                amplitude_generator.reset(new BsToPhiDileptonAmplitudes<tag::BFS2004>(p, o));
            }
            else if ("GvDV2020" == tag)
            {
                amplitude_generator.reset(new BsToPhiDileptonAmplitudes<tag::GvDV2020>(p, o));
            }
            else if ("Naive" == tag)
            {
                amplitude_generator.reset(new BsToPhiDileptonAmplitudes<tag::Naive>(p, o));
            }
            else
            {
                throw InternalError("BsToPhiDilepton: Unknown tag or no valid tag specified (tag = '" + tag + "')!");
            }

            u.uses(*amplitude_generator);
        }

        ~Implementation()
        {
        }

        inline std::array<double, 12> angular_coefficients_array(const BsToPhiDilepton::Amplitudes & A, const double & s) const
        {
            // cf. [BHvD2010], p. 26, eqs. (A1)-(A11)
            // cf. [BHvD2012], app B, eqs. (B1)-(B12)
            std::array<double, 12> result;

            double z = 4.0 * power_of<2>(m_l()) / s;
            double y = m_l / std::sqrt(s);
            double beta2 = 1.0 - z;
            double beta = std::sqrt(beta2);

            // j1s
            result[0] = 3.0 / 4.0 * (
                  (2.0 + beta2) / 4.0 * (norm(A.a_perp_left) + norm(A.a_perp_right) + norm(A.a_para_left) + norm(A.a_para_right))
                  + z * real(A.a_perp_left * conj(A.a_perp_right) + A.a_para_left * conj(A.a_para_right))
                  + 4.0 * beta2 * (norm(A.a_long_perp) + norm(A.a_long_para))
                  + 4.0 * (4.0 - 3.0 * beta2) * (norm(A.a_time_perp) + norm(A.a_time_para))
                  + 8.0 * std::sqrt(2.0) * y * real(
                       (A.a_para_left + A.a_para_right) * conj(A.a_time_para)
                     + (A.a_perp_left + A.a_perp_right) * conj(A.a_time_perp)
                  )
               );
            // j1c
            result[1] = 3.0 / 4.0 * (
                  norm(A.a_long_left) + norm(A.a_long_right)
                  + z * (norm(A.a_time) + 2.0 * real(A.a_long_left * conj(A.a_long_right)))
                  + beta2 * norm(A.a_scal)
                  + 8.0 * (2.0 - beta2) * norm(A.a_time_long)
                  + 8.0 * beta2 * norm(A.a_para_perp)
                  + 16.0 * y * real((A.a_long_left + A.a_long_right) * conj(A.a_time_long))
               );
            // j2s
            result[2] = 3.0 * beta2 / 16.0 * (
                  norm(A.a_perp_left) + norm(A.a_perp_right) + norm(A.a_para_left) + norm(A.a_para_right)
                  - 16.0 * (norm(A.a_time_perp) + norm(A.a_time_para) + norm(A.a_long_perp) + norm(A.a_long_para))
               );
            // j2c
            result[3] = -3.0 * beta2 / 4.0 * (
                  norm(A.a_long_left) + norm(A.a_long_right)
                  - 8.0 * (norm(A.a_time_long) + norm(A.a_para_perp))
               );
            // j3
            result[4] = 3.0 / 8.0 * beta2 * (
                  norm(A.a_perp_left) + norm(A.a_perp_right) - norm(A.a_para_left) - norm(A.a_para_right)
                  + 16.0 * (norm(A.a_time_para) - norm(A.a_time_perp) + norm(A.a_long_para) - norm(A.a_long_perp))
               );
            // j4
            result[5] = 3.0 / (4.0 * std::sqrt(2.0)) * beta2 * real(
                  A.a_long_left * conj(A.a_para_left) + A.a_long_right * conj(A.a_para_right)
                  - 8.0 * std::sqrt(2.0) * (A.a_time_long * conj(A.a_time_para) + A.a_para_perp * conj(A.a_long_para))
               );
            // j5
            result[6] = 3.0 * std::sqrt(2.0) / 4.0 * beta * real(
                  A.a_long_left * conj(A.a_perp_left) - A.a_long_right * conj(A.a_perp_right)
                  - 2.0 * std::sqrt(2.0) * A.a_time_para * conj(A.a_scal)
                  - y * (
                     (A.a_para_left + A.a_para_right) * conj(A.a_scal)
                     + 4.0 * std::sqrt(2.0) * A.a_long_para * conj(A.a_time)
                     - 4.0 * std::sqrt(2.0) * (A.a_long_left - A.a_long_right) * conj(A.a_time_perp)
                     - 4.0 * (A.a_perp_left - A.a_perp_right) * conj(A.a_time_long)
                  )
               );
            // j6s
            result[7] = 3.0 / 2.0 * beta * real(
                  A.a_para_left * conj(A.a_perp_left) - A.a_para_right * conj(A.a_perp_right)
                  + 4.0 * std::sqrt(2.0) * y * (
                       (A.a_perp_left - A.a_perp_right) * conj(A.a_time_para)
                     + (A.a_para_left - A.a_para_right) * conj(A.a_time_perp)
                  )
               );
            // j6c
            result[8] = 3.0 * beta * real(
                  2.0 * A.a_time_long * conj(A.a_scal)
                  + y * (
                     (A.a_long_left + A.a_long_right) * conj(A.a_scal)
                     + 4.0 * A.a_para_perp * conj(A.a_time)
                  )
               );
            // j7
            result[9] = 3.0 * std::sqrt(2.0) / 4.0 * beta * imag(
                  A.a_long_left * conj(A.a_para_left) - A.a_long_right * conj(A.a_para_right)
                  + 2.0 * std::sqrt(2.0) * A.a_time_perp * conj(A.a_scal)
                  + y * (
                     (A.a_perp_left + A.a_perp_right) * conj(A.a_scal)
                     + 4.0 * std::sqrt(2.0) * A.a_long_perp * conj(A.a_time)
                     + 4.0 * std::sqrt(2.0) * (A.a_long_left - A.a_long_right) * conj(A.a_time_para)
                     - 4.0 * (A.a_para_left - A.a_para_right) * conj(A.a_time_long)
                  )
               );
            // j8
            result[10] = 3.0 / 4.0 / std::sqrt(2.0) * beta2 * imag(
                  A.a_long_left * conj(A.a_perp_left) + A.a_long_right * conj(A.a_perp_right)
               );
            // j9
            result[11] = 3.0 / 4.0 * beta2 * imag(
                  conj(A.a_para_left) * A.a_perp_left + conj(A.a_para_right) * A.a_perp_right
               );

            return result;
        }

        inline std::array<double, 12> differential_angular_coefficients_array(const double & s) const
        {
            return angular_coefficients_array(amplitude_generator->amplitudes(s), s);
        }

        inline BsToPhiDilepton::AngularCoefficients differential_angular_coefficients(const double & s) const
        {
            return BsToPhiDilepton::AngularCoefficients(differential_angular_coefficients_array(s));
        }

        BsToPhiDilepton::AngularCoefficients integrated_angular_coefficients(const double & s_min, const double & s_max) const
        {
            std::function<std::array<double, 12> (const double &)> integrand =
                    std::bind(&Implementation<BsToPhiDilepton>::differential_angular_coefficients_array, this, std::placeholders::_1);
            std::array<double, 12> integrated_angular_coefficients_array = integrate1D(integrand, 64, s_min, s_max);

            return BsToPhiDilepton::AngularCoefficients(integrated_angular_coefficients_array);
        }

        inline double decay_width(const BsToPhiDilepton::AngularCoefficients & a_c)
        {
            // cf. [BHvD2010], p. 6, eq. (2.7)
            return 2.0 * a_c.j1s + a_c.j1c - 1.0 / 3.0 * (2.0 * a_c.j2s + a_c.j2c);
        }

        inline double beta_l(const double & s) const
        {
            return std::sqrt(1.0 - 4.0 * m_l * m_l / s);
        }

        double a_fb_zero_crossing() const
        {
            // We trust QCDF results in a validity range from 0.5 GeV^2 < s < 6.0 GeV^2
            static const double min_result = 0.5;
            static const double max_result = 7.0;

            // Use calT_perp / xi_perp = C_7 as start point.
            // Use hard coded values for mu=4.2 GeV and M_B = 5.2795 GeV here.
            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), amplitude_generator->lepton_flavor, amplitude_generator->cp_conjugate);
            const double start = -2.0 * model->m_b_msbar(4.2) * 5.2795 * real(wc.c7() / wc.c9());

            double result = start;
            // clamp result to QCDF validity region
            result = std::max(min_result, result);
            result = std::min(max_result, result);

            // perform a couple of Newton-Raphson steps
            for (unsigned i = 0 ; i < 100 ; ++i)
            {
                double xplus = result * 1.03;
                double xminus = result * 0.97;

                auto a_c_central = differential_angular_coefficients(result);
                double f = a_c_central.j6s + 0.5 * a_c_central.j6c;
                auto a_c_minus   = differential_angular_coefficients(xminus);
                double f_xminus = a_c_minus.j6s + 0.5 * a_c_minus.j6c;
                auto a_c_plus    = differential_angular_coefficients(xplus);
                double f_xplus = a_c_plus.j6s + 0.5 * a_c_plus.j6c;

                double fprime = (f_xplus - f_xminus) / (xplus - xminus);

                if (std::abs(f / fprime) < 1e-8)
                    break;

                result = result - f / fprime;
                // clamp result to QCDF validity region
                result = std::max(min_result, result);
                result = std::min(max_result, result);
            }

            return result;
        }
    };

    BsToPhiDilepton::BsToPhiDilepton(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BsToPhiDilepton>(new Implementation<BsToPhiDilepton>(parameters, options, *this))
    {
    }

    BsToPhiDilepton::~BsToPhiDilepton()
    {
    }

    const std::vector<OptionSpecification>
    Implementation<BsToPhiDilepton>::options
    {
        Model::option_specification(),
        { "l"_ok, { "e", "mu", "tau" }, "mu" }
    };

    double
    BsToPhiDilepton::a_fb_zero_crossing() const
    {
        return _imp->a_fb_zero_crossing();
    }

    double
    BsToPhiDilepton::decay_width(const double & s, const double & c_theta_l, const double & c_theta_k, const double & phi) const
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
        double Gamma = _imp->decay_width(_imp->integrated_angular_coefficients(1.00, 6.00));

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

    double
    BsToPhiDilepton::decay_width_LHCb(const double & s, const double & c_theta_l_LHCb, const double & c_theta_k_LHCb, const double & phi_LHCb) const
    {
        // compute d^4 Gamma, cf. [BHvD2010], p. 5, Eq. (2.6)
        // using the angular convention of the LHCb experiment

        return BsToPhiDilepton::decay_width(s, -c_theta_l_LHCb, +c_theta_k_LHCb, -phi_LHCb);
    }


    double
    BsToPhiDilepton::differential_decay_width(const double & s) const
    {
        return _imp->decay_width(_imp->differential_angular_coefficients(s));
    }

    double
    BsToPhiDilepton::differential_branching_ratio(const double & s) const
    {
        return differential_decay_width(s) * _imp->tau() / _imp->hbar();
    }

    double
    BsToPhiDilepton::differential_forward_backward_asymmetry(const double & s) const
    {
        // cf. [BHvD2010], p. 6, eq. (2.8)
        // cf. [BHvD2012], eq. (A7)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return (a_c.j6s + 0.5 * a_c.j6c) / _imp->decay_width(a_c);
    }

    double
    BsToPhiDilepton::differential_longitudinal_polarisation(const double & s) const
    {
        // cf. [BHvD2012], eq. (A9)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return (a_c.j1c - a_c.j2c / 3.0) / _imp->decay_width(a_c);
    }

    double
    BsToPhiDilepton::differential_transversal_polarisation(const double & s) const
    {
        // cf. [BHvD2012], eq. (A10)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 2.0 * (a_c.j1s - a_c.j2s / 3.0) / _imp->decay_width(a_c);
    }

    double
    BsToPhiDilepton::differential_transverse_asymmetry_2(const double & s) const
    {
        // cf. [BHvD2010], p. 6, eq. (2.10)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 0.5 * a_c.j3 / a_c.j2s;
    }

    double
    BsToPhiDilepton::differential_transverse_asymmetry_3(const double & s) const
    {
        // cf. [BHvD2010], p. 6, eq. (2.11)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return std::sqrt((4.0 * power_of<2>(a_c.j4) + power_of<2>(_imp->beta_l(s) * a_c.j7)) / (-2.0 * a_c.j2c * (2.0 * a_c.j2s + a_c.j3)));
    }

    double
    BsToPhiDilepton::differential_transverse_asymmetry_4(const double & s) const
    {
        // cf. [BHvD2010], p. 6, eq. (2.12)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return std::sqrt((power_of<2>(_imp->beta_l(s) * a_c.j5) + 4.0 * power_of<2>(a_c.j8)) / (4.0 * power_of<2>(a_c.j4) + power_of<2>(_imp->beta_l(s) * a_c.j7)));
    }

    double
    BsToPhiDilepton::differential_transverse_asymmetry_5(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);

        // cf. [BS2011], eq. (34), p. 9 for the massless case
        return std::sqrt(16.0 * power_of<2>(a_c.j2s) - power_of<2>(a_c.j6s) - 4.0 * (power_of<2>(a_c.j3) + power_of<2>(a_c.j9)))
            / 8.0 / a_c.j2s;
    }

    double
    BsToPhiDilepton::differential_transverse_asymmetry_re(const double & s) const
    {
        // cf. [BS2011], eq. (38), p. 10
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 0.25 * _imp->beta_l(s) * a_c.j6s / a_c.j2s;
    }

    double
    BsToPhiDilepton::differential_transverse_asymmetry_im(const double & s) const
    {
        // cf. [BS2011], eq. (30), p. 8
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 0.5 * a_c.j9 / a_c.j2s;
    }

    double
    BsToPhiDilepton::differential_h_1(const double & s) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.13)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return sqrt(2.0) * a_c.j4 / sqrt(-a_c.j2c * (2.0 * a_c.j2s - a_c.j3));
    }

    double
    BsToPhiDilepton::differential_h_2(const double & s) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.14)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return _imp->beta_l(s) * a_c.j5 / sqrt(-2.0 * a_c.j2c * (2.0 * a_c.j2s + a_c.j3));
    }

    double
    BsToPhiDilepton::differential_h_3(const double & s) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.15)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return _imp->beta_l(s) * a_c.j6s / (2.0 * sqrt(power_of<2>(2.0 * a_c.j2s) - power_of<2>(a_c.j3)));
    }

    double
    BsToPhiDilepton::differential_h_4(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return sqrt(2.0) * a_c.j8 / sqrt(-a_c.j2c * (2.0 * a_c.j2s + a_c.j3));
    }

    double
    BsToPhiDilepton::differential_h_5(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return -a_c.j9 / sqrt(power_of<2>(2.0 * a_c.j2s) + power_of<2>(a_c.j3));
    }

    // differential angular coefficients
    double
    BsToPhiDilepton::differential_j_1c(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j1c;
    }

    double
    BsToPhiDilepton::differential_j_1s(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j1s;
    }

    double
    BsToPhiDilepton::differential_j_2c(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j2c;
    }

    double
    BsToPhiDilepton::differential_j_2s(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j2s;
    }

    double
    BsToPhiDilepton::differential_j_3(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j3;
    }

    double
    BsToPhiDilepton::differential_j_4(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j4;
    }

    double
    BsToPhiDilepton::differential_j_5(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j5;
    }

    double
    BsToPhiDilepton::differential_j_6c(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j6c;
    }

    double
    BsToPhiDilepton::differential_j_6s(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j6s;
    }

    double
    BsToPhiDilepton::differential_j_7(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j7;
    }

    double
    BsToPhiDilepton::differential_j_8(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j8;
    }

    double
    BsToPhiDilepton::differential_j_9(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j9;
    }

    double
    BsToPhiDilepton::integrated_decay_width(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return _imp->decay_width(a_c);
    }

    double
    BsToPhiDilepton::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return integrated_decay_width(s_min, s_max) * _imp->tau() / _imp->hbar();
    }


    double
    BsToPhiDilepton::integrated_unnormalized_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        // Convert from asymmetry in the decay width to asymmetry in the BR
        // cf. [PDG2008] : Gamma = hbar / tau_B, pp. 5, 79
        // TODO: remove fixed value of tau_B and use the one from parameter.cc
        static const double Gamma(6.58211899e-22 * 1e-3 / 1.53e-12);

        // cf. [BHvD2010], eq. (2.8), p. 6
        // cf. [BHvD2012], eq. (A7)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);

        return (a_c.j6s + 0.5 * a_c.j6c) / Gamma;
     }

    double
    BsToPhiDilepton::integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.8), p. 6
        // cf. [BHvD2012], eq. (A7)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return (a_c.j6s + 0.5 * a_c.j6c) / _imp->decay_width(a_c);
    }

    double
    BsToPhiDilepton::integrated_longitudinal_polarisation(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2012], eq. (A9)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return (a_c.j1c - a_c.j2c / 3.0) / _imp->decay_width(a_c);
    }

    double
    BsToPhiDilepton::integrated_transversal_polarisation(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2012], eq. (A10)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 2.0 * (a_c.j1s - a_c.j2s / 3.0) / _imp->decay_width(a_c);
    }

    double
    BsToPhiDilepton::integrated_transverse_asymmetry_2(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.10), p. 6
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 0.5 * a_c.j3 / a_c.j2s;
    }

    double
    BsToPhiDilepton::integrated_transverse_asymmetry_3(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.11), p. 6
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);

        return sqrt((4.0 * power_of<2>(a_c.j4) + power_of<2>(a_c.j7)) / (-2.0 * a_c.j2c * (2.0 * a_c.j2s + a_c.j3)));
    }

    double
    BsToPhiDilepton::integrated_transverse_asymmetry_4(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.12), p. 6
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);

        return sqrt((power_of<2>(a_c.j5) + 4.0 * power_of<2>(a_c.j8)) / (4.0 * power_of<2>(a_c.j4) + power_of<2>(a_c.j7)));
    }

    double
    BsToPhiDilepton::integrated_transverse_asymmetry_5(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);

        // cf. [BS2011], eq. (34), p. 9 for the massless case
        return std::sqrt(16.0 * power_of<2>(a_c.j2s) - power_of<2>(a_c.j6s) - 4.0 * (power_of<2>(a_c.j3) + power_of<2>(a_c.j9)))
            / 8.0 / a_c.j2s;
    }

    double
    BsToPhiDilepton::integrated_transverse_asymmetry_re(const double & s_min, const double & s_max) const
    {
        // cf. [BS2011], eq. (38), p. 10
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 0.25 * a_c.j6s / a_c.j2s;
    }

    double
    BsToPhiDilepton::integrated_transverse_asymmetry_im(const double & s_min, const double & s_max) const
    {
        // cf. [BS2011], eq. (30), p. 8
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 0.5 * a_c.j9 / a_c.j2s;
    }

    double
    BsToPhiDilepton::integrated_h_1(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.13)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return sqrt(2.0) * a_c.j4 / sqrt(-a_c.j2c * (2.0 * a_c.j2s - a_c.j3));
    }

    double
    BsToPhiDilepton::integrated_h_2(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.14)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return  a_c.j5 / sqrt(-2.0 * a_c.j2c * (2.0 * a_c.j2s + a_c.j3));
    }

    double
    BsToPhiDilepton::integrated_h_3(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.15)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j6s / (2.0 * sqrt(power_of<2>(2.0 * a_c.j2s) - power_of<2>(a_c.j3)));
    }

    double
    BsToPhiDilepton::integrated_h_4(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return sqrt(2.0) * a_c.j8 / sqrt(-a_c.j2c * (2.0 * a_c.j2s + a_c.j3));
    }

    double
    BsToPhiDilepton::integrated_h_5(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return -a_c.j9 / sqrt(power_of<2>(2.0 * a_c.j2s) + power_of<2>(a_c.j3));
    }

    // integrated angular coefficients
    double
    BsToPhiDilepton::integrated_j_1c(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j1c;
    }

    double
    BsToPhiDilepton::integrated_j_1s(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j1s;
    }

    double
    BsToPhiDilepton::integrated_j_2c(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j2c;
    }

    double
    BsToPhiDilepton::integrated_j_2s(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j2s;
    }

    double
    BsToPhiDilepton::integrated_j_3(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j3;
    }

    double
    BsToPhiDilepton::integrated_j_4(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j4;
    }

    double
    BsToPhiDilepton::integrated_j_5(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j5;
    }

    double
    BsToPhiDilepton::integrated_j_6c(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j6c;
    }

    double
    BsToPhiDilepton::integrated_j_6s(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j6s;
    }

    double
    BsToPhiDilepton::integrated_j_7(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j7;
    }

    double
    BsToPhiDilepton::integrated_j_8(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j8;
    }

    double
    BsToPhiDilepton::integrated_j_9(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j9;
    }

    /*!
     * Probes of symmetry relations in the large-energy limit (q^2 << m_b^2)
     */
    double
    BsToPhiDilepton::differential_symrel_le_a1v(const double & q2) const
    {
        const auto & ag = *_imp->amplitude_generator;
        return power_of<2>(ag.m_B() + ag.m_V()) / (2.0 * ag.m_B() * ag.energy(q2))
            * ag.form_factors->a_1(q2) / ag.form_factors->v(q2);
    }

    double
    BsToPhiDilepton::differential_symrel_le_t1v(const double & q2) const
    {
        const auto & ag = *_imp->amplitude_generator;
        return (ag.m_B() + ag.m_V()) / ag.m_B()
            * ag.form_factors->t_1(q2) / ag.form_factors->v(q2);
    }

    double
    BsToPhiDilepton::differential_symrel_le_t2v(const double & q2) const
    {
        const auto & ag = *_imp->amplitude_generator;
        return (ag.m_B() + ag.m_V()) / (2.0 * ag.energy(q2))
            * ag.form_factors->t_2(q2) / ag.form_factors->v(q2);
    }

    double
    BsToPhiDilepton::real_C9_perp(const double & s) const
    {
        return _imp->amplitude_generator->real_C9_perp(s);
    }
    double
    BsToPhiDilepton::real_C9_para(const double & s) const
    {
        return _imp->amplitude_generator->real_C9_para(s);
    }
    double
    BsToPhiDilepton::imag_C9_perp(const double & s) const
    {
        return _imp->amplitude_generator->imag_C9_perp(s);
    }
    double
    BsToPhiDilepton::imag_C9_para(const double & s) const
    {
        return _imp->amplitude_generator->imag_C9_para(s);
    }

    const std::string
    BsToPhiDilepton::description = "\
The decay Bsbar->phi(-> Kbar K) l^+ l^-, with l=e,mu,tau \
a charged lepton. Various theory models can be selected using \
the 'tag' option";

    const std::string
    BsToPhiDilepton::kinematics_description_s = "\
The invariant mass of the charged lepton pair in GeV^2.";

    const std::string
    BsToPhiDilepton::kinematics_description_c_theta_l = "\
The cosine of the negatively-charged lepton l^-'s helicity angle theta_l in the l^+l^- rest frame.";

    const std::string
    BsToPhiDilepton::kinematics_description_c_theta_k = "\
The cosine of the Kbar's helicity angle theta_k in the Kbar-K rest frame.";

    const std::string
    BsToPhiDilepton::kinematics_description_phi = "\
The azimuthal angle between the Kbar-K plane and the l^+l^- plane.";

    BsToPhiDilepton::Amplitudes
    BsToPhiDilepton::amplitudes(const double & q2) const
    {
        return _imp->amplitude_generator->amplitudes(q2);
    }

    double
    BsToPhiDilepton::m_l() const
    {
        return _imp->m_l;
    }

    double
    BsToPhiDilepton::phiBs() const
    {
        return _imp->phiBs;
    }

    const std::set<ReferenceName>
    BsToPhiDilepton::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    BsToPhiDilepton::begin_options()
    {
        return Implementation<BsToPhiDilepton>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BsToPhiDilepton::end_options()
    {
        return Implementation<BsToPhiDilepton>::options.cend();
    }


    BsToPhiDileptonAndConjugate::BsToPhiDileptonAndConjugate(const Parameters & parameters, const Options & options):
        bstophidilepton(parameters, options + Options{{"cp-conjugate"_ok, "true"}}),
        bstophidilepton_conjugate(parameters, options + Options{{"cp-conjugate"_ok, "false"}})
    {
    }

    BsToPhiDileptonAndConjugate::~BsToPhiDileptonAndConjugate() = default;

    struct BsToPhiDileptonAndConjugate::AngularhCoefficients
    {
        double h1s, h1c;
        double h2s, h2c;
        double h3;
        double h4;
        double h5;
        double h6s, h6c;
        double h7;
        double h8;
        double h9;

        AngularhCoefficients()
        {
        }

        AngularhCoefficients(const std::array<double, 12> & a) :
            h1s(a[0]),
            h1c(a[1]),
            h2s(a[2]),
            h2c(a[3]),
            h3(a[4]),
            h4(a[5]),
            h5(a[6]),
            h6s(a[7]),
            h6c(a[8]),
            h7(a[9]),
            h8(a[10]),
            h9(a[11])
        {
        }
    };


    inline std::array<double, 12>
    BsToPhiDileptonAndConjugate::angular_h_coefficients_array(const BsToPhiDilepton::Amplitudes & A,
                                                              const BsToPhiDilepton::Amplitudes & Atilda,
                                                              const double & s) const
    {
        // cf. [DV2015], eqs. (117)-(128)
        std::array<double, 12> result;

        double m_l = bstophidilepton.m_l();
        double phiBs = bstophidilepton.phiBs();

        double z = 4.0 * power_of<2>(m_l) / s;
        double y = m_l / std::sqrt(s);
        double beta2 = 1.0 - z;
        double beta = std::sqrt(beta2);

        const complex<double> expiphi(cos(phiBs), sin(phiBs));

        // h1s
        result[0] = 3.0 / 4.0 * (
            (2.0 + beta2) / 2.0 * real( expiphi * (
                Atilda.a_perp_left * conj(A.a_perp_left) + Atilda.a_perp_right * conj(A.a_perp_right)
                + Atilda.a_para_left * conj(A.a_para_left) + Atilda.a_para_right * conj(A.a_para_right)
            ))
            + z * real( expiphi * (Atilda.a_perp_left * conj(A.a_perp_right) + Atilda.a_para_left * conj(A.a_para_right))
                - conj(expiphi) * (A.a_perp_left * conj(Atilda.a_perp_right) + A.a_para_left * conj(Atilda.a_para_right))
                )
        );
        // h1c
        result[1] = 3.0 / 2.0 * ( real(
            expiphi * (
                Atilda.a_long_left * conj(A.a_long_left) + Atilda.a_long_right * conj(A.a_long_right)
            ))
            + z * (
                real( expiphi * (Atilda.a_time * conj(A.a_time)))
                + real( expiphi * (Atilda.a_long_left * conj(A.a_long_right)) + conj(expiphi) * (A.a_long_left * conj(Atilda.a_long_right)) )
            + beta2 * real( expiphi * Atilda.a_scal * conj(A.a_scal)))
        );
        // h2s
        result[2] = 3.0 * beta2 / 8.0 * real( expiphi * (
            Atilda.a_perp_left * conj(A.a_perp_left) + Atilda.a_perp_right * conj(A.a_perp_right)
            + Atilda.a_para_left * conj(A.a_para_left) + Atilda.a_para_right * conj(A.a_para_right)
        ));
        // h2c
        result[3] = -3.0 * beta2 / 2.0 * real( expiphi * (
            Atilda.a_long_left * conj(A.a_long_left) + Atilda.a_long_right * conj(A.a_long_right)
        ));
        // h3
        result[4] = 3.0 / 4.0 * beta2 * real( expiphi * (
            Atilda.a_perp_left * conj(A.a_perp_left) + Atilda.a_perp_right * conj(A.a_perp_right)
            - Atilda.a_para_left * conj(A.a_para_left) - Atilda.a_para_right * conj(A.a_para_right)
        ));
        // h4
        result[5] = 3.0 / (4.0 * std::sqrt(2.0)) * beta2 * real(
            expiphi * (
                Atilda.a_long_left * conj(A.a_para_left) + Atilda.a_long_right * conj(A.a_para_right)
            ) + conj(expiphi) * (
                A.a_long_left * conj(Atilda.a_para_left) + A.a_long_right * conj(Atilda.a_para_right)
            ));
        // h5
        result[6] = 3.0 * std::sqrt(2.0) / 4.0 * beta * (real(
            expiphi * (
                Atilda.a_long_left * conj(A.a_perp_left) - Atilda.a_long_right * conj(A.a_perp_right)
            ) + conj(expiphi) * (
                A.a_long_left * conj(Atilda.a_perp_left) - A.a_long_right * conj(Atilda.a_perp_right)
            ))
            - y * real(
            expiphi * (
                Atilda.a_para_left * conj(A.a_scal) + Atilda.a_para_right * conj(A.a_scal)
            ) + conj(expiphi) * (
                A.a_para_left * conj(Atilda.a_scal) + A.a_para_right * conj(Atilda.a_scal)
            )));
        // h6s
        result[7] = 3.0 / 2.0 * beta * real(
            expiphi * (
                Atilda.a_para_left * conj(A.a_perp_left) - Atilda.a_para_right * conj(A.a_perp_right)
            ) + conj(expiphi) * (
                A.a_para_left * conj(Atilda.a_perp_left) - A.a_para_right * conj(Atilda.a_perp_right)
            ));
        // h6c
        result[8] = 3.0 * beta * y * real(
            expiphi * (
                Atilda.a_long_left * conj(A.a_scal) + Atilda.a_long_right * conj(A.a_scal)
            ) + conj(expiphi) * (
                A.a_long_left * conj(Atilda.a_scal) + A.a_long_right * conj(Atilda.a_scal)
            ));
        // h7
        result[9] = 3.0 * std::sqrt(2.0) / 4.0 * beta * (imag(
            expiphi * (
                Atilda.a_long_left * conj(A.a_para_left) - Atilda.a_long_right * conj(A.a_para_right)
            ) + conj(expiphi) * (
                A.a_long_left * conj(Atilda.a_para_left) - A.a_long_right * conj(Atilda.a_para_right)
            ))
            + y * imag(
            expiphi * (
                Atilda.a_perp_left * conj(A.a_scal) + Atilda.a_perp_right * conj(A.a_scal)
            ) + conj(expiphi) * (
                A.a_perp_left * conj(Atilda.a_scal) + A.a_perp_right * conj(Atilda.a_scal)
            )));
        // h8
        result[10] = 3.0 / 4.0 / std::sqrt(2.0) * beta2 * imag(
            expiphi * (
                Atilda.a_long_left * conj(A.a_perp_left) + Atilda.a_long_right * conj(A.a_perp_right)
            ) + conj(expiphi) * (
                A.a_long_left * conj(Atilda.a_perp_left) + A.a_long_right * conj(Atilda.a_perp_right)
            ));
        // h9
        result[11] = - 3.0 / 4.0 * beta2 * imag(
            expiphi * (
                Atilda.a_para_left * conj(A.a_perp_left) + Atilda.a_para_right * conj(A.a_perp_right)
            ) + conj(expiphi) * (
                A.a_para_left * conj(Atilda.a_perp_left) + A.a_para_right * conj(Atilda.a_perp_right)
            ));

        return result;
    }

    inline std::array<double, 12>
    BsToPhiDileptonAndConjugate::differential_angular_h_coefficients_array(const double & s) const
    {
        BsToPhiDilepton::Amplitudes A = bstophidilepton.amplitudes(s);
        BsToPhiDilepton::Amplitudes Atilda = bstophidilepton_conjugate.amplitudes(s);

        Atilda.a_perp_left *= -1;
        Atilda.a_perp_right *= -1;
        Atilda.a_scal *= -1;
        Atilda.a_time_perp *= -1; // TODO Check
        Atilda.a_long_perp *= -1; // TODO Check

        return angular_h_coefficients_array(A, Atilda, s);
    }

    BsToPhiDileptonAndConjugate::AngularhCoefficients
    BsToPhiDileptonAndConjugate::integrated_angular_h_coefficients(const double & s_min, const double & s_max) const
    {
        std::function<std::array<double, 12> (const double &)> integrand =
                std::bind(&BsToPhiDileptonAndConjugate::differential_angular_h_coefficients_array, this, std::placeholders::_1);
        std::array<double, 12> integrated_angular_h_coefficients_array = integrate1D(integrand, 64, s_min, s_max);

        return BsToPhiDileptonAndConjugate::AngularhCoefficients(integrated_angular_h_coefficients_array);
    }

    double
    BsToPhiDileptonAndConjugate::integrated_H_1c(const double & s_min, const double & s_max) const
    {
        AngularhCoefficients a_h_c = integrated_angular_h_coefficients(s_min, s_max);
        return a_h_c.h1c;
    }

    double
    BsToPhiDileptonAndConjugate::integrated_H_1s(const double & s_min, const double & s_max) const
    {
        AngularhCoefficients a_h_c = integrated_angular_h_coefficients(s_min, s_max);
        return a_h_c.h1s;
    }

    double
    BsToPhiDileptonAndConjugate::integrated_H_2c(const double & s_min, const double & s_max) const
    {
        AngularhCoefficients a_h_c = integrated_angular_h_coefficients(s_min, s_max);
        return a_h_c.h2c;
    }

    double
    BsToPhiDileptonAndConjugate::integrated_H_2s(const double & s_min, const double & s_max) const
    {
        AngularhCoefficients a_h_c = integrated_angular_h_coefficients(s_min, s_max);
        return a_h_c.h2s;
    }

    const std::set<ReferenceName>
    BsToPhiDileptonAndConjugate::references
    {
    };

    const std::vector<OptionSpecification>
    BsToPhiDileptonAndConjugate::options
    {
    };


    std::vector<OptionSpecification>::const_iterator
    BsToPhiDileptonAndConjugate::begin_options()
    {
        return options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    BsToPhiDileptonAndConjugate::end_options()
    {
        return options.cend();
    }

}
