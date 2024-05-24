/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2024 Danny van Dyk
 * Copyright (c) 2021 Méril Reboud
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
#include <eos/maths/power-of.hh>
#include <eos/s-decays/k-to-pi-ll-base.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>


namespace eos
{
    using std::abs;
    using std::norm;
    using std::sqrt;

    struct KToPiDilepton::AngularCoefficients
    {
        double a_l, b_l, c_l;

        AngularCoefficients()
        {
        }

        AngularCoefficients(const std::array<double, 3> & a) :
            a_l(a[0]),
            b_l(a[1]),
            c_l(a[2])
        {
        }
    };

    /*!
     * Implementation for the decay @f$\bar{B} \to \bar{K} \ell^+ \ell^-@f$.
     */
    template <>
    struct Implementation<KToPiDilepton>
    {
        std::shared_ptr<KToPiDilepton::AmplitudeGenerator> amplitude_generator;

        std::shared_ptr<Model> model;

        LeptonFlavorOption opt_l;
        QuarkFlavorOption opt_q;

        UsedParameter hbar;
        UsedParameter m_K;
        UsedParameter m_pi;
        UsedParameter m_l;
        UsedParameter tau;
        UsedParameter mu;

        static const std::vector<OptionSpecification> options;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "WET"), p, o)),
            opt_l(o, options, "l"),
            opt_q(o, options, "q"),
            hbar(p["QM::hbar"], u),
            m_K(p["mass::K^-"], u),
            m_pi(p["mass::pi^-"], u),
            m_l(p["mass::" + opt_l.str()], u),
            tau(p["life_time::K_u" + opt_q.str()], u),
            mu(p["ds" + opt_l.str() + opt_l.str() + "::mu"], u)
        {
            Context ctx("When constructing K->pill observables");

            std::string tag = o.get("tag", "");

            //if ("BFS2004" == tag)
            //{
            //    amplitude_generator.reset(new KToPiDileptonAmplitudes<tag::BFS2004>(p, o));
            //}
            //else
            {
                throw InternalError("KToPiDilepton: Unknown tag or no valid tag specified (tag = '" + tag + "')!");
            }

            u.uses(*amplitude_generator);
        }

        ~Implementation()
        {
        }

        inline std::array<double, 3> angular_coefficients_array(const KToPiDilepton::Amplitudes & A, const double & s) const
        {
            // cf. [BHP2007], Eq. (4.2) - (4.4)
            std::array<double, 3> result;

            // a_l
            result[0] = amplitude_generator->normalisation(s) * (
                s * (power_of<2>(beta_l(s)) * norm(A.F_S) + norm(A.F_P))
                + 0.25 * amplitude_generator->lambda(s) * (norm(A.F_A) + norm(A.F_V))
                + 2.0 * m_l * (m_K() * m_K() - m_pi() * m_pi() + s) * std::real(A.F_P * std::conj(A.F_A))
                + 4.0 * m_l * m_l * m_K() * m_K() * norm(A.F_A)
                );

            // b_l
            result[1] = 2.0 * amplitude_generator->normalisation(s) * (
                s * (power_of<2>(beta_l(s)) * std::real(A.F_S * std::conj(A.F_T))
                + std::real(A.F_P * std::conj(A.F_T5)))
                + m_l * (sqrt(amplitude_generator->lambda(s)) * beta_l(s) * std::real(A.F_S * std::conj(A.F_V))
                + (m_K() * m_K() - m_pi() * m_pi() + s) * std::real(A.F_T5 * std::conj(A.F_A)))
                );

            // c_l
            result[2] = amplitude_generator->normalisation(s) * (
                s * (power_of<2>(beta_l(s)) * norm(A.F_T) + norm(A.F_T5))
                - 0.25 * amplitude_generator->lambda(s) * power_of<2>(beta_l(s)) * (norm(A.F_A) + norm(A.F_V))
                + 2.0 * m_l * sqrt(amplitude_generator->lambda(s)) * beta_l(s) * std::real(A.F_T * std::conj(A.F_V))
                );

            return result;
        }

        inline std::array<double, 3> differential_angular_coefficients_array(const double & q2) const
        {
            return angular_coefficients_array(amplitude_generator->amplitudes(q2), q2);
        }

        inline KToPiDilepton::AngularCoefficients differential_angular_coefficients(const double & q2) const
        {
            return KToPiDilepton::AngularCoefficients(differential_angular_coefficients_array(q2));
        }

        // cf. [BHP2007], Eq. (4.8)
        inline double unnormalized_decay_width(const KToPiDilepton::AngularCoefficients & a) const
        {
            return 2.0 * (a.a_l + a.c_l / 3.0);
        }

        inline double differential_branching_ratio(const KToPiDilepton::AngularCoefficients & a) const
        {
            return unnormalized_decay_width(a) * tau() / hbar();
        }

        // cf. [BHP2007], Eq. (4.9)
        inline double differential_flat_term_numerator(const KToPiDilepton::AngularCoefficients & a) const
        {
            return 2.0 * (a.a_l + a.c_l);
        }

        inline double differential_forward_backward_asymmetry_numerator(const KToPiDilepton::AngularCoefficients & a) const
        {
            return a.b_l;
        }

        KToPiDilepton::AngularCoefficients integrated_angular_coefficients(const double & q2_min, const double & q2_max) const
        {
            std::function<std::array<double, 3> (const double &)> integrand =
                    std::bind(&Implementation<KToPiDilepton>::differential_angular_coefficients_array, this, std::placeholders::_1);
            std::array<double, 3> integrated_angular_coefficients_array = integrate1D(integrand, 64, q2_min, q2_max);

            return KToPiDilepton::AngularCoefficients(integrated_angular_coefficients_array);
        }

        inline double beta_l(const double & q2) const
        {
            return sqrt(1.0 - 4.0 * m_l * m_l / q2);
        }

    };

    KToPiDilepton::KToPiDilepton(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<KToPiDilepton>(new Implementation<KToPiDilepton>(parameters, options, *this))
    {
    }

    KToPiDilepton::~KToPiDilepton()
    {
    }

    const std::vector<OptionSpecification>
    Implementation<KToPiDilepton>::options
    {
        Model::option_specification(),
        {"l", { "e", "mu" }, "mu"},
        {"q", { "d", "u" },  "d"}
    };

    double
    KToPiDilepton::differential_branching_ratio(const double & s) const
    {
        return _imp->differential_branching_ratio(_imp->differential_angular_coefficients(s));
    }

    double
    KToPiDilepton::differential_flat_term(const double & s) const
    {
        AngularCoefficients a = _imp->differential_angular_coefficients(s);

        return _imp->differential_flat_term_numerator(a) / _imp->unnormalized_decay_width(a);
    }

    double
    KToPiDilepton::differential_forward_backward_asymmetry(const double & s) const
    {
        AngularCoefficients a = _imp->differential_angular_coefficients(s);

        return _imp->differential_forward_backward_asymmetry_numerator(a) / _imp->unnormalized_decay_width(a);
    }

    double
    KToPiDilepton::two_differential_decay_width(const double & s, const double & c_theta_l_LHCb) const
    {
        AngularCoefficients a = _imp->differential_angular_coefficients(s);

        // using the angular convention of the LHCb experiment
        const double c_theta_l = -c_theta_l_LHCb;

        // cf. [BHP2007], Eq. (4.1)
        return a.a_l + a.b_l * c_theta_l + a.c_l * c_theta_l * c_theta_l;
    }

    // Integrated Observables
    double
    KToPiDilepton::integrated_decay_width(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a = _imp->integrated_angular_coefficients(s_min, s_max);

        return _imp->unnormalized_decay_width(a);
    }

    double
    KToPiDilepton::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a = _imp->integrated_angular_coefficients(s_min, s_max);

        return _imp->differential_branching_ratio(a);
    }

    double
    KToPiDilepton::integrated_flat_term(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a = _imp->integrated_angular_coefficients(s_min, s_max);

        return _imp->differential_flat_term_numerator(a) / _imp->unnormalized_decay_width(a);
    }

    double
    KToPiDilepton::integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a = _imp->integrated_angular_coefficients(s_min, s_max);

        return _imp->differential_forward_backward_asymmetry_numerator(a) / _imp->unnormalized_decay_width(a);

    }


    const std::string
    KToPiDilepton::description = "\
The decay K->pi l^+ l^-, with l=e,mu a charged lepton.";

    const std::string
    KToPiDilepton::kinematics_description_q2 = "\
The invariant mass of the charged lepton pair in GeV^2.";

    const std::string
    KToPiDilepton::kinematics_description_c_theta_l = "\
The cosine of the lepton's helicity angle theta_l in the l^+l^- rest frame using the LHCb convention.";

    /*
     * For diagnostic purposes only!
     */
    KToPiDilepton::Amplitudes
    KToPiDilepton::amplitudes(const double & q2) const
    {
        return _imp->amplitude_generator->amplitudes(q2);
    }

    std::array<double, 3>
    KToPiDilepton::angular_coefficients(const double & q2) const
    {
        return _imp->angular_coefficients_array(_imp->amplitude_generator->amplitudes(q2), q2);
    }

    const std::set<ReferenceName>
    KToPiDilepton::references
    {
    };

    std::vector<OptionSpecification>::const_iterator
    KToPiDilepton::begin_options()
    {
        return Implementation<KToPiDilepton>::options.cbegin();
    }

    std::vector<OptionSpecification>::const_iterator
    KToPiDilepton::end_options()
    {
        return Implementation<KToPiDilepton>::options.cend();
    }
}
