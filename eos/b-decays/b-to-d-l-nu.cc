/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015, 2016, 2017 Danny van Dyk
 * Copyright (c) 2015 Marzia Bordone
 * Copyright (c) 2018, 2019 Ahmet Kokulu
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

#include <eos/form-factors/form-factors.hh>
#include <eos/b-decays/b-to-d-l-nu.hh>
#include <eos/utils/integrate.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/model.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/save.hh>

namespace eos
{
    template <>
    struct Implementation<BToDLeptonNeutrino>
    {
        std::shared_ptr<Model> model;

        std::shared_ptr<FormFactors<PToP>> form_factors;

        Parameters parameters;

        SwitchOption opt_q;

        UsedParameter m_B;

        UsedParameter mu;

        UsedParameter tau_B;

        UsedParameter m_D;

        SwitchOption opt_l;

        UsedParameter m_l;

        UsedParameter g_fermi;

        UsedParameter hbar;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "SM"), p, o)),
            parameters(p),
            opt_q(o, "q", { "u", "d" }, "d"),
            m_B(p["mass::B_" + opt_q.value()], u),
            mu(p["mu"], u),
            tau_B(p["life_time::B_" + opt_q.value()], u),
            m_D(p["mass::D_" + opt_q.value()], u),
            opt_l(o, "l", {"e", "mu", "tau"}, "mu"),
            m_l(p["mass::" + opt_l.value()], u),
            g_fermi(p["G_Fermi"], u),
            hbar(p["hbar"], u)
        {
            form_factors = FormFactorFactory<PToP>::create("B->D::" + o.get("form-factors", "BCL2008"), p, o);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            u.uses(*form_factors);
            u.uses(*model);
        }

        //* normalized(|Vcb|=1) two-fold-distribution, cf. [DSD2014], eq. (12), p. 6
        double normalized_two_differential_decay_width(const double & s, const double & c_theta_l) const
        {
            //  d^2 Gamma, cf. [DSD2014], p. 6, eq. (13)
            // Trigonometric identities
            double c_theta_l_2 = c_theta_l * c_theta_l;
            double s_theta_l_2 = 1.0 - c_theta_l_2;
            double s_theta_l = std::sqrt(s_theta_l_2);
            double c_2_theta_l = 2.0 * c_theta_l_2 - 1.0;
            double s_2_theta_l = 2.0 * s_theta_l * c_theta_l;

            double fp = form_factors->f_p(s);
            double f0 = form_factors->f_0(s);
            double fT = form_factors->f_t(s);
            // running quark masses
            double mbatmu = model->m_b_msbar(mu);
            double mcatmu = model->m_c_msbar(mu);

            double lam = lambda(m_B * m_B, m_D * m_D, s);
            double p = sqrt(lam) / (2.0 * m_B);
            // make sure we return NaN if s < m_l^2, v = lepton velocity in the dilepton rest frame
            double v = (1.0 - m_l * m_l / s);
            double v2 = v * v;
            double ml_hat = sqrt(1.0 - v);
            // universal electroweak correction, cf. [S1982]
            double etaEW = 1.0066;
            double nD = v2 * m_B * s * power_of<2>(g_fermi() * etaEW) / (256.0 * power_of<3>(M_PI * m_B));
            // NP contributions in EFT including tensor operator (cf. [DSD2014]). (in SM cvl=1, and all other couplings are zero)
            const WilsonCoefficients<BToC> wc = model->wilson_coefficients_b_to_c(opt_l.value(), false);
            const complex<double> vl = wc.cvl() - 1.0;
            const complex<double> vr = wc.cvr();
            const complex<double> sl = wc.csl();
            const complex<double> sr = wc.csr();
            const complex<double> gV = vr + vl;
            const complex<double> gA = vr - vl;
            const complex<double> gS = sr + sl;
            const complex<double> gP = sr - sl;
            const complex<double> tl = wc.ct();
            // helicity amplitudes, cf. [DSD2014] eqs. 13-14
            const complex<double> hh0 = 2.0 * m_B * p * fp * (1.0 + gV) / sqrt(s) ;
            const complex<double> hht = (1.0 + gV) * (m_B * m_B - m_D * m_D) * f0 / sqrt(s);
            const complex<double> hhS = - gS * (m_B * m_B - m_D * m_D) * f0 / (mbatmu - mcatmu);
            const complex<double> hhT = - 2.0 * m_B * p * fT * tl / (m_B + m_D);
            const complex<double> hhtS = hht - hhS / ml_hat;

            double result = 2.0 * nD * p * (
                                            std::norm(hh0) * s_theta_l_2 + (1.0 - v) * power_of<2>(std::abs(hh0) * c_theta_l - std::abs(hhtS))
                                            + 8.0 * ( ((2.0 - v) + v * c_2_theta_l) * std::norm(hhT) - ml_hat * std::real(hhT * (std::conj(hh0) - std::conj(hhtS) * c_theta_l)))
                                            );

            return result;
        }

        // normalized to V_cb = 1, obtained using cf. [DSD2014], eq. (12), agrees with Sakaki'13 et al cf. [STTW2013]
        double normalized_differential_decay_width(const double & s) const
        {
            double fp = form_factors->f_p(s);
            double f0 = form_factors->f_0(s);
            double fT = form_factors->f_t(s);
            // running quark masses
            double mbatmu = model->m_b_msbar(mu);
            double mcatmu = model->m_c_msbar(mu);

            double lam = lambda(m_B * m_B, m_D * m_D, s);
            double p = sqrt(lam) / (2.0 * m_B);
            // make sure we return NaN if s < m_l^2, v = lepton velocity in the dilepton rest frame
            double v = (1.0 - m_l * m_l / s);
            double v2 = v * v;
            double ml_hat = sqrt(1.0 - v);
            // universal electroweak correction, cf. [S1982]
            double etaEW = 1.0066;
            double nD = v2 * m_B * s * power_of<2>(g_fermi() * etaEW) / (256.0 * power_of<3>(M_PI * m_B));
            // NP contributions in EFT including tensor operator (cf. [DSD2014]). (in SM cvl=1, and all other couplings are zero)
            const WilsonCoefficients<BToC> wc = model->wilson_coefficients_b_to_c(opt_l.value(), false);
            const complex<double> vl = wc.cvl() - 1.0;
            const complex<double> vr = wc.cvr();
            const complex<double> sl = wc.csl();
            const complex<double> sr = wc.csr();
            const complex<double> gV = vr + vl;
            const complex<double> gA = vr - vl;
            const complex<double> gS = sr + sl;
            const complex<double> gP = sr - sl;
            const complex<double> tl = wc.ct();
            // helicity amplitudes, cf. [DSD2014] eqs. 13-14
            const complex<double> hh0 = 2.0 * m_B * p * fp * (1.0 + gV) / sqrt(s) ;
            const complex<double> hht = (1.0 + gV) * (m_B * m_B - m_D * m_D) * f0 / sqrt(s);
            const complex<double> hhS = - gS * (m_B * m_B - m_D * m_D) * f0 / (mbatmu - mcatmu);
            const complex<double> hhT = - 2.0 * m_B * p * fT * tl / (m_B + m_D);
            const complex<double> hhtS = hht - hhS / ml_hat;

            // normalized(|V_cb|=1) differential decay width
            return 4.0 / 3.0 * nD * p * ( std::norm(hh0) * (3.0 - v) + 3.0 * std::norm(hhtS) * (1.0 - v) + 16.0 * std::norm(hhT) * (3.0 - 2.0 * v) - 24.0 * ml_hat * std::real(hhT * std::conj(hh0)) );
        }

        // obtained using cf. [DSD2014], eq. (12), defined as int_1^0 d^2Gamma - int_0^-1 d^2Gamma
        double numerator_differential_a_fb_leptonic(const double & s) const
        {
            double fp = form_factors->f_p(s);
            double f0 = form_factors->f_0(s);
            double fT = form_factors->f_t(s);
            // running quark masses
            double mbatmu = model->m_b_msbar(mu);
            double mcatmu = model->m_c_msbar(mu);

            double lam = lambda(m_B * m_B, m_D * m_D, s);
            double p = sqrt(lam) / (2.0 * m_B);
            // make sure we return NaN if s < m_l^2, v = lepton velocity in the dilepton rest frame
            double v = (1.0 - m_l * m_l / s);
            double v2 = v * v;
            double ml_hat = sqrt(1.0 - v);
            // universal electroweak correction, cf. [S1982]
            double etaEW = 1.0066;
            double nD = v2 * m_B * s * power_of<2>(g_fermi() * etaEW) / (256.0 * power_of<3>(M_PI * m_B));
            // NP contributions in EFT including tensor operator (cf. [DSD2014]). (in SM cvl=1, and all other couplings are zero)
            const WilsonCoefficients<BToC> wc = model->wilson_coefficients_b_to_c(opt_l.value(), false);
            const complex<double> vl = wc.cvl() - 1.0;
            const complex<double> vr = wc.cvr();
            const complex<double> sl = wc.csl();
            const complex<double> sr = wc.csr();
            const complex<double> gV = vr + vl;
            const complex<double> gA = vr - vl;
            const complex<double> gS = sr + sl;
            const complex<double> gP = sr - sl;
            const complex<double> tl = wc.ct();
            // helicity amplitudes, cf. [DSD2014] eqs. 13-14
            const complex<double> hh0 = 2.0 * m_B * p * fp * (1.0 + gV) / sqrt(s) ;
            const complex<double> hht = (1.0 + gV) * (m_B * m_B - m_D * m_D) * f0 / sqrt(s);
            const complex<double> hhS = - gS * (m_B * m_B - m_D * m_D) * f0 / (mbatmu - mcatmu);
            const complex<double> hhT = - 2.0 * m_B * p * fT * tl / (m_B + m_D);
            const complex<double> hhtS = hht - hhS / ml_hat;

            return - 4.0 * nD * p * ( std::abs(hh0) * std::abs(hhtS) * (1.0 - v) - 4.0 * ml_hat * std::real(hhT * std::conj(hhtS)) );
        }

        // differential decay width
        double differential_decay_width(const double & s) const
        {
            return normalized_differential_decay_width(s) * std::norm(model->ckm_cb());
        }

        // differential branching_ratio
        double differential_branching_ratio(const double & s) const
        {
            return differential_decay_width(s) * tau_B / hbar;
        }

        // "normalized"(|Vcb|=1) differential branching_ratio
        double normalized_differential_branching_ratio(const double & s) const
        {
            return normalized_differential_decay_width(s) * tau_B / hbar;
        }

        double pdf_q2(const double & q2) const
        {
            const double q2_min = power_of<2>(m_l());
            const double q2_max = power_of<2>(m_B() - m_D());

            std::function<double (const double &)> f = std::bind(&Implementation<BToDLeptonNeutrino>::normalized_differential_branching_ratio, this, std::placeholders::_1);
            const double num   = normalized_differential_branching_ratio(q2);
            const double denom = integrate1D(f, 32, q2_min, q2_max);

            return num / denom;
        }

        double pdf_w(const double & w) const
        {
            const double m_B = this->m_B(), m_B2 = m_B * m_B;
            const double m_D = this->m_D(), m_D2 = m_D * m_D;
            const double q2  = m_B2 + m_D2 - 2.0 * m_B * m_D * w;

            return 2.0 * m_B * m_D * pdf_q2(q2);
        }

        double integrated_pdf_q2(const double & q2_min, const double & q2_max) const
        {
            const double q2_abs_min = power_of<2>(m_l());
            const double q2_abs_max = power_of<2>(m_B() - m_D());

            std::function<double (const double &)> f = std::bind(&Implementation<BToDLeptonNeutrino>::normalized_differential_branching_ratio, this, std::placeholders::_1);
            const double num   = integrate<GSL::QAGS>(f, q2_min,     q2_max);
            const double denom = integrate<GSL::QAGS>(f, q2_abs_min, q2_abs_max);

            return num / denom;
        }

        double integrated_pdf_w(const double & w_min, const double & w_max) const
        {
            const double m_B    = this->m_B(), m_B2 = m_B * m_B;
            const double m_D    = this->m_D(), m_D2 = m_D * m_D;
            const double q2_max = m_B2 + m_D2 - 2.0 * m_B * m_D * w_min;
            const double q2_min = m_B2 + m_D2 - 2.0 * m_B * m_D * w_max;

            return integrated_pdf_q2(q2_min, q2_max) / (w_max - w_min);
        }

        double lepton_polarization_numerator(const double & q2) const
        {
            const double m_l2 = m_l() * m_l();
            const double m_B  = this->m_B(), m_B2 = m_B * m_B;
            const double m_D  = this->m_D(), m_D2 = m_D * m_D;
            const double p_D  = sqrt(eos::lambda(m_B2, m_D2, q2)) / (2.0 * m_B);
            const double sqrt_q2  = sqrt(q2);

            const double f_p = form_factors->f_p(q2);
            const double f_0 = form_factors->f_0(q2);

            // cf. [CJLP2012]
            const double H_0 = 2.0 * m_B * p_D / sqrt_q2 * f_p;
            const double H_t = (m_B2 - m_D2) / sqrt_q2 * f_0;

            const double H_02 = H_0 * H_0;
            const double H_t2 = H_t * H_t;

            const double nf = p_D * q2 * pow(1.0 - m_l2 / q2, 2);

            // cf. [CJLP2012]], eq. (20), p. 13
            const double num  = H_02 * (1.0 - 0.5 * m_l2 / q2) - 3.0 / 2.0 * m_l2 / q2 * H_t2;

            return nf * num;
        }

        double lepton_polarization_denominator(const double & q2) const
        {
            const double m_l2 = m_l() * m_l();
            const double m_B  = this->m_B(), m_B2 = m_B * m_B;
            const double m_D  = this->m_D(), m_D2 = m_D * m_D;
            const double p_D  = sqrt(eos::lambda(m_B2, m_D2, q2)) / (2.0 * m_B);
            const double sqrt_q2  = sqrt(q2);

            const double f_p = form_factors->f_p(q2);
            const double f_0 = form_factors->f_0(q2);

            // cf. [CJLP2012]
            const double H_0 = 2.0 * m_B * p_D / sqrt_q2 * f_p;
            const double H_t = (m_B2 - m_D2) / sqrt_q2 * f_0;

            const double H_02 = H_0 * H_0;
            const double H_t2 = H_t * H_t;

            const double nf = p_D * q2 * pow(1.0 - m_l2 / q2, 2);

            // cf. [CJLP2012]], eq. (20), p. 13
            const double denom  = H_02 * (1.0 + 0.5 * m_l2 / q2) + 3.0 / 2.0 * m_l2 / q2 * H_t2;

            return nf * denom;
        }

        double lepton_polarization(const double & q2_min, const double & q2_max) const
        {
            std::function<double (const double &)> integrand_num   = std::bind(&Implementation<BToDLeptonNeutrino>::lepton_polarization_numerator,   this, std::placeholders::_1);
            std::function<double (const double &)> integrand_denom = std::bind(&Implementation<BToDLeptonNeutrino>::lepton_polarization_denominator, this, std::placeholders::_1);
            const auto num   = integrate<GSL::QAGS>(integrand_num,   q2_min, q2_max);
            const auto denom = integrate<GSL::QAGS>(integrand_denom, q2_min, q2_max);

            return num / denom;
        }
    };

    BToDLeptonNeutrino::BToDLeptonNeutrino(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToDLeptonNeutrino>(new Implementation<BToDLeptonNeutrino>(parameters, options, *this))
    {
    }

    BToDLeptonNeutrino::~BToDLeptonNeutrino()
    {
    }

    //* normalized(|Vcb|=1) two-fold-distribution, cf. [DSD2014], eq. (13), p. 6
    double
    BToDLeptonNeutrino::normalized_two_differential_decay_width(const double & s, const double & c_theta_l) const
    {
        return _imp->normalized_two_differential_decay_width(s, c_theta_l);
    }

    double
    BToDLeptonNeutrino::differential_branching_ratio(const double & s) const
    {
        return _imp->differential_branching_ratio(s);
    }

    double
    BToDLeptonNeutrino::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToDLeptonNeutrino>::differential_branching_ratio,
                _imp.get(), std::placeholders::_1);

        return integrate<GSL::QAGS>(f, s_min, s_max);
    }

    // normalized_differential_branching_ratio (|V_cb|=1)
    double
    BToDLeptonNeutrino::normalized_differential_branching_ratio(const double & s) const
    {
        return _imp->normalized_differential_branching_ratio(s);
    }

    // normalized(|Vcb|=1) integrated branching_ratio
    double
    BToDLeptonNeutrino::normalized_integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToDLeptonNeutrino>::normalized_differential_branching_ratio,
                                                             _imp.get(), std::placeholders::_1);
        
        return integrate<GSL::QAGS>(f, s_min, s_max);
    }

    double
    BToDLeptonNeutrino::differential_a_fb_leptonic(const double & s) const
    {
        return _imp->numerator_differential_a_fb_leptonic(s) / _imp->normalized_differential_decay_width(s);
    }

    double
    BToDLeptonNeutrino::integrated_lepton_polarization(const double & q2_min, const double & q2_max) const
    {
        return _imp->lepton_polarization(q2_min, q2_max);
    }

    double
    BToDLeptonNeutrino::integrated_a_fb_leptonic(const double & s_min, const double & s_max) const
    {
        double integrated_numerator;
        {
            std::function<double (const double &)> f = std::bind(&Implementation<BToDLeptonNeutrino>::numerator_differential_a_fb_leptonic,
                                                                 _imp.get(), std::placeholders::_1);
            integrated_numerator = integrate<GSL::QAGS>(f, s_min, s_max);
        }

        double integrated_denominator;
        {
            std::function<double (const double &)> f = std::bind(&Implementation<BToDLeptonNeutrino>::normalized_differential_decay_width,
                                                                 _imp.get(), std::placeholders::_1);
            integrated_denominator = integrate<GSL::QAGS>(f, s_min, s_max);
        }

        return integrated_numerator / integrated_denominator;
    }

    double
    BToDLeptonNeutrino::differential_r_d(const double & s) const
    {
        double br_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            Save<std::string> save_opt_l(_imp->opt_l._value, "mu");
            br_muons = _imp->differential_branching_ratio(s);
        }

        double br_taus;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::tau"]());
            Save<std::string> save_opt_l(_imp->opt_l._value, "tau");
            br_taus = _imp->differential_branching_ratio(s);
        }

        return br_taus / br_muons;
    }

    double
    BToDLeptonNeutrino::integrated_r_d(const double & s_min_mu, const double & s_min_tau, const double & s_max_mu, const double & s_max_tau) const
    {
        std::function<double (const double &)> f = std::bind(&Implementation<BToDLeptonNeutrino>::differential_branching_ratio,
                _imp.get(), std::placeholders::_1);

        double br_muons;
        {

            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            Save<std::string> save_opt_l(_imp->opt_l._value, "mu");

            br_muons = integrate<GSL::QAGS>(f, s_min_mu, s_max_mu);
        }

        double br_taus;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::tau"]());
            Save<std::string> save_opt_l(_imp->opt_l._value, "tau");

            br_taus = integrate<GSL::QAGS>(f, s_min_tau, s_max_tau);
        }

        return br_taus / br_muons;
    }

    double
    BToDLeptonNeutrino::differential_pdf_w(const double & w) const
    {
        return _imp->pdf_w(w);
    }

    double
    BToDLeptonNeutrino::integrated_pdf_w(const double & w_min, const double & w_max) const
    {
        return _imp->integrated_pdf_w(w_min, w_max);
    }

    const std::string
    BToDLeptonNeutrino::description = "\
    The decay B->D l nu, where l=e,mu,tau is a lepton.";

    const std::string
    BToDLeptonNeutrino::kinematics_description_s = "\
    The invariant mass of the l-nubar pair in GeV^2.";

    const std::string
    BToDLeptonNeutrino::kinematics_description_c_theta_l = "\
    The cosine of the polar angle theta_l between the charged lepton and the direction opposite to D meson in the l-nubar rest frame.";

}
