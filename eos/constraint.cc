/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2012, 2013, 2014 Danny van Dyk
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

#include <eos/constraint.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/log_likelihood.hh>
#include <eos/utils/observable_set.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <cmath>
#include <map>
#include <vector>

namespace eos
{
    UnknownConstraintError::UnknownConstraintError(const std::string & name) :
        Exception("Constraint '" + name + "' is unknown")
    {
    }

    struct GaussianConstraintTemplate
    {
        std::string observable;

        Kinematics kinematics;

        Options options;

        double central, sigma_hi_stat, sigma_lo_stat, sigma_hi_sys, sigma_lo_sys;

        unsigned number_of_observations;

        GaussianConstraintTemplate(const std::string & observable,
                const Kinematics & kinematics, const Options & options,
                const double & central,
                const double & sigma_hi_stat, const double & sigma_lo_stat,
                const double & sigma_hi_sys, const double & sigma_lo_sys,
                const unsigned & number_of_observations = 1u) :
            observable(observable),
            kinematics(kinematics),
            options(options),
            central(central),
            sigma_hi_stat(sigma_hi_stat),
            sigma_lo_stat(sigma_lo_stat),
            sigma_hi_sys(sigma_hi_sys),
            sigma_lo_sys(sigma_lo_sys),
            number_of_observations(number_of_observations)
        {
        }

        Constraint
        make(const std::string & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);

            ObservablePtr observable = Observable::make(this->observable, parameters, this->kinematics, this->options + options);
            if (! observable.get())
                throw InternalError("make_gaussian_constraint: " + name + ": '" + this->observable + "' is not a valid observable name");

            double min = 0.0, max = 0.0;
            if ("asymmetric+quadratic" == options.get("uncertainty", "asymmetric+quadratic"))
            {
                min = this->central - std::sqrt(power_of<2>(this->sigma_lo_stat) + power_of<2>(this->sigma_lo_sys));
                max = this->central + std::sqrt(power_of<2>(this->sigma_hi_stat) + power_of<2>(this->sigma_hi_sys));
            }

            LogLikelihoodBlockPtr block = LogLikelihoodBlock::Gaussian(cache, observable, min, this->central, max);

            return Constraint(name, { observable }, { block });
        }
    };

    struct AmorosoLimitConstraintTemplate
    {
        std::string observable;

        Kinematics kinematics;

        Options options;

        double physical_limit, upper_limit_90, upper_limit_95;

        double theta, alpha;

        Constraint
        make(const std::string & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);

            ObservablePtr observable = Observable::make(this->observable, parameters, this->kinematics, this->options + options);
            if (! observable.get())
                throw InternalError("make_amoroso_limit_constraint: " + name + ": '" + this->observable + "' is not a valid observable name");


            LogLikelihoodBlockPtr block = LogLikelihoodBlock::AmorosoLimit(cache, observable, physical_limit, upper_limit_90,
                                                                           upper_limit_95, theta, alpha);

            return Constraint(name, { observable }, { block });
        }
    };

    struct AmorosoModeConstraintTemplate
    {
        std::string observable;

        Kinematics kinematics;

        Options options;

        double physical_limit, mode, upper_limit_90, upper_limit_95;

        double theta, alpha, beta;

        Constraint
        make(const std::string & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);

            ObservablePtr observable = Observable::make(this->observable, parameters, this->kinematics, this->options + options);
            if (! observable.get())
                throw InternalError("make_amoroso_constraint: " + name + ": '" + this->observable + "' is not a valid observable name");


            LogLikelihoodBlockPtr block = LogLikelihoodBlock::AmorosoMode(cache, observable, physical_limit,
                                                                      mode, upper_limit_90, upper_limit_95,
                                                                      theta, alpha, beta);

            return Constraint(name, { observable }, { block });
        }
    };

    struct AmorosoTripleLimitConstraintTemplate
    {
        std::string observable;

        Kinematics kinematics;

        Options options;

        double physical_limit, upper_limit_10, upper_limit_50, upper_limit_90;

        double theta, alpha, beta;

        Constraint
        make(const std::string & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);

            ObservablePtr observable = Observable::make(this->observable, parameters, this->kinematics, this->options + options);
            if (! observable.get())
                throw InternalError("make_amoroso_constraint: " + name + ": '" + this->observable + "' is not a valid observable name");


            LogLikelihoodBlockPtr block = LogLikelihoodBlock::Amoroso(cache, observable, physical_limit,
                                                                      upper_limit_10, upper_limit_50, upper_limit_90,
                                                                      theta, alpha, beta);

            return Constraint(name, { observable }, { block });
        }
    };

    struct AmorosoConstraintTemplate
    {
        std::string observable;

        Kinematics kinematics;

        Options options;

        double physical_limit, theta, alpha, beta;

        Constraint
        make(const std::string & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);

            ObservablePtr observable = Observable::make(this->observable, parameters, this->kinematics, this->options + options);
            if (! observable.get())
                throw InternalError("make_amoroso_constraint: " + name + ": '" + this->observable + "' is not a valid observable name");

            LogLikelihoodBlockPtr block = LogLikelihoodBlock::Amoroso(cache, observable, physical_limit, theta, alpha, beta);

            return Constraint(name, { observable }, { block });
        }
    };

    template <size_t dim_>
    struct MultivariateGaussianConstraintTemplate
    {
        std::array<std::string, dim_> observables;

        std::array<Kinematics, dim_> kinematics;

        std::array<Options, dim_> options;

        std::array<double, dim_> means;

        std::array<double, dim_> sigma_stat_hi;
        std::array<double, dim_> sigma_stat_lo;

        std::array<double, dim_> sigma_sys;

        std::array<std::array<double, dim_>, dim_> correlation;

        unsigned number_of_observations;

        MultivariateGaussianConstraintTemplate(const std::array<std::string, dim_> & observables,
            const std::array<Kinematics, dim_> & kinematics,
            const std::array<Options, dim_> & options,
            const std::array<double, dim_> & means,
            const std::array<double, dim_> & sigma_stat_hi,
            const std::array<double, dim_> & sigma_stat_lo,
            const std::array<double, dim_> & sigma_sys,
            const std::array<std::array<double, dim_>, dim_> & correlation,
            const unsigned & number_of_observations = 1u) :
            observables(observables),
            kinematics(kinematics),
            options(options),
            means(means),
            sigma_stat_hi(sigma_stat_hi),
            sigma_stat_lo(sigma_stat_lo),
            sigma_sys(sigma_sys),
            correlation(correlation),
            number_of_observations(number_of_observations)
        {
        }

        Constraint
        make(const std::string & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);

            std::array<ObservablePtr, dim_> observables;
            for (auto i = 0u ; i < dim_ ; ++i)
            {
                observables[i] = Observable::make(this->observables[i], parameters, this->kinematics[i], this->options[i] + options);
                if (! observables[i].get())
                    throw InternalError("make_multivariate_gaussian_constraint<" + stringify(dim_) + ">: " + name + ": '" + this->observables[i] + "' is not a valid observable name");
            }

            std::array<double, dim_> variances;

            if ("symmetric+quadratic" == options.get("uncertainty", "symmetric+quadratic"))
            {
                for (auto i = 0u ; i < dim_ ; ++i)
                {
                    double combined_lo = power_of<2>(sigma_stat_lo[i]) + power_of<2>(sigma_sys[i]);
                    double combined_hi = power_of<2>(sigma_stat_hi[i]) + power_of<2>(sigma_sys[i]);

                    variances[i] = std::max(combined_lo, combined_hi);
                }
            }

            LogLikelihoodBlockPtr block = LogLikelihoodBlock::MultivariateGaussian(cache, observables, this->means, variances, this->correlation);

            return Constraint(name, std::vector<ObservablePtr>(observables.begin(), observables.end()), { block });
        }
    };

    namespace templates
    {
        ///@name 2000 Data
        ///@{
        /*
         * CLEO Collaboration
         *
         * Data taken from [CLEO:2000]
         */
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_gamma_BR_CLEO_2000
        {
            "B->K^*gamma::BRavg",
            Kinematics{ },
            Options{ { "q", "d" } },
            4.55e-5, 0.72e-5, 0.68e-5, 0.34e-5, 0.34e-5
        };
        static const GaussianConstraintTemplate Bplus_to_Kstarplus_gamma_BR_CLEO_2000
        {
            "B->K^*gamma::BRavg",
            Kinematics{ },
            Options{ { "q", "u" } },
            3.76e-5, 0.89e-5, 0.83e-5, 0.28e-5, 0.28e-5
        };
        ///@}

        ///@name 2004 Data
        ///@{
        /*
         * BaBar Collaboration
         *
         * Data taken from [BaBar:2004A]
         */
        static const GaussianConstraintTemplate Bmix_to_Xs_dilepton_BR_BaBar_2004A
        {
            "B->X_sll::BR@HLMW2005",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ },
            1.8e-6, 0.7e-6, 0.7e-6, 0.5e-6, 0.5e-6
        };

        /*
         * Belle Collaboration
         *
         * Data taken from [Belle:2004]
         */
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_gamma_BR_Belle_2004
        {
            "B->K^*gamma::BRavg",
            Kinematics{ },
            Options{ { "q", "d" } },
            4.01e-5, 0.21e-5, 0.21e-5, 0.17e-5, 0.17e-5
        };
        static const GaussianConstraintTemplate Bplus_to_Kstarplus_gamma_BR_Belle_2004
        {
            "B->K^*gamma::BRavg",
            Kinematics{ },
            Options{ { "q", "u" } },
            4.25e-5, 0.31e-5, 0.31e-5, 0.24e-5, 0.24e-5
        };
        ///@}

        ///@name 2005 Data
        ///@{
        static const GaussianConstraintTemplate Bmix_to_Xs_dilepton_BR_Belle_2005A
        {
            "B->X_sll::BR@HLMW2005",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ },
            1.493e-6, 0.504e-6, 0.504e-6, 0.411e-6, 0.321e-6
        };
        ///@}

        ///@name 2006 Data
        ///@{
        /*
         * Belle Collaboration
         *
         * Data taken from [Belle:2006]
         */
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_gamma_SKstargamma_Belle_2006
        {
            "B->K^*gamma::S_K^*gamma",
            Kinematics{ },
            Options{ { "q", "d" } },
            -0.32, +0.36, -0.33, +0.05, -0.05
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_gamma_CKstargamma_Belle_2006
        {
            "B->K^*gamma::C_K^*gamma",
            Kinematics{ },
            Options{ { "q", "d" } },
            +0.20, +0.24, -0.24, +0.05, -0.05
        };
        static const MultivariateGaussianConstraintTemplate<2> Bzero_to_Kstarzero_gamma_time_dependent_cp_asymmetries_Belle_2006
        {
            {{ "B->K^*gamma::S_K^*gamma", "B->K^*gamma::C_K^*gamma" }},
            {{ Kinematics{ }, Kinematics{ } }},
            {{ Options{ { "q", "d" } }, Options{ { "q", "d" } } }},
            {{ -0.32, +0.20 }},
            {{ +0.36, +0.24 }},
            {{ +0.33, +0.24 }},
            {{ +0.05, +0.05 }},
            {{
                {{ 1,  0.08 }},  /* Use correlation of the two measurements from */
                {{ 0.08,  1 }},  /* http://www.slac.stanford.edu/xorg/hfag/triangle/moriond2011/index.shtml#bsgamma */
            }}
        };
        ///@}

        ///@name 2008 Data
        ///@{
        /*
         * BaBar Collaboration
         *
         * Data taken from [BaBar:2008]
         */
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_gamma_SKstargamma_BaBar_2008
        {
            "B->K^*gamma::S_K^*gamma",
            Kinematics{ },
            Options{ { "q", "d" } },
            -0.03, +0.29, -0.29, +0.03, -0.03
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_gamma_CKstargamma_BaBar_2008
        {
            "B->K^*gamma::C_K^*gamma",
            Kinematics{ },
            Options{ { "q", "d" } },
            -0.14, +0.16, -0.16, +0.03, -0.03
        };
        static const MultivariateGaussianConstraintTemplate<2> Bzero_to_Kstarzero_gamma_time_dependent_cp_asymmetries_BaBar_2008
        {
            {{ "B->K^*gamma::S_K^*gamma", "B->K^*gamma::C_K^*gamma" }},
            {{ Kinematics{ }, Kinematics{ } }},
            {{ Options{ { "q", "d" } }, Options{ { "q", "d" } } }},
            {{ -0.03, -0.14 }},
            {{ +0.29, +0.16 }},
            {{ +0.29, +0.16 }},
            {{ +0.03, +0.03 }},
            {{
                {{ 1,  0.05 }},  // Use correlation of the two results from
                {{ 0.05,  1 }},  // http://www.slac.stanford.edu/xorg/hfag/triangle/moriond2011/index.shtml#bsgamma to calculate covariance matrix.
            }}
        };

        /*
         * Belle
         *
         * Data taken from [Belle:2008A]
         */
        static const GaussianConstraintTemplate B_to_Xs_gamma_E_1_1dot8_Belle_2008A
        {
            "B->X_sgamma::E_1(E_min)@NLO",
            Kinematics{ { "E_min", 1.8 } },
            Options{ },
            2.292, 0.027, 0.027, 0.033, 0.033
        };
        static const GaussianConstraintTemplate B_to_Xs_gamma_E_2_1dot8_Belle_2008A
        {
            "B->X_sgamma::E_2(E_min)@NLO",
            Kinematics{ { "E_min", 1.8 } },
            Options{ },
            0.0305, 0.079, 0.079, 0.099, 0.099
        };
        static const MultivariateGaussianConstraintTemplate<2> B_to_Xs_gamma_E_1_and_E_2_1dot8_Belle_2008A
        {
            {{ "B->X_sgamma::E_1(E_min)@NLO", "B->X_sgamma::E_2(E_min)@NLO" }},
            {{ Kinematics{ { "E_min", 1.8 } }, Kinematics{ { "E_min", 1.8 } } }},
            {{ Options{ }, Options{ } }},
            {{  2.292,  0.0305 }},
            {{  0.027,  0.0079 }},
            {{  0.027,  0.0079 }},
            {{  0.033,  0.0099 }},
            {{
                {{ +1.00,  -0.46 }},  // Use correlation of the two results from
                {{ -0.46,  +1.00 }},  // as given in Table VII in [Belle:2008]
            }}
        };
        ///@}

        ///@name 2009 Data
        ///@{
        /*
         * BaBar Collaboration
         *
         * Data taken from [BaBar:2009]
         */
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_gamma_BR_BaBar_2009
        {
            "B->K^*gamma::BRavg",
            Kinematics{ },
            Options{ { "q", "d" } },
            4.47e-5, +0.10e-5, -0.10e-5, +0.16e-5, -0.16e-5
        };
        static const GaussianConstraintTemplate Bplus_to_Kstarplus_gamma_BR_BaBar_2009
        {
            "B->K^*gamma::BRavg",
            Kinematics{ },
            Options{ { "q", "u" } },
            4.22e-5, +0.14e-5, -0.14e-5, +0.16e-5, -0.16e-5
        };

        /*
         * Belle Collaboration
         *
         * Data taken from [Belle:2009], Table I, p.5
         */
        // B^+ -> K^+ mu^+ mu^-
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_1_to_6_Belle_2009
        {
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            1.36e-7, +0.23e-7, -0.21e-7, +0.08e-7, -0.08e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_14dot18_to_16_Belle_2009
        {
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.38e-7, +0.19e-7, -0.12e-7, +0.02e-7, -0.02e-7
        };
        // BR in [16.00, 22.86]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_16_to_22dot86_Belle_2009
        {
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 22.86 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.98e-7, +0.20e-7, -0.19e-7, +0.06e-7, -0.06e-7
        };

        // A_FB in [1.0, 6.0]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_A_FB_1_to_6_Belle_2009
        {
            "B->Kll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            +0.04, +0.16, -0.13, +0.05, -0.05
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_Belle_2009
        {
            "B->Kll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            -0.04, +0.26, -0.32, +0.05, -0.05
        };
        // A_FB in [16.00, 22.86]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_A_FB_16_to_22dot86_Belle_2009
        {
            "B->Kll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 22.86 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            -0.02, +0.08, -0.11, +0.02, -0.02
        };

        // B^0 -> K^*0 mu^+ mu^-
        // BR in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_1_to_6_Belle_2009
        {
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.49e-7, +0.45e-7, -0.40e-7, +0.12e-7, -0.12e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_Belle_2009
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.05e-7, +0.29e-7, -0.26e-7, +0.08e-7, -0.08e-7
        };
        // BR in [16.00, 19.21]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_Belle_2009
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            2.04e-7, +0.27e-7, -0.24e-7, +0.16e-7, -0.16e-7
        };
        // A_FB in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_Belle_2009
        {
            "B->K^*ll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.26, +0.30, -0.27, +0.07, -0.07
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_Belle_2009
        {
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.70, +0.22, -0.16, +0.10, -0.10
        };
        // A_FB in [16.00, 19.21]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_Belle_2009
        {
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.66, +0.16, -0.11, +0.04, -0.04
        };
        // F_L in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_1_to_6_Belle_2009
        {
            "B->K^*ll::F_Lavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.67, +0.23, -0.23, +0.05, -0.05
        };
        // F_L in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_Belle_2009
        {
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.15, +0.27, -0.23, +0.07, -0.07
        };
        // F_L in [16.00, 19.21]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_Belle_2009
        {
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.12, +0.15, -0.13, +0.02, -0.02
        };

        /*
         * Belle Collaboration
         *
         * Data taken from [Belle:2009B], Tables I and II, pp. 7 and 9, respectively.
         */
        static const GaussianConstraintTemplate B_to_Xs_gamma_BR_1dot8_Belle_2009B
        {
            "B->X_sgamma::BR(E_min)@NLO",
            Kinematics{ { "E_min", 1.8 } },
            Options{ },
            3.36e-4, +0.13e-4, -0.13e-4, +0.25e-4, -0.25e-4
        };
        static const MultivariateGaussianConstraintTemplate<3> B_to_Xs_gamma_1dot8_Belle_2009B
        {
            {{ "B->X_sgamma::BR(E_min)@NLO", "B->X_sgamma::E_1(E_min)@NLO", "B->X_sgamma::E_2(E_min)@NLO" }},
            {{ Kinematics{ { "E_min", 1.8 } }, Kinematics{ { "E_min", 1.8 } }, Kinematics{ { "E_min", 1.8 } } }},
            {{ Options{ }, Options{ }, Options{ } }},
            {{ 3.36e-4, 2.294, 0.0370 }},
            {{ 0.13e-4, 0.011, 0.0029 }},
            {{ 0.13e-4, 0.011, 0.0029 }},
            {{ 0.25e-4, 0.028, 0.0081 }},
            {{
                 {{ 1.000, 0.670, 0.800 }},
                 {{ 0.670, 1.000, 0.780 }},
                 {{ 0.800, 0.780, 1.000 }},
            }}
        };
        ///@}

        ///@name 2010 Data
        ///@{
        /*
         * BaBar Collaboration
         *
         * Data taken from [BaBar:2010A], Table X, p. 33.
         */
        static const GaussianConstraintTemplate Bzero_to_pi_l_nu_BR_0_to_4_BaBar_2010A
        {
            "B->pilnu::BR",
            Kinematics{ { "s_min", 0.0 }, { "s_max", 4.0 } },
            Options{ { "q", "d" } },
            0.313e-4, +0.030e-4, -0.030e-4, +0.025e-4, -0.025e-4
        };
        static const GaussianConstraintTemplate Bzero_to_pi_l_nu_BR_4_to_8_BaBar_2010A
        {
            "B->pilnu::BR",
            Kinematics{ { "s_min", 4.0 }, { "s_max", 8.0 } },
            Options{ { "q", "d" } },
            0.329e-4, +0.018e-4, -0.018e-4, +0.016e-4, -0.016e-4
        };
        static const GaussianConstraintTemplate Bzero_to_pi_l_nu_BR_8_to_12_BaBar_2010A
        {
            "B->pilnu::BR",
            Kinematics{ { "s_min", 8.0 }, { "s_max", 12.0 } },
            Options{ { "q", "d" } },
            0.241e-4, +0.018e-4, -0.018e-4, +0.015e-4, -0.015e-4
        };

        /*
         * BaBar Collaboration
         *
         * Data taken from [BaBar:2010B], Table VIII, p. 14. We use the results
         * 'without FSR' (final state radiation). This data supercedes [BaBar:2010A].
         * We combine correlations between the statistic uncertainties and systematic uncertainties.
         */
        static const MultivariateGaussianConstraintTemplate<6> Bzero_to_pi_l_nu_BR_BaBar_2010B
        {
            {{ "B->pilnu::BR", "B->pilnu::BR", "B->pilnu::BR", "B->pilnu::BR", "B->pilnu::BR", "B->pilnu::BR" }},
            {{
                 Kinematics{ { "s_min",  0.0 }, { "s_max",  2.0 } },
                 Kinematics{ { "s_min",  2.0 }, { "s_max",  4.0 } },
                 Kinematics{ { "s_min",  4.0 }, { "s_max",  6.0 } },
                 Kinematics{ { "s_min",  6.0 }, { "s_max",  8.0 } },
                 Kinematics{ { "s_min",  8.0 }, { "s_max", 10.0 } },
                 Kinematics{ { "s_min", 10.0 }, { "s_max", 12.0 } }
            }},
            {{
                 Options{ { "q", "d" } },
                 Options{ { "q", "d" } },
                 Options{ { "q", "d" } },
                 Options{ { "q", "d" } },
                 Options{ { "q", "d" } },
                 Options{ { "q", "d" } }
            }},
            {{ 0.1280e-4, 0.1192e-4, 0.1446e-4, 0.1437e-4, 0.1525e-4, 0.1490e-4 }},
            {{ 0.0,       0.0,       0.0,       0.0,       0.0,       0.0       }},
            {{ 0.0,       0.0,       0.0,       0.0,       0.0,       0.0       }},
            {{ 0.0196e-4, 0.0135e-4, 0.0108e-4, 0.0105e-4, 0.0106e-4, 0.0111e-4 }},
            {{
                {{ +1.0000, -0.2465, +0.2417, +0.1055, +0.1695, +0.1839 }},
                {{ -0.2465, +1.0000, -0.2976, +0.0923, -0.0601, -0.0227 }},
                {{ +0.2417, -0.2976, +1.0000, -0.0192, +0.3264, +0.2382 }},
                {{ +0.1055, +0.0923, -0.0192, +1.0000, -0.0164, +0.2808 }},
                {{ +0.1695, -0.0601, +0.3264, -0.0164, +1.0000, +0.0202 }},
                {{ +0.1839, -0.0227, +0.2382, +0.2808, +0.0202, +1.0000 }}
            }}
        };

        /*
         * Belle Collaboration
         *
         * Data taken from [Belle:2010A], Table V, p. 14.
         * We combine correlations between the statistic uncertainties and systematic uncertainties.
         */
        static const MultivariateGaussianConstraintTemplate<6> Bzero_to_pi_l_nu_BR_Belle_2010A
        {
            {{ "B->pilnu::BR", "B->pilnu::BR", "B->pilnu::BR", "B->pilnu::BR", "B->pilnu::BR", "B->pilnu::BR" }},
            {{
                 Kinematics{ { "s_min",  0.0 }, { "s_max",  2.0 } },
                 Kinematics{ { "s_min",  2.0 }, { "s_max",  4.0 } },
                 Kinematics{ { "s_min",  4.0 }, { "s_max",  6.0 } },
                 Kinematics{ { "s_min",  6.0 }, { "s_max",  8.0 } },
                 Kinematics{ { "s_min",  8.0 }, { "s_max", 10.0 } },
                 Kinematics{ { "s_min", 10.0 }, { "s_max", 12.0 } }
            }},
            {{
                 Options{ { "q", "d" } },
                 Options{ { "q", "d" } },
                 Options{ { "q", "d" } },
                 Options{ { "q", "d" } },
                 Options{ { "q", "d" } },
                 Options{ { "q", "d" } }
            }},
            {{ 0.1173e-4, 0.1526e-4, 0.1213e-4, 0.1465e-4, 0.1473e-4, 0.1404e-4 }},
            {{ 0.0,       0.0,       0.0,       0.0,       0.0,       0.0       }},
            {{ 0.0,       0.0,       0.0,       0.0,       0.0,       0.0       }},
            {{ 0.0219e-4, 0.0103e-4, 0.0105e-4, 0.0102e-4, 0.0108e-4, 0.0124e-4 }},
            {{
                {{ +1.0000, -0.2638, +0.1400, -0.0635, +0.1670, +0.1084 }},
                {{ -0.2638, +1.0000, -0.1390, +0.3369, +0.0292, +0.0888 }},
                {{ +0.1400, -0.1390, +1.0000, -0.0484, +0.2842, +0.2288 }},
                {{ -0.0635, +0.3369, -0.0484, +1.0000, -0.1302, +0.1859 }},
                {{ +0.1670, +0.0292, +0.2842, -0.1302, +1.0000, +0.1476 }},
                {{ +0.1084, +0.0888, +0.2288, +0.1859, +0.1476, +1.0000 }}
            }}
        };
        ///@}

        ///@name 2011 Data
        ///@{
        // Use correlation of the results from Belle and BaBar on S_K and C_K
        // http://www.slac.stanford.edu/xorg/hfag/triangle/moriond2011/index.shtml#bsgamma
        static const MultivariateGaussianConstraintTemplate<2> Bzero_to_Kstarzero_gamma_time_dependent_cp_asymmetries_HFAG_2011
        {
            {{ "B->K^*gamma::S_K^*gamma", "B->K^*gamma::C_K^*gamma" }},
            {{ Kinematics{ }, Kinematics{ } }},
            {{ Options{ { "q", "d" } }, Options{ { "q", "d" } } }},
            {{ -0.16, -0.04 }},
            {{ +0.00, +0.00 }},
            {{ +0.00, +0.00 }},
            {{ +0.22, +0.14 }},
            {{
                {{ 1,  0.06 }},
                {{ 0.06,  1 }},
            }}
        };
        /*
         * CDF Collaboration
         *
         * Data taken from [CDF:2011A]
         */
        // B^0 -> K^0 mu^+ mu^-
        // BR in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kzero_dimuon_BR_1_to_6_CDF_2011
        {
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.98e-7, +0.614e-7, -0.614e-7, +0.074e-7, -0.074e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kzero_dimuon_BR_14dot18_to_16_CDF_2011
        {
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.726e-7, +0.257e-7, -0.257e-7, +0.054e-7, -0.054e-7
        };
        // BR in [16.00, 22.86]
        static const GaussianConstraintTemplate Bzero_to_Kzero_dimuon_BR_16_to_22dot86_CDF_2011
        {
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 22.86 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.214e-7, +0.182e-7, -0.182e-7, +0.016e-7, -0.016e-7
        };

        // B^+ -> K^+ mu^+ mu^-
        // BR in [1.0, 6.0]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_1_to_6_CDF_2011
        {
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            1.41e-7, +0.20e-7, -0.20e-7, +0.09e-7, -0.09e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_14dot18_to_16_CDF_2011
        {
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.53e-7, +0.10e-7, -0.19e-7, +0.03e-7, -0.03e-7
        };
        // BR in [16.00, 22.86]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_16_to_22dot86_CDF_2011
        {
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 22.86 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.48e-7, +0.11e-7, -0.11e-7, +0.03e-7, -0.03e-7
        };

        // B^0 -> K^*0 mu^+ mu^-
        // BR in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_1_to_6_CDF_2011
        {
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.42e-7, +0.41e-7, -0.4100e-7, +0.08e-7, -0.08e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_CDF_2011
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.34e-7, +0.26e-7, -0.26e-7, +0.08e-7, -0.08e-7
        };
        // BR in [16.00, 19.21]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_CDF_2011
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.97e-7, +0.26e-7, -0.26e-7, +0.06e-7, -0.06e-7
        };

        // B^+ -> K^*+ mu^+ mu^-
        // BR in [1.0, 6.0]
        static const GaussianConstraintTemplate Bplus_to_Kstarplus_dimuon_BR_1_to_6_CDF_2011
        {
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            2.57e-7, +1.61e-7, -1.61e-7, +0.22e-7, -0.22e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintTemplate Bplus_to_Kstarplus_dimuon_BR_14dot18_to_16_CDF_2011
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.52e-7, +0.61e-7, -0.61e-7, +0.04e-7, -0.04e-7
        };
        // BR in [16.00, 19.21]
        static const GaussianConstraintTemplate Bplus_to_Kstarplus_dimuon_BR_16_to_19dot21_CDF_2011
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            1.57e-7, +0.96e-7, -0.96e-7, +0.13e-7, -0.13e-7
        };

        /*
         * CDF Collaboration
         *
         * Data taken from [CDF:2011B]
         */
        // B^0 -> K^*0 mu^+ mu^-
        // A_FB in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_CDF_2011
        {
            "B->K^*ll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.36, +0.28, -0.46, +0.11, -0.11
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_CDF_2011
        {
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.40, +0.21, -0.18, +0.07, -0.07
        };
        // A_FB in [16.00, 19.21]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_CDF_2011
        {
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.66, +0.26, -0.18, +0.19, -0.19
        };
        // F_L in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_1_to_6_CDF_2011
        {
            "B->K^*ll::F_Lavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.60, +0.21, -0.23, +0.09, -0.09
        };
        // F_L in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_CDF_2011
        {
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.32, +0.14, -0.14, +0.03, -0.03
        };
        // F_L in [16.00, 19.21]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_CDF_2011
        {
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.16, +0.22, -0.18, +0.06, -0.06
        };
        // A_T^{2} in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_T_2_1_to_6_CDF_2011
        {
            "B->K^*ll::A_T^2avg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +1.6, +1.8, -1.9, +2.2, -2.2
        };
        // A_T^{2} in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_T_2_14dot18_to_16_CDF_2011
        {
            "B->K^*ll::A_T^2avg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.4, +0.8, -0.8, +0.2, -0.2
        };
        // A_T^{2} in [16.00, 19.21]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_T_2_16_to_19dot21_CDF_2011
        {
            "B->K^*ll::A_T^2avg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.9, +0.8, -0.8, +0.4, -0.4
        };
        // A_{im} in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_im_1_to_6_CDF_2011
        {
            "B->K^*ll::J_9normavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.02, +0.40, -0.40, +0.03, -0.03
        };
        // A_{im} in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_im_14dot18_to_16_CDF_2011
        {
            "B->K^*ll::J_9normavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.15, +0.25, -0.26, +0.01, -0.01
        };
        // A_{im} in [16.00, 19.21]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_im_16_to_19dot21_CDF_2011
        {
            "B->K^*ll::J_9normavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.30, +0.36, -0.35, +0.14, -0.14
        };

        // B^+ -> K^+ mu^+ mu^-
        // A_FB in [1.0, 6.0]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_A_FB_1_to_6_CDF_2011
        {
            "B->Kll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            -0.13, +0.09, -0.09, +0.02, -0.02
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_CDF_2011
        {
            "B->Kll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            +0.05, +0.11, -0.09, +0.03, -0.03
        };
        // A_FB in [16.00, 23]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_A_FB_16_to_22dot86_CDF_2011
        {
            "B->Kll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 22.86 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            -0.09, +0.13, -0.17, +0.03, -0.03
        };

        /*
         * CDF Collaboration
         *
         * Data taken from [CDF:2011C]
         */
        // limit on BR B^0_s -> mu^+ mu^-
        static const AmorosoLimitConstraintTemplate Bzero_to_dimuon_CDF_2011
        {
            "B_q->ll::BR@Untagged",
            Kinematics{ }, // kinematics are ignored
            Options{ { "q", "s"  }, { "l", "mu" } },
            0.0, 3.5e-8, 4.0e-8,
            3.60911035e-8, 0.30785263
        };

        /*
         * LHCb Collaboration
         *
         * Data taken from LHCb:2011B
         */
        // BR in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_1_to_6_LHCb_2011
        {
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            2.1000e-7, +0.3000e-7, -0.3000e-7, +0.1500e-7, -0.1500e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_LHCb_2011
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.1466e-7, +0.2002e-7, -0.2002e-7, +0.0910e-7, -0.0910e-7
        };
        // BR in [16.00, 19.00] (in the preliminary results, LHCb only integrated up to exactly 19.00 GeV^2!)
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_LHCb_2011
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.5000e-7, +0.2400e-7, -0.2400e-7, +0.1500e-7, -0.1500e-7
        };
        // A_FB in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_LHCb_2011
        {
            "B->K^*ll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.06, +0.14, -0.13, +0.04, -0.04
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_LHCb_2011
        {
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.47, +0.08, -0.06, +0.03, -0.03
        };
        // A_FB in [16.00, 19.00] (in the preliminary results, LHCb only integrated up to exactly 19.00 GeV^2!)
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_LHCb_2011
        {
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.16, +0.13, -0.11, +0.06, -0.06
        };
        // F_L in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_1_to_6_LHCb_2011
        {
            "B->K^*ll::F_Lavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.55, +0.10, -0.10, +0.03, -0.03
        };
        // F_L in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_LHCb_2011
        {
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.37, +0.09, -0.09, +0.05, -0.05
        };
        // F_L in [16.00, 19.00] (in the preliminary results, LHCb only integrated up to exactly 19.00 GeV^2!)
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_LHCb_2011
        {
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.26, +0.10, -0.08, +0.03, -0.03
        };


        /*
         * LHCb and CMS Collaboration
         *
         * Data taken from [LHCb:2011A]
         */
        // limit on BR B^0_s -> mu^+ mu^-
        static const AmorosoLimitConstraintTemplate Bzero_to_dimuon_LHCb_CMS_2011
        {
            "B_q->ll::BR@Untagged",
            Kinematics{ },
            Options{ { "q", "s"  }, { "l", "mu" } },
            0.0, 0.9e-8, 1.08e-8,
            0.74377978e-8, 0.53538044
        };
        // use the data from the Bayes-Heinrich method
        // the mode is not at zero, but around 3.1e-9
        static const AmorosoTripleLimitConstraintTemplate Bzero_to_dimuon_LHCb_CMS_2011_Bayes
        {
            "B_q->ll::BR@Untagged",
            Kinematics{ },
            Options{ { "q", "s"  }, { "l", "mu" } },
            0.0, 0.132749474699e-8, 0.446663009589e-8, 0.932149816388e-8,
            6.4184393253e-09, 8.1583565997e-01, 1.8230347158
        };
        ///@}

        ///@name 2012 Data
        ///@{
        /*
         * BaBar
         *
         * Data taken from the talks [BaBar:2012A]
         */

        // B^+ -> K^+ mu^+ mu^-

        // BR in [1.0, 6.0]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_1_to_6_BaBar_2012
        {
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            1.36e-7, +0.27e-7, -0.24e-7, +0.03e-7, -0.03e-7
        };
        // BR in [14.21, 16.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_14dot21_to_16_BaBar_2012
        {
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.21 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.49e-7, +0.15e-7, -0.14e-7, +0.02e-7, -0.02e-7
        };
        // BR in [16.00, 22.86]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_16_to_22dot86_BaBar_2012
        {
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 22.86 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.67e-7, +0.23e-7, -0.21e-7, +0.05e-7, -0.05e-7
        };

        // B^0 -> K^*0 mu^+ mu^-

        // BR in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_1_to_6_BaBar_2012
        {
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            2.05e-7, +0.53e-7, -0.48e-7, +0.07e-7, -0.07e-7
        };
        // BR in [14.21, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_14dot21_to_16_BaBar_2012
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.21 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.46e-7, +0.41e-7, -0.36e-7, +0.06e-7, -0.06e-7
        };
        // BR in [16.00, 22.86]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_BaBar_2012
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.02e-7, +0.47e-7, -0.42e-7, +0.06e-7, -0.06e-7
        };

        /*
         * BaBar
         *
         * Data taken from the talks [BaBar:2012B]
         */

        // B^0 -> K^*0 mu^+ mu^-

        // A_FB in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_BaBar_2012
        {
            "B->K^*ll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.17, +0.14, -0.12, +0.07, -0.07
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_BaBar_2012
        {
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.34, +0.15, -0.08, +0.07, -0.07
        };
        // A_FB in [16.00, 19.21]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_BaBar_2012
        {
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.34, +0.21, -0.19, +0.07, -0.07
        };
        // F_L in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_1_to_6_BaBar_2012
        {
            "B->K^*ll::F_Lavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.25, +0.09, -0.08, +0.03, -0.03
        };
        // F_L in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_BaBar_2012
        {
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.43, +0.10, -0.13, +0.09, -0.09
        };
        // F_L in [16.00, 19.21]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_BaBar_2012
        {
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.55, +0.15, -0.17, +0.03, -0.03
        };

        /*
         * BaBar
         *
         * Data taken from [BaBar:2012C]
         */
        static const GaussianConstraintTemplate B_to_Xs_gamma_BR_1dot8_BaBar_2012C
        {
            "B->X_sgamma::BR(E_min)@NLO",
            Kinematics{ { "E_min", 1.8 } },
            Options{ { "q", "d" } },
            +3.21e-4, +0.15e-4, -0.15e-4, +0.29e-4, -0.29e-4
        };
        static const GaussianConstraintTemplate B_to_Xs_gamma_E_1_1dot8_BaBar_2012C
        {
            "B->X_sgamma::E_1(E_min)@NLO",
            Kinematics{ { "E_min", 1.8 } },
            Options{ { "q", "d" } },
            +2.267, +0.019, -0.019, +0.032, -0.032
        };
        static const GaussianConstraintTemplate B_to_Xs_gamma_E_2_1dot8_BaBar_2012C
        {
            "B->X_sgamma::E_2(E_min)@NLO",
            Kinematics{ { "E_min", 1.8 } },
            Options{ { "q", "d"  } },
            +4.84e-2, +5.3e-3, -5.3e-3, +7.7e-3, -7.7e-3
        };
        static const MultivariateGaussianConstraintTemplate<2> B_to_Xs_gamma_E_1_and_E_2_1dot8_BaBar_2012C
        {
            {{ "B->X_sgamma::E_1(E_min)@NLO", "B->X_sgamma::E_2(E_min)@NLO" }},
            {{ Kinematics{ { "E_min", 1.8 } }, Kinematics{ { "E_min", 1.8 } } }},
            {{ Options{ { "q", "d"  } }, Options{ { "q", "d"  } } }},
            {{ 2.267, 0.0484 }},
            {{ 0.019, 0.0053 }},
            {{ 0.019, 0.0053 }},
            {{ 0.032, 0.0077 }},
            {{
                 {{ +1.00, -0.88 }},
                 {{ -0.88, +1.00 }}
            }}
        };

        /*
         * BaBar Collaboration
         *
         * Data taken from [BaBar:2012D], Table VIII, p. 14. We use the results
         * 'without FSR' (final state radiation).
         * We combine correlations between the statistic uncertainties and systematic uncertainties.
         */
        static const MultivariateGaussianConstraintTemplate<6> Bzero_to_pi_l_nu_BR_BaBar_2012D
        {
            {{ "B->pilnu::BR", "B->pilnu::BR", "B->pilnu::BR", "B->pilnu::BR", "B->pilnu::BR", "B->pilnu::BR" }},
            {{
                 Kinematics{ { "s_min",  0.0 }, { "s_max",  2.0 } },
                 Kinematics{ { "s_min",  2.0 }, { "s_max",  4.0 } },
                 Kinematics{ { "s_min",  4.0 }, { "s_max",  6.0 } },
                 Kinematics{ { "s_min",  6.0 }, { "s_max",  8.0 } },
                 Kinematics{ { "s_min",  8.0 }, { "s_max", 10.0 } },
                 Kinematics{ { "s_min", 10.0 }, { "s_max", 12.0 } }
            }},
            {{
                 Options{ { "q", "d" } },
                 Options{ { "q", "d" } },
                 Options{ { "q", "d" } },
                 Options{ { "q", "d" } },
                 Options{ { "q", "d" } },
                 Options{ { "q", "d" } }
            }},
            {{ 0.1225e-4, 0.1277e-4, 0.1274e-4, 0.1498e-4, 0.1405e-4, 0.1617e-4 }},
            {{ 0.0,       0.0,       0.0,       0.0,       0.0,       0.0       }},
            {{ 0.0,       0.0,       0.0,       0.0,       0.0,       0.0       }},
            {{ 0.0182e-4, 0.0128e-4, 0.0109e-4, 0.0103e-4, 0.0115e-4, 0.0104e-4 }},
            {{
                {{ +1.0000, -0.0759, +0.1679, +0.1432, +0.1831, +0.1471 }},
                {{ -0.0759, +1.0000, -0.1473, +0.2174, +0.0718, +0.0975 }},
                {{ +0.1679, -0.1473, +1.0000, -0.0889, +0.2250, +0.1076 }},
                {{ +0.1432, +0.2174, -0.0889, +1.0000, +0.0218, +0.2772 }},
                {{ +0.1831, +0.0718, +0.2250, +0.0218, +1.0000, -0.0425 }},
                {{ +0.1471, +0.0975, +0.1076, +0.2772, -0.0425, +1.0000 }},
            }}
        };

        /*
         * CDF Collaboration
         *
         * Data taken from [CDF:2012A]
         */
        // B^0 -> K^0 mu^+ mu^-
        // BR in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kzero_dimuon_BR_1_to_6_CDF_2012
        {
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.11e-7, +0.52e-7, -0.52e-7, +0.14e-7, -0.14e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kzero_dimuon_BR_14dot18_to_16_CDF_2012
        {
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.43e-7, +0.18e-7, -0.18e-7, +0.04e-7, -0.04e-7
        };
        // BR in [16.00, 22.86]
        static const GaussianConstraintTemplate Bzero_to_Kzero_dimuon_BR_16_to_22dot86_CDF_2012
        {
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 22.86 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.26e-7, +0.15e-7, -0.15e-7, +0.03e-7, -0.03e-7
        };

        // B^+ -> K^+ mu^+ mu^-
        // BR in [1.0, 6.0]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_1_to_6_CDF_2012
        {
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            1.35e-7, +0.18e-7, -0.18e-7, +0.08e-7, -0.08e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_14dot18_to_16_CDF_2012
        {
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.42e-7, +0.08e-7, -0.08e-7, +0.02e-7, -0.02e-7
        };
        // BR in [16.00, 22.86]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_16_to_22dot86_CDF_2012
        {
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 22.86 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.50e-7, +0.10e-7, -0.10e-7, +0.02e-7, -0.02e-7
        };

        // B^0 -> K^*0 mu^+ mu^-
        // BR in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_1_to_6_CDF_2012
        {
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.97e-7, +0.49e-7, -0.49e-7, +0.14e-7, -0.14e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_CDF_2012
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.35e-7, +0.21e-7, -0.21e-7, +0.08e-7, -0.08e-7
        };
        // BR in [16.00, 19.21]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_CDF_2012
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.98e-7, +0.22e-7, -0.22e-7, +0.06e-7, -0.06e-7
        };

        // B^+ -> K^*+ mu^+ mu^-
        // BR in [1.0, 6.0]
        static const GaussianConstraintTemplate Bplus_to_Kstarplus_dimuon_BR_1_to_6_CDF_2012
        {
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            3.56e-7, +1.38e-7, -1.38e-7, +0.43e-7, -0.43e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintTemplate Bplus_to_Kstarplus_dimuon_BR_14dot18_to_16_CDF_2012
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            1.02e-7, +0.58e-7, -0.58e-7, +0.13e-7, -0.13e-7
        };
        // BR in [16.00, 19.21]
        static const GaussianConstraintTemplate Bplus_to_Kstarplus_dimuon_BR_16_to_19dot21_CDF_2012
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            1.68e-7, +0.74e-7, -0.74e-7, +0.19e-7, -0.19e-7
        };

        // B^0 -> K^*0 mu^+ mu^-
        // A_FB in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_CDF_2012
        {
            "B->K^*ll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.29, +0.21, -0.25, +0.06, -0.06
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_CDF_2012
        {
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.49, +0.09, -0.10, +0.07, -0.07
        };
        // A_FB in [16.00, 19.21]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_CDF_2012
        {
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.42, +0.23, -0.22, +0.09, -0.09
        };
        // F_L in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_1_to_6_CDF_2012
        {
            "B->K^*ll::F_Lavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.78, +0.13, -0.15, +0.08, -0.08
        };
        // F_L in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_CDF_2012
        {
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.45, +0.12, -0.12, +0.04, -0.04
        };
        // F_L in [16.00, 19.21]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_CDF_2012
        {
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.09, +0.14, -0.12, +0.08, -0.08
        };
        // A_T^{2} in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_T_2_1_to_6_CDF_2012
        {
            "B->K^*ll::A_T^2avg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.45, +2.24, -2.22, +0.76, -0.76
        };
        // A_T^{2} in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_T_2_14dot18_to_16_CDF_2012
        {
            "B->K^*ll::A_T^2avg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.15, +0.72, -0.72, +0.14, -0.14
        };
        // A_T^{2} in [16.00, 19.21]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_T_2_16_to_19dot21_CDF_2012
        {
            "B->K^*ll::A_T^2avg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.62, +0.56, -0.53, +0.13, -0.13
        };
        // A_{im} in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_im_1_to_6_CDF_2012
        {
            "B->K^*ll::J_9normavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.51, +0.28, -0.29, +0.15, -0.15
        };
        // A_{im} in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_im_14dot18_to_16_CDF_2012
        {
            "B->K^*ll::J_9normavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.16, +0.21, -0.22, +0.03, -0.03
        };
        // A_{im} in [16.00, 19.21]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_im_16_to_19dot21_CDF_2012
        {
            "B->K^*ll::J_9normavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.02, +0.26, -0.27, +0.04, -0.04
        };

        // B^+ -> K^+ mu^+ mu^-
        // A_FB in [1.0, 6.0]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_A_FB_1_to_6_CDF_2012
        {
            "B->Kll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            -0.13, +0.10, -0.11, +0.02, -0.02
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_CDF_2012
        {
            "B->Kll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            +0.07, +0.08, -0.08, +0.01, -0.01
        };
        // A_FB in [16.00, 23]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_A_FB_16_to_22dot86_CDF_2012
        {
            "B->Kll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 22.86 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            -0.05, +0.10, -0.18, +0.05, -0.05
        };

        /*
         * LHCb Collaboration
         *
         * Data taken from LHCb:2012A
         */
        // BR in [1.0, 6.0], multiply with bin width 5 GeV^2
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_1_to_6_LHCb_2012
        {
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.42e-7 * 5, +0.04e-7 * 5, -0.04e-7 * 5, +0.04e-7 * 5, -0.04e-7 * 5
        };
        // BR in [14.18, 16.00], multiply with bin width 1.82 GeV^2
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_LHCb_2012
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.59e-7 * 1.82, +0.07e-7 * 1.82, -0.07e-7 * 1.82, +0.04e-7 * 1.82, -0.04e-7 * 1.82
        };
        // BR in [16.00, 19.00], multiply with bin width 3 GeV^2 (in the preliminary results, LHCb only integrated up to exactly 19.00 GeV^2!)
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_16_to_19_LHCb_2012
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.44e-7 * 3, +0.05e-7 * 3, -0.05e-7 * 3, +0.03e-7 * 3, -0.03e-7 * 3
        };
        // A_FB in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_LHCb_2012
        {
            "B->K^*ll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.18, +0.06, -0.06, +0.02, -0.01
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_LHCb_2012
        {
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.49, +0.06, -0.04, +0.05, -0.02
        };
        // A_FB in [16.00, 19.00] (in the preliminary results, LHCb only integrated up to exactly 19.00 GeV^2!)
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_16_to_19_LHCb_2012
        {
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.30, +0.07, -0.07, +0.01, -0.04
        };
        // F_L in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_1_to_6_LHCb_2012
        {
            "B->K^*ll::F_Lavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.66, +0.06, -0.06, +0.04, -0.03
        };
        // F_L in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_LHCb_2012
        {
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.35, +0.07, -0.06, +0.07, -0.02
        };
        // F_L in [16.00, 19.00] (in the preliminary results, LHCb only integrated up to exactly 19.00 GeV^2!)
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_16_to_19_LHCb_2012
        {
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.37, +0.06, -0.07, +0.03, -0.04
        };
        // A_{im} in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_im_1_to_6_LHCb_2012
        {
            "B->K^*ll::J_9normavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.07, +0.07, -0.07, +0.02, -0.01
        };
        // A_{im} in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_im_14dot18_to_16_LHCb_2012
        {
            "B->K^*ll::J_9normavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.01, +0.08, -0.07, +0.04, -0.02
        };
        // A_{im} in [16.00, 19.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_im_16_to_19_LHCb_2012
        {
            "B->K^*ll::J_9normavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.06, +0.09, -0.10, +0.03, -0.05
        };
        // S_3 in [1.0, 6.0], LHCb gives 2 * S_3
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_S_3_1_to_6_LHCb_2012
        {
            "B->K^*ll::J_3normavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.10 / 2, +0.15 / 2, -0.16 / 2, +0.02 / 2, -0.01 / 2
        };
        // S_3 in [14.18, 16.00], LHCb gives 2 * S_3
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_S_3_14dot18_to_16_LHCb_2012
        {
            "B->K^*ll::J_3normavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.04 / 2, +0.15 / 2, -0.19 / 2, +0.04 / 2, -0.02 / 2
        };
        // S_3 in [16.00, 19.00], LHCb gives 2 * S_3
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_S_3_16_to_19_LHCb_2012
        {
            "B->K^*ll::J_3normavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.47 / 2, +0.21 / 2, -0.10 / 2, +0.03 / 2, -0.05 / 2
        };

        /*
         * LHCb Collaboration
         *
         * Data taken from LHCb:2012C
         */
        // BR in [1.0, 6.0]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_1_to_6_LHCb_2012
        {
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            2.41e-8 * 5, +0.17e-8 * 5, -0.17e-8 * 5, +0.14e-8 * 5, -0.14e-8 * 5
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_14dot18_to_16_LHCb_2012
        {
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            2.07e-8 * 1.82, +0.20e-8 * 1.82, -0.20e-8 * 1.82, +0.08e-8 * 1.82, -0.08e-8 * 1.82
        };
        // BR in [16.00, 18.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_16_to_18_LHCb_2012
        {
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 18.00} },
            Options{ { "q", "u"  }, { "l", "mu" } },
            1.77e-8 * 2, +0.18e-8 * 2, -0.18e-8 * 2, +0.09e-8 * 2, -0.09e-8 * 2
        };
        // BR in [18.00, 22.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_18_to_22_LHCb_2012
        {
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 18.00 }, { "s_max", 22.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.78e-8 * 4, +0.10e-8 * 4, -0.10e-8 * 4, +0.04e-8 * 4, -0.04e-8 * 4
        };

        /* sign flipped wrt to table 1 in LHCb:2012C! */

        // A_FB in [1.0, 6.0]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_A_FB_1_to_6_LHCb_2012
        {
            "B->Kll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            -0.02, +0.03, -0.05, +0.01, -0.02
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_LHCb_2012
        {
            "B->Kll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            +0.01, +0.06, -0.12, +0.01, -0.01
        };
        // A_FB in [16.00, 18.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_A_FB_16_to_18_LHCb_2012
        {
            "B->Kll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 18.00} },
            Options{ { "q", "u"  }, { "l", "mu" } },
            +0.09, +0.09, -0.07, +0.01, -0.02
        };
        // A_FB in [18.00, 22.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_A_FB_18_to_22_LHCb_2012
        {
            "B->Kll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 18.00 }, { "s_max", 22.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            -0.02, +0.11, -0.11, +0.01, -0.01
        };

        // F_H in [1.0, 6.0]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_F_H_1_to_6_LHCb_2012
        {
            "B->Kll::F_Havg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.05, +0.08, -0.05, +0.04, -0.02
        };
        // F_H in [14.18, 16.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_F_H_14dot18_to_16_LHCb_2012
        {
            "B->Kll::F_Havg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.08, +0.28, -0.08, +0.02, -0.01
        };
        // F_H in [16.00, 18.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_F_H_16_to_18_LHCb_2012
        {
            "B->Kll::F_Havg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 18.00} },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.18, +0.22, -0.14, +0.01, -0.04
        };
        // F_H in [18.00, 22.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_F_H_18_to_22_LHCb_2012
        {
            "B->Kll::F_Havg@LowRecoil",
            Kinematics{ { "s_min", 18.00 }, { "s_max", 22.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.14, +0.31, -0.14, +0.01, -0.02
        };

        // Data taken from LHCb:2012B
        /* fit Amoroso to posterior supplied by Diego Martinez by fixing
         * a) cdf(0) = 0
         * b) cdf(x_10) = 0.1
         * c) cdf(x_50) = 0.5
         * d) cdf(x_90) = 0.9
         */
        static const AmorosoTripleLimitConstraintTemplate Bzero_to_dimuon_LHCb_2012
        {
            "B_q->ll::BR@Untagged",
            Kinematics(),
            Options{ { "q", "s"  }, { "l", "mu" } },
            0.0, 0.558367940293e-9, 2.03115589965e-9, 4.4528950788e-9,
            2.9708273062e-09, 8.2392613044e-01, 1.6993290032
        };

        // Data taken from LHCb:2012D
        /* fit Amoroso to result assuming
         * a) cdf(0) = 0
         * b) mode = 3.2
         * c) 90% in (1.3, 5.8)
         * d) pdf(1.3) = pdf(5.8) (smallest interval)
         */
        static const AmorosoModeConstraintTemplate Bzero_to_dimuon_LHCb_Nov_2012
        {
            "B_q->ll::BR@Untagged",
            Kinematics(),
            Options{ { "q", "s"  }, { "l", "mu" } },
            0.0, 3.2e-9, 5.479025195034372e-9, 6.110034104385014e-9,
            2.3625776605e-09, 2.2682156277, 1.7296007586
        };

        // Data taken from [LHCb:2012E]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_CP_1_to_6_LHCb_2012E
        {
            "B->K^*ll::A_CP@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.058, +0.064, -0.064, +0.009, -0.009
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_CP_14dot18_to_16_LHCb_2012E
        {
            "B->K^*ll::A_CP@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.201, +0.104, -0.104, +0.009, -0.009
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_CP_16_to_20_LHCb_2012E
        {
            "B->K^*ll::A_CP@LowRecoil",
            // [LHCb:2012E] gives unphysical upper kinematic range of 20.00GeV^2.
            // Reducing this to 19.81GeV^2.
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.81 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.089, +0.100, -0.100, +0.012, -0.012
        };

        /*
         * PDG
         *
         * Data taken from [PDG2012]
         */
        static const GaussianConstraintTemplate B_Bstar_mass_splitting_PDG_2012
        {
            "B::M_B^*-M_B",
            Kinematics{ },
            Options{ { "q", "d" } },
            0.04578, +0.00035, -0.00035, 0.0, 0.0
        };
        ///@}

        ///@name 2013 Data
        ///@{
        /*
         * ATLAS Collaboration
         *
         * Data taken from [ATLAS:2013A]
         */
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_ATLAS_2013A
        {
            "B->K^*ll::A_FB@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.07, +0.20, -0.20, +0.07, -0.07
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_ATLAS_2013A
        {
            "B->K^*ll::A_FB@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.48, +0.19, -0.19, +0.05, -0.05
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_16_to_19_ATLAS_2013A
        {
            "B->K^*ll::A_FB@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.16, +0.10, -0.10, +0.03, -0.03
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_1_to_6_ATLAS_2013A
        {
            "B->K^*ll::F_L@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.18, +0.15, -0.15, +0.03, -0.03
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_ATLAS_2013A
        {
            "B->K^*ll::F_L@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.28, +0.16, -0.16, +0.03, -0.03
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_16_to_19_ATLAS_2013A
        {
            "B->K^*ll::F_L@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.35, +0.08, -0.08, +0.02, -0.02
        };

        /*
         * CMS Collaboration
         *
         * Data taken from [CMS:2013A]
         */
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_1_to_6_CMS_2013A
        {
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            4.4e-8 * 5.0, +0.6e-8 * 5.0, -0.6e-8 * 5.0, +0.7e-8 * 5.0, -0.7e-8 * 5.0
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_CMS_2013A
        {
            "B->K^*ll::A_FB@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.07, +0.12, -0.12, +0.01, -0.01
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_1_to_6_CMS_2013A
        {
            "B->K^*ll::F_L@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.68, +0.10, -0.10, +0.02, -0.02
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_CMS_2013A
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            4.6e-8 * 1.82, +0.9e-8 * 1.82, -0.8e-8 * 1.82, +0.8e-8 * 1.82, -0.8e-8 * 1.82
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_CMS_2013A
        {
            "B->K^*ll::A_FB@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.29, +0.09, -0.09, +0.05, -0.05
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_CMS_2013A
        {
            "B->K^*ll::F_L@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.53, +0.12, -0.12, +0.03, -0.03
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_16_to_19_CMS_2013A
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            5.2e-8 * 3.0, +0.6e-8 * 3.0, -0.6e-8 * 3.0, +0.9e-8 * 3.0, -0.9e-8 * 3.0
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_16_to_19_CMS_2013A
        {
            "B->K^*ll::A_FB@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.41, +0.05, -0.05, +0.03, -0.03
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_16_to_19_CMS_2013A
        {
            "B->K^*ll::F_L@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.44, +0.07, -0.07, +0.03, -0.03
        };

        /*
         * Belle Collaboration
         *
         * Data taken from [Belle:2013A], Table XVI and XVII, pp. 37-38.
         */
        static const MultivariateGaussianConstraintTemplate<6> Bzero_to_pi_l_nu_BR_Belle_2013A
        {
            {{ "B->pilnu::BR", "B->pilnu::BR", "B->pilnu::BR", "B->pilnu::BR", "B->pilnu::BR", "B->pilnu::BR" }},
            {{
                 Kinematics{ { "s_min",  0.0 }, { "s_max",  2.0 } },
                 Kinematics{ { "s_min",  2.0 }, { "s_max",  4.0 } },
                 Kinematics{ { "s_min",  4.0 }, { "s_max",  6.0 } },
                 Kinematics{ { "s_min",  6.0 }, { "s_max",  8.0 } },
                 Kinematics{ { "s_min",  8.0 }, { "s_max", 10.0 } },
                 Kinematics{ { "s_min", 10.0 }, { "s_max", 12.0 } }
            }},
            {{
                 Options{ { "q", "d" } },
                 Options{ { "q", "d" } },
                 Options{ { "q", "d" } },
                 Options{ { "q", "d" } },
                 Options{ { "q", "d" } },
                 Options{ { "q", "d" } }
            }},
            {{ 0.195e-4, 0.106e-4, 0.151e-4, 0.097e-4, 0.078e-4, 0.159e-4 }},
            {{ 0.0,      0.0,      0.0,      0.0,      0.0,      0.0      }},
            {{ 0.0,      0.0,      0.0,      0.0,      0.0,      0.0      }},
            {{ 0.032e-4, 0.027e-4, 0.028e-4, 0.023e-4, 0.022e-4, 0.028e-4 }},
            {{
                {{ +1.000, -0.145, +0.010, -0.001, +0.000, +0.000 }},
                {{ -0.145, +1.000, -0.100, +0.008, +0.000, +0.001 }},
                {{ +0.010, -0.100, +1.000, -0.094, +0.003, -0.001 }},
                {{ -0.001, +0.008, -0.094, +1.000, -0.078, +0.005 }},
                {{ +0.000, +0.000, +0.003, -0.078, +1.000, -0.097 }},
                {{ +0.000, +0.001, -0.001, +0.005, -0.097, +1.000 }}
            }}
        };

        /*
         * CMS Collaboration
         *
         * Data taken from [CMS:2013B]
         */
        /* fit Amoroso to result assuming
         * a) cdf(0) = 0
         * b) mode = 3.0
         * c) 68% in (2.1, 4.0)
         * d) pdf(2.1) = pdf(4.0) (smallest interval)
         */
        static const AmorosoConstraintTemplate Bzero_to_dimuon_CMS_2013B
        {
            "B_q->ll::BR@Untagged",
            Kinematics(),
            Options{ { "q", "s"  }, { "l", "mu" } },
            0.0, 1.9859633460e-9, 2.7971996021, 2.0218845762
        };

        /*
         * LHCb Collaboration
         *
         * Data taken from [LHCb:2013B]
         */
        // BR
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_1_to_6_LHCb_2013B
        {
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            1.70e-7, +0.18e-7, -0.23e-7, +0.20e-7, -0.20e-7
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_LHCb_2013B
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            1.02e-7, +0.13e-7, -0.16e-7, +0.07e-7, -0.07e-7
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_16_to_19_LHCb_2013B
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.0 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            1.23e-7, +0.15e-7, -0.17e-7, +0.12e-7, -0.12e-7
        };
        // F_L
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_1_to_6_LHCb_2013B
        {
            "B->K^*ll::F_L@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.65, +0.08, -0.07, +0.03, -0.03
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_LHCb_2013B
        {
            "B->K^*ll::F_L@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.33, +0.08, -0.07, +0.02, -0.03
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_16_to_19_LHCb_2013B
        {
            "B->K^*ll::F_L@LowRecoil",
            Kinematics{ { "s_min", 16.0 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.38, +0.09, -0.07, +0.03, -0.03
        };
        // A_FB
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_LHCb_2013B
        {
            "B->K^*ll::A_FB@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.17, +0.06, -0.06, +0.01, -0.01
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_LHCb_2013B
        {
            "B->K^*ll::A_FB@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.51, +0.05, -0.07, +0.02, -0.02
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_16_to_19_LHCb_2013B
        {
            "B->K^*ll::A_FB@LowRecoil",
            Kinematics{ { "s_min", 16.0 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.30, +0.08, -0.08, +0.02, -0.01
        };
        // S_3
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_S_3_1_to_6_LHCb_2013B
        {
            "B->K^*ll::J_3normavg@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.03, +0.07, -0.07, +0.01, -0.01
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_S_3_14dot18_to_16_LHCb_2013B
        {
            "B->K^*ll::J_3normavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.03, +0.09, -0.10, +0.01, -0.01
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_S_3_16_to_19_LHCb_2013B
        {
            "B->K^*ll::J_3normavg@LowRecoil",
            Kinematics{ { "s_min", 16.0 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.22, +0.10, -0.09, +0.02, -0.01
        };
        // S_9
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_S_9_1_to_6_LHCb_2013B
        {
            "B->K^*ll::J_9normavg@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.07, +0.09, -0.08, +0.01, -0.01
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_S_9_14dot18_to_16_LHCb_2013B
        {
            "B->K^*ll::J_9normavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.00, +0.09, -0.08, +0.01, -0.01
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_S_9_16_to_19_LHCb_2013B
        {
            "B->K^*ll::J_9normavg@LowRecoil",
            Kinematics{ { "s_min", 16.0 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.06, +0.11, -0.10, +0.01, -0.01
        };
        // A_9
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_9_1_to_6_LHCb_2013B
        {
            "B->K^*ll::A_9@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.03, +0.08, -0.08, +0.01, -0.01
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_9_14dot18_to_16_LHCb_2013B
        {
            "B->K^*ll::A_9@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.06, +0.11, -0.08, +0.01, -0.01
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_9_16_to_19_LHCb_2013B
        {
            "B->K^*ll::A_9@LowRecoil",
            Kinematics{ { "s_min", 16.0 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.00, +0.11, -0.10, +0.01, -0.01
        };
        // A_T^2
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_T2_1_to_6_LHCb_2013B
        {
            "B->K^*ll::A_T^2@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.15, +0.39, -0.41, +0.03, -0.03
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_T2_14dot18_to_16_LHCb_2013B
        {
            "B->K^*ll::A_T^2@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.07, +0.26, -0.28, +0.02, -0.02
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_T2_16_to_19_LHCb_2013B
        {
            "B->K^*ll::A_T^2@LowRecoil",
            Kinematics{ { "s_min", 16.0 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.71, +0.35, -0.26, +0.06, -0.04
        };
        // A_T^re
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_Tre_1_to_6_LHCb_2013B
        {
            "B->K^*ll::A_T^re@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.66, +0.22, -0.24, +0.01, -0.04
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_Tre_14dot18_to_16_LHCb_2013B
        {
            "B->K^*ll::A_T^re@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -1.00, +0.05, -0.01, +0.02, -0.01
            // (unphysical) lower errors (both stat and syst) adjusted to -0.01 to work
            // around limitations in the asymmetric gaussian likelihood block.
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_Tre_16_to_19_LHCb_2013B
        {
            "B->K^*ll::A_T^re@LowRecoil",
            Kinematics{ { "s_min", 16.0 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.64, +0.15, -0.15, +0.02, -0.01
        };
        ///@}

        /*
         * LHCb Collaboration
         *
         * Data taken from [LHCb:2013C]
         */
        // P'_4, LHCb uses a different sign and a factor 1/2 relative to the theory prediction.
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_Pprime_4_1_to_6_LHCb_2013C
        {
            "B->K^*ll::P'_4@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.58, +0.32, -0.36, +0.06, -0.06
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_Pprime_4_14dot18_to_16_LHCb_2013C
        {
            "B->K^*ll::P'_4@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.18, +0.54, -0.70, +0.08, -0.08
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_Pprime_4_16_to_19_LHCb_2013C
        {
            "B->K^*ll::P'_4@LowRecoil",
            Kinematics{ { "s_min", 16.0 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.70, +0.44, -0.52, +0.06, -0.06
        };
        // P'_5
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_Pprime_5_1_to_6_LHCb_2013C
        {
            "B->K^*ll::P'_5@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.21, +0.20, -0.21, +0.03, -0.03
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_Pprime_5_14dot18_to_16_LHCb_2013C
        {
            "B->K^*ll::P'_5@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.79, +0.20, -0.13, +0.18, -0.18
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_Pprime_5_16_to_19_LHCb_2013C
        {
            "B->K^*ll::P'_5@LowRecoil",
            Kinematics{ { "s_min", 16.0 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.60, +0.19, -0.16, +0.09, -0.09
        };
        // P'_6, LHCb uses a different sign
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_Pprime_6_1_to_6_LHCb_2013C
        {
            "B->K^*ll::P'_6@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.18, +0.21, -0.21, +0.03, -0.03
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_Pprime_6_14dot18_to_16_LHCb_2013C
        {
            "B->K^*ll::P'_6@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.18, +0.24, -0.25, +0.03, -0.03
        };
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_Pprime_6_16_to_19_LHCb_2013C
        {
            "B->K^*ll::P'_6@LowRecoil",
            Kinematics{ { "s_min", 16.0 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.31, +0.37, -0.38, +0.03, -0.03
        };
        ///@}

        /*
         * LHCb Collaboration
         *
         * Data taken from [LHCb:2013D]
         */
        ///@name 2013
        ///@{

        /* fit Amoroso to result assuming
         * a) cdf(0) = 0
         * b) mode = 2.9
         * c) 68% in (1.90, 4.04)
         * d) pdf(1.90) = pdf(4.04) (smallest interval)
         */
        static const AmorosoConstraintTemplate Bzero_to_dimuon_LHCb_2013D
        {
            "B_q->ll::BR@Untagged",
            Kinematics(),
            Options{ { "q", "s"  }, { "l", "mu" } },
            0.0, 2.23067e-9, 2.19078, 2.00322
        };
        ///@}

        ///@name 2014
        ///@{
        /*
         * LHCb Collaboration
         *
         * Data taken from [LHCb:2014A]
         */
        // BR in [1.1, 2.0]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_1dot1_to_2_LHCb_2014
        {
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.1 }, { "s_max", 2.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            2.33e-8 * 0.9, +0.15e-8 * 0.9, -0.15e-8 * 0.9, +0.12e-8 * 0.9, -0.12e-8 * 0.9
        };
        // BR in [2.00, 3.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_2_to_3_LHCb_2014
        {
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 2.00 }, { "s_max", 3.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            2.82e-8 * 1.0, +0.16e-8 * 1.0, -0.16e-8 * 1.0, +0.14e-8 * 1.0, -0.14e-8 * 1.0
        };
        // BR in [3.00, 4.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_3_to_4_LHCb_2014
        {
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 3.00 }, { "s_max", 4.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            2.54e-8 * 1.0, +0.15e-8 * 1.0, -0.15e-8 * 1.0, +0.13e-8 * 1.0, -0.13e-8 * 1.0
        };
        // BR in [4.00, 5.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_4_to_5_LHCb_2014
        {
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 4.00 }, { "s_max", 5.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            2.21e-8 * 1.0, +0.14e-8 * 1.0, -0.14e-8 * 1.0, +0.11e-8 * 1.0, -0.11e-8 * 1.0
        };
        // BR in [5.00, 6.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_5_to_6_LHCb_2014
        {
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 5.00 }, { "s_max", 6.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            2.31e-8 * 1.0, +0.14e-8 * 1.0, -0.14e-8 * 1.0, +0.12e-8 * 1.0, -0.12e-8 * 1.0
        };
        // BR in [1.10, 6.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_1dot1_to_6_LHCb_2014
        {
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.10 }, { "s_max", 6.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            2.42e-8 * 4.9, +0.07e-8 * 4.9, -0.07e-8 * 4.9, +0.12e-8 * 4.9, -0.12e-8 * 4.9
        };

        // BR in [15.00, 22.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_BR_15_to_22_LHCb_2014
        {
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 15.00 }, { "s_max", 22.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            1.21e-8 * 7, +0.04e-8 * 7, -0.04e-8 * 7, +0.06e-8 * 7, -0.06e-8 * 7
        };


        /* Data taken from [LHCb:2014B] */
        // todo Fig. 4 shows best-fit point close to physical boundary
        // Gaussianization is a poor approximation
        // sign flipped wrt to table 1 in LHCb:2014A!
        // A_FB in [1.1, 6.0]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_A_FB_1dot1_to_6_LHCb_2014
        {
            "B->Kll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.1 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            -0.005, +0.015, -0.015, +0.01, -0.01
        };
        // A_FB in [15, 22]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_A_FB_15_to_22_LHCb_2014
        {
            "B->Kll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 15.00 }, { "s_max", 22.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.015, +0.015, -0.015, +0.01, -0.01
        };
        // F_H in [1.1, 6.0]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_F_H_1dot1_to_6_LHCb_2014
        {
            "B->Kll::F_Havg@LargeRecoil",
            Kinematics{ { "s_min", 1.1 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.03, +0.03, -0.03, +0.02, -0.02
        };
        // F_H in [15, 22]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_F_H_15_to_22_LHCb_2014
        {
            "B->Kll::F_Havg@LowRecoil",
            Kinematics{ { "s_min", 15.00 }, { "s_max", 22.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.035, +0.035, -0.035, +0.02, -0.02
        };
        ///@}

        /*
         * Theoretical Constraints from e.g. Lattice QCD.
         */

        ///@name 2013
        ///@{
        /*
         * Reproduced from [HPQCD:2013A].
         */
        static const MultivariateGaussianConstraintTemplate<3> B_to_K_fplus_17_to_23_HPQCD_2013A
        {
            {{ "B->K::f_+(s)", "B->K::f_+(s)", "B->K::f_+(s)" }},
            {{ Kinematics{ { "s", 17.0 } }, Kinematics{ { "s", 20.0 } }, Kinematics{ { "s", 23.0 } },}},
            {{ Options{ }, Options{ }, Options{ } }},
            {{ 1.07617,   1.50728,   2.34247   }},
            {{ 0.0265336, 0.0299249, 0.0696341 }},
            {{ 0.0265336, 0.0299249, 0.0696341 }},
            {{ +0.00, +0.00, +0.00 }}, // we assign no systematic uncertainty
            {{
                {{ 1.000000, 0.778675, 0.290551 }},
                {{ 0.778675, 1.000000, 0.708433 }},
                {{ 0.290551, 0.708433, 1.000000 }}
            }},
            0u
        };

        /*
         * Reproduced from [HPQCD:2013B].
         */
        static const MultivariateGaussianConstraintTemplate<2> B_to_Kstar_V_15_to_19dot21_HPQCD_2013B
        {
            {{ "B->K^*::V(s)", "B->K^*::V(s)" }},
            {{ Kinematics{ { "s", 15.0 } }, Kinematics{ { "s", 19.21 } } }},
            {{ Options{ }, Options{ } }},
            {{ 1.19122,   1.97758  }},
            {{ 0.0999166, 0.112599 }},
            {{ 0.0999166, 0.112599 }},
            {{ +0.00, +0.00 }}, // we assign no systematic uncertainty
            {{
                {{ 1.000000, 0.461948 }},
                {{ 0.461948, 1.000000 }},
            }},
            0u
        };
        static const MultivariateGaussianConstraintTemplate<2> B_to_Kstar_A1_15_to_19dot21_HPQCD_2013B
        {
            {{ "B->K^*::A_1(s)", "B->K^*::A_1(s)" }},
            {{ Kinematics{ { "s", 15.0 } }, Kinematics{ { "s", 19.21 } } }},
            {{ Options{ }, Options{ } }},
            {{ 0.515837,  0.644548  }},
            {{ 0.0293675, 0.0195013 }},
            {{ 0.0293675, 0.0195013 }},
            {{ +0.00, +0.00 }}, // we assign no systematic uncertainty
            {{
                {{ 1.000000, 0.52072  }},
                {{ 0.52072,  1.000000 }},
            }},
            0u
        };
        static const MultivariateGaussianConstraintTemplate<2> B_to_Kstar_A12_15_to_19dot21_HPQCD_2013B
        {
            {{ "B->K^*::A_12(s)", "B->K^*::A_12(s)" }},
            {{ Kinematics{ { "s", 15.0 } }, Kinematics{ { "s", 19.21 } } }},
            {{ Options{ }, Options{ } }},
            {{ 0.371041,  0.440076  }},
            {{ 0.0306946, 0.0273783 }},
            {{ 0.0306946, 0.0273783 }},
            {{ +0.00, +0.00 }}, // we assign no systematic uncertainty
            {{
                {{ 1.000000, 0.204495 }},
                {{ 0.204495, 1.000000 }},
            }},
        0u
        };
        ///@}

        ///@name 2014
        ///@{
        /*
         * From [IKMvD:2014].
         */
        static const MultivariateGaussianConstraintTemplate<6> B_to_pi_fp_IKMvD_2014
        {
            {{ "B->pi::f_+(s)", "B->pi::f_+'(s)", "B->pi::f_+''(s)", "B->pi::f_+(s)", "B->pi::f_+'(s)", "B->pi::f_+''(s)" }},
            {{ Kinematics{ { "s", 0.0 } }, Kinematics{ { "s", 0.0 } }, Kinematics{ { "s", 0.0 } }, Kinematics{ { "s", 10.0 } },Kinematics{ { "s", 10.0 } }, Kinematics{ { "s", 10.0 } } }},
            {{ Options{ }, Options{ }, Options{ }, Options{ }, Options{ }, Options{ } }},
            {{ 3.095910e-1, 1.5545000e-2, 1.242580e-3, 5.619880e-1, 4.033900e-2, 4.708910e-3 }},
            {{ 0.199207e-1, 0.0999788e-2, 0.108337e-3, 0.321127e-1, 0.237506e-2, 0.365924e-3 }},
            {{ 0.199207e-1, 0.0999788e-2, 0.108337e-3, 0.321127e-1, 0.237506e-2, 0.365924e-3 }},
            {{ 0.000000e-1, 0.0000000e-2, 0.000000e-3, 0.000000e-1, 0.000000e-2, 0.000000e-3 }},
            {{
                {{ 1.000, 0.735, 0.374, 0.925, 0.564, 0.313 }},
                {{ 0.735, 1.000, 0.867, 0.927, 0.863, 0.246 }},
                {{ 0.374, 0.867, 1.000, 0.682, 0.853, 0.221 }},
                {{ 0.925, 0.927, 0.682, 1.000, 0.814, 0.389 }},
                {{ 0.564, 0.863, 0.853, 0.814, 1.000, 0.647 }},
                {{ 0.313, 0.246, 0.221, 0.389, 0.647, 1.000 }}
            }},
        0u
        };
        ///@}

        ///@name 2014
        ///@{
        /*
         * From [BFvD2014], based on reproduction of data points from [DLMW2012] and [FY2011], and subleading Isgur-Wise functions
         * approximations.
         */
        static const MultivariateGaussianConstraintTemplate<4> LambdaB_to_Lambda_all_v_and_a_0_BFvD2014
        {
            {{ "Lambda_b->Lambda::f_perp^V(s)", "Lambda_b->Lambda::f_perp^A(s)", "Lambda_b->Lambda::f_long^V(s)", "Lambda_b->Lambda::f_long^A(s)" }},
            {{ Kinematics{ { "s", 0.0 } }, Kinematics{ { "s", 0.0 } }, Kinematics{ { "s", 0.0 } }, Kinematics{ { "s", 0.0 } } }},
            {{ Options{ }, Options{ }, Options{ }, Options{ } }},
            {{ 0.391, 0.389, 0.380, 0.380 }},
            {{ 0.226, 0.224, 0.221, 0.221 }},
            {{ 0.226, 0.224, 0.221, 0.221 }},
            {{ 0.000, 0.000, 0.000, 0.000 }},
            {{
                 {{ 1.000, 0.556, 0.773, 0.771 }},
                 {{ 0.556, 1.000, 0.773, 0.772 }},
                 {{ 0.773, 0.773, 1.000, 0.534 }},
                 {{ 0.771, 0.772, 0.534, 1.000 }}
            }}
        };
        static const MultivariateGaussianConstraintTemplate<2> LambdaB_to_Lambda_fperpV_13dot5_to_20dot3_BFvD2014
        {
            {{ "Lambda_b->Lambda::f_perp^V(s)", "Lambda_b->Lambda::f_perp^V(s)" }},
            {{ Kinematics{ { "s", 13.5 } }, Kinematics{ { "s", 20.5 } } }},
            {{ Options{ }, Options{ } }},
            {{ 0.73, 1.40 }},
            {{ 0.20, 0.20 }},
            {{ 0.20, 0.20 }},
            {{ 0.00, 0.00 }},
            {{
                 {{ 1.000, 0.000 }},
                 {{ 0.000, 1.000 }}
            }}
        };
        static const MultivariateGaussianConstraintTemplate<2> LambdaB_to_Lambda_fperpA_13dot5_to_20dot3_BFvD2014
        {
            {{ "Lambda_b->Lambda::f_perp^A(s)", "Lambda_b->Lambda::f_perp^A(s)" }},
            {{ Kinematics{ { "s", 13.5 } }, Kinematics{ { "s", 20.5 } } }},
            {{ Options{ }, Options{ } }},
            {{ 0.48, 0.84 }},
            {{ 0.19, 0.19 }},
            {{ 0.19, 0.19 }},
            {{ 0.00, 0.00 }},
            {{
                 {{ 1.000, 0.000 }},
                 {{ 0.000, 1.000 }}
            }}
        };
        static const MultivariateGaussianConstraintTemplate<2> LambdaB_to_Lambda_flongV_13dot5_to_20dot3_BFvD2014
        {
            {{ "Lambda_b->Lambda::f_long^V(s)", "Lambda_b->Lambda::f_long^V(s)" }},
            {{ Kinematics{ { "s", 13.5 } }, Kinematics{ { "s", 20.5 } } }},
            {{ Options{ }, Options{ } }},
            {{ 0.72, 1.39 }},
            {{ 0.21, 0.21 }},
            {{ 0.21, 0.21 }},
            {{ 0.00, 0.00 }},
            {{
                 {{ 1.000, 0.000 }},
                 {{ 0.000, 1.000 }}
            }}
        };
        static const MultivariateGaussianConstraintTemplate<2> LambdaB_to_Lambda_flongA_13dot5_to_20dot3_BFvD2014
        {
            {{ "Lambda_b->Lambda::f_long^A(s)", "Lambda_b->Lambda::f_long^A(s)" }},
            {{ Kinematics{ { "s", 13.5 } }, Kinematics{ { "s", 20.5 } } }},
            {{ Options{ }, Options{ } }},
            {{ 0.48, 0.85 }},
            {{ 0.19, 0.18 }},
            {{ 0.19, 0.18 }},
            {{ 0.00, 0.00 }},
            {{
                 {{ 1.000, 0.000 }},
                 {{ 0.000, 1.000 }}
            }}
        };
        ///@}
    }

    /* Constraint */
    template class WrappedForwardIterator<Constraint::BlockIteratorTag, LogLikelihoodBlockPtr>;
    template class WrappedForwardIterator<Constraint::ObservableIteratorTag, ObservablePtr>;

    template <>
    struct Implementation<Constraint>
    {
        std::string name;

        ObservableSet observables;

        std::vector<LogLikelihoodBlockPtr> blocks;

        Implementation(const std::string & name,
                const std::vector<ObservablePtr> & observables,
                const std::vector<LogLikelihoodBlockPtr> & blocks) :
            name(name),
            blocks(blocks)
        {
            for (auto o = observables.begin(), o_end = observables.end() ; o != o_end ; ++o)
            {
                this->observables.add(*o);
            }
        }
    };

    Constraint::Constraint(const std::string & name,
                const std::vector<ObservablePtr> & observables,
                const std::vector<LogLikelihoodBlockPtr> & blocks) :
        PrivateImplementationPattern<Constraint>(new Implementation<Constraint>(name, observables, blocks))
    {
    }

    Constraint::~Constraint()
    {
    }

    std::string
    Constraint::name() const
    {
        return _imp->name;
    }

    Constraint::BlockIterator
    Constraint::begin_blocks() const
    {
        return BlockIterator(_imp->blocks.begin());
    }

    Constraint::BlockIterator
    Constraint::end_blocks() const
    {
        return BlockIterator(_imp->blocks.end());
    }

    Constraint::ObservableIterator
    Constraint::begin_observables() const
    {
        return ObservableIterator(_imp->observables.begin());
    }

    Constraint::ObservableIterator
    Constraint::end_observables() const
    {
        return ObservableIterator(_imp->observables.end());
    }

    typedef std::function<Constraint (const std::string &, const Options & options)> ConstraintFactory;

    template <typename Factory_>
    ConstraintFactory make_factory(const Factory_ & f)
    {
        return std::bind(&Factory_::make, f, std::placeholders::_1, std::placeholders::_2);
    }

    /*
     * Adding a new constraint:
     * 1. Instantiate one of the existing ConstraintTemplate in namespace templates{...}
     * 2. Add an entry to the map in Constraint::make
     * 3. Add the constraint name to constraint_TEST.cc
     * 4. Run constraint_TEST and check text output for new constraint
     */
    Constraint
    Constraint::make(const std::string & name, const Options & options)
    {
        static const std::map<std::string, ConstraintFactory> factories
        {
            /* 2000 */
            // CLEO
            { "B^0->K^*0gamma::BR@CLEO-2000", make_factory(templates::Bzero_to_Kstarzero_gamma_BR_CLEO_2000) },
            { "B^+->K^*+gamma::BR@CLEO-2000", make_factory(templates::Bplus_to_Kstarplus_gamma_BR_CLEO_2000) },
            /* 2004 */
            // BaBar
            { "B->X_sll::BR[1.0,6.0]@BaBar-2004A", make_factory(templates::Bmix_to_Xs_dilepton_BR_BaBar_2004A) },
            // Belle
            { "B^0->K^*0gamma::BR@Belle-2004", make_factory(templates::Bzero_to_Kstarzero_gamma_BR_Belle_2004) },
            { "B^+->K^*+gamma::BR@Belle-2004", make_factory(templates::Bplus_to_Kstarplus_gamma_BR_Belle_2004) },
            /* 2005 */
            // Belle
            { "B->X_sll::BR[1.0,6.0]@Belle-2005A", make_factory(templates::Bmix_to_Xs_dilepton_BR_Belle_2005A) },
            /* 2006 */
            // Belle
            { "B^0->K^*0gamma::S_K@Belle-2006", make_factory(templates::Bzero_to_Kstarzero_gamma_SKstargamma_Belle_2006) },
            { "B^0->K^*0gamma::C_K@Belle-2006", make_factory(templates::Bzero_to_Kstarzero_gamma_CKstargamma_Belle_2006) },
            { "B^0->K^*0gamma::S_K+C_K@Belle-2006", make_factory(templates::Bzero_to_Kstarzero_gamma_time_dependent_cp_asymmetries_Belle_2006) },
            /* 2008 */
            // BaBar
            { "B^0->K^*0gamma::S_K@BaBar-2008", make_factory(templates::Bzero_to_Kstarzero_gamma_SKstargamma_BaBar_2008) },
            { "B^0->K^*0gamma::C_K@BaBar-2008", make_factory(templates::Bzero_to_Kstarzero_gamma_CKstargamma_BaBar_2008) },
            { "B^0->K^*0gamma::S_K+C_K@BaBar-2008", make_factory(templates::Bzero_to_Kstarzero_gamma_time_dependent_cp_asymmetries_BaBar_2008) },
            // Belle
            { "B->X_sgamma::E_1[1.8]@Belle-2008", make_factory(templates::B_to_Xs_gamma_E_1_1dot8_Belle_2008A) },
            { "B->X_sgamma::E_2[1.8]@Belle-2008", make_factory(templates::B_to_Xs_gamma_E_2_1dot8_Belle_2008A) },
            { "B->X_sgamma::E_1[1.8]+E_2[1.8]@Belle-2008", make_factory(templates::B_to_Xs_gamma_E_1_and_E_2_1dot8_Belle_2008A) },
            /* 2009 */
            // BaBar
            { "B^0->K^*0gamma::BR@BaBar-2009", make_factory(templates::Bzero_to_Kstarzero_gamma_BR_BaBar_2009) },
            { "B^+->K^*+gamma::BR@BaBar-2009", make_factory(templates::Bplus_to_Kstarplus_gamma_BR_BaBar_2009) },
            // Belle
            // B^+ -> K^+ mu^+ mu^-
            { "B^+->K^+mu^+mu^-::BR[1.00,6.00]@Belle-2009", make_factory(templates::Bplus_to_Kplus_dimuon_BR_1_to_6_Belle_2009) },
            { "B^+->K^+mu^+mu^-::BR[14.18,16.00]@Belle-2009", make_factory(templates::Bplus_to_Kplus_dimuon_BR_14dot18_to_16_Belle_2009) },
            { "B^+->K^+mu^+mu^-::BR[16.00,22.86]@Belle-2009", make_factory(templates::Bplus_to_Kplus_dimuon_BR_16_to_22dot86_Belle_2009) },

            { "B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@Belle-2009", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_1_to_6_Belle_2009) },
            { "B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@Belle-2009", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_Belle_2009) },
            { "B^+->K^+mu^+mu^-::A_FB[16.00,22.86]@Belle-2009", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_16_to_22dot86_Belle_2009) },

            // B^0 -> K^*0 mu^+ mu^-
            { "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@Belle-2009", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_1_to_6_Belle_2009) },
            { "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@Belle-2009", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_Belle_2009) },
            { "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@Belle-2009", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_Belle_2009) },
            { "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@Belle-2009", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_Belle_2009) },
            { "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@Belle-2009", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_Belle_2009) },
            { "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@Belle-2009", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_Belle_2009) },
            { "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@Belle-2009", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_1_to_6_Belle_2009) },
            { "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@Belle-2009", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_Belle_2009) },
            { "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@Belle-2009", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_Belle_2009) },
            // B->X_sgamma
            { "B->X_sgamma::BR[1.8]@Belle-2009B", make_factory(templates::B_to_Xs_gamma_BR_1dot8_Belle_2009B) },
            { "B->X_sgamma::BR[1.8]+E_1[1.8]+E_2[1.8]@Belle-2009B", make_factory(templates::B_to_Xs_gamma_1dot8_Belle_2009B) },

            /* 2010 */
            // BaBar
            { "B^0->pi^+lnu::BR[0.0,4.0]@BaBar-2010A", make_factory(templates::Bzero_to_pi_l_nu_BR_0_to_4_BaBar_2010A) },
            { "B^0->pi^+lnu::BR[4.0,8.0]@BaBar-2010A", make_factory(templates::Bzero_to_pi_l_nu_BR_4_to_8_BaBar_2010A) },
            { "B^0->pi^+lnu::BR[8.0,12.0]@BaBar-2010A", make_factory(templates::Bzero_to_pi_l_nu_BR_8_to_12_BaBar_2010A) },
            { "B^0->pi^+lnu::BR@BaBar-2010B", make_factory(templates::Bzero_to_pi_l_nu_BR_BaBar_2010B) },
            //Belle
            { "B^0->pi^+lnu::BR@Belle-2010A", make_factory(templates::Bzero_to_pi_l_nu_BR_Belle_2010A) },

            /* 2011 */
            // HFAG
            { "B^0->K^*0gamma::S_K+C_K@HFAG-2011", make_factory(templates::Bzero_to_Kstarzero_gamma_time_dependent_cp_asymmetries_HFAG_2011) },
            // CDF
            // B^0 -> K^*0 mu^+ mu^-
            { "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@CDF-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_1_to_6_CDF_2011) },
            { "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@CDF-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_CDF_2011) },
            { "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@CDF-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_CDF_2011) },
            // B^+ -> K^*+ mu^+ mu^-
            { "B^+->K^*+mu^+mu^-::BR[1.00,6.00]@CDF-2011", make_factory(templates::Bplus_to_Kstarplus_dimuon_BR_1_to_6_CDF_2011) },
            { "B^+->K^*+mu^+mu^-::BR[14.18,16.00]@CDF-2011", make_factory(templates::Bplus_to_Kstarplus_dimuon_BR_14dot18_to_16_CDF_2011) },
            { "B^+->K^*+mu^+mu^-::BR[16.00,19.21]@CDF-2011", make_factory(templates::Bplus_to_Kstarplus_dimuon_BR_16_to_19dot21_CDF_2011) },
            // B^0 -> K^*0 mu^+ mu^-
            { "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@CDF-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_CDF_2011) },
            { "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@CDF-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_CDF_2011) },
            { "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@CDF-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_CDF_2011) },
            { "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@CDF-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_1_to_6_CDF_2011) },
            { "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@CDF-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_CDF_2011) },
            { "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@CDF-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_CDF_2011) },
            { "B^0->K^*0mu^+mu^-::A_T^2[1.00,6.00]@CDF-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_T_2_1_to_6_CDF_2011) },
            { "B^0->K^*0mu^+mu^-::A_T^2[14.18,16.00]@CDF-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_T_2_14dot18_to_16_CDF_2011) },
            { "B^0->K^*0mu^+mu^-::A_T^2[16.00,19.21]@CDF-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_T_2_16_to_19dot21_CDF_2011) },
            { "B^0->K^*0mu^+mu^-::A_im[1.00,6.00]@CDF-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_im_1_to_6_CDF_2011) },
            { "B^0->K^*0mu^+mu^-::A_im[14.18,16.00]@CDF-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_im_14dot18_to_16_CDF_2011) },
            { "B^0->K^*0mu^+mu^-::A_im[16.00,19.21]@CDF-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_im_16_to_19dot21_CDF_2011) },
            // B^0 -> K^0 mu^+ mu^-
            { "B^0->K^0mu^+mu^-::BR[1.00,6.00]@CDF-2011", make_factory(templates::Bzero_to_Kzero_dimuon_BR_1_to_6_CDF_2011) },
            { "B^0->K^0mu^+mu^-::BR[14.18,16.00]@CDF-2011", make_factory(templates::Bzero_to_Kzero_dimuon_BR_14dot18_to_16_CDF_2011) },
            { "B^0->K^0mu^+mu^-::BR[16.00,22.86]@CDF-2011", make_factory(templates::Bzero_to_Kzero_dimuon_BR_16_to_22dot86_CDF_2011) },
            // B^+ -> K^+ mu^+ mu^-
            { "B^+->K^+mu^+mu^-::BR[1.00,6.00]@CDF-2011", make_factory(templates::Bplus_to_Kplus_dimuon_BR_1_to_6_CDF_2011) },
            { "B^+->K^+mu^+mu^-::BR[14.18,16.00]@CDF-2011", make_factory(templates::Bplus_to_Kplus_dimuon_BR_14dot18_to_16_CDF_2011) },
            { "B^+->K^+mu^+mu^-::BR[16.00,22.86]@CDF-2011", make_factory(templates::Bplus_to_Kplus_dimuon_BR_16_to_22dot86_CDF_2011) },

            { "B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@CDF-2011", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_1_to_6_CDF_2011) },
            { "B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@CDF-2011", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_CDF_2011) },
            { "B^+->K^+mu^+mu^-::A_FB[16.00,22.86]@CDF-2011", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_16_to_22dot86_CDF_2011) },

            // limit on B^0_s -> mu^+ mu^-
            { "B^0_s->mu^+mu^-::BR_limit@CDF-2011", make_factory(templates::Bzero_to_dimuon_CDF_2011) },
            // LHCb
            { "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@LHCb-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_1_to_6_LHCb_2011) },
            { "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@LHCb-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_LHCb_2011) },
            { "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@LHCb-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_LHCb_2011) },
            { "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@LHCb-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_LHCb_2011) },
            { "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@LHCb-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_LHCb_2011) },
            { "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@LHCb-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_LHCb_2011) },
            { "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@LHCb-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_1_to_6_LHCb_2011) },
            { "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@LHCb-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_LHCb_2011) },
            { "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@LHCb-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_LHCb_2011) },

            // limit on B^0_s -> mu^+ mu^-
            // LHCb + CMS
            { "B^0_s->mu^+mu^-::BR_limit@LHCb-CMS-2011", make_factory(templates::Bzero_to_dimuon_LHCb_CMS_2011) },
            { "B^0_s->mu^+mu^-::BR_limit@LHCb-CMS-2011-Bayes", make_factory(templates::Bzero_to_dimuon_LHCb_CMS_2011_Bayes) },

            /* 2012 */
            // BaBar
            // B^+ -> K^+ mu^+ mu^-
            { "B^+->K^+mu^+mu^-::BR[1.00,6.00]@BaBar-2012", make_factory(templates::Bplus_to_Kplus_dimuon_BR_1_to_6_BaBar_2012) },
            { "B^+->K^+mu^+mu^-::BR[14.21,16.00]@BaBar-2012", make_factory(templates::Bplus_to_Kplus_dimuon_BR_14dot21_to_16_BaBar_2012) },
            { "B^+->K^+mu^+mu^-::BR[16.00,22.86]@BaBar-2012", make_factory(templates::Bplus_to_Kplus_dimuon_BR_16_to_22dot86_BaBar_2012) },
            // B^0 -> K^*0 mu^+ mu^-
            { "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@BaBar-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_1_to_6_BaBar_2012) },
            { "B^0->K^*0mu^+mu^-::BR[14.21,16.00]@BaBar-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_14dot21_to_16_BaBar_2012) },
            { "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@BaBar-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_BaBar_2012) },
            { "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@BaBar-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_BaBar_2012) },
            { "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@BaBar-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_BaBar_2012) },
            { "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@BaBar-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_BaBar_2012) },
            { "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@BaBar-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_1_to_6_BaBar_2012) },
            { "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@BaBar-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_BaBar_2012) },
            { "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@BaBar-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_BaBar_2012) },
            // B -> X_s gamma
            { "B->X_sgamma::BR[1.8]@BaBar-2012", make_factory(templates::B_to_Xs_gamma_BR_1dot8_BaBar_2012C) },
            { "B->X_sgamma::E_1[1.8]@BaBar-2012", make_factory(templates::B_to_Xs_gamma_E_1_1dot8_BaBar_2012C) },
            { "B->X_sgamma::E_2[1.8]@BaBar-2012", make_factory(templates::B_to_Xs_gamma_E_2_1dot8_BaBar_2012C) },
            { "B->X_sgamma::E_1[1.8]+E_2[1.8]@BaBar-2012", make_factory(templates::B_to_Xs_gamma_E_1_and_E_2_1dot8_BaBar_2012C) },
            // B^0 -> pi^- l nu
            { "B^0->pi^+lnu::BR@BaBar-2012D", make_factory(templates::Bzero_to_pi_l_nu_BR_BaBar_2012D) },
            // CDF
            // B^0 -> K^*0 mu^+ mu^-
            { "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@CDF-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_1_to_6_CDF_2012) },
            { "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@CDF-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_CDF_2012) },
            { "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@CDF-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_CDF_2012) },
            // B^+ -> K^*+ mu^+ mu^-
            { "B^+->K^*+mu^+mu^-::BR[1.00,6.00]@CDF-2012", make_factory(templates::Bplus_to_Kstarplus_dimuon_BR_1_to_6_CDF_2012) },
            { "B^+->K^*+mu^+mu^-::BR[14.18,16.00]@CDF-2012", make_factory(templates::Bplus_to_Kstarplus_dimuon_BR_14dot18_to_16_CDF_2012) },
            { "B^+->K^*+mu^+mu^-::BR[16.00,19.21]@CDF-2012", make_factory(templates::Bplus_to_Kstarplus_dimuon_BR_16_to_19dot21_CDF_2012) },
            // B^0 -> K^*0 mu^+ mu^-
            { "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@CDF-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_CDF_2012) },
            { "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@CDF-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_CDF_2012) },
            { "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@CDF-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_CDF_2012) },
            { "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@CDF-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_1_to_6_CDF_2012) },
            { "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@CDF-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_CDF_2012) },
            { "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@CDF-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_CDF_2012) },
            { "B^0->K^*0mu^+mu^-::A_T^2[1.00,6.00]@CDF-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_T_2_1_to_6_CDF_2012) },
            { "B^0->K^*0mu^+mu^-::A_T^2[14.18,16.00]@CDF-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_T_2_14dot18_to_16_CDF_2012) },
            { "B^0->K^*0mu^+mu^-::A_T^2[16.00,19.21]@CDF-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_T_2_16_to_19dot21_CDF_2012) },
            { "B^0->K^*0mu^+mu^-::A_im[1.00,6.00]@CDF-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_im_1_to_6_CDF_2012) },
            { "B^0->K^*0mu^+mu^-::A_im[14.18,16.00]@CDF-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_im_14dot18_to_16_CDF_2012) },
            { "B^0->K^*0mu^+mu^-::A_im[16.00,19.21]@CDF-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_im_16_to_19dot21_CDF_2012) },
            // B^0 -> K^0 mu^+ mu^-
            { "B^0->K^0mu^+mu^-::BR[1.00,6.00]@CDF-2012", make_factory(templates::Bzero_to_Kzero_dimuon_BR_1_to_6_CDF_2012) },
            { "B^0->K^0mu^+mu^-::BR[14.18,16.00]@CDF-2012", make_factory(templates::Bzero_to_Kzero_dimuon_BR_14dot18_to_16_CDF_2012) },
            { "B^0->K^0mu^+mu^-::BR[16.00,22.86]@CDF-2012", make_factory(templates::Bzero_to_Kzero_dimuon_BR_16_to_22dot86_CDF_2012) },
            // B^+ -> K^+ mu^+ mu^-
            { "B^+->K^+mu^+mu^-::BR[1.00,6.00]@CDF-2012", make_factory(templates::Bplus_to_Kplus_dimuon_BR_1_to_6_CDF_2012) },
            { "B^+->K^+mu^+mu^-::BR[14.18,16.00]@CDF-2012", make_factory(templates::Bplus_to_Kplus_dimuon_BR_14dot18_to_16_CDF_2012) },
            { "B^+->K^+mu^+mu^-::BR[16.00,22.86]@CDF-2012", make_factory(templates::Bplus_to_Kplus_dimuon_BR_16_to_22dot86_CDF_2012) },

            { "B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@CDF-2012", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_1_to_6_CDF_2012) },
            { "B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@CDF-2012", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_CDF_2012) },
            { "B^+->K^+mu^+mu^-::A_FB[16.00,22.86]@CDF-2012", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_16_to_22dot86_CDF_2012) },

            // LHCb
            // B^0 -> K^*0 mu^+ mu^-
            { "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@LHCb-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_1_to_6_LHCb_2012) },
            { "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@LHCb-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_LHCb_2012) },
            { "B^0->K^*0mu^+mu^-::BR[16.00,19.00]@LHCb-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_16_to_19_LHCb_2012) },
            { "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@LHCb-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_LHCb_2012) },
            { "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@LHCb-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_LHCb_2012) },
            { "B^0->K^*0mu^+mu^-::A_FB[16.00,19.00]@LHCb-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_16_to_19_LHCb_2012) },
            { "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@LHCb-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_1_to_6_LHCb_2012) },
            { "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@LHCb-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_LHCb_2012) },
            { "B^0->K^*0mu^+mu^-::F_L[16.00,19.00]@LHCb-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_16_to_19_LHCb_2012) },
            { "B^0->K^*0mu^+mu^-::A_im[1.00,6.00]@LHCb-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_im_1_to_6_LHCb_2012) },
            { "B^0->K^*0mu^+mu^-::A_im[14.18,16.00]@LHCb-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_im_14dot18_to_16_LHCb_2012) },
            { "B^0->K^*0mu^+mu^-::A_im[16.00,19.00]@LHCb-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_im_16_to_19_LHCb_2012) },
            { "B^0->K^*0mu^+mu^-::S_3[1.00,6.00]@LHCb-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_S_3_1_to_6_LHCb_2012) },
            { "B^0->K^*0mu^+mu^-::S_3[14.18,16.00]@LHCb-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_S_3_14dot18_to_16_LHCb_2012) },
            { "B^0->K^*0mu^+mu^-::S_3[16.00,19.00]@LHCb-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_S_3_16_to_19_LHCb_2012) },
            { "B^0->K^*0mu^+mu^-::A_CP[1.00,6.00]@LHCb-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_CP_1_to_6_LHCb_2012E) },
            { "B^0->K^*0mu^+mu^-::A_CP[14.18,16.00]@LHCb-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_CP_14dot18_to_16_LHCb_2012E) },
            { "B^0->K^*0mu^+mu^-::A_CP[16.00,20.00]@LHCb-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_CP_16_to_20_LHCb_2012E) },
            // limit on B^0_s -> mu^+ mu^-
            { "B^0_s->mu^+mu^-::BR_limit@LHCb-2012", make_factory(templates::Bzero_to_dimuon_LHCb_2012) },
            // limit on B^0_s -> mu^+ mu^- of Nov 2012
            { "B^0_s->mu^+mu^-::BR_limit@LHCb-Nov-2012", make_factory(templates::Bzero_to_dimuon_LHCb_Nov_2012) },
            // B^+ -> K^+ mu^+ mu^-
            { "B^+->K^+mu^+mu^-::BR[1.00,6.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_BR_1_to_6_LHCb_2012) },
            { "B^+->K^+mu^+mu^-::BR[14.18,16.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_BR_14dot18_to_16_LHCb_2012) },
            { "B^+->K^+mu^+mu^-::BR[16.00,18.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_BR_16_to_18_LHCb_2012) },
            { "B^+->K^+mu^+mu^-::BR[18.00,22.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_BR_18_to_22_LHCb_2012) },

            { "B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_1_to_6_LHCb_2012) },
            { "B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_LHCb_2012) },
            { "B^+->K^+mu^+mu^-::A_FB[16.00,18.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_16_to_18_LHCb_2012) },
            { "B^+->K^+mu^+mu^-::A_FB[18.00,22.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_18_to_22_LHCb_2012) },

            { "B^+->K^+mu^+mu^-::F_H[1.00,6.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_F_H_1_to_6_LHCb_2012) },
            { "B^+->K^+mu^+mu^-::F_H[14.18,16.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_F_H_14dot18_to_16_LHCb_2012) },
            { "B^+->K^+mu^+mu^-::F_H[16.00,18.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_F_H_16_to_18_LHCb_2012) },
            { "B^+->K^+mu^+mu^-::F_H[18.00,22.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_F_H_18_to_22_LHCb_2012) },
            // PDG2012
            // B^0(*) Mass splitting
            { "B^0::M_B^*-M_B@PDG-2012", make_factory(templates::B_Bstar_mass_splitting_PDG_2012) },

            /* 2013 */
            // ATLAS
            { "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@ATLAS-2013A", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_ATLAS_2013A) },
            { "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@ATLAS-2013A", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_ATLAS_2013A) },
            { "B^0->K^*0mu^+mu^-::A_FB[16.00,19.00]@ATLAS-2013A", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_16_to_19_ATLAS_2013A) },
            { "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@ATLAS-2013A", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_1_to_6_ATLAS_2013A) },
            { "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@ATLAS-2013A", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_ATLAS_2013A) },
            { "B^0->K^*0mu^+mu^-::F_L[16.00,19.00]@ATLAS-2013A", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_16_to_19_ATLAS_2013A) },
            //Belle
            { "B^0->pi^+lnu::BR@Belle-2013A", make_factory(templates::Bzero_to_pi_l_nu_BR_Belle_2013A) },
            // CMS
            { "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@CMS-2013A", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_1_to_6_CMS_2013A) },
            { "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@CMS-2013A", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_CMS_2013A) },
            { "B^0->K^*0mu^+mu^-::BR[16.00,19.00]@CMS-2013A", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_16_to_19_CMS_2013A) },
            { "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@CMS-2013A", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_CMS_2013A) },
            { "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@CMS-2013A", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_CMS_2013A) },
            { "B^0->K^*0mu^+mu^-::A_FB[16.00,19.00]@CMS-2013A", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_16_to_19_CMS_2013A) },
            { "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@CMS-2013A", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_1_to_6_CMS_2013A) },
            { "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@CMS-2013A", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_CMS_2013A) },
            { "B^0->K^*0mu^+mu^-::F_L[16.00,19.00]@CMS-2013A", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_16_to_19_CMS_2013A) },
            { "B^0_s->mu^+mu^-::BR@CMS-2013B", make_factory(templates::Bzero_to_dimuon_CMS_2013B) },
            // LHCb
            { "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_1_to_6_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::BR[16.00,19.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_BR_16_to_19_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_1_to_6_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::F_L[16.00,19.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_F_L_16_to_19_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::A_FB[16.00,19.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_FB_16_to_19_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::S_3[1.00,6.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_S_3_1_to_6_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::S_3[14.18,16.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_S_3_14dot18_to_16_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::S_3[16.00,19.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_S_3_16_to_19_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::S_9[1.00,6.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_S_9_1_to_6_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::S_9[14.18,16.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_S_9_14dot18_to_16_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::S_9[16.00,19.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_S_9_16_to_19_LHCb_2013B) },
            // The following observables have not yet been implemented.
#if 0
            { "B^0->K^*0mu^+mu^-::A_9[1.00,6.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_9_1_to_6_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::A_9[14.18,16.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_9_14dot18_to_16_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::A_9[16.00,19.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_9_16_to_19_LHCb_2013B) },
#endif
            { "B^0->K^*0mu^+mu^-::A_T^2[1.00,6.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_T2_1_to_6_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::A_T^2[14.18,16.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_T2_14dot18_to_16_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::A_T^2[16.00,19.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_T2_16_to_19_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::A_T^re[1.00,6.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_Tre_1_to_6_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::A_T^re[14.18,16.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_Tre_14dot18_to_16_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::A_T^re[16.00,19.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_Tre_16_to_19_LHCb_2013B) },
            { "B^0->K^*0mu^+mu^-::P'_4[1.00,6.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_Pprime_4_1_to_6_LHCb_2013C) },
            { "B^0->K^*0mu^+mu^-::P'_4[14.18,16.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_Pprime_4_14dot18_to_16_LHCb_2013C) },
            { "B^0->K^*0mu^+mu^-::P'_4[16.00,19.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_Pprime_4_16_to_19_LHCb_2013C) },
            { "B^0->K^*0mu^+mu^-::P'_5[1.00,6.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_Pprime_5_1_to_6_LHCb_2013C) },
            { "B^0->K^*0mu^+mu^-::P'_5[14.18,16.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_Pprime_5_14dot18_to_16_LHCb_2013C) },
            { "B^0->K^*0mu^+mu^-::P'_5[16.00,19.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_Pprime_5_16_to_19_LHCb_2013C) },
            { "B^0->K^*0mu^+mu^-::P'_6[1.00,6.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_Pprime_6_1_to_6_LHCb_2013C) },
            { "B^0->K^*0mu^+mu^-::P'_6[14.18,16.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_Pprime_6_14dot18_to_16_LHCb_2013C) },
            { "B^0->K^*0mu^+mu^-::P'_6[16.00,19.00]@LHCb-2013", make_factory(templates::Bzero_to_Kstarzero_dimuon_Pprime_6_16_to_19_LHCb_2013C) },
            { "B^0_s->mu^+mu^-::BR@LHCb-2013D", make_factory(templates::Bzero_to_dimuon_LHCb_2013D) },

            /* 2014 */
            // LHCb
            { "B^+->K^+mu^+mu^-::BR[1.10,2.00]@LHCb-2014", make_factory(templates::Bplus_to_Kplus_dimuon_BR_1dot1_to_2_LHCb_2014) },
            { "B^+->K^+mu^+mu^-::BR[2.00,3.00]@LHCb-2014", make_factory(templates::Bplus_to_Kplus_dimuon_BR_2_to_3_LHCb_2014) },
            { "B^+->K^+mu^+mu^-::BR[3.00,4.00]@LHCb-2014", make_factory(templates::Bplus_to_Kplus_dimuon_BR_3_to_4_LHCb_2014) },
            { "B^+->K^+mu^+mu^-::BR[4.00,5.00]@LHCb-2014", make_factory(templates::Bplus_to_Kplus_dimuon_BR_4_to_5_LHCb_2014) },
            { "B^+->K^+mu^+mu^-::BR[5.00,6.00]@LHCb-2014", make_factory(templates::Bplus_to_Kplus_dimuon_BR_5_to_6_LHCb_2014) },
            { "B^+->K^+mu^+mu^-::BR[1.10,6.00]@LHCb-2014", make_factory(templates::Bplus_to_Kplus_dimuon_BR_1dot1_to_6_LHCb_2014) },
            { "B^+->K^+mu^+mu^-::BR[15.00,22.00]@LHCb-2014", make_factory(templates::Bplus_to_Kplus_dimuon_BR_15_to_22_LHCb_2014) },
            { "B^+->K^+mu^+mu^-::A_FB[1.10,6.00]@LHCb-2014", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_1dot1_to_6_LHCb_2014) },
            { "B^+->K^+mu^+mu^-::A_FB[15.00,22.00]@LHCb-2014", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_15_to_22_LHCb_2014) },
            { "B^+->K^+mu^+mu^-::F_H[1.10,6.00]@LHCb-2014", make_factory(templates::Bplus_to_Kplus_dimuon_F_H_1dot1_to_6_LHCb_2014) },
            { "B^+->K^+mu^+mu^-::F_H[15.00,22.00]@LHCb-2014", make_factory(templates::Bplus_to_Kplus_dimuon_F_H_15_to_22_LHCb_2014) },

            /* Theory Constraints */
            { "B->K::f_+@HPQCD-2013A", make_factory(templates::B_to_K_fplus_17_to_23_HPQCD_2013A) },

            { "B->K^*::V@HPQCD-2013B", make_factory(templates::B_to_Kstar_V_15_to_19dot21_HPQCD_2013B) },
            { "B->K^*::A_1@HPQCD-2013B", make_factory(templates::B_to_Kstar_A1_15_to_19dot21_HPQCD_2013B) },
            { "B->K^*::A_12@HPQCD-2013B", make_factory(templates::B_to_Kstar_A12_15_to_19dot21_HPQCD_2013B) },

            { "B->pi::f_+@IKMvD-2014", make_factory(templates::B_to_pi_fp_IKMvD_2014) },

            { "Lambda_b->Lambda::f_perp+long^V+A@BFvD2014", make_factory(templates::LambdaB_to_Lambda_all_v_and_a_0_BFvD2014) },
            { "Lambda_b->Lambda::f_perp^V@BFvD2014", make_factory(templates::LambdaB_to_Lambda_fperpV_13dot5_to_20dot3_BFvD2014) },
            { "Lambda_b->Lambda::f_perp^A@BFvD2014", make_factory(templates::LambdaB_to_Lambda_fperpA_13dot5_to_20dot3_BFvD2014) },
            { "Lambda_b->Lambda::f_long^V@BFvD2014", make_factory(templates::LambdaB_to_Lambda_flongV_13dot5_to_20dot3_BFvD2014) },
            { "Lambda_b->Lambda::f_long^A@BFvD2014", make_factory(templates::LambdaB_to_Lambda_flongA_13dot5_to_20dot3_BFvD2014) },
        };

        auto f = factories.find(name);
        if (f == factories.end())
            throw UnknownConstraintError(name);

        return f->second(f->first, options);
    }
}
