/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2013 Danny van Dyk
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

    struct AmorosoConstraintTemplate
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
        // Yields 0 without implementation of scalar operators
#if 0
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
#endif
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
        // Yields 0 without implementation of scalar operators
#if 0
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
            "B_q->ll::BR",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 22.86 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            -0.09, +0.13, -0.17, +0.03, -0.03
        };
#endif

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
        static const AmorosoConstraintTemplate Bzero_to_dimuon_LHCb_CMS_2011_Bayes
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
            -0.02, +0.18, -0.16, +0.07, -0.07
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_BaBar_2012
        {
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.31, +0.19, -0.11, +0.13, -0.13
        };
        // A_FB in [16.00, 19.21]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_BaBar_2012
        {
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.34, +0.26, -0.17, +0.08, -0.08
        };
        // F_L in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_1_to_6_BaBar_2012
        {
            "B->K^*ll::F_Lavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.47, +0.13, -0.13, +0.04, -0.04
        };
        // F_L in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_BaBar_2012
        {
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.42, +0.12, -0.16, +0.11, -0.11
        };
        // F_L in [16.00, 19.21]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_BaBar_2012
        {
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.47, +0.18, -0.20, +0.13, -0.13
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
        // Yields 0 without implementation of scalar operators
#if 0
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
            "B_q->ll::BR",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 22.86 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            -0.05, +0.10, -0.18, +0.05, -0.05
        };
#endif

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
        // sign flipped wrt to table 1 in LHCb:2012C!
#if 0
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
#endif
        // F_H in [1.0, 6.0]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_F_H_1_to_6_LHCb_2012
        {
            "B->Kll::F_H@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.05, +0.08, -0.05, +0.04, -0.02
        };
        // F_H in [14.18, 16.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_F_H_14dot18_to_16_LHCb_2012
        {
            "B->Kll::F_H@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.08, +0.28, -0.08, +0.02, -0.01
        };
        // F_H in [16.00, 18.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_F_H_16_to_18_LHCb_2012
        {
            "B->Kll::F_H@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 18.00} },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.18, +0.22, -0.14, +0.01, -0.04
        };
        // F_H in [18.00, 22.00]
        static const GaussianConstraintTemplate Bplus_to_Kplus_dimuon_F_H_18_to_22_LHCb_2012
        {
            "B->Kll::F_H@LowRecoil",
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
        static const AmorosoConstraintTemplate Bzero_to_dimuon_LHCb_2012
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
            // Belle
            { "B^0->K^*0gamma::BR@Belle-2004", make_factory(templates::Bzero_to_Kstarzero_gamma_BR_Belle_2004) },
            { "B^+->K^*+gamma::BR@Belle-2004", make_factory(templates::Bplus_to_Kstarplus_gamma_BR_Belle_2004) },
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
            /* 2009 */
            // BaBar
            { "B^0->K^*0gamma::BR@BaBar-2009", make_factory(templates::Bzero_to_Kstarzero_gamma_BR_BaBar_2009) },
            { "B^+->K^*+gamma::BR@BaBar-2009", make_factory(templates::Bplus_to_Kstarplus_gamma_BR_BaBar_2009) },
            // Belle
            // B^+ -> K^+ mu^+ mu^-
            { "B^+->K^+mu^+mu^-::BR[1.00,6.00]@Belle-2009", make_factory(templates::Bplus_to_Kplus_dimuon_BR_1_to_6_Belle_2009) },
            { "B^+->K^+mu^+mu^-::BR[14.18,16.00]@Belle-2009", make_factory(templates::Bplus_to_Kplus_dimuon_BR_14dot18_to_16_Belle_2009) },
            { "B^+->K^+mu^+mu^-::BR[16.00,22.86]@Belle-2009", make_factory(templates::Bplus_to_Kplus_dimuon_BR_16_to_22dot86_Belle_2009) },
            // Observable no yet implemented!
#if 0
            { "B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@Belle-2009", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_1_to_6_Belle_2009) },
            { "B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@Belle-2009", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_Belle_2009) },
            { "B^+->K^+mu^+mu^-::A_FB[16.00,22.86]@Belle-2009", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_16_to_22dot86_Belle_2009) },
#endif
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
            { "B^0->K^*0mu^+mu^-::A_T_2[1.00,6.00]@CDF-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_T_2_1_to_6_CDF_2011) },
            { "B^0->K^*0mu^+mu^-::A_T_2[14.18,16.00]@CDF-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_T_2_14dot18_to_16_CDF_2011) },
            { "B^0->K^*0mu^+mu^-::A_T_2[16.00,19.21]@CDF-2011", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_T_2_16_to_19dot21_CDF_2011) },
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
            // Observables not yet implemented!
#if 0
            { "B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@CDF-2011", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_1_to_6_CDF_2011) },
            { "B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@CDF-2011", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_CDF_2011) },
            { "B^+->K^+mu^+mu^-::A_FB[16.00,22.86]@CDF-2011", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_16_to_22dot86_CDF_2011) },
#endif
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
            { "B^0->K^*0mu^+mu^-::A_T_2[1.00,6.00]@CDF-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_T_2_1_to_6_CDF_2012) },
            { "B^0->K^*0mu^+mu^-::A_T_2[14.18,16.00]@CDF-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_T_2_14dot18_to_16_CDF_2012) },
            { "B^0->K^*0mu^+mu^-::A_T_2[16.00,19.21]@CDF-2012", make_factory(templates::Bzero_to_Kstarzero_dimuon_A_T_2_16_to_19dot21_CDF_2012) },
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
            // Observables not yet implemented!
#if 0
            { "B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@CDF-2012", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_1_to_6_CDF_2012) },
            { "B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@CDF-2012", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_CDF_2012) },
            { "B^+->K^+mu^+mu^-::A_FB[16.00,22.86]@CDF-2012", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_16_to_22dot86_CDF_2012) },
#endif
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
            // limit on B^0_s -> mu^+ mu^-
            { "B^0_s->mu^+mu^-::BR_limit@LHCb-2012", make_factory(templates::Bzero_to_dimuon_LHCb_2012) },
            // limit on B^0_s -> mu^+ mu^- of Nov 2012
            { "B^0_s->mu^+mu^-::BR_limit@LHCb-Nov-2012", make_factory(templates::Bzero_to_dimuon_LHCb_Nov_2012) },
            // B^+ -> K^+ mu^+ mu^-
            { "B^+->K^+mu^+mu^-::BR[1.00,6.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_BR_1_to_6_LHCb_2012) },
            { "B^+->K^+mu^+mu^-::BR[14.18,16.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_BR_14dot18_to_16_LHCb_2012) },
            { "B^+->K^+mu^+mu^-::BR[16.00,18.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_BR_16_to_18_LHCb_2012) },
            { "B^+->K^+mu^+mu^-::BR[18.00,22.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_BR_18_to_22_LHCb_2012) },
            // Observables not yet implemented!
#if 0
            { "B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_1_to_6_LHCb_2012) },
            { "B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_LHCb_2012) },
            { "B^+->K^+mu^+mu^-::A_FB[16.00,18.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_16_to_18_LHCb_2012) },
            { "B^+->K^+mu^+mu^-::A_FB[18.00,22.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_18_to_22_LHCb_2012) },
#endif
            { "B^+->K^+mu^+mu^-::F_H[1.00,6.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_F_H_1_to_6_LHCb_2012) },
            { "B^+->K^+mu^+mu^-::F_H[14.18,16.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_F_H_14dot18_to_16_LHCb_2012) },
            { "B^+->K^+mu^+mu^-::F_H[16.00,18.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_F_H_16_to_18_LHCb_2012) },
            { "B^+->K^+mu^+mu^-::F_H[18.00,22.00]@LHCb-2012", make_factory(templates::Bplus_to_Kplus_dimuon_F_H_18_to_22_LHCb_2012) },
            // PDG2012
            // B^0(*) Mass splitting
            { "B^0::M_B^*-M_B@PDG-2012", make_factory(templates::B_Bstar_mass_splitting_PDG_2012) },
        };

        auto f = factories.find(name);
        if (f == factories.end())
            throw UnknownConstraintError(name);

        return f->second(f->first, options);
    }
}
