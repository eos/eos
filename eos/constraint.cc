/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Danny van Dyk
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

    template <size_t dim_>
    struct MultivariateGaussianConstraintTemplate
    {
        std::array<std::string, dim_> observables;

        std::array<Kinematics, dim_> kinematics;

        std::array<Options, dim_> options;

        std::array<double, dim_> means;

        std::array<std::array<double, dim_>, dim_> covariance;

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

            LogLikelihoodBlockPtr block = LogLikelihoodBlock::MultivariateGaussian(cache, observables, this->means, this->covariance);

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
            {{
                {{ 0.1321,  0.0071 }},  /* Use correlation of the two results from */
                {{ 0.0071,  0.0601 }},  /* http://www.slac.stanford.edu/xorg/hfag/triangle/moriond2011/index.shtml#bsgamma
                                         * to calculate covariance matrix. Use the larger of the two statistical
                                         * uncertainties of S_K to form a combined, 2D block.
                                         */
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
            {{
                {{ 0.0850,  0.0024 }},  // Use correlation of the two results from
                {{ 0.0024,  0.0265 }},  // http://www.slac.stanford.edu/xorg/hfag/triangle/moriond2011/index.shtml#bsgamma to calculate covariance matrix.
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
        ///@}

        ///@name 2011 Data
        ///@{
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
        // BR in [16.00, 23.00]
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
        // BR in [16.00, 23.00]
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
        // A_T_2 in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_T_2_14dot18_to_16_CDF_2011
        {
            "B->K^*ll::A_T^2avg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.4, +0.8, -0.8, +0.2, -0.2
        };
        // A_T_2 in [16.00, 19.21]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_T_2_16_to_19dot21_CDF_2011
        {
            "B->K^*ll::A_T^2avg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.9, +0.8, -0.8, +0.4, -0.4
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
         * LHCb Collaboration
         *
         * Data taken from talk by M. Patel at EPS-HEP 2011
         */
        // BR in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_1_to_6_LHCb_2011
        {
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.9500e-7, +0.3000e-7, -0.3000e-7, +0.1000e-7, -0.1000e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_LHCb_2011
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.0738e-7, +0.1830e-7, -0.1830e-7, +0.0546e-7, -0.0546e-7
        };
        // BR in [16.00, 19.00] (in the preliminary results, LHCb only integrated up to exactly 19.00 GeV^2!)
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_LHCb_2011
        {
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.4400e-7, +0.2400e-7, -0.2400e-7, +0.0800e-7, -0.0800e-7
        };
        // A_FB in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_LHCb_2011
        {
            "B->K^*ll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.10, +0.14, -0.14, +0.05, -0.05
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_LHCb_2011
        {
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.50, +0.09, -0.06, +0.03, -0.03
        };
        // A_FB in [16.00, 19.00] (in the preliminary results, LHCb only integrated up to exactly 19.00 GeV^2!)
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_LHCb_2011
        {
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.10, +0.13, -0.13, +0.06, -0.06
        };
        // F_L in [1.0, 6.0]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_1_to_6_LHCb_2011
        {
            "B->K^*ll::F_Lavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.57, +0.11, -0.10, +0.03, -0.03
        };
        // F_L in [14.18, 16.00]
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_LHCb_2011
        {
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.33, +0.11, -0.08, +0.04, -0.04
        };
        // F_L in [16.00, 19.00] (in the preliminary results, LHCb only integrated up to exactly 19.00 GeV^2!)
        static const GaussianConstraintTemplate Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_LHCb_2011
        {
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.28, +0.10, -0.09, +0.04, -0.04
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
            /* 2011 */
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
            // B^0 -> K^0 mu^+ mu^-
            { "B^0->K^0mu^+mu^-::BR[1.00,6.00]@CDF-2011", make_factory(templates::Bzero_to_Kzero_dimuon_BR_1_to_6_CDF_2011) },
            { "B^0->K^0mu^+mu^-::BR[14.18,16.00]@CDF-2011", make_factory(templates::Bzero_to_Kzero_dimuon_BR_14dot18_to_16_CDF_2011) },
            { "B^0->K^0mu^+mu^-::BR[16.00,23.00]@CDF-2011", make_factory(templates::Bzero_to_Kzero_dimuon_BR_16_to_22dot86_CDF_2011) },
            // B^+ -> K^+ mu^+ mu^-
            { "B^+->K^+mu^+mu^-::BR[1.00,6.00]@CDF-2011", make_factory(templates::Bplus_to_Kplus_dimuon_BR_1_to_6_CDF_2011) },
            { "B^+->K^+mu^+mu^-::BR[14.18,16.00]@CDF-2011", make_factory(templates::Bplus_to_Kplus_dimuon_BR_14dot18_to_16_CDF_2011) },
            { "B^+->K^+mu^+mu^-::BR[16.00,23.00]@CDF-2011", make_factory(templates::Bplus_to_Kplus_dimuon_BR_16_to_22dot86_CDF_2011) },
            // Observable no yet implemented!
            { "B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@CDF-2011", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_1_to_6_CDF_2011) },
            { "B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@CDF-2011", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_CDF_2011) },
            { "B^+->K^+mu^+mu^-::A_FB[16.00,23.00]@CDF-2011", make_factory(templates::Bplus_to_Kplus_dimuon_A_FB_16_to_22dot86_CDF_2011) },
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
        };

        auto f = factories.find(name);
        if (f == factories.end())
            throw UnknownConstraintError(name);

        return f->second(f->first, options);
    }
}
