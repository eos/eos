/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011, 2012, 2013, 2014, 2015, 2016, 2017 Danny van Dyk
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
#include <eos/statistics/log-likelihood.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/observable_set.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qualified-name.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <cmath>
#include <map>
#include <vector>

namespace eos
{
    UnknownConstraintError::UnknownConstraintError(const QualifiedName & name) :
        Exception("Constraint '" + name.str() + "' is unknown")
    {
    }

    ConstraintEntry::~ConstraintEntry() = default;

    template <>
    struct WrappedForwardIteratorTraits<ConstraintEntry::ObservableNameIteratorTag>
    {
        typedef std::vector<QualifiedName>::const_iterator UnderlyingIterator;
    };
    template class WrappedForwardIterator<ConstraintEntry::ObservableNameIteratorTag, const QualifiedName>;

    class ConstraintEntryBase :
        public ConstraintEntry
    {
        protected:
            QualifiedName _name;

            std::vector<QualifiedName> _observable_names;

            ConstraintEntryBase(const QualifiedName & name,
                    std::vector<QualifiedName> && observable_names) :
                _name(name),
                _observable_names(std::move(observable_names))
            {
            }

        public:
            ConstraintEntryBase(const QualifiedName & name,
                    const QualifiedName & observable_name) :
                ConstraintEntryBase(name, std::vector<QualifiedName>{ observable_name })
            {
            }

            template <unsigned long n_>
            ConstraintEntryBase(const QualifiedName & name,
                    const std::array<QualifiedName, n_> & observable_names) :
                ConstraintEntryBase(name, std::vector<QualifiedName>(observable_names.begin(), observable_names.end()))
            {
            }

            ~ConstraintEntryBase() = default;

            virtual const QualifiedName & name() const { return _name; };

            virtual ConstraintEntry::ObservableNameIterator begin_observable_names() const
            {
                return _observable_names.begin();
            }

            virtual ConstraintEntry::ObservableNameIterator end_observable_names() const
            {
                return _observable_names.end();
            }
    };

    struct GaussianConstraintEntry :
        public ConstraintEntryBase
    {
        QualifiedName observable;

        Kinematics kinematics;

        Options options;

        double central, sigma_hi_stat, sigma_lo_stat, sigma_hi_sys, sigma_lo_sys;

        unsigned number_of_observations;

        GaussianConstraintEntry(const std::string & name,
                const QualifiedName & observable,
                const Kinematics & kinematics, const Options & options,
                const double & central,
                const double & sigma_hi_stat, const double & sigma_lo_stat,
                const double & sigma_hi_sys, const double & sigma_lo_sys,
                const unsigned & number_of_observations = 1u) :
            ConstraintEntryBase(name, observable),
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

        virtual ~GaussianConstraintEntry() = default;

        virtual const std::string & type() const
        {
            static const std::string type("Gaussian");

            return type;
        }

        virtual Constraint make(const QualifiedName & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);

            ObservablePtr observable = Observable::make(this->observable, parameters, this->kinematics, this->options + options);
            if (! observable.get())
                throw InternalError("make_gaussian_constraint: " + name.str() + ": '" + this->observable.str() + "' is not a valid observable name");

            double min = 0.0, max = 0.0;
            if ("asymmetric+quadratic" == options.get("uncertainty", "asymmetric+quadratic"))
            {
                min = this->central - std::sqrt(power_of<2>(this->sigma_lo_stat) + power_of<2>(this->sigma_lo_sys));
                max = this->central + std::sqrt(power_of<2>(this->sigma_hi_stat) + power_of<2>(this->sigma_hi_sys));
            }

            LogLikelihoodBlockPtr block = LogLikelihoodBlock::Gaussian(cache, observable, min, this->central, max);

            return Constraint(name, { observable }, { block });
        }

        virtual std::ostream & insert(std::ostream & os) const
        {
            os << _name.full() << ":" << std::endl;
            os << "    type: Gaussian" << std::endl;
            os << "    observable: " << observable << std::endl;

            return os;
        }
    };

    struct LogGammaConstraintEntry :
        public ConstraintEntryBase
    {
        QualifiedName observable;

        Kinematics kinematics;

        Options options;

        double central, sigma_hi, sigma_lo;

        unsigned number_of_observations;

        LogGammaConstraintEntry(const std::string & name,
                const QualifiedName & observable,
                const Kinematics & kinematics, const Options & options,
                const double & central,
                const double & sigma_hi, const double & sigma_lo,
                const unsigned & number_of_observations = 1u) :
            ConstraintEntryBase(name, observable),
            observable(observable),
            kinematics(kinematics),
            options(options),
            central(central),
            sigma_hi(sigma_hi),
            sigma_lo(sigma_lo),
            number_of_observations(number_of_observations)
        {
        }

        virtual ~LogGammaConstraintEntry() = default;

        virtual const std::string & type() const
        {
            static const std::string type("LogGamma");

            return type;
        }

        virtual Constraint make(const QualifiedName & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);

            ObservablePtr observable = Observable::make(this->observable, parameters, this->kinematics, this->options + options);
            if (! observable.get())
                throw InternalError("make_LogGamma_constraint: " + name.str() + ": '" + this->observable.str() + "' is not a valid observable name");

            double min = this->central - this->sigma_lo;
            double max = this->central + this->sigma_hi;

            LogLikelihoodBlockPtr block = LogLikelihoodBlock::LogGamma(cache, observable, min, this->central, max);

            return Constraint(name, { observable }, { block });
        }

        virtual std::ostream & insert(std::ostream & os) const
        {
            os << _name.full() << ":" << std::endl;
            os << "    type: LogGamma" << std::endl;
            os << "    observable: " << observable << std::endl;

            return os;
        }
    };

    struct AmorosoLimitConstraintEntry :
        public ConstraintEntryBase
    {
        QualifiedName observable;

        Kinematics kinematics;

        Options options;

        double physical_limit, upper_limit_90, upper_limit_95;

        double theta, alpha;

        AmorosoLimitConstraintEntry(const QualifiedName & name,
                const QualifiedName & observable,
                const Kinematics & kinematics,
                const Options & options,
                const double & physical_limit,
                const double & upper_limit_90,
                const double & upper_limit_95,
                const double & theta,
                const double & alpha) :
            ConstraintEntryBase(name, observable),
            observable(observable),
            kinematics(kinematics),
            options(options),
            physical_limit(physical_limit),
            upper_limit_90(upper_limit_90),
            upper_limit_95(upper_limit_95),
            theta(theta),
            alpha(alpha)
        {
        }

        virtual ~AmorosoLimitConstraintEntry() = default;

        virtual const std::string & type() const
        {
            static const std::string type("Amoroso (physical limit, 90% upper limit, 95% upper limit)");

            return type;
        }

        virtual Constraint make(const QualifiedName & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);

            ObservablePtr observable = Observable::make(this->observable, parameters, this->kinematics, this->options + options);
            if (! observable.get())
                throw InternalError("make_amoroso_limit_constraint: " + name.str() + ": '" + this->observable.str() + "' is not a valid observable name");


            LogLikelihoodBlockPtr block = LogLikelihoodBlock::AmorosoLimit(cache, observable, physical_limit, upper_limit_90,
                                                                           upper_limit_95, theta, alpha);

            return Constraint(name, { observable }, { block });
        }

        virtual std::ostream & insert(std::ostream & os) const
        {
            os << _name.full() << ":" << std::endl;
            os << "    type: AmorosoLimit" << std::endl;
            os << "    observable: " << observable << std::endl;

            return os;
        }
    };

    struct AmorosoModeConstraintEntry :
        public ConstraintEntryBase
    {
        QualifiedName observable;

        Kinematics kinematics;

        Options options;

        double physical_limit, mode, upper_limit_90, upper_limit_95;

        double theta, alpha, beta;

        AmorosoModeConstraintEntry(const QualifiedName & name,
                const QualifiedName & observable,
                const Kinematics & kinematics,
                const Options & options,
                const double & physical_limit,
                const double & mode,
                const double & upper_limit_90,
                const double & upper_limit_95,
                const double & theta,
                const double & alpha,
                const double & beta) :
            ConstraintEntryBase(name, observable),
            observable(observable),
            kinematics(kinematics),
            options(options),
            physical_limit(physical_limit),
            mode(mode),
            upper_limit_90(upper_limit_90),
            upper_limit_95(upper_limit_95),
            theta(theta),
            alpha(alpha),
            beta(beta)
        {
        }

        virtual ~AmorosoModeConstraintEntry() = default;

        virtual const std::string & type() const
        {
            static const std::string type("Amoroso (physical limit, mode, 90% upper limit, 95% upper limit)");

            return type;
        }

        virtual Constraint make(const QualifiedName & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);

            ObservablePtr observable = Observable::make(this->observable, parameters, this->kinematics, this->options + options);
            if (! observable.get())
                throw InternalError("make_amoroso_constraint: " + name.str() + ": '" + this->observable.str() + "' is not a valid observable name");


            LogLikelihoodBlockPtr block = LogLikelihoodBlock::AmorosoMode(cache, observable, physical_limit,
                                                                      mode, upper_limit_90, upper_limit_95,
                                                                      theta, alpha, beta);

            return Constraint(name, { observable }, { block });
        }

        virtual std::ostream & insert(std::ostream & os) const
        {
            os << _name.full() << ":" << std::endl;
            os << "    type: AmorosoMode" << std::endl;
            os << "    observable: " << observable << std::endl;

            return os;
        }
    };

    struct AmorosoTripleLimitConstraintEntry :
        public ConstraintEntryBase
    {
        QualifiedName observable;

        Kinematics kinematics;

        Options options;

        double physical_limit, upper_limit_10, upper_limit_50, upper_limit_90;

        double theta, alpha, beta;

        AmorosoTripleLimitConstraintEntry(const QualifiedName & name,
                const QualifiedName & observable,
                const Kinematics & kinematics,
                const Options & options,
                const double & physical_limit,
                const double & upper_limit_10,
                const double & upper_limit_50,
                const double & upper_limit_90,
                const double & theta,
                const double & alpha,
                const double & beta) :
            ConstraintEntryBase(name, observable),
            observable(observable),
            kinematics(kinematics),
            options(options),
            physical_limit(physical_limit),
            upper_limit_10(upper_limit_10),
            upper_limit_50(upper_limit_50),
            upper_limit_90(upper_limit_90),
            theta(theta),
            alpha(alpha),
            beta(beta)
        {
        }

        virtual ~AmorosoTripleLimitConstraintEntry() = default;

        virtual const std::string & type() const
        {
            static const std::string type("Amoroso (physical limit, 10% upper limit, 90% upper limit, 95% upper limit)");

            return type;
        }

        virtual Constraint make(const QualifiedName & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);

            ObservablePtr observable = Observable::make(this->observable, parameters, this->kinematics, this->options + options);
            if (! observable.get())
                throw InternalError("make_amoroso_constraint: " + name.str() + ": '" + this->observable.str() + "' is not a valid observable name");


            LogLikelihoodBlockPtr block = LogLikelihoodBlock::Amoroso(cache, observable, physical_limit,
                                                                      upper_limit_10, upper_limit_50, upper_limit_90,
                                                                      theta, alpha, beta);

            return Constraint(name, { observable }, { block });
        }

        virtual std::ostream & insert(std::ostream & os) const
        {
            os << _name.full() << ":" << std::endl;
            os << "    type: AmorosoTripleLimit" << std::endl;
            os << "    observable: " << observable << std::endl;

            return os;
        }
    };

    struct AmorosoConstraintEntry :
        public ConstraintEntryBase
    {
        QualifiedName observable;

        Kinematics kinematics;

        Options options;

        double physical_limit, theta, alpha, beta;

        AmorosoConstraintEntry(const QualifiedName & name,
                const QualifiedName & observable,
                const Kinematics & kinematics,
                const Options & options,
                const double & physical_limit,
                const double & theta,
                const double & alpha,
                const double & beta) :
            ConstraintEntryBase(name, observable),
            observable(observable),
            kinematics(kinematics),
            options(options),
            physical_limit(physical_limit),
            theta(theta),
            alpha(alpha),
            beta(beta)
        {
        }

        virtual ~AmorosoConstraintEntry() = default;

        virtual const std::string & type() const
        {
            static const std::string type("Amoroso (physical limit, theta, alpha, beta)");

            return type;
        }

        virtual Constraint make(const QualifiedName & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);

            ObservablePtr observable = Observable::make(this->observable, parameters, this->kinematics, this->options + options);
            if (! observable.get())
                throw InternalError("make_amoroso_constraint: " + name.str() + ": '" + this->observable.str() + "' is not a valid observable name");

            LogLikelihoodBlockPtr block = LogLikelihoodBlock::Amoroso(cache, observable, physical_limit, theta, alpha, beta);

            return Constraint(name, { observable }, { block });
        }

        virtual std::ostream & insert(std::ostream & os) const
        {
            os << _name.full() << ":" << std::endl;
            os << "    type: Amoroso" << std::endl;
            os << "    observable: " << observable << std::endl;

            return os;
        }
    };

    template <size_t dim_>
    struct MultivariateGaussianConstraintEntry :
        public ConstraintEntryBase
    {
        std::array<QualifiedName, dim_> observable_names;

        std::array<Kinematics, dim_> kinematics;

        std::array<Options, dim_> options;

        std::array<double, dim_> means;

        std::array<double, dim_> sigma_stat_hi;
        std::array<double, dim_> sigma_stat_lo;

        std::array<double, dim_> sigma_sys;

        std::array<std::array<double, dim_>, dim_> correlation;

        unsigned number_of_observations;

        MultivariateGaussianConstraintEntry(const QualifiedName & name,
                const std::array<QualifiedName, dim_> & observable_names,
                const std::array<Kinematics, dim_> & kinematics,
                const std::array<Options, dim_> & options,
                const std::array<double, dim_> & means,
                const std::array<double, dim_> & sigma_stat_hi,
                const std::array<double, dim_> & sigma_stat_lo,
                const std::array<double, dim_> & sigma_sys,
                const std::array<std::array<double, dim_>, dim_> & correlation,
                const unsigned number_of_observations = dim_) :
            ConstraintEntryBase(name, observable_names),
            observable_names(observable_names),
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

        virtual ~MultivariateGaussianConstraintEntry() = default;

        virtual const std::string & type() const
        {
            static const std::string type("MultivariateGaussian<" + stringify(dim_) + "> (using correlation matrix)");

            return type;
        }

        virtual Constraint make(const QualifiedName & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);

            std::array<ObservablePtr, dim_> observables;
            for (auto i = 0u ; i < dim_ ; ++i)
            {
                observables[i] = Observable::make(this->observable_names[i], parameters, this->kinematics[i], this->options[i] + options);
                if (! observables[i].get())
                    throw InternalError("make_multivariate_gaussian_constraint<" + stringify(dim_) + ">: " + name.str() + ": '" + this->observable_names[i].str() + "' is not a valid observable name");
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

            LogLikelihoodBlockPtr block = LogLikelihoodBlock::MultivariateGaussian(cache, observables, this->means, variances, this->correlation, number_of_observations);

            return Constraint(name, std::vector<ObservablePtr>(observables.begin(), observables.end()), { block });
        }

        virtual std::ostream & insert(std::ostream & os) const
        {
            os << _name.full() << ":" << std::endl;
            os << "    type: MultivariateGaussian<" << dim_ << ">" << std::endl;

            return os;
        }
    };

    template <size_t dim_>
    struct MultivariateGaussianCovarianceConstraintEntry :
        public ConstraintEntryBase
    {
        std::array<QualifiedName, dim_> observables;

        std::array<Kinematics, dim_> kinematics;

        std::array<Options, dim_> options;

        std::array<double, dim_> means;

        std::array<std::array<double, dim_>, dim_> covariance;

        unsigned number_of_observations;

        MultivariateGaussianCovarianceConstraintEntry(const QualifiedName & name,
                const std::array<QualifiedName, dim_> & observables,
                const std::array<Kinematics, dim_> & kinematics,
                const std::array<Options, dim_> & options,
                const std::array<double, dim_> & means,
                const std::array<std::array<double, dim_>, dim_> & covariance,
                const unsigned number_of_observations = dim_) :
            ConstraintEntryBase(name, observables),
            observables(observables),
            kinematics(kinematics),
            options(options),
            means(means),
            covariance(covariance),
            number_of_observations(number_of_observations)
        {
        }

        virtual ~MultivariateGaussianCovarianceConstraintEntry() = default;

        virtual const std::string & type() const
        {
            static const std::string type("MultivariateGaussian<" + stringify(dim_) + "> (using covariance matrix)");

            return type;
        }

        virtual Constraint make(const QualifiedName & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);

            std::array<ObservablePtr, dim_> observables;
            for (auto i = 0u ; i < dim_ ; ++i)
            {
                observables[i] = Observable::make(this->observables[i], parameters, this->kinematics[i], this->options[i] + options);
                if (! observables[i].get())
                    throw InternalError("make_multivariate_gaussian_covariance_constraint<" + stringify(dim_) + ">: " + name.str() + ": '" + this->observables[i].str() + "' is not a valid observable name");
            }

            auto block = LogLikelihoodBlock::MultivariateGaussian(cache, observables, means, covariance, number_of_observations);

            return Constraint(name, std::vector<ObservablePtr>(observables.begin(), observables.end()), { block });
        }

        virtual std::ostream & insert(std::ostream & os) const
        {
            os << _name.full() << ":" << std::endl;
            os << "    type: MultivariateGaussianCovariance<" << dim_ << ">" << std::endl;

            return os;
        }
    };

    namespace entries
    {
        ///@name 2000 Data
        ///@{
        /*
         * CLEO Collaboration
         *
         * Data taken from [CLEO:2000]
         */
        static const GaussianConstraintEntry Bzero_to_Kstarzero_gamma_BR_CLEO_2000
        {
            "B^0->K^*0gamma::BR@CLEO-2000",
            "B->K^*gamma::BRavg",
            Kinematics{ },
            Options{ { "q", "d" } },
            4.55e-5, 0.72e-5, 0.68e-5, 0.34e-5, 0.34e-5
        };
        static const GaussianConstraintEntry Bplus_to_Kstarplus_gamma_BR_CLEO_2000
        {
            "B^+->K^*+gamma::BR@CLEO-2000",
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
        static const GaussianConstraintEntry Bmix_to_Xs_dilepton_BR_BaBar_2004A
        {
            "B->X_sll::BR[1.0,6.0]@BaBar-2004A",
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
        static const GaussianConstraintEntry Bzero_to_Kstarzero_gamma_BR_Belle_2004
        {
            "B^0->K^*0gamma::BR@Belle-2004",
            "B->K^*gamma::BRavg",
            Kinematics{ },
            Options{ { "q", "d" } },
            4.01e-5, 0.21e-5, 0.21e-5, 0.17e-5, 0.17e-5
        };
        static const GaussianConstraintEntry Bplus_to_Kstarplus_gamma_BR_Belle_2004
        {
            "B^+->K^*+gamma::BR@Belle-2004",
            "B->K^*gamma::BRavg",
            Kinematics{ },
            Options{ { "q", "u" } },
            4.25e-5, 0.31e-5, 0.31e-5, 0.24e-5, 0.24e-5
        };
        ///@}

        ///@name 2005 Data
        ///@{
        static const GaussianConstraintEntry Bmix_to_Xs_dilepton_BR_Belle_2005A
        {
            "B->X_sll::BR[1.0,6.0]@Belle-2005A",
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
        static const GaussianConstraintEntry Bzero_to_Kstarzero_gamma_SKstargamma_Belle_2006
        {
            "B^0->K^*0gamma::S_K@Belle-2006",
            "B->K^*gamma::S_K^*gamma",
            Kinematics{ },
            Options{ { "q", "d" } },
            -0.32, +0.36, -0.33, +0.05, -0.05
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_gamma_CKstargamma_Belle_2006
        {
            "B^0->K^*0gamma::C_K@Belle-2006",
            "B->K^*gamma::C_K^*gamma",
            Kinematics{ },
            Options{ { "q", "d" } },
            +0.20, +0.24, -0.24, +0.05, -0.05
        };
        static const MultivariateGaussianConstraintEntry<2> Bzero_to_Kstarzero_gamma_time_dependent_cp_asymmetries_Belle_2006
        {
            "B^0->K^*0gamma::S_K+C_K@Belle-2006",
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
        static const GaussianConstraintEntry Bzero_to_Kstarzero_gamma_SKstargamma_BaBar_2008
        {
            "B^0->K^*0gamma::S_K@BaBar-2008",
            "B->K^*gamma::S_K^*gamma",
            Kinematics{ },
            Options{ { "q", "d" } },
            -0.03, +0.29, -0.29, +0.03, -0.03
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_gamma_CKstargamma_BaBar_2008
        {
            "B^0->K^*0gamma::C_K@BaBar-2008",
            "B->K^*gamma::C_K^*gamma",
            Kinematics{ },
            Options{ { "q", "d" } },
            -0.14, +0.16, -0.16, +0.03, -0.03
        };
        static const MultivariateGaussianConstraintEntry<2> Bzero_to_Kstarzero_gamma_time_dependent_cp_asymmetries_BaBar_2008
        {
            "B^0->K^*0gamma::S_K+C_K@BaBar-2008",
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
        static const GaussianConstraintEntry B_to_Xs_gamma_E_1_1dot8_Belle_2008A
        {
            "B->X_sgamma::E_1[1.8]@Belle-2008",
            "B->X_sgamma::E_1(E_min)@NLO",
            Kinematics{ { "E_min", 1.8 } },
            Options{ },
            2.292, 0.027, 0.027, 0.033, 0.033
        };
        static const GaussianConstraintEntry B_to_Xs_gamma_E_2_1dot8_Belle_2008A
        {
            "B->X_sgamma::E_2[1.8]@Belle-2008",
            "B->X_sgamma::E_2(E_min)@NLO",
            Kinematics{ { "E_min", 1.8 } },
            Options{ },
            0.0305, 0.079, 0.079, 0.099, 0.099
        };
        static const MultivariateGaussianConstraintEntry<2> B_to_Xs_gamma_E_1_and_E_2_1dot8_Belle_2008A
        {
            "B->X_sgamma::E_1[1.8]+E_2[1.8]@Belle-2008",
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
        static const GaussianConstraintEntry Bzero_to_Kstarzero_gamma_BR_BaBar_2009
        {
            "B^0->K^*0gamma::BR@BaBar-2009",
            "B->K^*gamma::BRavg",
            Kinematics{ },
            Options{ { "q", "d" } },
            4.47e-5, +0.10e-5, -0.10e-5, +0.16e-5, -0.16e-5
        };
        static const GaussianConstraintEntry Bplus_to_Kstarplus_gamma_BR_BaBar_2009
        {
            "B^+->K^*+gamma::BR@BaBar-2009",
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
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_1_to_6_Belle_2009
        {
            "B^+->K^+mu^+mu^-::BR[1.00,6.00]@Belle-2009",
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            1.36e-7, +0.23e-7, -0.21e-7, +0.08e-7, -0.08e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_14dot18_to_16_Belle_2009
        {
            "B^+->K^+mu^+mu^-::BR[14.18,16.00]@Belle-2009",
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.38e-7, +0.19e-7, -0.12e-7, +0.02e-7, -0.02e-7
        };
        // BR in [16.00, 22.86]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_16_to_22dot86_Belle_2009
        {
            "B^+->K^+mu^+mu^-::BR[16.00,22.86]@Belle-2009",
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 22.86 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.98e-7, +0.20e-7, -0.19e-7, +0.06e-7, -0.06e-7
        };

        // A_FB in [1.0, 6.0]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_A_FB_1_to_6_Belle_2009
        {
            "B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@Belle-2009",
            "B->Kll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            +0.04, +0.16, -0.13, +0.05, -0.05
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_Belle_2009
        {
            "B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@Belle-2009",
            "B->Kll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            -0.04, +0.26, -0.32, +0.05, -0.05
        };
        // A_FB in [16.00, 22.86]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_A_FB_16_to_22dot86_Belle_2009
        {
            "B^+->K^+mu^+mu^-::A_FB[16.00,22.86]@Belle-2009",
            "B->Kll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 22.86 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            -0.02, +0.08, -0.11, +0.02, -0.02
        };

        // B^0 -> K^*0 mu^+ mu^-
        // BR in [1.0, 6.0]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_1_to_6_Belle_2009
        {
            "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@Belle-2009",
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.49e-7, +0.45e-7, -0.40e-7, +0.12e-7, -0.12e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_Belle_2009
        {
            "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@Belle-2009",
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.05e-7, +0.29e-7, -0.26e-7, +0.08e-7, -0.08e-7
        };
        // BR in [16.00, 19.21]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_Belle_2009
        {
            "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@Belle-2009",
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            2.04e-7, +0.27e-7, -0.24e-7, +0.16e-7, -0.16e-7
        };
        // A_FB in [1.0, 6.0]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_Belle_2009
        {
            "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@Belle-2009",
            "B->K^*ll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.26, +0.30, -0.27, +0.07, -0.07
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_Belle_2009
        {
            "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@Belle-2009",
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.70, +0.22, -0.16, +0.10, -0.10
        };
        // A_FB in [16.00, 19.21]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_Belle_2009
        {
            "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@Belle-2009",
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.66, +0.16, -0.11, +0.04, -0.04
        };
        // F_L in [1.0, 6.0]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_1_to_6_Belle_2009
        {
            "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@Belle-2009",
            "B->K^*ll::F_Lavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.67, +0.23, -0.23, +0.05, -0.05
        };
        // F_L in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_Belle_2009
        {
            "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@Belle-2009",
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.15, +0.27, -0.23, +0.07, -0.07
        };
        // F_L in [16.00, 19.21]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_Belle_2009
        {
            "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@Belle-2009",
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
        static const GaussianConstraintEntry B_to_Xs_gamma_BR_1dot8_Belle_2009B
        {
            "B->X_sgamma::BR[1.8]@Belle-2009B",
            "B->X_sgamma::BR(E_min)@NLO",
            Kinematics{ { "E_min", 1.8 } },
            Options{ },
            3.36e-4, +0.13e-4, -0.13e-4, +0.25e-4, -0.25e-4
        };
        static const MultivariateGaussianConstraintEntry<3> B_to_Xs_gamma_1dot8_Belle_2009B
        {
            "B->X_sgamma::BR[1.8]+E_1[1.8]+E_2[1.8]@Belle-2009B",
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
        static const GaussianConstraintEntry Bzero_to_pi_l_nu_BR_0_to_4_BaBar_2010A
        {
            "B^0->pi^+lnu::BR[0.0,4.0]@BaBar-2010A",
            "B->pilnu::BR",
            Kinematics{ { "s_min", 0.0 }, { "s_max", 4.0 } },
            Options{ { "q", "d" } },
            0.313e-4, +0.030e-4, -0.030e-4, +0.025e-4, -0.025e-4
        };
        static const GaussianConstraintEntry Bzero_to_pi_l_nu_BR_4_to_8_BaBar_2010A
        {
            "B^0->pi^+lnu::BR[4.0,8.0]@BaBar-2010A",
            "B->pilnu::BR",
            Kinematics{ { "s_min", 4.0 }, { "s_max", 8.0 } },
            Options{ { "q", "d" } },
            0.329e-4, +0.018e-4, -0.018e-4, +0.016e-4, -0.016e-4
        };
        static const GaussianConstraintEntry Bzero_to_pi_l_nu_BR_8_to_12_BaBar_2010A
        {
            "B^0->pi^+lnu::BR[8.0,12.0]@BaBar-2010A",
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
        static const MultivariateGaussianConstraintEntry<6> Bzero_to_pi_l_nu_BR_BaBar_2010B
        {
            "B^0->pi^+lnu::BR@BaBar-2010B",
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
        static const MultivariateGaussianConstraintEntry<6> Bzero_to_pi_l_nu_BR_Belle_2010A
        {
            "B^0->pi^+lnu::BR@Belle-2010A",
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
        static const MultivariateGaussianConstraintEntry<2> Bzero_to_Kstarzero_gamma_time_dependent_cp_asymmetries_HFAG_2011
        {"B^0->K^*0gamma::S_K+C_K@HFAG-2011",
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
        static const GaussianConstraintEntry Bzero_to_Kzero_dimuon_BR_1_to_6_CDF_2011
        {
            "B^0->K^0mu^+mu^-::BR[1.00,6.00]@CDF-2011",
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.98e-7, +0.614e-7, -0.614e-7, +0.074e-7, -0.074e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kzero_dimuon_BR_14dot18_to_16_CDF_2011
        {
            "B^0->K^0mu^+mu^-::BR[14.18,16.00]@CDF-2011",
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.726e-7, +0.257e-7, -0.257e-7, +0.054e-7, -0.054e-7
        };
        // BR in [16.00, 22.86]
        static const GaussianConstraintEntry Bzero_to_Kzero_dimuon_BR_16_to_22dot86_CDF_2011
        {
            "B^0->K^0mu^+mu^-::BR[16.00,22.86]@CDF-2011",
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 22.86 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.214e-7, +0.182e-7, -0.182e-7, +0.016e-7, -0.016e-7
        };

        // B^+ -> K^+ mu^+ mu^-
        // BR in [1.0, 6.0]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_1_to_6_CDF_2011
        {
            "B^+->K^+mu^+mu^-::BR[1.00,6.00]@CDF-2011",
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            1.41e-7, +0.20e-7, -0.20e-7, +0.09e-7, -0.09e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_14dot18_to_16_CDF_2011
        {
            "B^+->K^+mu^+mu^-::BR[14.18,16.00]@CDF-2011",
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.53e-7, +0.10e-7, -0.19e-7, +0.03e-7, -0.03e-7
        };
        // BR in [16.00, 22.86]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_16_to_22dot86_CDF_2011
        {
            "B^+->K^+mu^+mu^-::BR[16.00,22.86]@CDF-2011",
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 22.86 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.48e-7, +0.11e-7, -0.11e-7, +0.03e-7, -0.03e-7
        };

        // B^0 -> K^*0 mu^+ mu^-
        // BR in [1.0, 6.0]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_1_to_6_CDF_2011
        {
            "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@CDF-2011",
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.42e-7, +0.41e-7, -0.4100e-7, +0.08e-7, -0.08e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_CDF_2011
        {
            "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@CDF-2011",
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.34e-7, +0.26e-7, -0.26e-7, +0.08e-7, -0.08e-7
        };
        // BR in [16.00, 19.21]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_CDF_2011
        {
            "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@CDF-2011",
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.97e-7, +0.26e-7, -0.26e-7, +0.06e-7, -0.06e-7
        };

        // B^+ -> K^*+ mu^+ mu^-
        // BR in [1.0, 6.0]
        static const GaussianConstraintEntry Bplus_to_Kstarplus_dimuon_BR_1_to_6_CDF_2011
        {
            "B^+->K^*+mu^+mu^-::BR[1.00,6.00]@CDF-2011",
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            2.57e-7, +1.61e-7, -1.61e-7, +0.22e-7, -0.22e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintEntry Bplus_to_Kstarplus_dimuon_BR_14dot18_to_16_CDF_2011
        {
            "B^+->K^*+mu^+mu^-::BR[14.18,16.00]@CDF-2011",
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.52e-7, +0.61e-7, -0.61e-7, +0.04e-7, -0.04e-7
        };
        // BR in [16.00, 19.21]
        static const GaussianConstraintEntry Bplus_to_Kstarplus_dimuon_BR_16_to_19dot21_CDF_2011
        {
            "B^+->K^*+mu^+mu^-::BR[16.00,19.21]@CDF-2011",
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
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_CDF_2011
        {
            "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@CDF-2011",
            "B->K^*ll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.36, +0.28, -0.46, +0.11, -0.11
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_CDF_2011
        {
            "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@CDF-2011",
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.40, +0.21, -0.18, +0.07, -0.07
        };
        // A_FB in [16.00, 19.21]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_CDF_2011
        {
            "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@CDF-2011",
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.66, +0.26, -0.18, +0.19, -0.19
        };
        // F_L in [1.0, 6.0]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_1_to_6_CDF_2011
        {
            "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@CDF-2011",
            "B->K^*ll::F_Lavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.60, +0.21, -0.23, +0.09, -0.09
        };
        // F_L in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_CDF_2011
        {
            "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@CDF-2011",
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.32, +0.14, -0.14, +0.03, -0.03
        };
        // F_L in [16.00, 19.21]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_CDF_2011
        {
            "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@CDF-2011",
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.16, +0.22, -0.18, +0.06, -0.06
        };
        // A_T^{2} in [1.0, 6.0]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_T_2_1_to_6_CDF_2011
        {
            "B^0->K^*0mu^+mu^-::A_T^2[1.00,6.00]@CDF-2011",
            "B->K^*ll::A_T^2avg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +1.6, +1.8, -1.9, +2.2, -2.2
        };
        // A_T^{2} in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_T_2_14dot18_to_16_CDF_2011
        {
            "B^0->K^*0mu^+mu^-::A_T^2[14.18,16.00]@CDF-2011",
            "B->K^*ll::A_T^2avg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.4, +0.8, -0.8, +0.2, -0.2
        };
        // A_T^{2} in [16.00, 19.21]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_T_2_16_to_19dot21_CDF_2011
        {
            "B^0->K^*0mu^+mu^-::A_T^2[16.00,19.21]@CDF-2011",
            "B->K^*ll::A_T^2avg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.9, +0.8, -0.8, +0.4, -0.4
        };
        // A_{im} in [1.0, 6.0]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_im_1_to_6_CDF_2011
        {
            "B^0->K^*0mu^+mu^-::A_im[1.00,6.00]@CDF-2011",
            "B->K^*ll::J_9normavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.02, +0.40, -0.40, +0.03, -0.03
        };
        // A_{im} in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_im_14dot18_to_16_CDF_2011
        {
            "B^0->K^*0mu^+mu^-::A_im[14.18,16.00]@CDF-2011",
            "B->K^*ll::J_9normavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.15, +0.25, -0.26, +0.01, -0.01
        };
        // A_{im} in [16.00, 19.21]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_im_16_to_19dot21_CDF_2011
        {
            "B^0->K^*0mu^+mu^-::A_im[16.00,19.21]@CDF-2011",
            "B->K^*ll::J_9normavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.30, +0.36, -0.35, +0.14, -0.14
        };

        // B^+ -> K^+ mu^+ mu^-
        // A_FB in [1.0, 6.0]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_A_FB_1_to_6_CDF_2011
        {
            "B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@CDF-2011",
            "B->Kll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            -0.13, +0.09, -0.09, +0.02, -0.02
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_CDF_2011
        {
            "B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@CDF-2011",
            "B->Kll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            +0.05, +0.11, -0.09, +0.03, -0.03
        };
        // A_FB in [16.00, 23]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_A_FB_16_to_22dot86_CDF_2011
        {
            "B^+->K^+mu^+mu^-::A_FB[16.00,22.86]@CDF-2011",
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
        static const AmorosoLimitConstraintEntry Bzero_to_dimuon_CDF_2011
        {
            "B^0_s->mu^+mu^-::BR_limit@CDF-2011",
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
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_1_to_6_LHCb_2011
        {
            "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@LHCb-2011",
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            2.1000e-7, +0.3000e-7, -0.3000e-7, +0.1500e-7, -0.1500e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_LHCb_2011
        {
            "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@LHCb-2011",
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.1466e-7, +0.2002e-7, -0.2002e-7, +0.0910e-7, -0.0910e-7
        };
        // BR in [16.00, 19.00] (in the preliminary results, LHCb only integrated up to exactly 19.00 GeV^2!)
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_LHCb_2011
        {
            "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@LHCb-2011",
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.5000e-7, +0.2400e-7, -0.2400e-7, +0.1500e-7, -0.1500e-7
        };
        // A_FB in [1.0, 6.0]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_LHCb_2011
        {
            "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@LHCb-2011",
            "B->K^*ll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.06, +0.14, -0.13, +0.04, -0.04
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_LHCb_2011
        {
            "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@LHCb-2011",
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.47, +0.08, -0.06, +0.03, -0.03
        };
        // A_FB in [16.00, 19.00] (in the preliminary results, LHCb only integrated up to exactly 19.00 GeV^2!)
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_LHCb_2011
        {
            "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@LHCb-2011",
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.16, +0.13, -0.11, +0.06, -0.06
        };
        // F_L in [1.0, 6.0]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_1_to_6_LHCb_2011
        {
            "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@LHCb-2011",
            "B->K^*ll::F_Lavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.55, +0.10, -0.10, +0.03, -0.03
        };
        // F_L in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_LHCb_2011
        {
            "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@LHCb-2011",
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.37, +0.09, -0.09, +0.05, -0.05
        };
        // F_L in [16.00, 19.00] (in the preliminary results, LHCb only integrated up to exactly 19.00 GeV^2!)
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_LHCb_2011
        {
            "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@LHCb-2011",
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
        static const AmorosoLimitConstraintEntry Bzero_to_dimuon_LHCb_CMS_2011
        {
            "B^0_s->mu^+mu^-::BR_limit@LHCb-CMS-2011",
            "B_q->ll::BR@Untagged",
            Kinematics{ },
            Options{ { "q", "s"  }, { "l", "mu" } },
            0.0, 0.9e-8, 1.08e-8,
            0.74377978e-8, 0.53538044
        };
        // use the data from the Bayes-Heinrich method
        // the mode is not at zero, but around 3.1e-9
        static const AmorosoTripleLimitConstraintEntry Bzero_to_dimuon_LHCb_CMS_2011_Bayes
        {
            "B^0_s->mu^+mu^-::BR_limit@LHCb-CMS-2011-Bayes",
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
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_1_to_6_BaBar_2012
        {
            "B^+->K^+mu^+mu^-::BR[1.00,6.00]@BaBar-2012",
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            1.36e-7, +0.27e-7, -0.24e-7, +0.03e-7, -0.03e-7
        };
        // BR in [14.21, 16.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_14dot21_to_16_BaBar_2012
        {
            "B^+->K^+mu^+mu^-::BR[14.21,16.00]@BaBar-2012",
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.21 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.49e-7, +0.15e-7, -0.14e-7, +0.02e-7, -0.02e-7
        };
        // BR in [16.00, 22.86]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_16_to_22dot86_BaBar_2012
        {
            "B^+->K^+mu^+mu^-::BR[16.00,22.86]@BaBar-2012",
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 22.86 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.67e-7, +0.23e-7, -0.21e-7, +0.05e-7, -0.05e-7
        };

        // B^0 -> K^*0 mu^+ mu^-

        // BR in [1.0, 6.0]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_1_to_6_BaBar_2012
        {
            "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@BaBar-2012",
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            2.05e-7, +0.53e-7, -0.48e-7, +0.07e-7, -0.07e-7
        };
        // BR in [14.21, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_14dot21_to_16_BaBar_2012
        {
            "B^0->K^*0mu^+mu^-::BR[14.21,16.00]@BaBar-2012",
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.21 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.46e-7, +0.41e-7, -0.36e-7, +0.06e-7, -0.06e-7
        };
        // BR in [16.00, 22.86]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_BaBar_2012
        {
            "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@BaBar-2012",
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
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_BaBar_2012
        {
            "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@BaBar-2012",
            "B->K^*ll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.17, +0.14, -0.12, +0.07, -0.07
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_BaBar_2012
        {
            "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@BaBar-2012",
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.34, +0.15, -0.08, +0.07, -0.07
        };
        // A_FB in [16.00, 19.21]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_BaBar_2012
        {
            "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@BaBar-2012",
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.34, +0.21, -0.19, +0.07, -0.07
        };
        // F_L in [1.0, 6.0]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_1_to_6_BaBar_2012
        {
            "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@BaBar-2012",
            "B->K^*ll::F_Lavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.25, +0.09, -0.08, +0.03, -0.03
        };
        // F_L in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_BaBar_2012
        {
            "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@BaBar-2012",
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.43, +0.10, -0.13, +0.09, -0.09
        };
        // F_L in [16.00, 19.21]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_BaBar_2012
        {
            "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@BaBar-2012",
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
        static const GaussianConstraintEntry B_to_Xs_gamma_BR_1dot8_BaBar_2012C
        {
            "B->X_sgamma::BR[1.8]@BaBar-2012",
            "B->X_sgamma::BR(E_min)@NLO",
            Kinematics{ { "E_min", 1.8 } },
            Options{ { "q", "d" } },
            +3.21e-4, +0.15e-4, -0.15e-4, +0.29e-4, -0.29e-4
        };
        static const GaussianConstraintEntry B_to_Xs_gamma_E_1_1dot8_BaBar_2012C
        {
            "B->X_sgamma::E_1[1.8]@BaBar-2012",
            "B->X_sgamma::E_1(E_min)@NLO",
            Kinematics{ { "E_min", 1.8 } },
            Options{ { "q", "d" } },
            +2.267, +0.019, -0.019, +0.032, -0.032
        };
        static const GaussianConstraintEntry B_to_Xs_gamma_E_2_1dot8_BaBar_2012C
        {
            "B->X_sgamma::E_2[1.8]@BaBar-2012",
            "B->X_sgamma::E_2(E_min)@NLO",
            Kinematics{ { "E_min", 1.8 } },
            Options{ { "q", "d"  } },
            +4.84e-2, +5.3e-3, -5.3e-3, +7.7e-3, -7.7e-3
        };
        static const MultivariateGaussianConstraintEntry<2> B_to_Xs_gamma_E_1_and_E_2_1dot8_BaBar_2012C
        {
            "B->X_sgamma::E_1[1.8]+E_2[1.8]@BaBar-2012",
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
        static const MultivariateGaussianConstraintEntry<6> Bzero_to_pi_l_nu_BR_BaBar_2012D
        {
            "B^0->pi^+lnu::BR@BaBar-2012D",
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
        static const GaussianConstraintEntry Bzero_to_Kzero_dimuon_BR_1_to_6_CDF_2012
        {
            "B^0->K^0mu^+mu^-::BR[1.00,6.00]@CDF-2012",
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.11e-7, +0.52e-7, -0.52e-7, +0.14e-7, -0.14e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kzero_dimuon_BR_14dot18_to_16_CDF_2012
        {
            "B^0->K^0mu^+mu^-::BR[14.18,16.00]@CDF-2012",
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.43e-7, +0.18e-7, -0.18e-7, +0.04e-7, -0.04e-7
        };
        // BR in [16.00, 22.86]
        static const GaussianConstraintEntry Bzero_to_Kzero_dimuon_BR_16_to_22dot86_CDF_2012
        {
            "B^0->K^0mu^+mu^-::BR[16.00,22.86]@CDF-2012",
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 22.86 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.26e-7, +0.15e-7, -0.15e-7, +0.03e-7, -0.03e-7
        };

        // B^+ -> K^+ mu^+ mu^-
        // BR in [1.0, 6.0]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_1_to_6_CDF_2012
        {
            "B^+->K^+mu^+mu^-::BR[1.00,6.00]@CDF-2012",
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            1.35e-7, +0.18e-7, -0.18e-7, +0.08e-7, -0.08e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_14dot18_to_16_CDF_2012
        {
            "B^+->K^+mu^+mu^-::BR[14.18,16.00]@CDF-2012",
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.42e-7, +0.08e-7, -0.08e-7, +0.02e-7, -0.02e-7
        };
        // BR in [16.00, 22.86]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_16_to_22dot86_CDF_2012
        {
            "B^+->K^+mu^+mu^-::BR[16.00,22.86]@CDF-2012",
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 22.86 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.50e-7, +0.10e-7, -0.10e-7, +0.02e-7, -0.02e-7
        };

        // B^0 -> K^*0 mu^+ mu^-
        // BR in [1.0, 6.0]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_1_to_6_CDF_2012
        {
            "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@CDF-2012",
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.97e-7, +0.49e-7, -0.49e-7, +0.14e-7, -0.14e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_CDF_2012
        {
            "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@CDF-2012",
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            1.35e-7, +0.21e-7, -0.21e-7, +0.08e-7, -0.08e-7
        };
        // BR in [16.00, 19.21]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_CDF_2012
        {
            "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@CDF-2012",
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.98e-7, +0.22e-7, -0.22e-7, +0.06e-7, -0.06e-7
        };

        // B^+ -> K^*+ mu^+ mu^-
        // BR in [1.0, 6.0]
        static const GaussianConstraintEntry Bplus_to_Kstarplus_dimuon_BR_1_to_6_CDF_2012
        {
            "B^+->K^*+mu^+mu^-::BR[1.00,6.00]@CDF-2012",
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            3.56e-7, +1.38e-7, -1.38e-7, +0.43e-7, -0.43e-7
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintEntry Bplus_to_Kstarplus_dimuon_BR_14dot18_to_16_CDF_2012
        {
            "B^+->K^*+mu^+mu^-::BR[14.18,16.00]@CDF-2012",
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            1.02e-7, +0.58e-7, -0.58e-7, +0.13e-7, -0.13e-7
        };
        // BR in [16.00, 19.21]
        static const GaussianConstraintEntry Bplus_to_Kstarplus_dimuon_BR_16_to_19dot21_CDF_2012
        {
            "B^+->K^*+mu^+mu^-::BR[16.00,19.21]@CDF-2012",
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            1.68e-7, +0.74e-7, -0.74e-7, +0.19e-7, -0.19e-7
        };

        // B^0 -> K^*0 mu^+ mu^-
        // A_FB in [1.0, 6.0]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_CDF_2012
        {
            "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@CDF-2012",
            "B->K^*ll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.29, +0.21, -0.25, +0.06, -0.06
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_CDF_2012
        {
            "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@CDF-2012",
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.49, +0.09, -0.10, +0.07, -0.07
        };
        // A_FB in [16.00, 19.21]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_CDF_2012
        {
            "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@CDF-2012",
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.42, +0.23, -0.22, +0.09, -0.09
        };
        // F_L in [1.0, 6.0]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_1_to_6_CDF_2012
        {
            "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@CDF-2012",
            "B->K^*ll::F_Lavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.78, +0.13, -0.15, +0.08, -0.08
        };
        // F_L in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_CDF_2012
        {
            "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@CDF-2012",
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.45, +0.12, -0.12, +0.04, -0.04
        };
        // F_L in [16.00, 19.21]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_CDF_2012
        {
            "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@CDF-2012",
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.09, +0.14, -0.12, +0.08, -0.08
        };
        // A_T^{2} in [1.0, 6.0]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_T_2_1_to_6_CDF_2012
        {
            "B^0->K^*0mu^+mu^-::A_T^2[1.00,6.00]@CDF-2012",
            "B->K^*ll::A_T^2avg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.45, +2.24, -2.22, +0.76, -0.76
        };
        // A_T^{2} in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_T_2_14dot18_to_16_CDF_2012
        {
            "B^0->K^*0mu^+mu^-::A_T^2[14.18,16.00]@CDF-2012",
            "B->K^*ll::A_T^2avg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.15, +0.72, -0.72, +0.14, -0.14
        };
        // A_T^{2} in [16.00, 19.21]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_T_2_16_to_19dot21_CDF_2012
        {
            "B^0->K^*0mu^+mu^-::A_T^2[16.00,19.21]@CDF-2012",
            "B->K^*ll::A_T^2avg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.62, +0.56, -0.53, +0.13, -0.13
        };
        // A_{im} in [1.0, 6.0]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_im_1_to_6_CDF_2012
        {
            "B^0->K^*0mu^+mu^-::A_im[1.00,6.00]@CDF-2012",
            "B->K^*ll::J_9normavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.51, +0.28, -0.29, +0.15, -0.15
        };
        // A_{im} in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_im_14dot18_to_16_CDF_2012
        {
            "B^0->K^*0mu^+mu^-::A_im[14.18,16.00]@CDF-2012",
            "B->K^*ll::J_9normavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.16, +0.21, -0.22, +0.03, -0.03
        };
        // A_{im} in [16.00, 19.21]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_im_16_to_19dot21_CDF_2012
        {
            "B^0->K^*0mu^+mu^-::A_im[16.00,19.21]@CDF-2012",
            "B->K^*ll::J_9normavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.21 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.02, +0.26, -0.27, +0.04, -0.04
        };

        // B^+ -> K^+ mu^+ mu^-
        // A_FB in [1.0, 6.0]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_A_FB_1_to_6_CDF_2012
        {
            "B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@CDF-2012",
            "B->Kll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            -0.13, +0.10, -0.11, +0.02, -0.02
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_CDF_2012
        {
            "B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@CDF-2012",
            "B->Kll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            +0.07, +0.08, -0.08, +0.01, -0.01
        };
        // A_FB in [16.00, 23]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_A_FB_16_to_22dot86_CDF_2012
        {
            "B^+->K^+mu^+mu^-::A_FB[16.00,22.86]@CDF-2012",
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
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_1_to_6_LHCb_2012
        {
            "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@LHCb-2012",
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.42e-7 * 5, +0.04e-7 * 5, -0.04e-7 * 5, +0.04e-7 * 5, -0.04e-7 * 5
        };
        // BR in [14.18, 16.00], multiply with bin width 1.82 GeV^2
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_LHCb_2012
        {
            "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@LHCb-2012",
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.59e-7 * 1.82, +0.07e-7 * 1.82, -0.07e-7 * 1.82, +0.04e-7 * 1.82, -0.04e-7 * 1.82
        };
        // BR in [16.00, 19.00], multiply with bin width 3 GeV^2 (in the preliminary results, LHCb only integrated up to exactly 19.00 GeV^2!)
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_16_to_19_LHCb_2012
        {
            "B^0->K^*0mu^+mu^-::BR[16.00,19.00]@LHCb-2012",
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.44e-7 * 3, +0.05e-7 * 3, -0.05e-7 * 3, +0.03e-7 * 3, -0.03e-7 * 3
        };
        // A_FB in [1.0, 6.0]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_LHCb_2012
        {
            "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@LHCb-2012",
            "B->K^*ll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.18, +0.06, -0.06, +0.02, -0.01
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_LHCb_2012
        {
            "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@LHCb-2012",
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.49, +0.06, -0.04, +0.05, -0.02
        };
        // A_FB in [16.00, 19.00] (in the preliminary results, LHCb only integrated up to exactly 19.00 GeV^2!)
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_16_to_19_LHCb_2012
        {
            "B^0->K^*0mu^+mu^-::A_FB[16.00,19.00]@LHCb-2012",
            "B->K^*ll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.30, +0.07, -0.07, +0.01, -0.04
        };
        // F_L in [1.0, 6.0]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_1_to_6_LHCb_2012
        {
            "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@LHCb-2012",
            "B->K^*ll::F_Lavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.66, +0.06, -0.06, +0.04, -0.03
        };
        // F_L in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_LHCb_2012
        {
            "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@LHCb-2012",
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.35, +0.07, -0.06, +0.07, -0.02
        };
        // F_L in [16.00, 19.00] (in the preliminary results, LHCb only integrated up to exactly 19.00 GeV^2!)
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_16_to_19_LHCb_2012
        {
            "B^0->K^*0mu^+mu^-::F_L[16.00,19.00]@LHCb-2012",
            "B->K^*ll::F_Lavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.37, +0.06, -0.07, +0.03, -0.04
        };
        // A_{im} in [1.0, 6.0]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_im_1_to_6_LHCb_2012
        {
            "B^0->K^*0mu^+mu^-::A_im[1.00,6.00]@LHCb-2012",
            "B->K^*ll::J_9normavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.07, +0.07, -0.07, +0.02, -0.01
        };
        // A_{im} in [14.18, 16.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_im_14dot18_to_16_LHCb_2012
        {
            "B^0->K^*0mu^+mu^-::A_im[14.18,16.00]@LHCb-2012",
            "B->K^*ll::J_9normavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            -0.01, +0.08, -0.07, +0.04, -0.02
        };
        // A_{im} in [16.00, 19.00]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_im_16_to_19_LHCb_2012
        {
            "B^0->K^*0mu^+mu^-::A_im[16.00,19.00]@LHCb-2012",
            "B->K^*ll::J_9normavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.06, +0.09, -0.10, +0.03, -0.05
        };
        // S_3 in [1.0, 6.0], LHCb gives 2 * S_3
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_S_3_1_to_6_LHCb_2012
        {
            "B^0->K^*0mu^+mu^-::S_3[1.00,6.00]@LHCb-2012",
            "B->K^*ll::J_3normavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.10 / 2, +0.15 / 2, -0.16 / 2, +0.02 / 2, -0.01 / 2
        };
        // S_3 in [14.18, 16.00], LHCb gives 2 * S_3
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_S_3_14dot18_to_16_LHCb_2012
        {
            "B^0->K^*0mu^+mu^-::S_3[14.18,16.00]@LHCb-2012",
            "B->K^*ll::J_3normavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d"  }, { "l", "mu" } },
            +0.04 / 2, +0.15 / 2, -0.19 / 2, +0.04 / 2, -0.02 / 2
        };
        // S_3 in [16.00, 19.00], LHCb gives 2 * S_3
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_S_3_16_to_19_LHCb_2012
        {
            "B^0->K^*0mu^+mu^-::S_3[16.00,19.00]@LHCb-2012",
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
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_1_to_6_LHCb_2012
        {
            "B^+->K^+mu^+mu^-::BR[1.00,6.00]@LHCb-2012",
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            2.41e-8 * 5, +0.17e-8 * 5, -0.17e-8 * 5, +0.14e-8 * 5, -0.14e-8 * 5
        };
        // BR in [14.18, 16.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_14dot18_to_16_LHCb_2012
        {
            "B^+->K^+mu^+mu^-::BR[14.18,16.00]@LHCb-2012",
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            2.07e-8 * 1.82, +0.20e-8 * 1.82, -0.20e-8 * 1.82, +0.08e-8 * 1.82, -0.08e-8 * 1.82
        };
        // BR in [16.00, 18.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_16_to_18_LHCb_2012
        {
            "B^+->K^+mu^+mu^-::BR[16.00,18.00]@LHCb-2012",
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 18.00} },
            Options{ { "q", "u"  }, { "l", "mu" } },
            1.77e-8 * 2, +0.18e-8 * 2, -0.18e-8 * 2, +0.09e-8 * 2, -0.09e-8 * 2
        };
        // BR in [18.00, 22.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_18_to_22_LHCb_2012
        {
            "B^+->K^+mu^+mu^-::BR[18.00,22.00]@LHCb-2012",
            "B->Kll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 18.00 }, { "s_max", 22.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.78e-8 * 4, +0.10e-8 * 4, -0.10e-8 * 4, +0.04e-8 * 4, -0.04e-8 * 4
        };

        /* sign flipped wrt to table 1 in LHCb:2012C! */

        // A_FB in [1.0, 6.0]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_A_FB_1_to_6_LHCb_2012
        {
            "B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@LHCb-2012",
            "B->Kll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            -0.02, +0.03, -0.05, +0.01, -0.02
        };
        // A_FB in [14.18, 16.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_LHCb_2012
        {
            "B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@LHCb-2012",
            "B->Kll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            +0.01, +0.06, -0.12, +0.01, -0.01
        };
        // A_FB in [16.00, 18.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_A_FB_16_to_18_LHCb_2012
        {
            "B^+->K^+mu^+mu^-::A_FB[16.00,18.00]@LHCb-2012",
            "B->Kll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 18.00} },
            Options{ { "q", "u"  }, { "l", "mu" } },
            +0.09, +0.09, -0.07, +0.01, -0.02
        };
        // A_FB in [18.00, 22.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_A_FB_18_to_22_LHCb_2012
        {
            "B^+->K^+mu^+mu^-::A_FB[18.00,22.00]@LHCb-2012",
            "B->Kll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 18.00 }, { "s_max", 22.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            -0.02, +0.11, -0.11, +0.01, -0.01
        };

        // F_H in [1.0, 6.0]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_F_H_1_to_6_LHCb_2012
        {
            "B^+->K^+mu^+mu^-::F_H[1.00,6.00]@LHCb-2012",
            "B->Kll::F_Havg@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.05, +0.08, -0.05, +0.04, -0.02
        };
        // F_H in [14.18, 16.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_F_H_14dot18_to_16_LHCb_2012
        {
            "B^+->K^+mu^+mu^-::F_H[14.18,16.00]@LHCb-2012",
            "B->Kll::F_Havg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.08, +0.28, -0.08, +0.02, -0.01
        };
        // F_H in [16.00, 18.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_F_H_16_to_18_LHCb_2012
        {
            "B^+->K^+mu^+mu^-::F_H[16.00,18.00]@LHCb-2012",
            "B->Kll::F_Havg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 18.00} },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.18, +0.22, -0.14, +0.01, -0.04
        };
        // F_H in [18.00, 22.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_F_H_18_to_22_LHCb_2012
        {
            "B^+->K^+mu^+mu^-::F_H[18.00,22.00]@LHCb-2012",
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
        static const AmorosoTripleLimitConstraintEntry Bzero_to_dimuon_LHCb_2012
        {
            "B^0_s->mu^+mu^-::BR_limit@LHCb-2012",
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
        static const AmorosoModeConstraintEntry Bzero_to_dimuon_LHCb_Nov_2012
        {
            "B^0_s->mu^+mu^-::BR_limit@LHCb-Nov-2012",
            "B_q->ll::BR@Untagged",
            Kinematics(),
            Options{ { "q", "s"  }, { "l", "mu" } },
            0.0, 3.2e-9, 5.479025195034372e-9, 6.110034104385014e-9,
            2.3625776605e-09, 2.2682156277, 1.7296007586
        };

        // Data taken from [LHCb:2012E]
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_CP_1_to_6_LHCb_2012E
        {
            "B^0->K^*0mu^+mu^-::A_CP[1.00,6.00]@LHCb-2012",
            "B->K^*ll::A_CP@LargeRecoil",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.058, +0.064, -0.064, +0.009, -0.009
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_CP_14dot18_to_16_LHCb_2012E
        {
            "B^0->K^*0mu^+mu^-::A_CP[14.18,16.00]@LHCb-2012",
            "B->K^*ll::A_CP@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.201, +0.104, -0.104, +0.009, -0.009
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_CP_16_to_20_LHCb_2012E
        {
            "B^0->K^*0mu^+mu^-::A_CP[16.00,20.00]@LHCb-2012",
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
        static const GaussianConstraintEntry B_Bstar_mass_splitting_PDG_2012
        {
            "B^0::M_B^*-M_B@PDG-2012",
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
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_ATLAS_2013A
        {
            "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@ATLAS-2013A",
            "B->K^*ll::A_FB@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.07, +0.20, -0.20, +0.07, -0.07
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_ATLAS_2013A
        {
            "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@ATLAS-2013A",
            "B->K^*ll::A_FB@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.48, +0.19, -0.19, +0.05, -0.05
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_16_to_19_ATLAS_2013A
        {
            "B^0->K^*0mu^+mu^-::A_FB[16.00,19.00]@ATLAS-2013A",
            "B->K^*ll::A_FB@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.16, +0.10, -0.10, +0.03, -0.03
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_1_to_6_ATLAS_2013A
        {
            "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@ATLAS-2013A",
            "B->K^*ll::F_L@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.18, +0.15, -0.15, +0.03, -0.03
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_ATLAS_2013A
        {
            "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@ATLAS-2013A",
            "B->K^*ll::F_L@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.28, +0.16, -0.16, +0.03, -0.03
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_16_to_19_ATLAS_2013A
        {
            "B^0->K^*0mu^+mu^-::F_L[16.00,19.00]@ATLAS-2013A",
            "B->K^*ll::F_L@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.35, +0.08, -0.08, +0.02, -0.02
        };

        /*
         * BaBar Collaboration
         *
         * Data taken from [BaBar:2013A], Table I, p. 6.
         */
        static const GaussianConstraintEntry Bmix_to_Xs_dilepton_BR_BaBar_2013A
        {
            "B->X_sll::BR[1.0,6.0]@BaBar-2013A",
            "B->X_sll::BR@HLMW2005",
            Kinematics{ { "s_min", 1.0 }, { "s_max", 6.0 } },
            Options{ },
            1.60e-6, 0.41e-6, 0.39e-6, 0.25e-6, 0.22e-6
        };

        /*
         * Belle Collaboration
         *
         * Data taken from [Belle:2013A], Table XVI and XVII, pp. 37-38.
         */
        static const MultivariateGaussianConstraintEntry<6> Bzero_to_pi_l_nu_BR_Belle_2013A
        {
            "B^0->pi^+lnu::BR@Belle-2013A",
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
         * Data taken from [CMS:2013A]
         */
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_1_to_6_CMS_2013A
        {
            "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@CMS-2013A",
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            4.4e-8 * 5.0, +0.6e-8 * 5.0, -0.6e-8 * 5.0, +0.7e-8 * 5.0, -0.7e-8 * 5.0
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_CMS_2013A
        {
            "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@CMS-2013A",
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            4.6e-8 * 1.82, +0.9e-8 * 1.82, -0.8e-8 * 1.82, +0.8e-8 * 1.82, -0.8e-8 * 1.82
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_16_to_19_CMS_2013A
        {
            "B^0->K^*0mu^+mu^-::BR[16.00,19.00]@CMS-2013A",
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            5.2e-8 * 3.0, +0.6e-8 * 3.0, -0.6e-8 * 3.0, +0.9e-8 * 3.0, -0.9e-8 * 3.0
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_CMS_2013A
        {
            "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@CMS-2013A",
            "B->K^*ll::A_FB@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.07, +0.12, -0.12, +0.01, -0.01
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_CMS_2013A
        {
            "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@CMS-2013A",
            "B->K^*ll::A_FB@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.29, +0.09, -0.09, +0.05, -0.05
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_16_to_19_CMS_2013A
        {
            "B^0->K^*0mu^+mu^-::A_FB[16.00,19.00]@CMS-2013A",
            "B->K^*ll::A_FB@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.41, +0.05, -0.05, +0.03, -0.03
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_1_to_6_CMS_2013A
        {
            "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@CMS-2013A",
            "B->K^*ll::F_L@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.68, +0.10, -0.10, +0.02, -0.02
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_CMS_2013A
        {
            "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@CMS-2013A",
            "B->K^*ll::F_L@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.53, +0.12, -0.12, +0.03, -0.03
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_16_to_19_CMS_2013A
        {
            "B^0->K^*0mu^+mu^-::F_L[16.00,19.00]@CMS-2013A",
            "B->K^*ll::F_L@LowRecoil",
            Kinematics{ { "s_min", 16.00 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.44, +0.07, -0.07, +0.03, -0.03
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
        static const AmorosoConstraintEntry Bzero_to_dimuon_CMS_2013B
        {
            "B^0_s->mu^+mu^-::BR@CMS-2013B",
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
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_1_to_6_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@LHCb-2013",
            "B->K^*ll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            1.70e-7, +0.18e-7, -0.23e-7, +0.20e-7, -0.20e-7
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@LHCb-2013",
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            1.02e-7, +0.13e-7, -0.16e-7, +0.07e-7, -0.07e-7
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_BR_16_to_19_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::BR[16.00,19.00]@LHCb-2013",
            "B->K^*ll::BRavg@LowRecoil",
            Kinematics{ { "s_min", 16.0 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            1.23e-7, +0.15e-7, -0.17e-7, +0.12e-7, -0.12e-7
        };
        // F_L
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_1_to_6_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@LHCb-2013",
            "B->K^*ll::F_L@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.65, +0.08, -0.07, +0.03, -0.03
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@LHCb-2013",
            "B->K^*ll::F_L@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.33, +0.08, -0.07, +0.02, -0.03
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_F_L_16_to_19_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::F_L[16.00,19.00]@LHCb-2013",
            "B->K^*ll::F_L@LowRecoil",
            Kinematics{ { "s_min", 16.0 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.38, +0.09, -0.07, +0.03, -0.03
        };
        // A_FB
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@LHCb-2013",
            "B->K^*ll::A_FB@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.17, +0.06, -0.06, +0.01, -0.01
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@LHCb-2013",
            "B->K^*ll::A_FB@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.51, +0.05, -0.07, +0.02, -0.02
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_FB_16_to_19_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::A_FB[16.00,19.00]@LHCb-2013",
            "B->K^*ll::A_FB@LowRecoil",
            Kinematics{ { "s_min", 16.0 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.30, +0.08, -0.08, +0.02, -0.01
        };
        // S_3
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_S_3_1_to_6_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::S_3[1.00,6.00]@LHCb-2013",
            "B->K^*ll::J_3normavg@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.03, +0.07, -0.07, +0.01, -0.01
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_S_3_14dot18_to_16_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::S_3[14.18,16.00]@LHCb-2013",
            "B->K^*ll::J_3normavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.03, +0.09, -0.10, +0.01, -0.01
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_S_3_16_to_19_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::S_3[16.00,19.00]@LHCb-2013",
            "B->K^*ll::J_3normavg@LowRecoil",
            Kinematics{ { "s_min", 16.0 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.22, +0.10, -0.09, +0.02, -0.01
        };
        // S_9
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_S_9_1_to_6_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::S_9[1.00,6.00]@LHCb-2013",
            "B->K^*ll::J_9normavg@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.07, +0.09, -0.08, +0.01, -0.01
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_S_9_14dot18_to_16_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::S_9[14.18,16.00]@LHCb-2013",
            "B->K^*ll::J_9normavg@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.00, +0.09, -0.08, +0.01, -0.01
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_S_9_16_to_19_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::S_9[16.00,19.00]@LHCb-2013",
            "B->K^*ll::J_9normavg@LowRecoil",
            Kinematics{ { "s_min", 16.0 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.06, +0.11, -0.10, +0.01, -0.01
        };
#if 0
        // A_9
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_9_1_to_6_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::A_9[1.00,6.00]@LHCb-2013",
            "B->K^*ll::A_9@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.03, +0.08, -0.08, +0.01, -0.01
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_9_14dot18_to_16_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::A_9[14.18,16.00]@LHCb-2013",
            "B->K^*ll::A_9@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.06, +0.11, -0.08, +0.01, -0.01
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_9_16_to_19_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::A_9[16.00,19.00]@LHCb-2013",
            "B->K^*ll::A_9@LowRecoil",
            Kinematics{ { "s_min", 16.0 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.00, +0.11, -0.10, +0.01, -0.01
        };
#endif
        // A_T^2
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_T2_1_to_6_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::A_T^2[1.00,6.00]@LHCb-2013",
            "B->K^*ll::A_T^2@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.15, +0.39, -0.41, +0.03, -0.03
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_T2_14dot18_to_16_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::A_T^2[14.18,16.00]@LHCb-2013",
            "B->K^*ll::A_T^2@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.07, +0.26, -0.28, +0.02, -0.02
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_T2_16_to_19_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::A_T^2[16.00,19.00]@LHCb-2013",
            "B->K^*ll::A_T^2@LowRecoil",
            Kinematics{ { "s_min", 16.0 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.71, +0.35, -0.26, +0.06, -0.04
        };
        // A_T^re
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_Tre_1_to_6_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::A_T^re[1.00,6.00]@LHCb-2013",
            "B->K^*ll::A_T^re@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.66, +0.22, -0.24, +0.01, -0.04
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_Tre_14dot18_to_16_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::A_T^re[14.18,16.00]@LHCb-2013",
            "B->K^*ll::A_T^re@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -1.00, +0.05, -0.01, +0.02, -0.01
            // (unphysical) lower errors (both stat and syst) adjusted to -0.01 to work
            // around limitations in the asymmetric gaussian likelihood block.
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_A_Tre_16_to_19_LHCb_2013B
        {
            "B^0->K^*0mu^+mu^-::A_T^re[16.00,19.00]@LHCb-2013",
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
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_Pprime_4_1_to_6_LHCb_2013C
        {
            "B^0->K^*0mu^+mu^-::P'_4[1.00,6.00]@LHCb-2013",
            "B->K^*ll::P'_4@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.58, +0.32, -0.36, +0.06, -0.06
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_Pprime_4_14dot18_to_16_LHCb_2013C
        {
            "B^0->K^*0mu^+mu^-::P'_4[16.00,19.00]@LHCb-2013",
            "B->K^*ll::P'_4@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.18, +0.54, -0.70, +0.08, -0.08
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_Pprime_4_16_to_19_LHCb_2013C
        {
            "B^0->K^*0mu^+mu^-::P'_4[16.00,19.00]@LHCb-2013",
            "B->K^*ll::P'_4@LowRecoil",
            Kinematics{ { "s_min", 16.0 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.70, +0.44, -0.52, +0.06, -0.06
        };
        // P'_5
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_Pprime_5_1_to_6_LHCb_2013C
        {
            "B^0->K^*0mu^+mu^-::P'_5[1.00,6.00]@LHCb-2013",
            "B->K^*ll::P'_5@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.21, +0.20, -0.21, +0.03, -0.03
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_Pprime_5_14dot18_to_16_LHCb_2013C
        {
            "B^0->K^*0mu^+mu^-::P'_5[14.18,16.00]@LHCb-2013",
            "B->K^*ll::P'_5@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.79, +0.20, -0.13, +0.18, -0.18
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_Pprime_5_16_to_19_LHCb_2013C
        {
            "B^0->K^*0mu^+mu^-::P'_5[16.00,19.00]@LHCb-2013",
            "B->K^*ll::P'_5@LowRecoil",
            Kinematics{ { "s_min", 16.0 }, { "s_max", 19.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            -0.60, +0.19, -0.16, +0.09, -0.09
        };
        // P'_6, LHCb uses a different sign
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_Pprime_6_1_to_6_LHCb_2013C
        {
            "B^0->K^*0mu^+mu^-::P'_6[1.00,6.00]@LHCb-2013",
            "B->K^*ll::P'_6@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.18, +0.21, -0.21, +0.03, -0.03
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_Pprime_6_14dot18_to_16_LHCb_2013C
        {
            "B^0->K^*0mu^+mu^-::P'_6[14.18,16.00]@LHCb-2013",
            "B->K^*ll::P'_6@LowRecoil",
            Kinematics{ { "s_min", 14.18 }, { "s_max", 16.00 } },
            Options{ { "q", "d" }, { "l", "mu" } },
            +0.18, +0.24, -0.25, +0.03, -0.03
        };
        static const GaussianConstraintEntry Bzero_to_Kstarzero_dimuon_Pprime_6_16_to_19_LHCb_2013C
        {
            "B^0->K^*0mu^+mu^-::P'_6[16.00,19.00]@LHCb-2013",
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
        static const AmorosoConstraintEntry Bzero_to_dimuon_LHCb_2013D
        {
            "B^0_s->mu^+mu^-::BR@LHCb-2013D",
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
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_1dot1_to_2_LHCb_2014
        {
            "B^+->K^+mu^+mu^-::BR[1.10,2.00]@LHCb-2014",
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.1 }, { "s_max", 2.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            2.33e-8 * 0.9, +0.15e-8 * 0.9, -0.15e-8 * 0.9, +0.12e-8 * 0.9, -0.12e-8 * 0.9
        };
        // BR in [2.00, 3.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_2_to_3_LHCb_2014
        {
            "B^+->K^+mu^+mu^-::BR[2.00,3.00]@LHCb-2014",
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 2.00 }, { "s_max", 3.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            2.82e-8 * 1.0, +0.16e-8 * 1.0, -0.16e-8 * 1.0, +0.14e-8 * 1.0, -0.14e-8 * 1.0
        };
        // BR in [3.00, 4.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_3_to_4_LHCb_2014
        {
            "B^+->K^+mu^+mu^-::BR[3.00,4.00]@LHCb-2014",
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 3.00 }, { "s_max", 4.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            2.54e-8 * 1.0, +0.15e-8 * 1.0, -0.15e-8 * 1.0, +0.13e-8 * 1.0, -0.13e-8 * 1.0
        };
        // BR in [4.00, 5.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_4_to_5_LHCb_2014
        {
            "B^+->K^+mu^+mu^-::BR[4.00,5.00]@LHCb-2014",
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 4.00 }, { "s_max", 5.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            2.21e-8 * 1.0, +0.14e-8 * 1.0, -0.14e-8 * 1.0, +0.11e-8 * 1.0, -0.11e-8 * 1.0
        };
        // BR in [5.00, 6.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_5_to_6_LHCb_2014
        {
            "B^+->K^+mu^+mu^-::BR[5.00,6.00]@LHCb-2014",
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 5.00 }, { "s_max", 6.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            2.31e-8 * 1.0, +0.14e-8 * 1.0, -0.14e-8 * 1.0, +0.12e-8 * 1.0, -0.12e-8 * 1.0
        };
        // BR in [1.10, 6.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_1dot1_to_6_LHCb_2014
        {
            "B^+->K^+mu^+mu^-::BR[1.10,6.00]@LHCb-2014",
            "B->Kll::BRavg@LargeRecoil",
            Kinematics{ { "s_min", 1.10 }, { "s_max", 6.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            2.42e-8 * 4.9, +0.07e-8 * 4.9, -0.07e-8 * 4.9, +0.12e-8 * 4.9, -0.12e-8 * 4.9
        };

        // BR in [15.00, 22.00]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_BR_15_to_22_LHCb_2014
        {
            "B^+->K^+mu^+mu^-::BR[15.00,22.00]@LHCb-2014",
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
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_A_FB_1dot1_to_6_LHCb_2014
        {
            "B^+->K^+mu^+mu^-::A_FB[1.10,6.00]@LHCb-2014",
            "B->Kll::A_FBavg@LargeRecoil",
            Kinematics{ { "s_min", 1.1 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            -0.005, +0.015, -0.015, +0.01, -0.01
        };
        // A_FB in [15, 22]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_A_FB_15_to_22_LHCb_2014
        {
            "B^+->K^+mu^+mu^-::A_FB[15.00,22.00]@LHCb-2014",
            "B->Kll::A_FBavg@LowRecoil",
            Kinematics{ { "s_min", 15.00 }, { "s_max", 22.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.015, +0.015, -0.015, +0.01, -0.01
        };
        // F_H in [1.1, 6.0]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_F_H_1dot1_to_6_LHCb_2014
        {
            "B^+->K^+mu^+mu^-::F_H[1.10,6.00]@LHCb-2014",
            "B->Kll::F_Havg@LargeRecoil",
            Kinematics{ { "s_min", 1.1 }, { "s_max", 6.0 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.03, +0.03, -0.03, +0.02, -0.02
        };
        // F_H in [15, 22]
        static const GaussianConstraintEntry Bplus_to_Kplus_dimuon_F_H_15_to_22_LHCb_2014
        {
            "B^+->K^+mu^+mu^-::F_H[15.00,22.00]@LHCb-2014",
            "B->Kll::F_Havg@LowRecoil",
            Kinematics{ { "s_min", 15.00 }, { "s_max", 22.00 } },
            Options{ { "q", "u"  }, { "l", "mu" } },
            0.035, +0.035, -0.035, +0.02, -0.02
        };

        /*
         * Data taken from [LHCb:2014C]
         */
        static const GaussianConstraintEntry Bplus_to_Kplus_dilepton_r_k_1_to_6_LHCb_2014C
        {
            "B^+->K^+l^+l^-::R_K[1.00,6.00]@LHCb-2014",
            "B->Kll::R_K@LargeRecoil",
            Kinematics{ { "s_min", 1.00 }, { "s_max", 6.00 } },
            Options{ },
            0.745, 0.090, 0.074, 0.036, 0.036
        };

        /*
         * CMS + LHCb Collaborations
         *
         * Data taken from [CMS-LHCb:2014A]
         */

        /* fit Amoroso to result assuming
         * a) cdf(0) = 0
         * b) mode = 2.8
         * c) 68% in (2.2, 3.5)
         * d) pdf(2.2) = pdf(3.5) (smallest interval)
         */
        static const AmorosoConstraintEntry B_s_to_dimuon_CMS_LHCb_2014A
        {
            "B^0_s->mu^+mu^-::BR@CMS-LHCb-2014",
            "B_q->ll::BR@Untagged",
            Kinematics(),
            Options{ { "q", "s"  }, { "l", "mu" } },
            0.0, 0.154243e-9, 19.4025, 1.0048
        };

        static const AmorosoConstraintEntry Bzero_to_dimuon_CMS_LHCb_2014A
        {
            "B^0_d->mu^+mu^-::BR@CMS-LHCb-2014",
            "B_q->ll::BR@Untagged",
            Kinematics(),
            Options{ { "q", "d"  }, { "l", "mu" } },
            0.0, 3.09469e-10, 2.08667, 1.98708
        };
        ///@}

        ///@name 2015
        ///@{
        /*
         * LHCb Collaboration
         *
         * Data taken from [LHCb:2015A], table 6, 50 and appendix F.
         * We combine correlations between the statistic uncertainties and systematic uncertainties.
         */
        static const MultivariateGaussianConstraintEntry<5> Bzero_to_Kstarzero_dimuon_aobs_moments_1dot1_to_2_LHCb_2015A
        {
            "B^0->K^*0mu^+mu^-::AngularObservables[1.10,2.00]@LHCb-2015A",
            {{ "B->K^*ll::F_L@LargeRecoil", "B->K^*ll::S_3@LargeRecoil", "B->K^*ll::S_4@LargeRecoil","B->K^*ll::S_5@LargeRecoil", "B->K^*ll::A_FB@LargeRecoil" }},
            {{
                 Kinematics{ { "s_min", 1.1 }, { "s_max", 2.0 } },
                 Kinematics{ { "s_min", 1.1 }, { "s_max", 2.0 } },
                 Kinematics{ { "s_min", 1.1 }, { "s_max", 2.0 } },
                 Kinematics{ { "s_min", 1.1 }, { "s_max", 2.0 } },
                 Kinematics{ { "s_min", 1.1 }, { "s_max", 2.0 } }
            }},
            {{
                 Options{ { "q", "d" }, { "l", "mu" } },
                 Options{ { "q", "d" }, { "l", "mu" } },
                 Options{ { "q", "d" }, { "l", "mu" } },
                 Options{ { "q", "d" }, { "l", "mu" } },
                 Options{ { "q", "d" }, { "l", "mu" } }
            }},
            {{  0.768,  0.065,  0.127,  0.286, -0.333 }},
            {{  0.0,    0.0,    0.0,    0.0,    0.0   }},
            {{  0.0,    0.0,    0.0,    0.0,    0.0,  }},
            {{  0.138,  0.132,  0.187,  0.170,  0.123 }},
            {{
                {{ +1.00, +0.02, +0.05, +0.12, +0.01 }},
                {{ +0.02, +1.00, +0.00, +0.05, +0.05 }},
                {{ +0.05, +0.00, +1.00, -0.00, -0.03 }},
                {{ +0.12, +0.05, -0.00, +1.00, -0.06 }},
                {{ +0.01, +0.05, -0.03, -0.06, +1.00 }},
            }}
        };
        static const MultivariateGaussianConstraintEntry<5> Bzero_to_Kstarzero_dimuon_aobs_moments_2_to_3_LHCb_2015A
        {
            "B^0->K^*0mu^+mu^-::AngularObservables[2.00,3.00]@LHCb-2015A",
            {{ "B->K^*ll::F_L@LargeRecoil", "B->K^*ll::S_3@LargeRecoil", "B->K^*ll::S_4@LargeRecoil","B->K^*ll::S_5@LargeRecoil", "B->K^*ll::A_FB@LargeRecoil" }},
            {{
                 Kinematics{ { "s_min", 2.0 }, { "s_max", 3.0 } },
                 Kinematics{ { "s_min", 2.0 }, { "s_max", 3.0 } },
                 Kinematics{ { "s_min", 2.0 }, { "s_max", 3.0 } },
                 Kinematics{ { "s_min", 2.0 }, { "s_max", 3.0 } },
                 Kinematics{ { "s_min", 2.0 }, { "s_max", 3.0 } }
            }},
            {{
                 Options{ { "q", "d" }, { "l", "mu" } },
                 Options{ { "q", "d" }, { "l", "mu" } },
                 Options{ { "q", "d" }, { "l", "mu" } },
                 Options{ { "q", "d" }, { "l", "mu" } },
                 Options{ { "q", "d" }, { "l", "mu" } }
            }},
            {{  0.690,  0.006, -0.339,  0.206, -0.158 }},
            {{  0.0,    0.0,    0.0,    0.0,    0.0   }},
            {{  0.0,    0.0,    0.0,    0.0,    0.0,  }},
            {{  0.100,  0.100,  0.134,  0.123,  0.085 }},
            {{
                {{ +1.00, -0.11, -0.14, -0.05, +0.12 }},
                {{ -0.11, +1.00, -0.10, +0.07, +0.03 }},
                {{ -0.14, -0.10, +1.00, -0.08, -0.00 }},
                {{ -0.05, +0.07, -0.08, +1.00, -0.10 }},
                {{ +0.12, +0.03, -0.00, -0.10, +1.00 }},
            }}
        };
        static const MultivariateGaussianConstraintEntry<5> Bzero_to_Kstarzero_dimuon_aobs_moments_3_to_4_LHCb_2015A
        {
            "B^0->K^*0mu^+mu^-::AngularObservables[3.00,4.00]@LHCb-2015A",
            {{ "B->K^*ll::F_L@LargeRecoil", "B->K^*ll::S_3@LargeRecoil", "B->K^*ll::S_4@LargeRecoil","B->K^*ll::S_5@LargeRecoil", "B->K^*ll::A_FB@LargeRecoil" }},
            {{
                 Kinematics{ { "s_min", 3.0 }, { "s_max", 4.0 } },
                 Kinematics{ { "s_min", 3.0 }, { "s_max", 4.0 } },
                 Kinematics{ { "s_min", 3.0 }, { "s_max", 4.0 } },
                 Kinematics{ { "s_min", 3.0 }, { "s_max", 4.0 } },
                 Kinematics{ { "s_min", 3.0 }, { "s_max", 4.0 } }
            }},
            {{
                 Options{ { "q", "d" }, { "l", "mu" } },
                 Options{ { "q", "d" }, { "l", "mu" } },
                 Options{ { "q", "d" }, { "l", "mu" } },
                 Options{ { "q", "d" }, { "l", "mu" } },
                 Options{ { "q", "d" }, { "l", "mu" } }
            }},
            {{  0.873,  0.078, -0.046, -0.110, -0.041 }},
            {{  0.0,    0.0,    0.0,    0.0,    0.0   }},
            {{  0.0,    0.0,    0.0,    0.0,    0.0,  }},
            {{  0.132,  0.127,  0.200,  0.166,  0.091 }},
            {{
                {{ +1.00, +0.04, -0.00, -0.05, -0.02 }},
                {{ +0.04, +1.00, +0.00, -0.05, +0.02 }},
                {{ -0.00, +0.00, +1.00, +0.17, +0.05 }},
                {{ -0.05, -0.05, +0.17, +1.00, +0.01 }},
                {{ -0.02, +0.02, +0.05, +0.01, +1.00 }},
            }}
        };
        static const MultivariateGaussianConstraintEntry<5> Bzero_to_Kstarzero_dimuon_aobs_moments_4_to_5_LHCb_2015A
        {
            "B^0->K^*0mu^+mu^-::AngularObservables[4.00,5.00]@LHCb-2015A",
            {{ "B->K^*ll::F_L@LargeRecoil", "B->K^*ll::S_3@LargeRecoil", "B->K^*ll::S_4@LargeRecoil","B->K^*ll::S_5@LargeRecoil", "B->K^*ll::A_FB@LargeRecoil" }},
            {{
                 Kinematics{ { "s_min", 4.0 }, { "s_max", 5.0 } },
                 Kinematics{ { "s_min", 4.0 }, { "s_max", 5.0 } },
                 Kinematics{ { "s_min", 4.0 }, { "s_max", 5.0 } },
                 Kinematics{ { "s_min", 4.0 }, { "s_max", 5.0 } },
                 Kinematics{ { "s_min", 4.0 }, { "s_max", 5.0 } }
            }},
            {{
                 Options{ { "q", "d" }, { "l", "mu" } },
                 Options{ { "q", "d" }, { "l", "mu" } },
                 Options{ { "q", "d" }, { "l", "mu" } },
                 Options{ { "q", "d" }, { "l", "mu" } },
                 Options{ { "q", "d" }, { "l", "mu" } }
            }},
            {{  0.899,  0.200, -0.148, -0.306, +0.052 }},
            {{  0.0,    0.0,    0.0,    0.0,    0.0   }},
            {{  0.0,    0.0,    0.0,    0.0,    0.0,  }},
            {{  0.107,  0.099,  0.161,  0.140,  0.080 }},
            {{
                {{ +1.00, -0.02, +0.02, -0.10, -0.03 }},
                {{ -0.02, +1.00, -0.10, -0.11, +0.03 }},
                {{ +0.02, -0.10, +1.00, +0.15, -0.03 }},
                {{ -0.10, -0.11, +0.15, +1.00, -0.03 }},
                {{ -0.03, +0.03, -0.03, -0.03, +1.00 }},
            }}
        };
        static const MultivariateGaussianConstraintEntry<5> Bzero_to_Kstarzero_dimuon_aobs_moments_5_to_6_LHCb_2015A
        {
            "B^0->K^*0mu^+mu^-::AngularObservables[5.00,6.00]@LHCb-2015A",
            {{ "B->K^*ll::F_L@LargeRecoil", "B->K^*ll::S_3@LargeRecoil", "B->K^*ll::S_4@LargeRecoil","B->K^*ll::S_5@LargeRecoil", "B->K^*ll::A_FB@LargeRecoil" }},
            {{
                 Kinematics{ { "s_min", 5.0 }, { "s_max", 6.0 } },
                 Kinematics{ { "s_min", 5.0 }, { "s_max", 6.0 } },
                 Kinematics{ { "s_min", 5.0 }, { "s_max", 6.0 } },
                 Kinematics{ { "s_min", 5.0 }, { "s_max", 6.0 } },
                 Kinematics{ { "s_min", 5.0 }, { "s_max", 6.0 } }
            }},
            {{
                 Options{ { "q", "d" }, { "l", "mu" } },
                 Options{ { "q", "d" }, { "l", "mu" } },
                 Options{ { "q", "d" }, { "l", "mu" } },
                 Options{ { "q", "d" }, { "l", "mu" } },
                 Options{ { "q", "d" }, { "l", "mu" } }
            }},
            {{  0.644, -0.122, -0.273, -0.095, +0.057 }},
            {{  0.0,    0.0,    0.0,    0.0,    0.0   }},
            {{  0.0,    0.0,    0.0,    0.0,    0.0,  }},
            {{  0.128,  0.123,  0.185,  0.140,  0.135 }},
            {{
                {{ +1.00, +0.01, -0.02, -0.02, -0.02 }},
                {{ +0.01, +1.00, -0.03, -0.07, -0.10 }},
                {{ -0.02, -0.03, +1.00, +0.09, +0.03 }},
                {{ -0.02, -0.07, +0.09, +1.00, -0.07 }},
                {{ -0.02, -0.10, +0.03, -0.07, +1.00 }},
            }}
        };

        /*
         * LHCb Collaboration
         *
         * Data taken from [LHCb:2015B], table 6, 50 and appendix F.
         * We combine correlations between the statistic uncertainties and systematic uncertainties.
         */
        static const GaussianConstraintEntry Lambdab_to_Lambda_dimuon_br_15_to_20_LHCb_2015B
        {
            "Lambda_b->Lambdamu^+mu^-::BR[15.0,20.0]@LHCb-2015B",
            "Lambda_b->Lambdall::BR@LowRecoil",
            Kinematics{ { "s_min", 15.0 }, { "s_max", 20.00 } },
            Options{ { "l", "mu" } },
            1.20e-7 * 5.0, +0.09e-7 * 5.0, -0.09e-7 * 5.0, +0.25e-7 * 5.0, -0.25e-7 * 5.0
        };
        static const GaussianConstraintEntry Lambdab_to_Lambda_dimuon_f_0_15_to_20_LHCb_2015B
        {
            "Lambda_b->Lambdamu^+mu^-::F_0[15.0,20.0]@LHCb-2015B",
            "Lambda_b->Lambdall::F_0@LowRecoil",
            Kinematics{ { "s_min", 15.0 }, { "s_max", 20.00 } },
            Options{ { "l", "mu" } },
            0.61, +0.11, -0.14, +0.03, -0.03
        };
        static const GaussianConstraintEntry Lambdab_to_Lambda_dimuon_a_fb_l_15_to_20_LHCb_2015B
        {
            "Lambda_b->Lambdamu^+mu^-::A_FB^l[15.0,20.0]@LHCb-2015B",
            "Lambda_b->Lambdall::A_FB^l@LowRecoil",
            Kinematics{ { "s_min", 15.0 }, { "s_max", 20.00 } },
            Options{ { "l", "mu" } },
            -0.05, +0.09, -0.09, +0.03, -0.03
        };
        static const GaussianConstraintEntry Lambdab_to_Lambda_dimuon_a_fb_h_15_to_20_LHCb_2015B
        {
            "Lambda_b->Lambdamu^+mu^-::A_FB^h[15.0,20.0]@LHCb-2015B",
            "Lambda_b->Lambdall::A_FB^h@LowRecoil",
            Kinematics{ { "s_min", 15.0 }, { "s_max", 20.00 } },
            Options{ { "l", "mu" } },
            -0.29, +0.07, -0.07, +0.03, -0.03
        };
        ///@}
        ///@}

        /*
         * Theoretical Constraints from e.g. Lattice QCD.
         */

        ///@name 2013
        ///@{
        /*
         * Reproduced from [HPQCD:2013A].
         */
        static const MultivariateGaussianConstraintEntry<9> B_to_K_fzero_fplus_ftensor_17_to_23_HPQCD_2013A
        {
            "B->K::f_0+f_++f_T@HPQCD-2013A",
            {{ "B->K::f_0(s)", "B->K::f_0(s)", "B->K::f_0(s)",
               "B->K::f_+(s)", "B->K::f_+(s)", "B->K::f_+(s)",
               "B->K::f_T(s)", "B->K::f_T(s)", "B->K::f_T(s)" }},
            {{ Kinematics{ { "s", 17.0 } }, Kinematics{ { "s", 20.0 } }, Kinematics{ { "s", 23.0 } },
               Kinematics{ { "s", 17.0 } }, Kinematics{ { "s", 20.0 } }, Kinematics{ { "s", 23.0 } },
               Kinematics{ { "s", 17.0 } }, Kinematics{ { "s", 20.0 } }, Kinematics{ { "s", 23.0 } } }},
            {{ Options{ { "form-factors", "KMPW2010" } }, Options{ { "form-factors", "KMPW2010" } }, Options{ { "form-factors", "KMPW2010" } },
               Options{ { "form-factors", "KMPW2010" } }, Options{ { "form-factors", "KMPW2010" } }, Options{ { "form-factors", "KMPW2010" } },
               Options{ { "form-factors", "KMPW2010" } }, Options{ { "form-factors", "KMPW2010" } }, Options{ { "form-factors", "KMPW2010" } } }},
            {{ 0.615697,  0.722728,  0.86973,   1.13053,   1.62747,   2.67559,  1.01869,   1.47371, 2.42477  }},
            {{ 0.0423107, 0.0464784, 0.0545373, 0.0730107, 0.0973055, 0.186149, 0.0842842, 0.11937, 0.249658 }},
            {{ 0.0423107, 0.0464784, 0.0545373, 0.0730107, 0.0973055, 0.186149, 0.0842842, 0.11937, 0.249658 }},
            {{ +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00, +0.00 }}, // we assign no systematic uncertainty
            {{
                {{ 1.,        0.929602,  0.727828,  0.576305,  0.477212, 0.230181,  0.128042,  0.0998671, 0.0477771 }},
                {{ 0.929602,  1.,        0.872508,  0.405678,  0.451924, 0.243271,  0.0586115, 0.0933393, 0.0749594 }},
                {{ 0.727828,  0.872508,  1.,        0.280099,  0.34324,  0.317962,  0.0103097, 0.046188,  0.0580487 }},
                {{ 0.576305,  0.405678,  0.280099,  1.,        0.782097, 0.298457,  0.363865,  0.212712,  0.0508819 }},
                {{ 0.477212,  0.451924,  0.34324,   0.782097,  1.,       0.712792,  0.222467,  0.262759,  0.197032  }},
                {{ 0.230181,  0.243271,  0.317962,  0.298457,  0.712792, 1.,       -0.0249157, 0.140985,  0.220042  }},
                {{ 0.128042,  0.0586115, 0.0103097, 0.363865,  0.222467,-0.0249157, 1.,        0.785177,  0.457065  }},
                {{ 0.0998671, 0.0933393, 0.046188,  0.212712,  0.262759, 0.140985,  0.785177,  1.,        0.853052  }},
                {{ 0.0477771, 0.0749594, 0.0580487, 0.0508819, 0.197032, 0.220042,  0.457065,  0.853052,  1.        }}
            }},
            0u
        };

        /*
         * Reproduced from [MILC:2013A].
         */
        // B -> K^*
        static const MultivariateGaussianConstraintEntry<2> B_to_Kstar_V_15_to_19dot21_MILC_2013A
        {
            "B->K^*::V@MILC-2013A",
            {{ "B->K^*::V(s)", "B->K^*::V(s)" }},
            {{ Kinematics{ { "s", 15.0 } }, Kinematics{ { "s", 19.21 } } }},
            {{ Options{ { "form-factors", "KMPW2010" } }, Options{ { "form-factors", "KMPW2010" } } }},
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
        static const MultivariateGaussianConstraintEntry<2> B_to_Kstar_A1_15_to_19dot21_MILC_2013A
        {
            "B->K^*::A_1@MILC-2013A",
            {{ "B->K^*::A_1(s)", "B->K^*::A_1(s)" }},
            {{ Kinematics{ { "s", 15.0 } }, Kinematics{ { "s", 19.21 } } }},
            {{ Options{ { "form-factors", "KMPW2010" } }, Options{ { "form-factors", "KMPW2010" } } }},
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
        static const MultivariateGaussianConstraintEntry<2> B_to_Kstar_A12_15_to_19dot21_MILC_2013A
        {
            "B->K^*::A_12@MILC-2013A",
            {{ "B->K^*::A_12(s)", "B->K^*::A_12(s)" }},
            {{ Kinematics{ { "s", 15.0 } }, Kinematics{ { "s", 19.21 } } }},
            {{ Options{ { "form-factors", "KMPW2010" } }, Options{ { "form-factors", "KMPW2010" } } }},
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
        // B_s -> K^*
        static const MultivariateGaussianConstraintEntry<2> Bs_to_Kstar_V_15_to_19dot21_MILC_2013A
        {
            "B_s->K^*::V@MILC-2013A",
            {{ "B_s->K^*::V(s)", "B_s->K^*::V(s)" }},
            {{ Kinematics{ { "s", 15.0 } }, Kinematics{ { "s", 19.21 } } }},
            {{ Options{ { "form-factors", "FMvD2015" } }, Options{ { "form-factors", "FMvD2015" } } }},
            {{ 0.872, 1.722 }},
            {{ 0.066, 0.066 }},
            {{ 0.066, 0.066 }},
            {{ 0.000, 0.000 }}, // we assign no systematic uncertainty
            {{
                {{ 1.000, 0.271 }},
                {{ 0.271, 1.000 }},
            }},
        0u
        };
        static const MultivariateGaussianConstraintEntry<2> Bs_to_Kstar_A1_15_to_19dot21_MILC_2013A
        {
            "B_s->K^*::A_1@MILC-2013A",
            {{ "B_s->K^*::A_1(s)", "B_s->K^*::A_1(s)" }},
            {{ Kinematics{ { "s", 15.0 } }, Kinematics{ { "s", 19.21 } } }},
            {{ Options{ { "form-factors", "FMvD2015" } }, Options{ { "form-factors", "FMvD2015" } } }},
            {{ 0.427, 0.548 }},
            {{ 0.015, 0.015 }},
            {{ 0.015, 0.015 }},
            {{ 0.000, 0.000 }}, // we assign no systematic uncertainty
            {{
                {{ 1.000, 0.305 }},
                {{ 0.305, 1.000 }},
            }},
        0u
        };
        static const MultivariateGaussianConstraintEntry<2> Bs_to_Kstar_A12_15_to_19dot21_MILC_2013A
        {
            "B_s->K^*::A_12@MILC-2013A",
            {{ "B_s->K^*::A_12(s)", "B_s->K^*::A_12(s)" }},
            {{ Kinematics{ { "s", 15.0 } }, Kinematics{ { "s", 19.21 } } }},
            {{ Options{ { "form-factors", "FMvD2015" } }, Options{ { "form-factors", "FMvD2015" } } }},
            {{ 0.342, 0.408 }},
            {{ 0.016, 0.016 }},
            {{ 0.016, 0.016 }},
            {{ 0.000, 0.000 }}, // we assign no systematic uncertainty
            {{
                {{ 1.000, 0.334 }},
                {{ 0.334, 1.000 }},
            }},
        0u
        };
        ///@}

        ///@name 2014
        ///@{
        /*
         * From [IKMvD:2014].
         */
        static const MultivariateGaussianConstraintEntry<6> B_to_pi_fp_IKMvD_2014
        {
            "B->pi::f_+@IKMvD-2014",
            {{ "B->pi::f_+(s)", "B->pi::f_+'(s)", "B->pi::f_+''(s)", "B->pi::f_+(s)", "B->pi::f_+'(s)", "B->pi::f_+''(s)" }},
            {{ Kinematics{ { "s", 0.0 } }, Kinematics{ { "s", 0.0 } }, Kinematics{ { "s", 0.0 } }, Kinematics{ { "s", 10.0 } },Kinematics{ { "s", 10.0 } }, Kinematics{ { "s", 10.0 } } }},
            {{ Options{ { "form-factors", "BCL2008" } }, Options{ { "form-factors", "BCL2008" } }, Options{ { "form-factors", "BCL2008" } }, Options{ { "form-factors", "BCL2008" } }, Options{ { "form-factors", "BCL2008" } }, Options{ { "form-factors", "BCL2008" } } }},
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

        /*
         * From [BFvD2014], based on reproduction of data points from [DLMW2012] and [FY2011], and subleading Isgur-Wise functions
         * approximations.
         */
        static const MultivariateGaussianConstraintEntry<4> LambdaB_to_Lambda_all_v_and_a_0_BFvD2014
        {
            "Lambda_b->Lambda::f_perp+long^V+A@BFvD2014",
            {{ "Lambda_b->Lambda::f_perp^V(s)", "Lambda_b->Lambda::f_perp^A(s)", "Lambda_b->Lambda::f_long^V(s)", "Lambda_b->Lambda::f_long^A(s)" }},
            {{ Kinematics{ { "s", 0.0 } }, Kinematics{ { "s", 0.0 } }, Kinematics{ { "s", 0.0 } }, Kinematics{ { "s", 0.0 } } }},
            {{ Options{ { "form-factors", "BFvD2014" } }, Options{ { "form-factors", "BFvD2014" } }, Options{ { "form-factors", "BFvD2014" } }, Options{ { "form-factors", "BFvD2014" } } }},
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
        static const MultivariateGaussianConstraintEntry<2> LambdaB_to_Lambda_fperpV_13dot5_to_20dot3_BFvD2014
        {
            "Lambda_b->Lambda::f_perp^V@BFvD2014",
            {{ "Lambda_b->Lambda::f_perp^V(s)", "Lambda_b->Lambda::f_perp^V(s)" }},
            {{ Kinematics{ { "s", 13.5 } }, Kinematics{ { "s", 20.5 } } }},
            {{ Options{ { "form-factors", "BFvD2014" } }, Options{ { "form-factors", "BFvD2014" } } }},
            {{ 0.73, 1.40 }},
            {{ 0.20, 0.20 }},
            {{ 0.20, 0.20 }},
            {{ 0.00, 0.00 }},
            {{
                 {{ 1.000, 0.000 }},
                 {{ 0.000, 1.000 }}
            }}
        };
        static const MultivariateGaussianConstraintEntry<2> LambdaB_to_Lambda_fperpA_13dot5_to_20dot3_BFvD2014
        {
            "Lambda_b->Lambda::f_perp^A@BFvD2014",
            {{ "Lambda_b->Lambda::f_perp^A(s)", "Lambda_b->Lambda::f_perp^A(s)" }},
            {{ Kinematics{ { "s", 13.5 } }, Kinematics{ { "s", 20.5 } } }},
            {{ Options{ { "form-factors", "BFvD2014" } }, Options{ { "form-factors", "BFvD2014" } } }},
            {{ 0.48, 0.84 }},
            {{ 0.19, 0.19 }},
            {{ 0.19, 0.19 }},
            {{ 0.00, 0.00 }},
            {{
                 {{ 1.000, 0.000 }},
                 {{ 0.000, 1.000 }}
            }}
        };
        static const MultivariateGaussianConstraintEntry<2> LambdaB_to_Lambda_flongV_13dot5_to_20dot3_BFvD2014
        {
            "Lambda_b->Lambda::f_long^V@BFvD2014",
            {{ "Lambda_b->Lambda::f_long^V(s)", "Lambda_b->Lambda::f_long^V(s)" }},
            {{ Kinematics{ { "s", 13.5 } }, Kinematics{ { "s", 20.5 } } }},
            {{ Options{ { "form-factors", "BFvD2014" } }, Options{ { "form-factors", "BFvD2014" } } }},
            {{ 0.72, 1.39 }},
            {{ 0.21, 0.21 }},
            {{ 0.21, 0.21 }},
            {{ 0.00, 0.00 }},
            {{
                 {{ 1.000, 0.000 }},
                 {{ 0.000, 1.000 }}
            }}
        };
        static const MultivariateGaussianConstraintEntry<2> LambdaB_to_Lambda_flongA_13dot5_to_20dot3_BFvD2014
        {
            "Lambda_b->Lambda::f_long^A@BFvD2014",
            {{ "Lambda_b->Lambda::f_long^A(s)", "Lambda_b->Lambda::f_long^A(s)" }},
            {{ Kinematics{ { "s", 13.5 } }, Kinematics{ { "s", 20.5 } } }},
            {{ Options{ { "form-factors", "BFvD2014" } }, Options{ { "form-factors", "BFvD2014" } } }},
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

        ///@name 2015
        ///@{
        /*
         * B -> K* LCSR form factors from [BSZ2015]
         *
         * Central values, uncertainties and correlation matrix as extracted
         * from  the ancillary file http://arxiv.org/src/1503.05534v2/anc/BKstar_LCSR.json.
         */
        static const MultivariateGaussianConstraintEntry<11> B_to_Kstar_V_A0_A1_A2_0dot1_to_12dot1_BSZ_2015_lcsr_parameters
        {
            "B->K^*::V+A_0+A_1+A_2[LCSR]@BSZ2015",
            {{
                "B->K^*::alpha^A0_0@BSZ2015",  "B->K^*::alpha^A0_1@BSZ2015",  "B->K^*::alpha^A0_2@BSZ2015",
                "B->K^*::alpha^A1_0@BSZ2015",  "B->K^*::alpha^A1_1@BSZ2015",  "B->K^*::alpha^A1_2@BSZ2015",
                                               "B->K^*::alpha^A12_1@BSZ2015", "B->K^*::alpha^A12_2@BSZ2015",
                "B->K^*::alpha^V_0@BSZ2015",   "B->K^*::alpha^V_1@BSZ2015",   "B->K^*::alpha^V_2@BSZ2015"
            }},
            {{ Kinematics{{ }}, Kinematics{{ }}, Kinematics{{ }},
               Kinematics{{ }}, Kinematics{{ }}, Kinematics{{ }},
               /* A12_0 */      Kinematics{{ }}, Kinematics{{ }},
               Kinematics{{ }}, Kinematics{{ }}, Kinematics{{ }}
            }},
            {{
               Options{}, Options{}, Options{},
               Options{}, Options{}, Options{},
               /*A12_0*/  Options{}, Options{},
               Options{}, Options{}, Options{}
            }},
            {{
                0.355851, -1.04363,   1.12403,
                0.269264,  0.304578, -0.10662,
                           0.601902,  0.117626,
                0.341428, -1.04834,   2.37143
            }},
            {{
                0.0462313, 0.268857, 1.35295,
                0.0293981, 0.189926, 0.479212,
                           0.203756, 0.839054,
                0.0360552, 0.240855, 1.38825
            }},
            {{
                0.0462313, 0.268857, 1.35295,
                0.0293981, 0.189926, 0.479212,
                           0.203756, 0.839054,
                0.0360552, 0.240855, 1.38825
            }},
            {{
                 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0,
                      0.0, 0.0,
                 0.0, 0.0, 0.0,
            }}, // we assign no systematic uncertainty
            {{
                 {{1., 0.317804, -0.357411, 0.07425, -0.0776887, 0.0472153, 0.546571, 0.35966, 0.0487672, -0.0813098, 0.124826}},
                 {{0.317804, 1., -0.453809, 0.0123412, 0.179589, 0.143343, 0.813819, 0.802623, -0.0419978, 0.254118, -0.0000121958}},
                 {{-0.357411, -0.453809, 1., -0.0405148, -0.0558461, 0.00366388, -0.435518, -0.280068, -0.0268353, -0.101528, 0.235726}},
                 {{0.07425, 0.0123412, -0.0405148, 1., 0.466474, 0.280936, -0.17383, -0.271301, 0.919697, 0.467005, -0.485362}},
                 {{-0.0776887, 0.179589, -0.0558461, 0.466474, 1., 0.717856, 0.123579, -0.156642, 0.454858, 0.965108, -0.526988}},
                 {{0.0472153, 0.143343, 0.00366388, 0.280936, 0.717856, 1., 0.121826, 0.00220348, 0.251687, 0.752819, -0.0991215}},
                 {{0.546571, 0.813819, -0.435518, -0.17383, 0.123579, 0.121826, 1., 0.830963, -0.207207, 0.208704, -0.0424022}},
                 {{0.35966, 0.802623, -0.280068, -0.271301, -0.156642, 0.00220348, 0.830963, 1., -0.311073, -0.0367557, 0.116866}},
                 {{0.0487672, -0.0419978, -0.0268353, 0.919697, 0.454858, 0.251687, -0.207207, -0.311073, 1., 0.500818, -0.569767}},
                 {{-0.0813098, 0.254118, -0.101528, 0.467005, 0.965108, 0.752819, 0.208704, -0.0367557, 0.500818, 1., -0.582944}},
                 {{0.124826, -0.0000121958, 0.235726, -0.485362, -0.526988, -0.0991215, -0.0424022, 0.116866, -0.569767, -0.582944, 1.}}
            }},
            0u,
        };

        /*
         * B -> K* LCSR form factors from [BSZ2015]
         *
         * Central values, uncertainties and correlation matrix as extracted
         * from  the ancillary file http://arxiv.org/src/1503.05534v2/anc/BKstar_LCSR-Lattice.json.
         */
        static const MultivariateGaussianConstraintEntry<11> B_to_Kstar_V_A0_A1_A2_0dot1_to_12dot1_BSZ_2015_lcsrlattice_parameters
        {
            "B->K^*::V+A_0+A_1+A_2[LCSR+Lattice]@BSZ2015",
            {{
                "B->K^*::alpha^A0_0@BSZ2015",  "B->K^*::alpha^A0_1@BSZ2015",  "B->K^*::alpha^A0_2@BSZ2015",
                "B->K^*::alpha^A1_0@BSZ2015",  "B->K^*::alpha^A1_1@BSZ2015",  "B->K^*::alpha^A1_2@BSZ2015",
                                               "B->K^*::alpha^A12_1@BSZ2015", "B->K^*::alpha^A12_2@BSZ2015",
                "B->K^*::alpha^V_0@BSZ2015",   "B->K^*::alpha^V_1@BSZ2015",   "B->K^*::alpha^V_2@BSZ2015"
            }},
            {{ Kinematics{{ }}, Kinematics{{ }}, Kinematics{{ }},
               Kinematics{{ }}, Kinematics{{ }}, Kinematics{{ }},
               /* A12_0 */      Kinematics{{ }}, Kinematics{{ }},
               Kinematics{{ }}, Kinematics{{ }}, Kinematics{{ }}
            }},
            {{
               Options{}, Options{}, Options{},
               Options{}, Options{}, Options{},
               /*A12_0*/  Options{}, Options{},
               Options{}, Options{}, Options{}
            }},
            {{
                0.369196, -1.36584,  0.128191,
                0.29725,   0.392378, 1.18916,
                /*A_12*/   0.533638, 0.483166,
                0.376313, -1.16597, 2.42443
            }},
            {{
                0.0289419, 0.256849, 1.63438,
                0.026356,  0.187894, 1.02531,
                           0.128772, 0.656273,
                0.0332944, 0.261268, 1.53102
            }},
            {{
                0.0289419, 0.256849, 1.63438,
                0.026356,  0.187894, 1.02531,
                           0.128772, 0.656273,
                0.0332944, 0.261268, 1.53102
            }},
            {{
                 0.0, 0.0, 0.0,
                 0.0, 0.0, 0.0,
                      0.0, 0.0,
                 0.0, 0.0, 0.0,
            }}, // we assign no systematic uncertainty
            {{
                {{1., 0.633689, 0.0575305, 0.0615373, 0.169044, 0.240647, 0.609495, 0.0355029, 0.0393316, 0.178168, 0.195421}},
                {{0.633689, 1., 0.486314,  0.300568, 0.472167, 0.321569, 0.756564, 0.46913, 0.280475, 0.553833, 0.0752289}},
                {{0.0575305, 0.486314, 1., 0.515475, 0.743538, 0.661689, 0.326536, 0.665149, 0.545574, 0.731129, 0.00860747}},
                {{0.0615373, 0.300568, 0.515475, 1., 0.702494, 0.153604, 0.0105278, 0.205325, 0.922982, 0.731071, -0.35833}},
                {{0.169044, 0.472167, 0.743538, 0.702494, 1., 0.747682, 0.284455, 0.480939, 0.678427, 0.862475, -0.0633373}},
                {{0.240647, 0.321569, 0.661689, 0.153604, 0.747682, 1., 0.316222, 0.535592, 0.194572, 0.517205, 0.253395}},
                {{0.609495, 0.756564, 0.326536, 0.0105278, 0.284455, 0.316222, 1., 0.700096, -0.0140921, 0.350603, 0.234661}},
                {{0.0355029, 0.46913, 0.665149, 0.205325, 0.480939, 0.535592, 0.700096, 1., 0.196602, 0.485841, 0.139392}},
                {{0.0393316, 0.280475, 0.545574, 0.922982, 0.678427, 0.194572, -0.0140921, 0.196602, 1., 0.757379, -0.397005}},
                {{0.178168, 0.553833, 0.731129, 0.731071, 0.862475, 0.517205, 0.350603, 0.485841, 0.757379, 1., 0.0346143}},
                {{0.195421, 0.0752289, 0.00860747, -0.35833, -0.0633373, 0.253395, 0.234661, 0.139392, -0.397005, 0.0346143, 1.}}
            }},
            0u,
        };
        /*
         * B -> K* lattice form factors from [HLMW2015]
         *
         * We got input files from Matthew Wingate with form factors extrapolated to the physical limit:
         * a) the correlation between (V, A0, A1, A12) and (T1, T2, T23) is neglected;
         *
         * b) we add 5 % systematic uncertainty relative to the mean. This induces extra
         * positive correlation with the two set of FF, we assume maximum correlation
         * rho=1, so the covariance matrix is updated as \Sigma_{ij} \to \Sigma_{ij} +
         * 0.05 * \sigma_i \sigma_j;
         *
         * c) Mathew's data includes extrapolation of the q^2 value to the physical
         * value. We neglect any uncertainty from this procedure.
         */
        static const MultivariateGaussianCovarianceConstraintEntry<48> B_to_Kstar_V_A0_A1_A12_11dot9_to_17dot8_HLMW_2015
        {
            "B->K^*::V+A_0+A_1+A_12@HLMW2015",
            {{ "B->K^*::V(s)", "B->K^*::V(s)", "B->K^*::V(s)", "B->K^*::V(s)", "B->K^*::V(s)", "B->K^*::V(s)",
               "B->K^*::V(s)", "B->K^*::V(s)", "B->K^*::V(s)", "B->K^*::V(s)", "B->K^*::V(s)", "B->K^*::V(s)",
               "B->K^*::A_1(s)", "B->K^*::A_1(s)", "B->K^*::A_1(s)", "B->K^*::A_1(s)", "B->K^*::A_1(s)", "B->K^*::A_1(s)",
               "B->K^*::A_1(s)", "B->K^*::A_1(s)", "B->K^*::A_1(s)", "B->K^*::A_1(s)", "B->K^*::A_1(s)", "B->K^*::A_1(s)",
               "B->K^*::A_12(s)", "B->K^*::A_12(s)", "B->K^*::A_12(s)", "B->K^*::A_12(s)", "B->K^*::A_12(s)", "B->K^*::A_12(s)",
               "B->K^*::A_12(s)", "B->K^*::A_12(s)", "B->K^*::A_12(s)", "B->K^*::A_12(s)", "B->K^*::A_12(s)", "B->K^*::A_12(s)",
               "B->K^*::A_0(s)", "B->K^*::A_0(s)", "B->K^*::A_0(s)", "B->K^*::A_0(s)", "B->K^*::A_0(s)", "B->K^*::A_0(s)",
               "B->K^*::A_0(s)", "B->K^*::A_0(s)", "B->K^*::A_0(s)", "B->K^*::A_0(s)", "B->K^*::A_0(s)", "B->K^*::A_0(s)", }},
            {{
               Kinematics{{ "s", 16.5179098782 }}, Kinematics{{ "s", 14.9296311557 }}, Kinematics{{ "s", 13.9898947981 }},
               Kinematics{{ "s", 12.9663104802 }}, Kinematics{{ "s", 16.2551930084 }}, Kinematics{{ "s", 14.9875855712 }},
               Kinematics{{ "s", 13.9849812172 }}, Kinematics{{ "s", 12.8588109331 }}, Kinematics{{ "s", 15.6490320061 }},
               Kinematics{{ "s", 14.2961458423 }}, Kinematics{{ "s", 13.2122218257 }}, Kinematics{{ "s", 11.9180885363 }},
               Kinematics{{ "s", 17.7671223270 }},  Kinematics{{ "s", 16.5179098782 }}, Kinematics{{ "s", 14.9296311557 }},
               Kinematics{{ "s", 12.9663104802 }}, Kinematics{{ "s", 17.7899990577 }}, Kinematics{{ "s", 16.2551930084 }},
               Kinematics{{ "s", 14.9875855712 }}, Kinematics{{ "s", 12.8588109331 }}, Kinematics{{ "s", 17.0322061886 }},
               Kinematics{{ "s", 15.6490320061 }}, Kinematics{{ "s", 14.2961458423 }}, Kinematics{{ "s", 11.9180885363 }},
               Kinematics{{ "s", 16.5179098782 }}, Kinematics{{ "s", 14.9296311557 }}, Kinematics{{ "s", 13.9898947981 }},
               Kinematics{{ "s", 12.9663104802 }}, Kinematics{{ "s", 16.2551930084 }}, Kinematics{{ "s", 14.9875855712 }},
               Kinematics{{ "s", 13.9849812172 }}, Kinematics{{ "s", 12.8588109331 }}, Kinematics{{ "s", 15.6490320061 }},
               Kinematics{{ "s", 14.2961458423 }}, Kinematics{{ "s", 13.2122218257 }}, Kinematics{{ "s", 11.9180885363 }},
               Kinematics{{ "s", 16.5179098782 }}, Kinematics{{ "s", 14.9296311557 }}, Kinematics{{ "s", 13.9898947981 }},
               Kinematics{{ "s", 12.9663104802 }}, Kinematics{{ "s", 16.2551930084 }}, Kinematics{{ "s", 14.9875855712 }},
               Kinematics{{ "s", 13.9849812172 }}, Kinematics{{ "s", 12.8588109331 }}, Kinematics{{ "s", 15.6490320061 }},
               Kinematics{{ "s", 14.2961458423 }}, Kinematics{{ "s", 13.2122218257 }}, Kinematics{{ "s", 11.9180885363 }},  }},
         // no options; default ctor just fine
         {{
              Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } },
              Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } },
              Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } },
              Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } },
              Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } },
              Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } },
              Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } },
              Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } },
              Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } },
              Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } },
              Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } },
              Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } },
              Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } },
              Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } },
              Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } },
              Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } }, Options{ { "form-factors", "BSZ2015" } }
         }},
         {{
              1.41654756,  1.19402316,  1.06139874,  0.98029398,  1.41456527,  1.38416036,  1.3568178 ,  1.0031494 ,
              1.30654947,  1.10920364,  1.03058181,  0.72171754,  0.56552245,  0.4972255 ,  0.46514096,  0.45514081,
              0.54684044,  0.52813813,  0.51173916,  0.40282807,  0.54058157,  0.51625489,  0.47855314,  0.41713015,
              0.37980242,  0.34367319,  0.34415972,  0.33971277,  0.34482022,  0.42729289,  0.28990601,  0.31845106,
              0.37294807,  0.33592982,  0.34252166,  0.34954018,  1.40033868,  1.22824918,  1.17492672,  0.95528844,
              1.21238441,  1.30580617,  0.98476889,  0.84125097,  1.35984599,  1.02937473,  1.0462483 ,  0.89997334 }},
         {{
       {{ 1.16627106e-02,   7.49153547e-03,   6.45838440e-03,   8.15007038e-03,   6.55522244e-03,   5.21724930e-03,
          5.81927477e-03,   5.32914011e-03,   1.10130135e-02,   9.79601309e-03,   8.96424430e-03,   6.79804829e-03,
          2.00272362e-03,   1.76085892e-03,   1.64723575e-03,   1.61182151e-03,   1.93656372e-03,   1.87033196e-03,
          1.81225715e-03,   1.42656280e-03,   1.91439876e-03,   1.82824903e-03,   1.69473322e-03,   1.47721175e-03,
          1.34502050e-03,   1.21707356e-03,   1.21879652e-03,   1.20304825e-03,   1.22113562e-03,   1.51320174e-03,
          1.02666411e-03,   1.12775268e-03,   1.32074669e-03,   1.18965143e-03,   1.21299554e-03,   1.23785073e-03,
          4.95911584e-03,   4.34968347e-03,   4.16084895e-03,   3.38302876e-03,   4.29350045e-03,   4.62434137e-03,
          3.48742991e-03,   2.97918002e-03,   4.81571629e-03,   3.64539565e-03,   3.70515119e-03,   3.18713761e-03}},
       {{ 7.49153547e-03,   7.96299365e-03,   7.12459730e-03,   6.27992088e-03,   6.19790293e-03,   5.74700092e-03,
          6.68600916e-03,   5.61712859e-03,   8.94093087e-03,   7.42951189e-03,   7.66363426e-03,   5.49036526e-03,
          1.68811726e-03,   1.48424690e-03,   1.38847271e-03,   1.35862166e-03,   1.63235036e-03,   1.57652290e-03,
          1.52757102e-03,   1.20246511e-03,   1.61366728e-03,   1.54105075e-03,   1.42850884e-03,   1.24515765e-03,
          1.13373223e-03,   1.02588438e-03,   1.02733668e-03,   1.01406229e-03,   1.02930833e-03,   1.27549400e-03,
          8.65386208e-04,   9.50594845e-04,   1.11327157e-03,   1.00276997e-03,   1.02244697e-03,   1.04339768e-03,
          4.18009202e-03,   3.66639492e-03,   3.50722428e-03,   2.85159129e-03,   3.61903765e-03,   3.89790701e-03,
          2.93959213e-03,   2.51118284e-03,   4.05921900e-03,   3.07274315e-03,   3.12311174e-03,   2.68647253e-03}},
       {{ 6.45838440e-03,   7.12459730e-03,   8.81059561e-03,   5.75039374e-03,   5.38962278e-03,   5.27100095e-03,
          6.50638686e-03,   5.51867632e-03,   8.04887956e-03,   7.59889340e-03,   6.48615802e-03,   6.18843541e-03,
          1.50061204e-03,   1.31938629e-03,   1.23425008e-03,   1.20771470e-03,   1.45103937e-03,   1.40141287e-03,
          1.35789825e-03,   1.06890301e-03,   1.43443149e-03,   1.36988073e-03,   1.26983926e-03,   1.10685354e-03,
          1.00780453e-03,   9.11935732e-04,   9.13226720e-04,   9.01426766e-04,   9.14979375e-04,   1.13382032e-03,
          7.69264668e-04,   8.45008877e-04,   9.89616517e-04,   8.91388721e-04,   9.08880133e-04,   9.27503766e-04,
          3.71579425e-03,   3.25915533e-03,   3.11766434e-03,   2.53485485e-03,   3.21705820e-03,   3.46495255e-03,
          2.61308112e-03,   2.23225678e-03,   3.60834703e-03,   2.73144259e-03,   2.77621655e-03,   2.38807642e-03}},
       {{ 8.15007038e-03,   6.27992088e-03,   5.75039374e-03,   2.28639916e-02,   5.82493292e-03,   3.53803685e-03,
          3.76667462e-03,   7.43132038e-03,   8.00208318e-03,   4.02834482e-03,   6.64685465e-03,   4.47052703e-03,
          1.38594564e-03,   1.21856792e-03,   1.13993723e-03,   1.11542949e-03,   1.34016098e-03,   1.29432659e-03,
          1.25413705e-03,   9.87224834e-04,   1.32482215e-03,   1.26520392e-03,   1.17280692e-03,   1.02227545e-03,
          9.30795081e-04,   8.42251910e-04,   8.43444250e-04,   8.32545967e-04,   8.45062980e-04,   1.04718161e-03,
          7.10482783e-04,   7.80439143e-04,   9.13996867e-04,   8.23274959e-04,   8.39429799e-04,   8.56630343e-04,
          3.43185895e-03,   3.01011322e-03,   2.87943399e-03,   2.34115877e-03,   2.97123286e-03,   3.20018484e-03,
          2.41340754e-03,   2.06168316e-03,   3.33262210e-03,   2.52272463e-03,   2.56407728e-03,   2.20559614e-03}},
       {{ 6.55522244e-03,   6.19790293e-03,   5.38962278e-03,   5.82493292e-03,   1.51327700e-02,   1.18107593e-02,
          8.69859305e-03,   8.74551276e-03,   9.35557411e-03,   7.51380725e-03,   8.38428241e-03,   7.78812179e-03,
          1.99992104e-03,   1.75839480e-03,   1.64493063e-03,   1.60956595e-03,   1.93385372e-03,   1.86771465e-03,
          1.80972111e-03,   1.42456649e-03,   1.91171978e-03,   1.82569061e-03,   1.69236164e-03,   1.47514456e-03,
          1.34313830e-03,   1.21537041e-03,   1.21709095e-03,   1.20136472e-03,   1.21942678e-03,   1.51108419e-03,
          1.02522741e-03,   1.12617452e-03,   1.31889845e-03,   1.18798665e-03,   1.21129810e-03,   1.23611850e-03,
          4.95217613e-03,   4.34359658e-03,   4.15502633e-03,   3.37829460e-03,   4.28749219e-03,   4.61787013e-03,
          3.48254965e-03,   2.97501100e-03,   4.80897725e-03,   3.64029434e-03,   3.69996625e-03,   3.18267758e-03}},
       {{ 5.21724930e-03,   5.74700092e-03,   5.27100095e-03,   3.53803685e-03,   1.18107593e-02,   1.96905962e-02,
          1.16325353e-02,   6.37505422e-03,   8.31954749e-03,   5.66319186e-03,   5.93629621e-03,   5.08714414e-03,
          1.95693440e-03,   1.72059957e-03,   1.60957422e-03,   1.57496967e-03,   1.89228715e-03,   1.82756968e-03,
          1.77082266e-03,   1.39394662e-03,   1.87062896e-03,   1.78644891e-03,   1.65598574e-03,   1.44343756e-03,
          1.31426866e-03,   1.18924703e-03,   1.19093060e-03,   1.17554238e-03,   1.19321622e-03,   1.47860469e-03,
          1.00319100e-03,   1.10196833e-03,   1.29054983e-03,   1.16245186e-03,   1.18526225e-03,   1.20954916e-03,
          4.84573323e-03,   4.25023459e-03,   4.06571749e-03,   3.30568097e-03,   4.19533611e-03,   4.51861286e-03,
          3.40769515e-03,   2.91106561e-03,   4.70561229e-03,   3.56204924e-03,   3.62043856e-03,   3.11426858e-03}},
       {{ 5.81927477e-03,   6.68600916e-03,   6.50638686e-03,   3.76667462e-03,   8.69859305e-03,   1.16325353e-02,
          2.47754589e-02,   9.00544376e-03,   6.56037367e-03,   4.54395419e-03,   6.04464713e-03,   5.92041532e-03,
          1.91827732e-03,   1.68661102e-03,   1.57777885e-03,   1.54385787e-03,   1.85490710e-03,   1.79146805e-03,
          1.73584200e-03,   1.36641073e-03,   1.83367674e-03,   1.75115957e-03,   1.62327356e-03,   1.41492403e-03,
          1.28830672e-03,   1.16575476e-03,   1.16740507e-03,   1.15232084e-03,   1.16964554e-03,   1.44939648e-03,
          9.83374068e-04,   1.08020016e-03,   1.26505644e-03,   1.13948890e-03,   1.16184870e-03,   1.18565585e-03,
          4.75001110e-03,   4.16627588e-03,   3.98540371e-03,   3.24038088e-03,   4.11246186e-03,   4.42935263e-03,
          3.34037987e-03,   2.85356071e-03,   4.61265809e-03,   3.49168488e-03,   3.54892077e-03,   3.05274962e-03}},
       {{ 5.32914011e-03,   5.61712859e-03,   5.51867632e-03,   7.43132038e-03,   8.74551276e-03,   6.37505422e-03,
          9.00544376e-03,   2.76107295e-02,   7.78861971e-03,   8.55657526e-03,   5.49022011e-03,   2.73002423e-03,
          1.41825877e-03,   1.24697866e-03,   1.16651470e-03,   1.14143557e-03,   1.37140664e-03,   1.32450363e-03,
          1.28337708e-03,   1.01024184e-03,   1.35571020e-03,   1.29470197e-03,   1.20015075e-03,   1.04610966e-03,
          9.52496438e-04,   8.61888895e-04,   8.63109035e-04,   8.51956660e-04,   8.64765504e-04,   1.07159651e-03,
          7.27047589e-04,   7.98634972e-04,   9.35306575e-04,   8.42469500e-04,   8.59000987e-04,   8.76602559e-04,
          3.51187226e-03,   3.08029358e-03,   2.94656759e-03,   2.39574256e-03,   3.04050674e-03,   3.27479670e-03,
          2.46967580e-03,   2.10975101e-03,   3.41032172e-03,   2.58154161e-03,   2.62385838e-03,   2.25701930e-03}},
       {{ 1.10130135e-02,   8.94093087e-03,   8.04887956e-03,   8.00208318e-03,   9.35557411e-03,   8.31954749e-03,
          6.56037367e-03,   7.78861971e-03,   2.61871870e-02,   2.16239884e-02,   1.90189039e-02,   1.72815590e-02,
          1.84720764e-03,   1.62412428e-03,   1.51932420e-03,   1.48665995e-03,   1.78618521e-03,   1.72509649e-03,
          1.67153132e-03,   1.31578700e-03,   1.76574141e-03,   1.68628139e-03,   1.56313339e-03,   1.36250294e-03,
          1.24057664e-03,   1.12256507e-03,   1.12415424e-03,   1.10962885e-03,   1.12631170e-03,   1.39569823e-03,
          9.46941342e-04,   1.04018015e-03,   1.21818775e-03,   1.09727233e-03,   1.11880372e-03,   1.14172884e-03,
          4.57402938e-03,   4.01192079e-03,   3.83774971e-03,   3.12032900e-03,   3.96010051e-03,   4.26525089e-03,
          3.21662316e-03,   2.74784001e-03,   4.44176513e-03,   3.36232251e-03,   3.41743789e-03,   2.93964923e-03}},
       {{ 9.79601309e-03,   7.42951189e-03,   7.59889340e-03,   4.02834482e-03,   7.51380725e-03,   5.66319186e-03,
          4.54395419e-03,   8.55657526e-03,   2.16239884e-02,   3.01068213e-02,   1.71198873e-02,   1.49421281e-02,
          1.56819890e-03,   1.37881083e-03,   1.28984012e-03,   1.26210960e-03,   1.51639350e-03,   1.46453184e-03,
          1.41905734e-03,   1.11704590e-03,   1.49903761e-03,   1.43157952e-03,   1.32703222e-03,   1.15670570e-03,
          1.05319558e-03,   9.53008889e-04,   9.54358022e-04,   9.42026604e-04,   9.56189616e-04,   1.18488706e-03,
          8.03911987e-04,   8.83067680e-04,   1.03418838e-03,   9.31536450e-04,   9.49815667e-04,   9.69278099e-04,
          3.88315188e-03,   3.40594615e-03,   3.25808248e-03,   2.64902352e-03,   3.36195299e-03,   3.62101238e-03,
          2.73077307e-03,   2.33279658e-03,   3.77086528e-03,   2.85446548e-03,   2.90125604e-03,   2.49563426e-03}},
       {{ 8.96424430e-03,   7.66363426e-03,   6.48615802e-03,   6.64685465e-03,   8.38428241e-03,   5.93629621e-03,
          6.04464713e-03,   5.49022011e-03,   1.90189039e-02,   1.71198873e-02,   2.01116134e-02,   1.32479150e-02,
          1.45704289e-03,   1.28107889e-03,   1.19841455e-03,   1.17264960e-03,   1.40890953e-03,   1.36072389e-03,
          1.31847268e-03,   1.03786821e-03,   1.39278384e-03,   1.33010727e-03,   1.23297042e-03,   1.07471687e-03,
          9.78543681e-04,   8.85458358e-04,   8.86711864e-04,   8.75254513e-04,   8.88413632e-04,   1.10090069e-03,
          7.46929642e-04,   8.20474675e-04,   9.60883739e-04,   8.65507914e-04,   8.82491476e-04,   9.00574386e-04,
          3.60790894e-03,   3.16452818e-03,   3.02714528e-03,   2.46125723e-03,   3.12365331e-03,   3.36435023e-03,
          2.53721226e-03,   2.16744487e-03,   3.50358136e-03,   2.65213719e-03,   2.69561117e-03,   2.31874040e-03}},
       {{ 6.79804829e-03,   5.49036526e-03,   6.18843541e-03,   4.47052703e-03,   7.78812179e-03,   5.08714414e-03,
          5.92041532e-03,   2.73002423e-03,   1.72815590e-02,   1.49421281e-02,   1.32479150e-02,   3.91671013e-02,
          1.02036868e-03,   8.97140909e-04,   8.39250981e-04,   8.21207760e-04,   9.86660837e-04,   9.52916383e-04,
          9.23327819e-04,   7.26820207e-04,   9.75368000e-04,   9.31475529e-04,   8.63450493e-04,   7.52625366e-04,
          6.85275178e-04,   6.20087428e-04,   6.20965259e-04,   6.12941664e-04,   6.22157009e-04,   7.70961924e-04,
          5.23075621e-04,   5.74579285e-04,   6.72907902e-04,   6.06116110e-04,   6.18009717e-04,   6.30673198e-04,
          2.52662246e-03,   2.21612244e-03,   2.11991305e-03,   1.72362105e-03,   2.18749773e-03,   2.35605804e-03,
          1.77681244e-03,   1.51786394e-03,   2.45356175e-03,   1.85729449e-03,   1.88773936e-03,   1.62381637e-03}},
       {{ 2.00272362e-03,   1.68811726e-03,   1.50061204e-03,   1.38594564e-03,   1.99992104e-03,   1.95693440e-03,
          1.91827732e-03,   1.41825877e-03,   1.84720764e-03,   1.56819890e-03,   1.45704289e-03,   1.02036868e-03,
          1.05033096e-03,   9.47651027e-04,   8.32132830e-04,   7.60606608e-04,   8.33584901e-04,   7.74835597e-04,
          7.93930031e-04,   6.63219191e-04,   9.71372459e-04,   9.19579521e-04,   8.22721346e-04,   8.03966624e-04,
          5.36966995e-04,   4.85887267e-04,   4.86575116e-04,   4.80287999e-04,   4.87508946e-04,   6.04109300e-04,
          4.09870887e-04,   4.50228058e-04,   5.27276262e-04,   4.74939641e-04,   4.84259217e-04,   4.94182050e-04,
          1.97980740e-03,   1.73650622e-03,   1.66111860e-03,   1.35059265e-03,   1.71407651e-03,   1.84615677e-03,
          1.39227228e-03,   1.18936577e-03,   1.92255859e-03,   1.45533630e-03,   1.47919225e-03,   1.27238783e-03}},
       {{ 1.76085892e-03,   1.48424690e-03,   1.31938629e-03,   1.21856792e-03,   1.75839480e-03,   1.72059957e-03,
          1.68661102e-03,   1.24697866e-03,   1.62412428e-03,   1.37881083e-03,   1.28107889e-03,   8.97140909e-04,
          9.47651027e-04,   1.45432384e-03,   1.08504506e-03,   9.71088387e-04,   7.60269530e-04,   6.88964942e-04,
          8.56791039e-04,   5.68351501e-04,   8.76328386e-04,   7.97295160e-04,   6.49094182e-04,   6.56763984e-04,
          4.72118626e-04,   4.27207688e-04,   4.27812467e-04,   4.22284632e-04,   4.28633521e-04,   5.31152296e-04,
          3.60371645e-04,   3.95854966e-04,   4.63598222e-04,   4.17582184e-04,   4.25776254e-04,   4.34500728e-04,
          1.74071024e-03,   1.52679203e-03,   1.46050881e-03,   1.18748442e-03,   1.50707111e-03,   1.62320031e-03,
          1.22413050e-03,   1.04572858e-03,   1.69037525e-03,   1.27957841e-03,   1.30055333e-03,   1.11872424e-03}},
       {{ 1.64723575e-03,   1.38847271e-03,   1.23425008e-03,   1.13993723e-03,   1.64493063e-03,   1.60957422e-03,
          1.57777885e-03,   1.16651470e-03,   1.51932420e-03,   1.28984012e-03,   1.19841455e-03,   8.39250981e-04,
          8.32132830e-04,   1.08504506e-03,   1.69953051e-03,   9.07113975e-04,   7.18386601e-04,   6.86671338e-04,
          6.40331006e-04,   5.78498872e-04,   7.82851370e-04,   7.25792847e-04,   6.33879008e-04,   4.89937091e-04,
          4.41654166e-04,   3.99641202e-04,   4.00206957e-04,   3.95035816e-04,   4.00975030e-04,   4.96878563e-04,
          3.37117897e-04,   3.70311582e-04,   4.33683559e-04,   3.90636804e-04,   3.98302135e-04,   4.06463643e-04,
          1.62838721e-03,   1.42827253e-03,   1.36626637e-03,   1.11085946e-03,   1.40982413e-03,   1.51845986e-03,
          1.14514087e-03,   9.78250717e-04,   1.58130019e-03,   1.19701089e-03,   1.21663236e-03,   1.04653617e-03}},
       {{ 1.61182151e-03,   1.35862166e-03,   1.20771470e-03,   1.11542949e-03,   1.60956595e-03,   1.57496967e-03,
          1.54385787e-03,   1.14143557e-03,   1.48665995e-03,   1.26210960e-03,   1.17264960e-03,   8.21207760e-04,
          7.60606608e-04,   9.71088387e-04,   9.07113975e-04,   3.49877827e-03,   5.53739752e-04,   4.57241916e-04,
          6.20183364e-04,   6.98431016e-04,   6.08299927e-04,   7.19901822e-04,   7.10451341e-04,   5.35225527e-04,
          4.32158957e-04,   3.91049238e-04,   3.91602829e-04,   3.86542864e-04,   3.92354389e-04,   4.86196074e-04,
          3.29870134e-04,   3.62350181e-04,   4.24359712e-04,   3.82238427e-04,   3.89738959e-04,   3.97725002e-04,
          1.59337819e-03,   1.39756582e-03,   1.33689274e-03,   1.08697688e-03,   1.37951405e-03,   1.48581419e-03,
          1.12052127e-03,   9.57219114e-04,   1.54730351e-03,   1.17127611e-03,   1.19047574e-03,   1.02403649e-03}},
       {{ 1.93656372e-03,   1.63235036e-03,   1.45103937e-03,   1.34016098e-03,   1.93385372e-03,   1.89228715e-03,
          1.85490710e-03,   1.37140664e-03,   1.78618521e-03,   1.51639350e-03,   1.40890953e-03,   9.86660837e-04,
          8.33584901e-04,   7.60269530e-04,   7.18386601e-04,   5.53739752e-04,   1.52849980e-03,   9.83620792e-04,
          7.98642156e-04,   4.33500450e-04,   8.90977070e-04,   9.11169739e-04,   8.80551105e-04,   1.02269121e-03,
          5.19228311e-04,   4.69835999e-04,   4.70501125e-04,   4.64421702e-04,   4.71404106e-04,   5.84152571e-04,
          3.96330817e-04,   4.35354791e-04,   5.09857710e-04,   4.59250027e-04,   4.68261731e-04,   4.77856764e-04,
          1.91440454e-03,   1.67914080e-03,   1.60624361e-03,   1.30597587e-03,   1.65745205e-03,   1.78516905e-03,
          1.34627862e-03,   1.15007512e-03,   1.85904694e-03,   1.40725932e-03,   1.43032719e-03,   1.23035454e-03}},
       {{ 1.87033196e-03,   1.57652290e-03,   1.40141287e-03,   1.29432659e-03,   1.86771465e-03,   1.82756968e-03,
          1.79146805e-03,   1.32450363e-03,   1.72509649e-03,   1.46453184e-03,   1.36072389e-03,   9.52916383e-04,
          7.74835597e-04,   6.88964942e-04,   6.86671338e-04,   4.57241916e-04,   9.83620792e-04,   1.60625073e-03,
          1.10058093e-03,   8.03549827e-04,   8.74166680e-04,   8.15572225e-04,   9.57805085e-04,   8.75448896e-04,
          5.01470359e-04,   4.53767297e-04,   4.54409675e-04,   4.48538173e-04,   4.55281774e-04,   5.64174167e-04,
          3.82776041e-04,   4.20465369e-04,   4.92420239e-04,   4.43543372e-04,   4.52246870e-04,   4.61513747e-04,
          1.84893063e-03,   1.62171308e-03,   1.55130901e-03,   1.26131063e-03,   1.60076610e-03,   1.72411508e-03,
          1.30023500e-03,   1.11074179e-03,   1.79546630e-03,   1.35913012e-03,   1.38140906e-03,   1.18827560e-03}},
       {{ 1.81225715e-03,   1.52757102e-03,   1.35789825e-03,   1.25413705e-03,   1.80972111e-03,   1.77082266e-03,
          1.73584200e-03,   1.28337708e-03,   1.67153132e-03,   1.41905734e-03,   1.31847268e-03,   9.23327819e-04,
          7.93930031e-04,   8.56791039e-04,   6.40331006e-04,   6.20183364e-04,   7.98642156e-04,   1.10058093e-03,
          1.41464015e-03,   9.53790256e-04,   8.56684261e-04,   7.89427904e-04,   8.03500320e-04,   7.76700961e-04,
          4.85899435e-04,   4.39677579e-04,   4.40300011e-04,   4.34610822e-04,   4.41145031e-04,   5.46656257e-04,
          3.70890639e-04,   4.07409694e-04,   4.77130327e-04,   4.29771113e-04,   4.38204363e-04,   4.47183497e-04,
          1.79152035e-03,   1.57135802e-03,   1.50314004e-03,   1.22214626e-03,   1.55106145e-03,   1.67058038e-03,
          1.25986201e-03,   1.07625266e-03,   1.73971611e-03,   1.31692840e-03,   1.33851556e-03,   1.15137901e-03}},
       {{ 1.42656280e-03,   1.20246511e-03,   1.06890301e-03,   9.87224834e-04,   1.42456649e-03,   1.39394662e-03,
          1.36641073e-03,   1.01024184e-03,   1.31578700e-03,   1.11704590e-03,   1.03786821e-03,   7.26820207e-04,
          6.63219191e-04,   5.68351501e-04,   5.78498872e-04,   6.98431016e-04,   4.33500450e-04,   8.03549827e-04,
          9.53790256e-04,   3.43421858e-03,   6.68183621e-04,   6.96162747e-04,   7.73513378e-04,   4.92728641e-04,
          3.82487694e-04,   3.46103023e-04,   3.46592986e-04,   3.42114600e-04,   3.47258163e-04,   4.30313920e-04,
          2.91955691e-04,   3.20702563e-04,   3.75584875e-04,   3.38304904e-04,   3.44943344e-04,   3.52011491e-04,
          1.41023931e-03,   1.23693312e-03,   1.18323366e-03,   9.62042492e-04,   1.22095618e-03,   1.31503845e-03,
          9.91731373e-04,   8.47198758e-04,   1.36946033e-03,   1.03665259e-03,   1.05364545e-03,   9.06336312e-04}},
       {{ 1.91439876e-03,   1.61366728e-03,   1.43443149e-03,   1.32482215e-03,   1.91171978e-03,   1.87062896e-03,
          1.83367674e-03,   1.35571020e-03,   1.76574141e-03,   1.49903761e-03,   1.39278384e-03,   9.75368000e-04,
          9.71372459e-04,   8.76328386e-04,   7.82851370e-04,   6.08299927e-04,   8.90977070e-04,   8.74166680e-04,
          8.56684261e-04,   6.68183621e-04,   1.39854787e-03,   1.25147913e-03,   1.17263191e-03,   9.83628094e-04,
          5.13285477e-04,   4.64458486e-04,   4.65116000e-04,   4.59106159e-04,   4.66008645e-04,   5.77466647e-04,
          3.91794608e-04,   4.30371933e-04,   5.04022129e-04,   4.53993676e-04,   4.62902237e-04,   4.72387449e-04,
          1.89249320e-03,   1.65992218e-03,   1.58785933e-03,   1.29102831e-03,   1.63848167e-03,   1.76473687e-03,
          1.33086977e-03,   1.13691192e-03,   1.83776920e-03,   1.39115252e-03,   1.41395637e-03,   1.21627251e-03}},
       {{ 1.82824903e-03,   1.54105075e-03,   1.36988073e-03,   1.26520392e-03,   1.82569061e-03,   1.78644891e-03,
          1.75115957e-03,   1.29470197e-03,   1.68628139e-03,   1.43157952e-03,   1.33010727e-03,   9.31475529e-04,
          9.19579521e-04,   7.97295160e-04,   7.25792847e-04,   7.19901822e-04,   9.11169739e-04,   8.15572225e-04,
          7.89427904e-04,   6.96162747e-04,   1.25147913e-03,   1.40899169e-03,   1.23860234e-03,   1.13685249e-03,
          4.90187152e-04,   4.43557420e-04,   4.44185346e-04,   4.38445953e-04,   4.45037822e-04,   5.51480109e-04,
          3.74163485e-04,   4.11004794e-04,   4.81340662e-04,   4.33563537e-04,   4.42071204e-04,   4.51129573e-04,
          1.80732924e-03,   1.58522413e-03,   1.51640418e-03,   1.23293083e-03,   1.56474846e-03,   1.68532207e-03,
          1.27097939e-03,   1.08574982e-03,   1.75506787e-03,   1.32854935e-03,   1.35032701e-03,   1.16153911e-03}},
       {{ 1.69473322e-03,   1.42850884e-03,   1.26983926e-03,   1.17280692e-03,   1.69236164e-03,   1.65598574e-03,
          1.62327356e-03,   1.20015075e-03,   1.56313339e-03,   1.32703222e-03,   1.23297042e-03,   8.63450493e-04,
          8.22721346e-04,   6.49094182e-04,   6.33879008e-04,   7.10451341e-04,   8.80551105e-04,   9.57805085e-04,
          8.03500320e-04,   7.73513378e-04,   1.17263191e-03,   1.23860234e-03,   1.78635029e-03,   1.14936947e-03,
          4.54389112e-04,   4.11164718e-04,   4.11746786e-04,   4.06426538e-04,   4.12537007e-04,   5.11205885e-04,
          3.46838576e-04,   3.80989388e-04,   4.46188675e-04,   4.01900681e-04,   4.09787039e-04,   4.18183882e-04,
          1.67534119e-03,   1.46945627e-03,   1.40566219e-03,   1.14289071e-03,   1.45047593e-03,   1.56224412e-03,
          1.17816062e-03,   1.00645824e-03,   1.62689643e-03,   1.23152628e-03,   1.25171353e-03,   1.07671268e-03}},
       {{ 1.47721175e-03,   1.24515765e-03,   1.10685354e-03,   1.02227545e-03,   1.47514456e-03,   1.44343756e-03,
          1.41492403e-03,   1.04610966e-03,   1.36250294e-03,   1.15670570e-03,   1.07471687e-03,   7.52625366e-04,
          8.03966624e-04,   6.56763984e-04,   4.89937091e-04,   5.35225527e-04,   1.02269121e-03,   8.75448896e-04,
          7.76700961e-04,   4.92728641e-04,   9.83628094e-04,   1.13685249e-03,   1.14936947e-03,   3.45963121e-03,
          3.96067608e-04,   3.58391128e-04,   3.58898487e-04,   3.54261100e-04,   3.59587281e-04,   4.45591865e-04,
          3.02321340e-04,   3.32088846e-04,   3.88919709e-04,   3.50316144e-04,   3.57190276e-04,   3.64509372e-04,
          1.46030871e-03,   1.28084942e-03,   1.22524340e-03,   9.96199026e-04,   1.26430523e-03,   1.36172782e-03,
          1.02694199e-03,   8.77277860e-04,   1.41808191e-03,   1.07345809e-03,   1.09105428e-03,   9.38515043e-04}},
       {{ 1.34502050e-03,   1.13373223e-03,   1.00780453e-03,   9.30795081e-04,   1.34313830e-03,   1.31426866e-03,
          1.28830672e-03,   9.52496438e-04,   1.24057664e-03,   1.05319558e-03,   9.78543681e-04,   6.85275178e-04,
          5.36966995e-04,   4.72118626e-04,   4.41654166e-04,   4.32158957e-04,   5.19228311e-04,   5.01470359e-04,
          4.85899435e-04,   3.82487694e-04,   5.13285477e-04,   4.90187152e-04,   4.54389112e-04,   3.96067608e-04,
          7.04849673e-04,   5.46701282e-04,   4.93484802e-04,   7.44829962e-04,   4.37867700e-04,   4.67171881e-04,
          4.56147849e-04,   4.49669426e-04,   7.34567494e-04,   6.12454427e-04,   5.08546766e-04,   6.28601201e-04,
          1.56385469e-03,   1.31293893e-03,   1.28654836e-03,   1.43809711e-03,   1.01443327e-03,   1.27128326e-03,
          8.35335045e-04,   9.15299594e-04,   1.22312381e-03,   8.60802781e-04,   7.63265976e-04,   7.94814777e-04}},
       {{ 1.21707356e-03,   1.02588438e-03,   9.11935732e-04,   8.42251910e-04,   1.21537041e-03,   1.18924703e-03,
          1.16575476e-03,   8.61888895e-04,   1.12256507e-03,   9.53008889e-04,   8.85458358e-04,   6.20087428e-04,
          4.85887267e-04,   4.27207688e-04,   3.99641202e-04,   3.91049238e-04,   4.69835999e-04,   4.53767297e-04,
          4.39677579e-04,   3.46103023e-04,   4.64458486e-04,   4.43557420e-04,   4.11164718e-04,   3.58391128e-04,
          5.46701282e-04,   8.45756086e-04,   5.34228970e-04,   7.91804920e-04,   5.71888753e-04,   5.68230054e-04,
          7.81061250e-04,   3.75531134e-04,   6.07606788e-04,   6.24590634e-04,   5.17497842e-04,   5.33033829e-04,
          1.35645413e-03,   1.66752429e-03,   1.07724973e-03,   1.45066027e-03,   1.31063879e-03,   1.34265288e-03,
          1.31170005e-03,   7.47138541e-04,   1.01714463e-03,   8.94399194e-04,   1.14275093e-03,   8.08201181e-04}},
       {{ 1.21879652e-03,   1.02733668e-03,   9.13226720e-04,   8.43444250e-04,   1.21709095e-03,   1.19093060e-03,
          1.16740507e-03,   8.63109035e-04,   1.12415424e-03,   9.54358022e-04,   8.86711864e-04,   6.20965259e-04,
          4.86575116e-04,   4.27812467e-04,   4.00206957e-04,   3.91602829e-04,   4.70501125e-04,   4.54409675e-04,
          4.40300011e-04,   3.46592986e-04,   4.65116000e-04,   4.44185346e-04,   4.11746786e-04,   3.58898487e-04,
          4.93484802e-04,   5.34228970e-04,   8.35492795e-04,   7.00333288e-04,   3.37784113e-04,   4.12025170e-04,
         -4.24580452e-04,   4.08990045e-04,   6.05495682e-04,   5.88293491e-04,   5.48859822e-04,   5.31118908e-04,
          1.15629002e-03,   1.20158071e-03,   1.87883558e-03,   1.14047363e-03,   1.17349678e-03,   1.00235517e-03,
          2.03437868e-04,   1.26090814e-03,   1.07363243e-03,   9.16279529e-04,   6.88048208e-04,   1.31847608e-03}},
       {{ 1.20304825e-03,   1.01406229e-03,   9.01426766e-04,   8.32545967e-04,   1.20136472e-03,   1.17554238e-03,
          1.15232084e-03,   8.51956660e-04,   1.10962885e-03,   9.42026604e-04,   8.75254513e-04,   6.12941664e-04,
          4.80287999e-04,   4.22284632e-04,   3.95035816e-04,   3.86542864e-04,   4.64421702e-04,   4.48538173e-04,
          4.34610822e-04,   3.42114600e-04,   4.59106159e-04,   4.38445953e-04,   4.06426538e-04,   3.54261100e-04,
          7.44829962e-04,   7.91804920e-04,   7.00333288e-04,   7.56789522e-03,   1.99374349e-04,   1.03150868e-03,
          1.83817346e-04,   6.37348484e-04,   4.21077641e-04,   4.72122327e-04,   2.31408775e-04,   1.24400834e-03,
          1.42279050e-03,   1.43298233e-03,   9.30704025e-04,   8.05921777e-03,   5.45398848e-04,   2.57677852e-03,
         -1.06622020e-04,   7.16501311e-05,   1.49589387e-03,   4.98951942e-04,   1.65300146e-03,   2.28772211e-03}},
       {{ 1.22113562e-03,   1.02930833e-03,   9.14979375e-04,   8.45062980e-04,   1.21942678e-03,   1.19321622e-03,
          1.16964554e-03,   8.64765504e-04,   1.12631170e-03,   9.56189616e-04,   8.88413632e-04,   6.22157009e-04,
          4.87508946e-04,   4.28633521e-04,   4.00975030e-04,   3.92354389e-04,   4.71404106e-04,   4.55281774e-04,
          4.41145031e-04,   3.47258163e-04,   4.66008645e-04,   4.45037822e-04,   4.12537007e-04,   3.59587281e-04,
          4.37867700e-04,   5.71888753e-04,   3.37784113e-04,   1.99374349e-04,   3.18161174e-03,   6.75273607e-04,
          1.07490238e-03,   6.44888647e-04,   6.10649172e-04,   5.77466040e-04,   4.80122632e-04,   6.11664360e-04,
          1.66383290e-03,   1.36807255e-03,   9.62439621e-04,   5.04849264e-04,   3.52613692e-03,   1.00496214e-03,
          1.18660022e-03,   7.42587266e-04,   9.40734522e-04,   9.26615671e-04,   1.36626987e-03,   7.18390101e-04}},
       {{ 1.51320174e-03,   1.27549400e-03,   1.13382032e-03,   1.04718161e-03,   1.51108419e-03,   1.47860469e-03,
          1.44939648e-03,   1.07159651e-03,   1.39569823e-03,   1.18488706e-03,   1.10090069e-03,   7.70961924e-04,
          6.04109300e-04,   5.31152296e-04,   4.96878563e-04,   4.86196074e-04,   5.84152571e-04,   5.64174167e-04,
          5.46656257e-04,   4.30313920e-04,   5.77466647e-04,   5.51480109e-04,   5.11205885e-04,   4.45591865e-04,
          4.67171881e-04,   5.68230054e-04,   4.12025170e-04,   1.03150868e-03,   6.75273607e-04,   4.80366347e-03,
          4.95696885e-03,   4.77526132e-05,   8.65745168e-04,   4.86454913e-04,   5.29006027e-04,   7.52931937e-04,
          6.28924249e-04,   1.13444548e-03,   5.91434309e-04,   1.04509283e-03,   6.79229881e-04,   6.56975882e-03,
          5.62161868e-03,   6.12289899e-04,   1.83546934e-03,   1.26296729e-03,   1.81250600e-03,   1.20100596e-03}},
       {{ 1.02666411e-03,   8.65386208e-04,   7.69264668e-04,   7.10482783e-04,   1.02522741e-03,   1.00319100e-03,
          9.83374068e-04,   7.27047589e-04,   9.46941342e-04,   8.03911987e-04,   7.46929642e-04,   5.23075621e-04,
          4.09870887e-04,   3.60371645e-04,   3.37117897e-04,   3.29870134e-04,   3.96330817e-04,   3.82776041e-04,
          3.70890639e-04,   2.91955691e-04,   3.91794608e-04,   3.74163485e-04,   3.46838576e-04,   3.02321340e-04,
          4.56147849e-04,   7.81061250e-04,  -4.24580452e-04,   1.83817346e-04,   1.07490238e-03,   4.95696885e-03,
          3.08127198e-02,  -3.20993863e-04,   9.73979683e-04,   4.45762161e-04,  -1.54865753e-04,   1.50384140e-03,
         -2.61163315e-04,   2.06161181e-04,  -1.71824128e-03,  -4.58779393e-04,  -3.79452796e-03,   7.03343102e-03,
          2.95689467e-02,  -1.23587563e-03,   2.58422584e-03,   2.77929929e-03,  -3.78601882e-04,   3.06037890e-03}},
       {{ 1.12775268e-03,   9.50594845e-04,   8.45008877e-04,   7.80439143e-04,   1.12617452e-03,   1.10196833e-03,
          1.08020016e-03,   7.98634972e-04,   1.04018015e-03,   8.83067680e-04,   8.20474675e-04,   5.74579285e-04,
          4.50228058e-04,   3.95854966e-04,   3.70311582e-04,   3.62350181e-04,   4.35354791e-04,   4.20465369e-04,
          4.07409694e-04,   3.20702563e-04,   4.30371933e-04,   4.11004794e-04,   3.80989388e-04,   3.32088846e-04,
          4.49669426e-04,   3.75531134e-04,   4.08990045e-04,   6.37348484e-04,   6.44888647e-04,   4.77526132e-05,
         -3.20993863e-04,   2.95120706e-03,   3.81616148e-04,   4.80852027e-04,   3.54022259e-04,   4.34691956e-04,
          1.24061065e-03,   1.13303966e-03,   1.13317572e-03,   2.66518032e-04,   2.47238976e-03,   8.96510580e-04,
          1.98897459e-04,   4.11588312e-03,   6.46648779e-04,   7.96914278e-04,   9.29845625e-04,   5.56817267e-04}},
       {{ 1.32074669e-03,   1.11327157e-03,   9.89616517e-04,   9.13996867e-04,   1.31889845e-03,   1.29054983e-03,
          1.26505644e-03,   9.35306575e-04,   1.21818775e-03,   1.03418838e-03,   9.60883739e-04,   6.72907902e-04,
          5.27276262e-04,   4.63598222e-04,   4.33683559e-04,   4.24359712e-04,   5.09857710e-04,   4.92420239e-04,
          4.77130327e-04,   3.75584875e-04,   5.04022129e-04,   4.81340662e-04,   4.46188675e-04,   3.88919709e-04,
          7.34567494e-04,   6.07606788e-04,   6.05495682e-04,   4.21077641e-04,   6.10649172e-04,   8.65745168e-04,
          9.73979683e-04,   3.81616148e-04,   1.85006292e-03,   1.12631785e-03,   9.72520588e-04,   1.14094680e-03,
          1.19282083e-03,   1.15436836e-03,   8.69858885e-04,   8.29136183e-04,   1.03012551e-03,   1.57146619e-03,
          1.14946730e-03,   7.13462790e-04,   2.12278595e-03,   1.06578187e-03,   1.00354814e-03,   1.13488677e-03}},
       {{ 1.18965143e-03,   1.00276997e-03,   8.91388721e-04,   8.23274959e-04,   1.18798665e-03,   1.16245186e-03,
          1.13948890e-03,   8.42469500e-04,   1.09727233e-03,   9.31536450e-04,   8.65507914e-04,   6.06116110e-04,
          4.74939641e-04,   4.17582184e-04,   3.90636804e-04,   3.82238427e-04,   4.59250027e-04,   4.43543372e-04,
          4.29771113e-04,   3.38304904e-04,   4.53993676e-04,   4.33563537e-04,   4.01900681e-04,   3.50316144e-04,
          6.12454427e-04,   6.24590634e-04,   5.88293491e-04,   4.72122327e-04,   5.77466040e-04,   4.86454913e-04,
          4.45762161e-04,   4.80852027e-04,   1.12631785e-03,   1.56085654e-03,   1.08972498e-03,   9.82908081e-04,
          1.19713615e-03,   1.19499304e-03,   9.53489597e-04,   7.28034578e-04,   1.09610665e-03,   8.00113601e-04,
          6.70931636e-04,   3.94186445e-04,   1.24725275e-03,   1.71791996e-03,   1.03758754e-03,   1.20993724e-03}},
       {{ 1.21299554e-03,   1.02244697e-03,   9.08880133e-04,   8.39429799e-04,   1.21129810e-03,   1.18526225e-03,
          1.16184870e-03,   8.59000987e-04,   1.11880372e-03,   9.49815667e-04,   8.82491476e-04,   6.18009717e-04,
          4.84259217e-04,   4.25776254e-04,   3.98302135e-04,   3.89738959e-04,   4.68261731e-04,   4.52246870e-04,
          4.38204363e-04,   3.44943344e-04,   4.62902237e-04,   4.42071204e-04,   4.09787039e-04,   3.57190276e-04,
          5.08546766e-04,   5.17497842e-04,   5.48859822e-04,   2.31408775e-04,   4.80122632e-04,   5.29006027e-04,
         -1.54865753e-04,   3.54022259e-04,   9.72520588e-04,   1.08972498e-03,   2.34400043e-03,   1.07495317e-03,
          1.17876923e-03,   1.19304944e-03,   8.85787609e-04,   3.62201277e-04,   1.02750321e-03,   1.43589426e-03,
          8.26537984e-04,   7.89878079e-04,   1.07513394e-03,   1.25030310e-03,   2.48522661e-03,   1.89995014e-03}},
       {{ 1.23785073e-03,   1.04339768e-03,   9.27503766e-04,   8.56630343e-04,   1.23611850e-03,   1.20954916e-03,
          1.18565585e-03,   8.76602559e-04,   1.14172884e-03,   9.69278099e-04,   9.00574386e-04,   6.30673198e-04,
          4.94182050e-04,   4.34500728e-04,   4.06463643e-04,   3.97725002e-04,   4.77856764e-04,   4.61513747e-04,
          4.47183497e-04,   3.52011491e-04,   4.72387449e-04,   4.51129573e-04,   4.18183882e-04,   3.64509372e-04,
          6.28601201e-04,   5.33033829e-04,   5.31118908e-04,   1.24400834e-03,   6.11664360e-04,   7.52931937e-04,
          1.50384140e-03,   4.34691956e-04,   1.14094680e-03,   9.82908081e-04,   1.07495317e-03,   2.40248634e-03,
          1.69973282e-03,   1.30324740e-03,   7.58106794e-04,   1.83244405e-03,   1.38243050e-03,   1.99675742e-03,
          1.79287163e-03,   4.14454010e-04,   1.42580254e-03,   1.15308962e-03,   1.40323524e-03,   4.02199278e-03}},
       {{ 4.95911584e-03,   4.18009202e-03,   3.71579425e-03,   3.43185895e-03,   4.95217613e-03,   4.84573323e-03,
          4.75001110e-03,   3.51187226e-03,   4.57402938e-03,   3.88315188e-03,   3.60790894e-03,   2.52662246e-03,
          1.97980740e-03,   1.74071024e-03,   1.62838721e-03,   1.59337819e-03,   1.91440454e-03,   1.84893063e-03,
          1.79152035e-03,   1.41023931e-03,   1.89249320e-03,   1.80732924e-03,   1.67534119e-03,   1.46030871e-03,
          1.56385469e-03,   1.35645413e-03,   1.15629002e-03,   1.42279050e-03,   1.66383290e-03,   6.28924249e-04,
         -2.61163315e-04,   1.24061065e-03,   1.19282083e-03,   1.19713615e-03,   1.17876923e-03,   1.69973282e-03,
          9.17190457e-03,   7.10852472e-03,   5.85389292e-03,   6.09036969e-03,   6.56707634e-03,   4.92078547e-03,
          3.58487836e-03,   3.70572531e-03,   7.50604852e-03,   6.19473144e-03,   6.18040130e-03,   5.75711022e-03}},
       {{ 4.34968347e-03,   3.66639492e-03,   3.25915533e-03,   3.01011322e-03,   4.34359658e-03,   4.25023459e-03,
          4.16627588e-03,   3.08029358e-03,   4.01192079e-03,   3.40594615e-03,   3.16452818e-03,   2.21612244e-03,
          1.73650622e-03,   1.52679203e-03,   1.42827253e-03,   1.39756582e-03,   1.67914080e-03,   1.62171308e-03,
          1.57135802e-03,   1.23693312e-03,   1.65992218e-03,   1.58522413e-03,   1.46945627e-03,   1.28084942e-03,
          1.31293893e-03,   1.66752429e-03,   1.20158071e-03,   1.43298233e-03,   1.36807255e-03,   1.13444548e-03,
          2.06161181e-04,   1.13303966e-03,   1.15436836e-03,   1.19499304e-03,   1.19304944e-03,   1.30324740e-03,
          7.10852472e-03,   7.71641829e-03,   5.46786661e-03,   5.31871486e-03,   5.35749409e-03,   4.93872224e-03,
          3.11714586e-03,   3.90659631e-03,   6.49534653e-03,   5.40851364e-03,   6.36079949e-03,   4.92747310e-03}},
       {{ 4.16084895e-03,   3.50722428e-03,   3.11766434e-03,   2.87943399e-03,   4.15502633e-03,   4.06571749e-03,
          3.98540371e-03,   2.94656759e-03,   3.83774971e-03,   3.25808248e-03,   3.02714528e-03,   2.11991305e-03,
          1.66111860e-03,   1.46050881e-03,   1.36626637e-03,   1.33689274e-03,   1.60624361e-03,   1.55130901e-03,
          1.50314004e-03,   1.18323366e-03,   1.58785933e-03,   1.51640418e-03,   1.40566219e-03,   1.22524340e-03,
          1.28654836e-03,   1.07724973e-03,   1.87883558e-03,   9.30704025e-04,   9.62439621e-04,   5.91434309e-04,
         -1.71824128e-03,   1.13317572e-03,   8.69858885e-04,   9.53489597e-04,   8.85787609e-04,   7.58106794e-04,
          5.85389292e-03,   5.46786661e-03,   1.25358200e-02,   4.08022244e-03,   3.94710135e-03,   4.12160187e-03,
          1.14233505e-03,   5.54794485e-03,   5.50647287e-03,   5.43680078e-03,   4.29615181e-03,   5.33418622e-03}},
       {{ 3.38302876e-03,   2.85159129e-03,   2.53485485e-03,   2.34115877e-03,   3.37829460e-03,   3.30568097e-03,
          3.24038088e-03,   2.39574256e-03,   3.12032900e-03,   2.64902352e-03,   2.46125723e-03,   1.72362105e-03,
          1.35059265e-03,   1.18748442e-03,   1.11085946e-03,   1.08697688e-03,   1.30597587e-03,   1.26131063e-03,
          1.22214626e-03,   9.62042492e-04,   1.29102831e-03,   1.23293083e-03,   1.14289071e-03,   9.96199026e-04,
          1.43809711e-03,   1.45066027e-03,   1.14047363e-03,   8.05921777e-03,   5.04849264e-04,   1.04509283e-03,
         -4.58779393e-04,   2.66518032e-04,   8.29136183e-04,   7.28034578e-04,   3.62201277e-04,   1.83244405e-03,
          6.09036969e-03,   5.31871486e-03,   4.08022244e-03,   1.38208790e-02,   4.20450879e-03,   4.95541051e-03,
          1.18456306e-03,   2.04428001e-03,   5.98432041e-03,   3.70430611e-03,   4.89779393e-03,   4.82114306e-03}},
       {{ 4.29350045e-03,   3.61903765e-03,   3.21705820e-03,   2.97123286e-03,   4.28749219e-03,   4.19533611e-03,
          4.11246186e-03,   3.04050674e-03,   3.96010051e-03,   3.36195299e-03,   3.12365331e-03,   2.18749773e-03,
          1.71407651e-03,   1.50707111e-03,   1.40982413e-03,   1.37951405e-03,   1.65745205e-03,   1.60076610e-03,
          1.55106145e-03,   1.22095618e-03,   1.63848167e-03,   1.56474846e-03,   1.45047593e-03,   1.26430523e-03,
          1.01443327e-03,   1.31063879e-03,   1.17349678e-03,   5.45398848e-04,   3.52613692e-03,   6.79229881e-04,
         -3.79452796e-03,   2.47238976e-03,   1.03012551e-03,   1.09610665e-03,   1.02750321e-03,   1.38243050e-03,
          6.56707634e-03,   5.35749409e-03,   3.94710135e-03,   4.20450879e-03,   2.70032750e-02,   5.58448977e-03,
         -1.02964508e-03,   4.42564003e-03,   6.93922375e-03,   7.01259759e-03,   7.37404471e-03,   3.21683952e-03}},
       {{ 4.62434137e-03,   3.89790701e-03,   3.46495255e-03,   3.20018484e-03,   4.61787013e-03,   4.51861286e-03,
          4.42935263e-03,   3.27479670e-03,   4.26525089e-03,   3.62101238e-03,   3.36435023e-03,   2.35605804e-03,
          1.84615677e-03,   1.62320031e-03,   1.51845986e-03,   1.48581419e-03,   1.78516905e-03,   1.72411508e-03,
          1.67058038e-03,   1.31503845e-03,   1.76473687e-03,   1.68532207e-03,   1.56224412e-03,   1.36172782e-03,
          1.27128326e-03,   1.34265288e-03,   1.00235517e-03,   2.57677852e-03,   1.00496214e-03,   6.56975882e-03,
          7.03343102e-03,   8.96510580e-04,   1.57146619e-03,   8.00113601e-04,   1.43589426e-03,   1.99675742e-03,
          4.92078547e-03,   4.93872224e-03,   4.12160187e-03,   4.95541051e-03,   5.58448977e-03,   1.71561029e-02,
          1.32668164e-02,   5.03812979e-03,   7.60764388e-03,   5.92736822e-03,   7.07916854e-03,   6.12496655e-03}},
       {{ 3.48742991e-03,   2.93959213e-03,   2.61308112e-03,   2.41340754e-03,   3.48254965e-03,   3.40769515e-03,
          3.34037987e-03,   2.46967580e-03,   3.21662316e-03,   2.73077307e-03,   2.53721226e-03,   1.77681244e-03,
          1.39227228e-03,   1.22413050e-03,   1.14514087e-03,   1.12052127e-03,   1.34627862e-03,   1.30023500e-03,
          1.25986201e-03,   9.91731373e-04,   1.33086977e-03,   1.27097939e-03,   1.17816062e-03,   1.02694199e-03,
          8.35335045e-04,   1.31170005e-03,   2.03437868e-04,  -1.06622020e-04,   1.18660022e-03,   5.62161868e-03,
          2.95689467e-02,   1.98897459e-04,   1.14946730e-03,   6.70931636e-04,   8.26537984e-04,   1.79287163e-03,
          3.58487836e-03,   3.11714586e-03,   1.14233505e-03,   1.18456306e-03,  -1.02964508e-03,   1.32668164e-02,
          4.19731825e-02,   1.86568323e-03,   6.89543333e-03,   5.23019625e-03,   3.70609172e-03,   7.84528369e-03}},
       {{ 2.97918002e-03,   2.51118284e-03,   2.23225678e-03,   2.06168316e-03,   2.97501100e-03,   2.91106561e-03,
          2.85356071e-03,   2.10975101e-03,   2.74784001e-03,   2.33279658e-03,   2.16744487e-03,   1.51786394e-03,
          1.18936577e-03,   1.04572858e-03,   9.78250717e-04,   9.57219114e-04,   1.15007512e-03,   1.11074179e-03,
          1.07625266e-03,   8.47198758e-04,   1.13691192e-03,   1.08574982e-03,   1.00645824e-03,   8.77277860e-04,
          9.15299594e-04,   7.47138541e-04,   1.26090814e-03,   7.16501311e-05,   7.42587266e-04,   6.12289899e-04,
         -1.23587563e-03,   4.11588312e-03,   7.13462790e-04,   3.94186445e-04,   7.89878079e-04,   4.14454010e-04,
          3.70572531e-03,   3.90659631e-03,   5.54794485e-03,   2.04428001e-03,   4.42564003e-03,   5.03812979e-03,
          1.86568323e-03,   2.86333952e-02,   3.93070528e-03,   3.53116766e-03,   6.11464535e-03,   3.22423463e-03}},
       {{ 4.81571629e-03,   4.05921900e-03,   3.60834703e-03,   3.33262210e-03,   4.80897725e-03,   4.70561229e-03,
          4.61265809e-03,   3.41032172e-03,   4.44176513e-03,   3.77086528e-03,   3.50358136e-03,   2.45356175e-03,
          1.92255859e-03,   1.69037525e-03,   1.58130019e-03,   1.54730351e-03,   1.85904694e-03,   1.79546630e-03,
          1.73971611e-03,   1.36946033e-03,   1.83776920e-03,   1.75506787e-03,   1.62689643e-03,   1.41808191e-03,
          1.22312381e-03,   1.01714463e-03,   1.07363243e-03,   1.49589387e-03,   9.40734522e-04,   1.83546934e-03,
          2.58422584e-03,   6.46648779e-04,   2.12278595e-03,   1.24725275e-03,   1.07513394e-03,   1.42580254e-03,
          7.50604852e-03,   6.49534653e-03,   5.50647287e-03,   5.98432041e-03,   6.93922375e-03,   7.60764388e-03,
          6.89543333e-03,   3.93070528e-03,   1.58276633e-02,   1.06037071e-02,   1.01667597e-02,   9.57967990e-03}},
       {{ 3.64539565e-03,   3.07274315e-03,   2.73144259e-03,   2.52272463e-03,   3.64029434e-03,   3.56204924e-03,
          3.49168488e-03,   2.58154161e-03,   3.36232251e-03,   2.85446548e-03,   2.65213719e-03,   1.85729449e-03,
          1.45533630e-03,   1.27957841e-03,   1.19701089e-03,   1.17127611e-03,   1.40725932e-03,   1.35913012e-03,
          1.31692840e-03,   1.03665259e-03,   1.39115252e-03,   1.32854935e-03,   1.23152628e-03,   1.07345809e-03,
          8.60802781e-04,   8.94399194e-04,   9.16279529e-04,   4.98951942e-04,   9.26615671e-04,   1.26296729e-03,
          2.77929929e-03,   7.96914278e-04,   1.06578187e-03,   1.71791996e-03,   1.25030310e-03,   1.15308962e-03,
          6.19473144e-03,   5.40851364e-03,   5.43680078e-03,   3.70430611e-03,   7.01259759e-03,   5.92736822e-03,
          5.23019625e-03,   3.53116766e-03,   1.06037071e-02,   1.37781543e-02,   9.61972241e-03,   8.49261230e-03}},
       {{ 3.70515119e-03,   3.12311174e-03,   2.77621655e-03,   2.56407728e-03,   3.69996625e-03,   3.62043856e-03,
          3.54892077e-03,   2.62385838e-03,   3.41743789e-03,   2.90125604e-03,   2.69561117e-03,   1.88773936e-03,
          1.47919225e-03,   1.30055333e-03,   1.21663236e-03,   1.19047574e-03,   1.43032719e-03,   1.38140906e-03,
          1.33851556e-03,   1.05364545e-03,   1.41395637e-03,   1.35032701e-03,   1.25171353e-03,   1.09105428e-03,
          7.63265976e-04,   1.14275093e-03,   6.88048208e-04,   1.65300146e-03,   1.36626987e-03,   1.81250600e-03,
         -3.78601882e-04,   9.29845625e-04,   1.00354814e-03,   1.03758754e-03,   2.48522661e-03,   1.40323524e-03,
          6.18040130e-03,   6.36079949e-03,   4.29615181e-03,   4.89779393e-03,   7.37404471e-03,   7.07916854e-03,
          3.70609172e-03,   6.11464535e-03,   1.01667597e-02,   9.61972241e-03,   1.45128431e-02,   8.92758070e-03}},
       {{ 3.18713761e-03,   2.68647253e-03,   2.38807642e-03,   2.20559614e-03,   3.18267758e-03,   3.11426858e-03,
          3.05274962e-03,   2.25701930e-03,   2.93964923e-03,   2.49563426e-03,   2.31874040e-03,   1.62381637e-03,
          1.27238783e-03,   1.11872424e-03,   1.04653617e-03,   1.02403649e-03,   1.23035454e-03,   1.18827560e-03,
          1.15137901e-03,   9.06336312e-04,   1.21627251e-03,   1.16153911e-03,   1.07671268e-03,   9.38515043e-04,
          7.94814777e-04,   8.08201181e-04,   1.31847608e-03,   2.28772211e-03,   7.18390101e-04,   1.20100596e-03,
          3.06037890e-03,   5.56817267e-04,   1.13488677e-03,   1.20993724e-03,   1.89995014e-03,   4.02199278e-03,
          5.75711022e-03,   4.92747310e-03,   5.33418622e-03,   4.82114306e-03,   3.21683952e-03,   6.12496655e-03,
          7.84528369e-03,   3.22423463e-03,   9.57967990e-03,   8.49261230e-03,   8.92758070e-03,   2.18743795e-02}}}},
            0u
        };

        /*
         * B -> D lattice form factors from [HPQCD2015A]
         */
        static const MultivariateGaussianCovarianceConstraintEntry<7> B_to_D_f_plus_f_zero_0_to_11dot62_HPQCD_2015A
        {
            "B->D::f_++f_0@HPQCD2015A",
            {{ "B->D::f_+(s)", "B->D::f_+(s)", "B->D::f_+(s)", "B->D::f_+(s)", "B->D::f_0(s)", "B->D::f_0(s)", "B->D::f_0(s)" }},
            {{
                 Kinematics{ { "s",  0.00 } },
                 Kinematics{ { "s",  4.00 } },
                 Kinematics{ { "s",  8.00 } },
                 Kinematics{ { "s", 11.62 } },
                 Kinematics{ { "s",  4.00 } },
                 Kinematics{ { "s",  8.00 } },
                 Kinematics{ { "s", 11.62 } }
            }},
            {{
                 Options{ { "form-factors", "BCL2008" } },
                 Options{ { "form-factors", "BCL2008" } },
                 Options{ { "form-factors", "BCL2008" } },
                 Options{ { "form-factors", "BCL2008" } },
                 Options{ { "form-factors", "BCL2008" } },
                 Options{ { "form-factors", "BCL2008" } },
                 Options{ { "form-factors", "BCL2008" } }
            }},
            {{ 0.665,          0.798,          0.972,         1.177,           0.729,          0.810,          0.901          }},
            {{
                 {{ 1.128e-3, 1.042e-3, 9.230e-4, 7.727e-4, 1.093e-3, 1.063e-3, 1.045e-3 }},
                 {{ 1.042e-3, 1.079e-3, 1.108e-3, 1.123e-3, 1.026e-3, 1.017e-3, 1.021e-3 }},
                 {{ 9.230e-4, 1.108e-3, 1.331e-3, 1.576e-3, 9.307e-4, 9.511e-4, 9.865e-4 }},
                 {{ 7.727e-4, 1.123e-3, 1.576e-3, 2.112e-3, 8.108e-4, 8.681e-4, 9.425e-4 }},
                 {{ 1.093e-3, 1.026e-3, 9.307e-4, 8.108e-4, 1.126e-3, 1.165e-3, 1.210e-3 }},
                 {{ 1.063e-3, 1.017e-3, 9.511e-4, 8.681e-4, 1.165e-3, 1.283e-3, 1.410e-3 }},
                 {{ 1.045e-3, 1.021e-3, 9.865e-4, 9.425e-4, 1.210e-3, 1.410e-3, 1.635e-3 }}
            }},
            0u
        };

        /*
         * B -> pi lattice form factors from [FNALMILC2015A]
         */
        static const MultivariateGaussianCovarianceConstraintEntry<6> B_to_pi_f_plus_f_zero_18_to_26_FNALMILC_2015A
        {
            "B->pi::f_++f_0@FNALMILC2015A",
            {{ "B->pi::f_+(s)", "B->pi::f_+(s)", "B->pi::f_+(s)", "B->pi::f_0(s)", "B->pi::f_0(s)", "B->pi::f_0(s)" }},
            {{
                 Kinematics{ { "s", 18.00 } },
                 Kinematics{ { "s", 22.00 } },
                 Kinematics{ { "s", 26.00 } },
                 Kinematics{ { "s", 18.00 } },
                 Kinematics{ { "s", 22.00 } },
                 Kinematics{ { "s", 26.00 } }
            }},
            {{
                 Options{ { "form-factors", "BCL2008" } },
                 Options{ { "form-factors", "BCL2008" } },
                 Options{ { "form-factors", "BCL2008" } },
                 Options{ { "form-factors", "BCL2008" } },
                 Options{ { "form-factors", "BCL2008" } },
                 Options{ { "form-factors", "BCL2008" } }
            }},
            {{ 1.016,           1.971,           6.443,           0.417,           0.609,           0.961           }},
            {{
                 {{ 3.492e-3, 1.642e-3, 1.648e-3, 1.067e-3, 2.904e-4, 1.096e-4 }},
                 {{ 1.997e-3, 3.371e-3, 6.193e-3, 2.123e-4, 2.167e-4, 1.294e-4 }},
                 {{ 1.648e-3, 6.193e-3, 7.419e-2, 2.064e-3, 1.139e-3, 1.346e-3 }},
                 {{ 1.067e-3, 2.123e-4, 2.064e-3, 8.478e-4, 4.266e-4, 3.150e-4 }},
                 {{ 2.904e-4, 2.167e-4, 1.139e-3, 4.266e-4, 3.923e-4, 4.009e-4 }},
                 {{ 1.096e-4, 1.294e-4, 1.346e-3, 3.150e-4, 4.009e-4, 6.467e-4 }},
            }},
            0u
        };

        /*
         * Lambda_b -> Lambda lattice form factors from [DM2016]
         *
         * We got input files from Stefan Meinel on the total covariance for the
         * form factors including systematic effects from a higher-order fit.
         * This works well, as long as they are used only at large q^2 values.
         *
         * The covariance matrix covers both the statistical and systematic errors.
         */
        static const MultivariateGaussianCovarianceConstraintEntry<22> Lambda_b_to_Lambda_v_a_t_t5_parameters_DM_2016
        {
            "Lambda_b->Lambda::f_perp+long^V+A+T+T5@DM2016",
            {{
                "Lambda_b->Lambda::a_0_long^V@DM2016",  "Lambda_b->Lambda::a_1_long^V@DM2016",  "Lambda_b->Lambda::a_2_long^V@DM2016",
                "Lambda_b->Lambda::a_0_perp^V@DM2016",  "Lambda_b->Lambda::a_1_perp^V@DM2016",  "Lambda_b->Lambda::a_2_perp^V@DM2016",
                "Lambda_b->Lambda::a_0_long^A@DM2016",  "Lambda_b->Lambda::a_1_long^A@DM2016",  "Lambda_b->Lambda::a_2_long^A@DM2016",
                /* perp^A replaced through FF rel. */   "Lambda_b->Lambda::a_1_perp^A@DM2016",  "Lambda_b->Lambda::a_2_perp^A@DM2016",
                "Lambda_b->Lambda::a_0_long^T@DM2016",  "Lambda_b->Lambda::a_1_long^T@DM2016",  "Lambda_b->Lambda::a_2_long^T@DM2016",
                "Lambda_b->Lambda::a_0_perp^T@DM2016",  "Lambda_b->Lambda::a_1_perp^T@DM2016",  "Lambda_b->Lambda::a_2_perp^T@DM2016",
                "Lambda_b->Lambda::a_0_long^T5@DM2016", "Lambda_b->Lambda::a_1_long^T5@DM2016", "Lambda_b->Lambda::a_2_long^T5@DM2016",
                /* perp^T5 replaced through FF rel. */  "Lambda_b->Lambda::a_1_perp^T5@DM2016", "Lambda_b->Lambda::a_2_perp^T5@DM2016",
            }},
            {{
                Kinematics{},                           Kinematics{},                           Kinematics{},
                Kinematics{},                           Kinematics{},                           Kinematics{},
                Kinematics{},                           Kinematics{},                           Kinematics{},
                /* replaced */                          Kinematics{},                           Kinematics{},
                Kinematics{},                           Kinematics{},                           Kinematics{},
                Kinematics{},                           Kinematics{},                           Kinematics{},
                Kinematics{},                           Kinematics{},                           Kinematics{},
                /* replaced */                          Kinematics{},                           Kinematics{},
            }},
            // no options; default ctor just fine
            {{
                Options{},                              Options{},                              Options{},
                Options{},                              Options{},                              Options{},
                Options{},                              Options{},                              Options{},
                /* replaced */                          Options{},                              Options{},
                Options{},                              Options{},                              Options{},
                Options{},                              Options{},                              Options{},
                Options{},                              Options{},                              Options{},
                /* replaced */                          Options{},                              Options{},
            }},
            // mean value
            {{
                0.42208407714270,                      -1.13856687165290,                       0.0,
                0.51820088929706,                      -1.34951141147510,                       0.0,
                0.35627001295791,                      -1.06124674284310,                       0.0,
                /* replaced */                         -1.13570008046520,                       0.0,
                0.49599719500981,                      -1.12751367398150,                       0.0,
                0.38763266069053,                      -0.96229696628603,                       0.0,
                0.34026877180649,                      -0.76969429814750,                       0.0,
                /* replaced */                         -0.80082045409023,                       0.0
            }},
            // covariance matrix
            {{
                {{ 0.000749594187095, -0.006164273736522, 0.011372890644823, 0.000670223129197, -0.004243056424882, 0.001311870031246, 0.000264886750161, -0.002053357392716, -0.000398016204121, -0.002115221757594, 0.001820772880284, 0.000499996242162, -0.003191406233943, -0.001640271158008, 0.000383714387856, -0.002901531299397, 0.000472774777764, 0.000241060019024, -0.001700916517912, -0.000697409076608, -0.001581214048131, -0.000775332925801 }},
                {{ -0.006164273736522, 0.094116175041483, -0.371146407307241, -0.004710830597348, 0.042219205628345, -0.019976547027344, -0.0019281304725, 0.018735543980848, 0.014036368652166, 0.019201191872417, -0.020106226373253, -0.003835093433078, 0.031430053093567, 0.007276614667715, -0.002993782257011, 0.027808187373224, -0.015267877320676, -0.001716349735084, 0.014137211484609, 0.017181121597857, 0.013650147261026, 0.01326332411015 }},
                {{ 0.011372890644823, -0.371146407307241, 3.22993771179459, 0.004972780726833, -0.072506782444349, 0.261221262531002, -0.000215467856408, 0.00838207516522, 0.1945442807215, 0.009519525198421, 0.181693364650388, 0.003361077035058, -0.049693892474132, 0.105918716815944, 0.002547113451561, -0.042437381280063, 0.177689190424785, -0.000254410809612, 0.00281801275487, 0.197330994815814, 0.001615329143183, 0.159483826612556 }},
                {{ 0.000670223129197, -0.004710830597348, 0.004972780726833, 0.001244975293461, -0.010103508140104, 0.009465313700532, 0.00031710571775, -0.002370606620954, -0.002073093853597, -0.002455839770168, 0.002035470147828, 0.000674156597717, -0.004529175663623, -0.002211827382624, 0.000502830828959, -0.00395072788092, 0.000807832408955, 0.000286488048264, -0.001981930307026, -0.001683025856184, -0.001827629003901, -0.001755499943138 }},
                {{ -0.004243056424882, 0.042219205628345, -0.072506782444349, -0.010103508140104, 0.16237356863464, -0.393141309425425, -0.002102727635876, 0.020507351111945, 0.034673239099397, 0.021108524650448, -0.020472349960809, -0.004840912074375, 0.044918377921202, 0.00493978888909, -0.003614639566594, 0.036227885425357, -0.022107436919726, -0.001846776263455, 0.015230274625778, 0.033150706453419, 0.014969174695481, 0.024981996683362 }},
                {{ 0.001311870031246, -0.019976547027344, 0.261221262531002, 0.009465313700532, -0.393141309425425, 2.6882275522931, -0.000209850641227, 0.000573909527041, 0.314953624848913, 0.003987601147349, 0.095782792052462, 0.000919027439209, -0.017126393605156, 0.193442612467125, 0.000759784430623, -0.018270569805431, 0.170856919496887, -0.000120107739701, 0.000253200262659, 0.07238090091922, 0.000126527927104, 0.086301413616358 }},
                {{ 0.000264886750161, -0.0019281304725, -0.000215467856408, 0.00031710571775, -0.002102727635876, -0.000209850641227, 0.000421991521668, -0.003979786111674, 0.014411659822313, -0.003877135897698, 0.008461260962696, 0.000291906584519, -0.001881247064266, -0.0014014576582, 0.000227998420793, -0.001642649592213, -0.000390194788151, 0.000206965783915, -0.001540812967208, -0.001282055263547, -0.001478113690651, -0.00100909599116 }},
                {{ -0.002053357392716, 0.018735543980848, 0.00838207516522, -0.002370606620954, 0.020507351111945, 0.000573909527041, -0.003979786111674, 0.08362176103837, -0.56600519794316, 0.055036258016742, -0.138287010845185, -0.002143620553598, 0.017487112571067, 0.008695701716403, -0.001736129561578, 0.015050849220897, 0.001930445648249, -0.001746797606144, 0.020411788809007, 0.006368955993306, 0.019360782343715, -0.001941016331335 }},
                {{ -0.000398016204121, 0.014036368652166, 0.1945442807215, -0.002073093853597, 0.034673239099397, 0.314953624848913, 0.014411659822313, -0.56600519794316, 7.34730246672972, -0.240144672696889, 1.34324198503075, -0.002103935457479, 0.030878873940266, 0.204760275234375, -0.000955655887978, 0.009344977098389, 0.171318684690066, 0.001800545651617, -0.062693187561659, 0.565109329672463, -0.054982167413635, 0.573291529393943 }},
                {{ -0.002115221757594, 0.019201191872417, 0.009519525198421, -0.002455839770168, 0.021108524650448, 0.003987601147349, -0.003877135897698, 0.055036258016742, -0.240144672696889, 0.087151626684731, -0.43623988664169, -0.002237558080252, 0.018243425717256, 0.011287492152364, -0.00179377376775, 0.015468573621476, 0.003619796423987, -0.001712096300904, 0.020580755619435, 0.00243286851162, 0.019511123699285, -0.003792141576066 }},
                {{ 0.001820772880284, -0.020106226373253, 0.181693364650388, 0.002035470147828, -0.020472349960809, 0.095782792052462, 0.008461260962696, -0.138287010845185, 1.34324198503075, -0.43623988664169, 6.06210830469137, 0.001354763557436, -0.01667194422627, 0.079580777434124, 0.001213509391727, -0.016206827037401, 0.077336984425539, 0.000934270925133, -0.036656076080858, 0.554264828901894, -0.032688115526928, 0.489039800962392 }},
                {{ 0.000499996242162, -0.003835093433078, 0.003361077035058, 0.000674156597717, -0.004840912074375, 0.000919027439209, 0.000291906584519, -0.002143620553598, -0.002103935457479, -0.002237558080252, 0.001354763557436, 0.001787395798472, -0.010435717161899, 0.003619965157987, 0.001054263300679, -0.005109429595942, 0.002104183253969, 0.000674377659471, -0.002897656132027, 0.003227125938417, -0.002707762322253, 0.001867376623833 }},
                {{ -0.003191406233943, 0.031430053093567, -0.049693892474132, -0.004529175663623, 0.044918377921202, -0.017126393605156, -0.001881247064266, 0.017487112571067, 0.030878873940266, 0.018243425717256, -0.01667194422627, -0.010435717161899, 0.15975337705658, -0.3512302745458, -0.004258118714408, 0.035002257332446, -0.021197955591748, -0.002478298810949, 0.01546678047878, 0.024862517939355, 0.015504170032356, 0.01906297476362 }},
                {{ -0.001640271158008, 0.007276614667715, 0.105918716815944, -0.002211827382624, 0.00493978888909, 0.193442612467125, -0.0014014576582, 0.008695701716403, 0.204760275234375, 0.011287492152364, 0.079580777434124,
0.003619965157987, -0.3512302745458, 2.8088303844212, -0.002137591323261, 0.002336940401947, 0.170948236544902, -0.001684398181157, 0.010911341318195, 0.063831820629273, 0.009468286998763, 0.07495950955294 }},
                {{ 0.000383714387856, -0.002993782257011, 0.002547113451561, 0.000502830828959, -0.003614639566594, 0.000759784430623, 0.000227998420793, -0.001736129561578, -0.000955655887978, -0.00179377376775, 0.001213509391727, 0.001054263300679, -0.004258118714408, -0.002137591323261, 0.000976602997817, -0.005798480068419, 0.005542408630269, 0.000531506161576, -0.002337876066771, 0.002913031047101, -0.002181920159564, 0.001975057127825 }},
                {{ -0.002901531299397, 0.027808187373224, -0.042437381280063, -0.00395072788092, 0.036227885425357, -0.018270569805431, -0.001642649592213, 0.015050849220897, 0.009344977098389, 0.015468573621476, -0.016206827037401, -0.005109429595942, 0.035002257332446, 0.002336940401947, -0.005798480068419, 0.076534052830695, -0.189181102504368, -0.002250134070338, 0.014176317637204, 0.011809686997914, 0.013599017609615, 0.006770890198818 }},
                {{ 0.000472774777764, -0.015267877320676, 0.177689190424785, 0.000807832408955, -0.022107436919726, 0.170856919496887, -0.000390194788151, 0.001930445648249, 0.171318684690066, 0.003619796423987, 0.077336984425539, 0.002104183253969, -0.021197955591748, 0.170948236544902, 0.005542408630269, -0.189181102504368, 1.34695571318776, 0.000435170076675, 0.00050048269771, 0.068672132008896, 0.000015364979713, 0.070620062283742 }},
                {{ 0.000241060019024, -0.001716349735084, -0.000254410809612, 0.000286488048264, -0.001846776263455, -0.000120107739701, 0.000206965783915, -0.001746797606144, 0.001800545651617, -0.001712096300904, 0.000934270925133, 0.000674377659471, -0.002478298810949, -0.001684398181157, 0.000531506161576, -0.002250134070338, 0.000435170076675, 0.000613916259219, -0.003480209951551, 0.008227405441161, -0.003287132276617, 0.007217165089323 }},
                {{ -0.001700916517912, 0.014137211484609, 0.00281801275487, -0.001981930307026, 0.015230274625778, 0.000253200262659, -0.001540812967208, 0.020411788809007, -0.062693187561659, 0.020580755619435, -0.036656076080858, -0.002897656132027, 0.01546678047878, 0.010911341318195, -0.002337876066771, 0.014176317637204, 0.00050048269771, -0.003480209951551, 0.062308556818102, -0.343929817844085, 0.04036330201319, -0.125254141399107 }},
                {{ -0.000697409076608, 0.017181121597857, 0.197330994815814, -0.001683025856184, 0.033150706453419, 0.07238090091922, -0.001282055263547, 0.006368955993306, 0.565109329672463, 0.00243286851162, 0.554264828901894, 0.003227125938417, 0.024862517939355, 0.063831820629273, 0.002913031047101, 0.011809686997914, 0.068672132008896, 0.008227405441161, -0.343929817844085, 8.86837050853116, -0.110102609787211, 3.19403648331734 }},
                {{ -0.001581214048131, 0.013650147261026, 0.001615329143183, -0.001827629003901, 0.014969174695481, 0.000126527927104, -0.001478113690651, 0.019360782343715, -0.054982167413635, 0.019511123699285, -0.032688115526928, -0.002707762322253, 0.015504170032356, 0.009468286998763, -0.002181920159564, 0.013599017609615, 0.000015364979713, -0.003287132276617, 0.04036330201319, -0.110102609787211, 0.051437631908735, -0.291661725840191 }},
                {{ -0.000775332925801, 0.01326332411015, 0.159483826612556, -0.001755499943138, 0.024981996683362, 0.086301413616358, -0.00100909599116, -0.001941016331335, 0.573291529393943, -0.003792141576066, 0.489039800962392, 0.001867376623833, 0.01906297476362, 0.07495950955294, 0.001975057127825, 0.006770890198818, 0.070620062283742, 0.007217165089323, -0.125254141399107, 3.19403648331734, -0.291661725840191, 6.14163645986141 }},
            }},
            0u
        };
        ///@}
    }

    /* Constraint */
    template <>
    struct WrappedForwardIteratorTraits<Constraint::BlockIteratorTag>
    {
        typedef std::vector<LogLikelihoodBlockPtr>::iterator UnderlyingIterator;
    };
    template class WrappedForwardIterator<Constraint::BlockIteratorTag, LogLikelihoodBlockPtr>;

    template <>
    struct WrappedForwardIteratorTraits<Constraint::ObservableIteratorTag>
    {
        typedef ObservableSet::Iterator UnderlyingIterator;
    };
    template class WrappedForwardIterator<Constraint::ObservableIteratorTag, ObservablePtr>;

    template <>
    struct Implementation<Constraint>
    {
        QualifiedName name;

        ObservableSet observables;

        std::vector<LogLikelihoodBlockPtr> blocks;

        Implementation(const QualifiedName & name,
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

    Constraint::Constraint(const QualifiedName & name,
                const std::vector<ObservablePtr> & observables,
                const std::vector<LogLikelihoodBlockPtr> & blocks) :
        PrivateImplementationPattern<Constraint>(new Implementation<Constraint>(name, observables, blocks))
    {
    }

    Constraint::~Constraint()
    {
    }

    const QualifiedName &
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

    typedef std::function<Constraint (const QualifiedName &, const Options & options)> ConstraintFactory;

    template <typename Factory_>
    ConstraintFactory make_factory(const Factory_ & f)
    {
        return std::bind(&Factory_::make, f, std::placeholders::_1, std::placeholders::_2);
    }

    const std::map<QualifiedName, const ConstraintEntry *> &
    make_constraint_entries()
    {
        using ValueType = std::map<QualifiedName, const ConstraintEntry *>::value_type;

        static const std::map<QualifiedName, const ConstraintEntry *> constraint_entries =
        {
            /* 2000 */
            // CLEO
            ValueType{ "B^0->K^*0gamma::BR@CLEO-2000", &entries::Bzero_to_Kstarzero_gamma_BR_CLEO_2000 },
            ValueType{ "B^+->K^*+gamma::BR@CLEO-2000", &entries::Bplus_to_Kstarplus_gamma_BR_CLEO_2000 },

            /* 2004 */
            // BaBar
            ValueType{ "B->X_sll::BR[1.0,6.0]@BaBar-2004A", &entries::Bmix_to_Xs_dilepton_BR_BaBar_2004A },
            // Belle
            ValueType{ "B^0->K^*0gamma::BR@Belle-2004", &entries::Bzero_to_Kstarzero_gamma_BR_Belle_2004 },
            ValueType{ "B^+->K^*+gamma::BR@Belle-2004", &entries::Bplus_to_Kstarplus_gamma_BR_Belle_2004 },

            /* 2005 */
            // Belle
            ValueType{ "B->X_sll::BR[1.0,6.0]@Belle-2005A", &entries::Bmix_to_Xs_dilepton_BR_Belle_2005A },

            /* 2006 */
            // Belle
            ValueType{ "B^0->K^*0gamma::S_K@Belle-2006", &entries::Bzero_to_Kstarzero_gamma_SKstargamma_Belle_2006 },
            ValueType{ "B^0->K^*0gamma::C_K@Belle-2006", &entries::Bzero_to_Kstarzero_gamma_CKstargamma_Belle_2006 },
            ValueType{ "B^0->K^*0gamma::S_K+C_K@Belle-2006", &entries::Bzero_to_Kstarzero_gamma_time_dependent_cp_asymmetries_Belle_2006 },

            /* 2008 */
            // BaBar
            ValueType{ "B^0->K^*0gamma::S_K@BaBar-2008", &entries::Bzero_to_Kstarzero_gamma_SKstargamma_BaBar_2008 },
            ValueType{ "B^0->K^*0gamma::C_K@BaBar-2008", &entries::Bzero_to_Kstarzero_gamma_CKstargamma_BaBar_2008 },
            ValueType{ "B^0->K^*0gamma::S_K+C_K@BaBar-2008", &entries::Bzero_to_Kstarzero_gamma_time_dependent_cp_asymmetries_BaBar_2008 },
            // Belle
            ValueType{ "B->X_sgamma::E_1[1.8]@Belle-2008", &entries::B_to_Xs_gamma_E_1_1dot8_Belle_2008A },
            ValueType{ "B->X_sgamma::E_2[1.8]@Belle-2008", &entries::B_to_Xs_gamma_E_2_1dot8_Belle_2008A },
            ValueType{ "B->X_sgamma::E_1[1.8]+E_2[1.8]@Belle-2008", &entries::B_to_Xs_gamma_E_1_and_E_2_1dot8_Belle_2008A },

            /* 2009 */
            // BaBar
            ValueType{ "B^0->K^*0gamma::BR@BaBar-2009", &entries::Bzero_to_Kstarzero_gamma_BR_BaBar_2009 },
            ValueType{ "B^+->K^*+gamma::BR@BaBar-2009", &entries::Bplus_to_Kstarplus_gamma_BR_BaBar_2009 },
            // Belle
            // B^+ -> K^+ mu^+ mu^-
            ValueType{ "B^+->K^+mu^+mu^-::BR[1.00,6.00]@Belle-2009", &entries::Bplus_to_Kplus_dimuon_BR_1_to_6_Belle_2009 },
            ValueType{ "B^+->K^+mu^+mu^-::BR[14.18,16.00]@Belle-2009", &entries::Bplus_to_Kplus_dimuon_BR_14dot18_to_16_Belle_2009 },
            ValueType{ "B^+->K^+mu^+mu^-::BR[16.00,22.86]@Belle-2009", &entries::Bplus_to_Kplus_dimuon_BR_16_to_22dot86_Belle_2009 },

            ValueType{ "B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@Belle-2009", &entries::Bplus_to_Kplus_dimuon_A_FB_1_to_6_Belle_2009 },
            ValueType{ "B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@Belle-2009", &entries::Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_Belle_2009 },
            ValueType{ "B^+->K^+mu^+mu^-::A_FB[16.00,22.86]@Belle-2009", &entries::Bplus_to_Kplus_dimuon_A_FB_16_to_22dot86_Belle_2009 },

            // B^0 -> K^*0 mu^+ mu^-
            ValueType{ "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@Belle-2009", &entries::Bzero_to_Kstarzero_dimuon_BR_1_to_6_Belle_2009 },
            ValueType{ "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@Belle-2009", &entries::Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_Belle_2009 },
            ValueType{ "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@Belle-2009", &entries::Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_Belle_2009 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@Belle-2009", &entries::Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_Belle_2009 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@Belle-2009", &entries::Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_Belle_2009 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@Belle-2009", &entries::Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_Belle_2009 },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@Belle-2009", &entries::Bzero_to_Kstarzero_dimuon_F_L_1_to_6_Belle_2009 },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@Belle-2009", &entries::Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_Belle_2009 },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@Belle-2009", &entries::Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_Belle_2009 },
            // B->X_sgamma
            ValueType{ "B->X_sgamma::BR[1.8]@Belle-2009B", &entries::B_to_Xs_gamma_BR_1dot8_Belle_2009B },
            ValueType{ "B->X_sgamma::BR[1.8]+E_1[1.8]+E_2[1.8]@Belle-2009B", &entries::B_to_Xs_gamma_1dot8_Belle_2009B },

            /* 2010 */
            // BaBar
            ValueType{ "B^0->pi^+lnu::BR[0.0,4.0]@BaBar-2010A", &entries::Bzero_to_pi_l_nu_BR_0_to_4_BaBar_2010A },
            ValueType{ "B^0->pi^+lnu::BR[4.0,8.0]@BaBar-2010A", &entries::Bzero_to_pi_l_nu_BR_4_to_8_BaBar_2010A },
            ValueType{ "B^0->pi^+lnu::BR[8.0,12.0]@BaBar-2010A", &entries::Bzero_to_pi_l_nu_BR_8_to_12_BaBar_2010A },
            ValueType{ "B^0->pi^+lnu::BR@BaBar-2010B", &entries::Bzero_to_pi_l_nu_BR_BaBar_2010B },
            //Belle
            ValueType{ "B^0->pi^+lnu::BR@Belle-2010A", &entries::Bzero_to_pi_l_nu_BR_Belle_2010A },

            /* 2011 */
            // HFAG
            ValueType{ "B^0->K^*0gamma::S_K+C_K@HFAG-2011", &entries::Bzero_to_Kstarzero_gamma_time_dependent_cp_asymmetries_HFAG_2011 },
            // CDF
            // B^0 -> K^0 mu^+ mu^-
            ValueType{ "B^0->K^0mu^+mu^-::BR[1.00,6.00]@CDF-2011", &entries::Bzero_to_Kzero_dimuon_BR_1_to_6_CDF_2011 },
            ValueType{ "B^0->K^0mu^+mu^-::BR[14.18,16.00]@CDF-2011", &entries::Bzero_to_Kzero_dimuon_BR_14dot18_to_16_CDF_2011 },
            ValueType{ "B^0->K^0mu^+mu^-::BR[16.00,22.86]@CDF-2011", &entries::Bzero_to_Kzero_dimuon_BR_16_to_22dot86_CDF_2011 },
            // B^+ -> K^+ mu^+ mu^-
            ValueType{ "B^+->K^+mu^+mu^-::BR[1.00,6.00]@CDF-2011", &entries::Bplus_to_Kplus_dimuon_BR_1_to_6_CDF_2011 },
            ValueType{ "B^+->K^+mu^+mu^-::BR[14.18,16.00]@CDF-2011", &entries::Bplus_to_Kplus_dimuon_BR_14dot18_to_16_CDF_2011 },
            ValueType{ "B^+->K^+mu^+mu^-::BR[16.00,22.86]@CDF-2011", &entries::Bplus_to_Kplus_dimuon_BR_16_to_22dot86_CDF_2011 },
            // B^0 -> K^*0 mu^+ mu^-
            ValueType{ "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@CDF-2011", &entries::Bzero_to_Kstarzero_dimuon_BR_1_to_6_CDF_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@CDF-2011", &entries::Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_CDF_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@CDF-2011", &entries::Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_CDF_2011 },
            // B^+ -> K^*+ mu^+ mu^-
            ValueType{ "B^+->K^*+mu^+mu^-::BR[1.00,6.00]@CDF-2011", &entries::Bplus_to_Kstarplus_dimuon_BR_1_to_6_CDF_2011 },
            ValueType{ "B^+->K^*+mu^+mu^-::BR[14.18,16.00]@CDF-2011", &entries::Bplus_to_Kstarplus_dimuon_BR_14dot18_to_16_CDF_2011 },
            ValueType{ "B^+->K^*+mu^+mu^-::BR[16.00,19.21]@CDF-2011", &entries::Bplus_to_Kstarplus_dimuon_BR_16_to_19dot21_CDF_2011 },
            // B^0 -> K^*0 mu^+ mu^-
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@CDF-2011", &entries::Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_CDF_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@CDF-2011", &entries::Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_CDF_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@CDF-2011", &entries::Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_CDF_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@CDF-2011", &entries::Bzero_to_Kstarzero_dimuon_F_L_1_to_6_CDF_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@CDF-2011", &entries::Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_CDF_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@CDF-2011", &entries::Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_CDF_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_T^2[1.00,6.00]@CDF-2011", &entries::Bzero_to_Kstarzero_dimuon_A_T_2_1_to_6_CDF_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_T^2[14.18,16.00]@CDF-2011", &entries::Bzero_to_Kstarzero_dimuon_A_T_2_14dot18_to_16_CDF_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_T^2[16.00,19.21]@CDF-2011", &entries::Bzero_to_Kstarzero_dimuon_A_T_2_16_to_19dot21_CDF_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_im[1.00,6.00]@CDF-2011", &entries::Bzero_to_Kstarzero_dimuon_A_im_1_to_6_CDF_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_im[14.18,16.00]@CDF-2011", &entries::Bzero_to_Kstarzero_dimuon_A_im_14dot18_to_16_CDF_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_im[16.00,19.21]@CDF-2011", &entries::Bzero_to_Kstarzero_dimuon_A_im_16_to_19dot21_CDF_2011 },
            ValueType{ "B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@CDF-2011", &entries::Bplus_to_Kplus_dimuon_A_FB_1_to_6_CDF_2011 },
            ValueType{ "B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@CDF-2011", &entries::Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_CDF_2011 },
            ValueType{ "B^+->K^+mu^+mu^-::A_FB[16.00,22.86]@CDF-2011", &entries::Bplus_to_Kplus_dimuon_A_FB_16_to_22dot86_CDF_2011 },

            // limit on B^0_s -> mu^+ mu^-
            ValueType{ "B^0_s->mu^+mu^-::BR_limit@CDF-2011", &entries::Bzero_to_dimuon_CDF_2011 },
            // LHCb
            ValueType{ "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@LHCb-2011", &entries::Bzero_to_Kstarzero_dimuon_BR_1_to_6_LHCb_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@LHCb-2011", &entries::Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_LHCb_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@LHCb-2011", &entries::Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_LHCb_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@LHCb-2011", &entries::Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_LHCb_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@LHCb-2011", &entries::Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_LHCb_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@LHCb-2011", &entries::Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_LHCb_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@LHCb-2011", &entries::Bzero_to_Kstarzero_dimuon_F_L_1_to_6_LHCb_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@LHCb-2011", &entries::Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_LHCb_2011 },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@LHCb-2011", &entries::Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_LHCb_2011 },

            // limit on B^0_s -> mu^+ mu^-
            // LHCb + CMS
            ValueType{ "B^0_s->mu^+mu^-::BR_limit@LHCb-CMS-2011", &entries::Bzero_to_dimuon_LHCb_CMS_2011 },
            ValueType{ "B^0_s->mu^+mu^-::BR_limit@LHCb-CMS-2011-Bayes", &entries::Bzero_to_dimuon_LHCb_CMS_2011_Bayes },

            /* 2012 */
            // BaBar
            // B^+ -> K^+ mu^+ mu^-
            ValueType{ "B^+->K^+mu^+mu^-::BR[1.00,6.00]@BaBar-2012", &entries::Bplus_to_Kplus_dimuon_BR_1_to_6_BaBar_2012 },
            ValueType{ "B^+->K^+mu^+mu^-::BR[14.21,16.00]@BaBar-2012", &entries::Bplus_to_Kplus_dimuon_BR_14dot21_to_16_BaBar_2012 },
            ValueType{ "B^+->K^+mu^+mu^-::BR[16.00,22.86]@BaBar-2012", &entries::Bplus_to_Kplus_dimuon_BR_16_to_22dot86_BaBar_2012 },
            // B^0 -> K^*0 mu^+ mu^-
            ValueType{ "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@BaBar-2012", &entries::Bzero_to_Kstarzero_dimuon_BR_1_to_6_BaBar_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::BR[14.21,16.00]@BaBar-2012", &entries::Bzero_to_Kstarzero_dimuon_BR_14dot21_to_16_BaBar_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@BaBar-2012", &entries::Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_BaBar_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@BaBar-2012", &entries::Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_BaBar_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@BaBar-2012", &entries::Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_BaBar_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@BaBar-2012", &entries::Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_BaBar_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@BaBar-2012", &entries::Bzero_to_Kstarzero_dimuon_F_L_1_to_6_BaBar_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@BaBar-2012", &entries::Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_BaBar_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@BaBar-2012", &entries::Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_BaBar_2012 },
            // B -> X_s gamma
            ValueType{ "B->X_sgamma::BR[1.8]@BaBar-2012", &entries::B_to_Xs_gamma_BR_1dot8_BaBar_2012C },
            ValueType{ "B->X_sgamma::E_1[1.8]@BaBar-2012", &entries::B_to_Xs_gamma_E_1_1dot8_BaBar_2012C },
            ValueType{ "B->X_sgamma::E_2[1.8]@BaBar-2012", &entries::B_to_Xs_gamma_E_2_1dot8_BaBar_2012C },
            ValueType{ "B->X_sgamma::E_1[1.8]+E_2[1.8]@BaBar-2012", &entries::B_to_Xs_gamma_E_1_and_E_2_1dot8_BaBar_2012C },
            // B^0 -> pi^- l nu
            ValueType{ "B^0->pi^+lnu::BR@BaBar-2012D", &entries::Bzero_to_pi_l_nu_BR_BaBar_2012D },
            // CDF
            // B^0 -> K^0 mu^+ mu^-
            ValueType{ "B^0->K^0mu^+mu^-::BR[1.00,6.00]@CDF-2012", &entries::Bzero_to_Kzero_dimuon_BR_1_to_6_CDF_2012 },
            ValueType{ "B^0->K^0mu^+mu^-::BR[14.18,16.00]@CDF-2012", &entries::Bzero_to_Kzero_dimuon_BR_14dot18_to_16_CDF_2012 },
            ValueType{ "B^0->K^0mu^+mu^-::BR[16.00,22.86]@CDF-2012", &entries::Bzero_to_Kzero_dimuon_BR_16_to_22dot86_CDF_2012 },
            // B^+ -> K^+ mu^+ mu^-
            ValueType{ "B^+->K^+mu^+mu^-::BR[1.00,6.00]@CDF-2012", &entries::Bplus_to_Kplus_dimuon_BR_1_to_6_CDF_2012 },
            ValueType{ "B^+->K^+mu^+mu^-::BR[14.18,16.00]@CDF-2012", &entries::Bplus_to_Kplus_dimuon_BR_14dot18_to_16_CDF_2012 },
            ValueType{ "B^+->K^+mu^+mu^-::BR[16.00,22.86]@CDF-2012", &entries::Bplus_to_Kplus_dimuon_BR_16_to_22dot86_CDF_2012 },
            // B^0 -> K^*0 mu^+ mu^-
            ValueType{ "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@CDF-2012", &entries::Bzero_to_Kstarzero_dimuon_BR_1_to_6_CDF_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@CDF-2012", &entries::Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_CDF_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::BR[16.00,19.21]@CDF-2012", &entries::Bzero_to_Kstarzero_dimuon_BR_16_to_19dot21_CDF_2012 },
            // B^+ -> K^*+ mu^+ mu^-
            ValueType{ "B^+->K^*+mu^+mu^-::BR[1.00,6.00]@CDF-2012", &entries::Bplus_to_Kstarplus_dimuon_BR_1_to_6_CDF_2012 },
            ValueType{ "B^+->K^*+mu^+mu^-::BR[14.18,16.00]@CDF-2012", &entries::Bplus_to_Kstarplus_dimuon_BR_14dot18_to_16_CDF_2012 },
            ValueType{ "B^+->K^*+mu^+mu^-::BR[16.00,19.21]@CDF-2012", &entries::Bplus_to_Kstarplus_dimuon_BR_16_to_19dot21_CDF_2012 },
            // B^0 -> K^*0 mu^+ mu^-
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@CDF-2012", &entries::Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_CDF_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@CDF-2012", &entries::Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_CDF_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[16.00,19.21]@CDF-2012", &entries::Bzero_to_Kstarzero_dimuon_A_FB_16_to_19dot21_CDF_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@CDF-2012", &entries::Bzero_to_Kstarzero_dimuon_F_L_1_to_6_CDF_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@CDF-2012", &entries::Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_CDF_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[16.00,19.21]@CDF-2012", &entries::Bzero_to_Kstarzero_dimuon_F_L_16_to_19dot21_CDF_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_T^2[1.00,6.00]@CDF-2012", &entries::Bzero_to_Kstarzero_dimuon_A_T_2_1_to_6_CDF_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_T^2[14.18,16.00]@CDF-2012", &entries::Bzero_to_Kstarzero_dimuon_A_T_2_14dot18_to_16_CDF_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_T^2[16.00,19.21]@CDF-2012", &entries::Bzero_to_Kstarzero_dimuon_A_T_2_16_to_19dot21_CDF_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_im[1.00,6.00]@CDF-2012", &entries::Bzero_to_Kstarzero_dimuon_A_im_1_to_6_CDF_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_im[14.18,16.00]@CDF-2012", &entries::Bzero_to_Kstarzero_dimuon_A_im_14dot18_to_16_CDF_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_im[16.00,19.21]@CDF-2012", &entries::Bzero_to_Kstarzero_dimuon_A_im_16_to_19dot21_CDF_2012 },
            // B^+ -> K^+ mu^+ mu^-
            ValueType{ "B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@CDF-2012", &entries::Bplus_to_Kplus_dimuon_A_FB_1_to_6_CDF_2012 },
            ValueType{ "B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@CDF-2012", &entries::Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_CDF_2012 },
            ValueType{ "B^+->K^+mu^+mu^-::A_FB[16.00,22.86]@CDF-2012", &entries::Bplus_to_Kplus_dimuon_A_FB_16_to_22dot86_CDF_2012 },

            // LHCb
            // B^0 -> K^*0 mu^+ mu^-
            ValueType{ "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@LHCb-2012", &entries::Bzero_to_Kstarzero_dimuon_BR_1_to_6_LHCb_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@LHCb-2012", &entries::Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_LHCb_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::BR[16.00,19.00]@LHCb-2012", &entries::Bzero_to_Kstarzero_dimuon_BR_16_to_19_LHCb_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@LHCb-2012", &entries::Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_LHCb_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@LHCb-2012", &entries::Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_LHCb_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[16.00,19.00]@LHCb-2012", &entries::Bzero_to_Kstarzero_dimuon_A_FB_16_to_19_LHCb_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@LHCb-2012", &entries::Bzero_to_Kstarzero_dimuon_F_L_1_to_6_LHCb_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@LHCb-2012", &entries::Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_LHCb_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[16.00,19.00]@LHCb-2012", &entries::Bzero_to_Kstarzero_dimuon_F_L_16_to_19_LHCb_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_im[1.00,6.00]@LHCb-2012", &entries::Bzero_to_Kstarzero_dimuon_A_im_1_to_6_LHCb_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_im[14.18,16.00]@LHCb-2012", &entries::Bzero_to_Kstarzero_dimuon_A_im_14dot18_to_16_LHCb_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::A_im[16.00,19.00]@LHCb-2012", &entries::Bzero_to_Kstarzero_dimuon_A_im_16_to_19_LHCb_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::S_3[1.00,6.00]@LHCb-2012", &entries::Bzero_to_Kstarzero_dimuon_S_3_1_to_6_LHCb_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::S_3[14.18,16.00]@LHCb-2012", &entries::Bzero_to_Kstarzero_dimuon_S_3_14dot18_to_16_LHCb_2012 },
            ValueType{ "B^0->K^*0mu^+mu^-::S_3[16.00,19.00]@LHCb-2012", &entries::Bzero_to_Kstarzero_dimuon_S_3_16_to_19_LHCb_2012 },
            // B^+ -> K^+ mu^+ mu^-
            ValueType{ "B^+->K^+mu^+mu^-::BR[1.00,6.00]@LHCb-2012", &entries::Bplus_to_Kplus_dimuon_BR_1_to_6_LHCb_2012 },
            ValueType{ "B^+->K^+mu^+mu^-::BR[14.18,16.00]@LHCb-2012", &entries::Bplus_to_Kplus_dimuon_BR_14dot18_to_16_LHCb_2012 },
            ValueType{ "B^+->K^+mu^+mu^-::BR[16.00,18.00]@LHCb-2012", &entries::Bplus_to_Kplus_dimuon_BR_16_to_18_LHCb_2012 },
            ValueType{ "B^+->K^+mu^+mu^-::BR[18.00,22.00]@LHCb-2012", &entries::Bplus_to_Kplus_dimuon_BR_18_to_22_LHCb_2012 },

            ValueType{ "B^+->K^+mu^+mu^-::A_FB[1.00,6.00]@LHCb-2012", &entries::Bplus_to_Kplus_dimuon_A_FB_1_to_6_LHCb_2012 },
            ValueType{ "B^+->K^+mu^+mu^-::A_FB[14.18,16.00]@LHCb-2012", &entries::Bplus_to_Kplus_dimuon_A_FB_14dot18_to_16_LHCb_2012 },
            ValueType{ "B^+->K^+mu^+mu^-::A_FB[16.00,18.00]@LHCb-2012", &entries::Bplus_to_Kplus_dimuon_A_FB_16_to_18_LHCb_2012 },
            ValueType{ "B^+->K^+mu^+mu^-::A_FB[18.00,22.00]@LHCb-2012", &entries::Bplus_to_Kplus_dimuon_A_FB_18_to_22_LHCb_2012 },

            ValueType{ "B^+->K^+mu^+mu^-::F_H[1.00,6.00]@LHCb-2012", &entries::Bplus_to_Kplus_dimuon_F_H_1_to_6_LHCb_2012 },
            ValueType{ "B^+->K^+mu^+mu^-::F_H[14.18,16.00]@LHCb-2012", &entries::Bplus_to_Kplus_dimuon_F_H_14dot18_to_16_LHCb_2012 },
            ValueType{ "B^+->K^+mu^+mu^-::F_H[16.00,18.00]@LHCb-2012", &entries::Bplus_to_Kplus_dimuon_F_H_16_to_18_LHCb_2012 },
            ValueType{ "B^+->K^+mu^+mu^-::F_H[18.00,22.00]@LHCb-2012", &entries::Bplus_to_Kplus_dimuon_F_H_18_to_22_LHCb_2012 },
            // limit on B^0_s -> mu^+ mu^-
            ValueType{ "B^0_s->mu^+mu^-::BR_limit@LHCb-2012", &entries::Bzero_to_dimuon_LHCb_2012 },
            // limit on B^0_s -> mu^+ mu^- of Nov 2012
            ValueType{ "B^0_s->mu^+mu^-::BR_limit@LHCb-Nov-2012", &entries::Bzero_to_dimuon_LHCb_Nov_2012 },
            // B^0 -> K^*0 mu^+ mu^-
            ValueType{ "B^0->K^*0mu^+mu^-::A_CP[1.00,6.00]@LHCb-2012", &entries::Bzero_to_Kstarzero_dimuon_A_CP_1_to_6_LHCb_2012E },
            ValueType{ "B^0->K^*0mu^+mu^-::A_CP[14.18,16.00]@LHCb-2012", &entries::Bzero_to_Kstarzero_dimuon_A_CP_14dot18_to_16_LHCb_2012E },
            ValueType{ "B^0->K^*0mu^+mu^-::A_CP[16.00,20.00]@LHCb-2012", &entries::Bzero_to_Kstarzero_dimuon_A_CP_16_to_20_LHCb_2012E },
            // PDG2012
            // B^0(*) Mass splitting
            ValueType{ "B^0::M_B^*-M_B@PDG-2012", &entries::B_Bstar_mass_splitting_PDG_2012 },

            /* 2013 */
            // ATLAS
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@ATLAS-2013A", &entries::Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_ATLAS_2013A },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@ATLAS-2013A", &entries::Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_ATLAS_2013A },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[16.00,19.00]@ATLAS-2013A", &entries::Bzero_to_Kstarzero_dimuon_A_FB_16_to_19_ATLAS_2013A },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@ATLAS-2013A", &entries::Bzero_to_Kstarzero_dimuon_F_L_1_to_6_ATLAS_2013A },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@ATLAS-2013A", &entries::Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_ATLAS_2013A },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[16.00,19.00]@ATLAS-2013A", &entries::Bzero_to_Kstarzero_dimuon_F_L_16_to_19_ATLAS_2013A },
            // BaBar
            ValueType{ "B->X_sll::BR[1.0,6.0]@BaBar-2013A", &entries::Bmix_to_Xs_dilepton_BR_BaBar_2013A },
            //Belle
            ValueType{ "B^0->pi^+lnu::BR@Belle-2013A", &entries::Bzero_to_pi_l_nu_BR_Belle_2013A },
            // CMS
            ValueType{ "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@CMS-2013A", &entries::Bzero_to_Kstarzero_dimuon_BR_1_to_6_CMS_2013A },
            ValueType{ "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@CMS-2013A", &entries::Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_CMS_2013A },
            ValueType{ "B^0->K^*0mu^+mu^-::BR[16.00,19.00]@CMS-2013A", &entries::Bzero_to_Kstarzero_dimuon_BR_16_to_19_CMS_2013A },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@CMS-2013A", &entries::Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_CMS_2013A },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@CMS-2013A", &entries::Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_CMS_2013A },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[16.00,19.00]@CMS-2013A", &entries::Bzero_to_Kstarzero_dimuon_A_FB_16_to_19_CMS_2013A },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@CMS-2013A", &entries::Bzero_to_Kstarzero_dimuon_F_L_1_to_6_CMS_2013A },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@CMS-2013A", &entries::Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_CMS_2013A },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[16.00,19.00]@CMS-2013A", &entries::Bzero_to_Kstarzero_dimuon_F_L_16_to_19_CMS_2013A },
            ValueType{ "B^0_s->mu^+mu^-::BR@CMS-2013B", &entries::Bzero_to_dimuon_CMS_2013B },
            // LHCb
            ValueType{ "B^0->K^*0mu^+mu^-::BR[1.00,6.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_BR_1_to_6_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::BR[14.18,16.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_BR_14dot18_to_16_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::BR[16.00,19.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_BR_16_to_19_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[1.00,6.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_F_L_1_to_6_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[14.18,16.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_F_L_14dot18_to_16_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::F_L[16.00,19.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_F_L_16_to_19_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[1.00,6.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_A_FB_1_to_6_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[14.18,16.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_A_FB_14dot18_to_16_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::A_FB[16.00,19.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_A_FB_16_to_19_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::S_3[1.00,6.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_S_3_1_to_6_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::S_3[14.18,16.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_S_3_14dot18_to_16_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::S_3[16.00,19.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_S_3_16_to_19_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::S_9[1.00,6.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_S_9_1_to_6_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::S_9[14.18,16.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_S_9_14dot18_to_16_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::S_9[16.00,19.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_S_9_16_to_19_LHCb_2013B },
            // The following observables have not yet been implemented.
#if 0
            ValueType{ "B^0->K^*0mu^+mu^-::A_9[1.00,6.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_A_9_1_to_6_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::A_9[14.18,16.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_A_9_14dot18_to_16_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::A_9[16.00,19.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_A_9_16_to_19_LHCb_2013B },
#endif
            ValueType{ "B^0->K^*0mu^+mu^-::A_T^2[1.00,6.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_A_T2_1_to_6_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::A_T^2[14.18,16.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_A_T2_14dot18_to_16_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::A_T^2[16.00,19.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_A_T2_16_to_19_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::A_T^re[1.00,6.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_A_Tre_1_to_6_LHCb_2013B },
            // The following constraint does not conform to a gaussian likelihood.
            ValueType{ "B^0->K^*0mu^+mu^-::A_T^re[14.18,16.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_A_Tre_14dot18_to_16_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::A_T^re[16.00,19.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_A_Tre_16_to_19_LHCb_2013B },
            ValueType{ "B^0->K^*0mu^+mu^-::P'_4[1.00,6.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_Pprime_4_1_to_6_LHCb_2013C },
            ValueType{ "B^0->K^*0mu^+mu^-::P'_4[14.18,16.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_Pprime_4_14dot18_to_16_LHCb_2013C },
            ValueType{ "B^0->K^*0mu^+mu^-::P'_4[16.00,19.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_Pprime_4_16_to_19_LHCb_2013C },
            ValueType{ "B^0->K^*0mu^+mu^-::P'_5[1.00,6.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_Pprime_5_1_to_6_LHCb_2013C },
            ValueType{ "B^0->K^*0mu^+mu^-::P'_5[14.18,16.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_Pprime_5_14dot18_to_16_LHCb_2013C },
            ValueType{ "B^0->K^*0mu^+mu^-::P'_5[16.00,19.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_Pprime_5_16_to_19_LHCb_2013C },
            ValueType{ "B^0->K^*0mu^+mu^-::P'_6[1.00,6.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_Pprime_6_1_to_6_LHCb_2013C },
            ValueType{ "B^0->K^*0mu^+mu^-::P'_6[14.18,16.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_Pprime_6_14dot18_to_16_LHCb_2013C },
            ValueType{ "B^0->K^*0mu^+mu^-::P'_6[16.00,19.00]@LHCb-2013", &entries::Bzero_to_Kstarzero_dimuon_Pprime_6_16_to_19_LHCb_2013C },
            ValueType{ "B^0_s->mu^+mu^-::BR@LHCb-2013D", &entries::Bzero_to_dimuon_LHCb_2013D },

            /* 2014 */
            // LHCb
            ValueType{ "B^+->K^+mu^+mu^-::BR[1.10,2.00]@LHCb-2014", &entries::Bplus_to_Kplus_dimuon_BR_1dot1_to_2_LHCb_2014 },
            ValueType{ "B^+->K^+mu^+mu^-::BR[2.00,3.00]@LHCb-2014", &entries::Bplus_to_Kplus_dimuon_BR_2_to_3_LHCb_2014 },
            ValueType{ "B^+->K^+mu^+mu^-::BR[3.00,4.00]@LHCb-2014", &entries::Bplus_to_Kplus_dimuon_BR_3_to_4_LHCb_2014 },
            ValueType{ "B^+->K^+mu^+mu^-::BR[4.00,5.00]@LHCb-2014", &entries::Bplus_to_Kplus_dimuon_BR_4_to_5_LHCb_2014 },
            ValueType{ "B^+->K^+mu^+mu^-::BR[5.00,6.00]@LHCb-2014", &entries::Bplus_to_Kplus_dimuon_BR_5_to_6_LHCb_2014 },
            ValueType{ "B^+->K^+mu^+mu^-::BR[1.10,6.00]@LHCb-2014", &entries::Bplus_to_Kplus_dimuon_BR_1dot1_to_6_LHCb_2014 },
            ValueType{ "B^+->K^+mu^+mu^-::BR[15.00,22.00]@LHCb-2014", &entries::Bplus_to_Kplus_dimuon_BR_15_to_22_LHCb_2014 },
            ValueType{ "B^+->K^+mu^+mu^-::A_FB[1.10,6.00]@LHCb-2014", &entries::Bplus_to_Kplus_dimuon_A_FB_1dot1_to_6_LHCb_2014 },
            ValueType{ "B^+->K^+mu^+mu^-::A_FB[15.00,22.00]@LHCb-2014", &entries::Bplus_to_Kplus_dimuon_A_FB_15_to_22_LHCb_2014 },
            ValueType{ "B^+->K^+mu^+mu^-::F_H[1.10,6.00]@LHCb-2014", &entries::Bplus_to_Kplus_dimuon_F_H_1dot1_to_6_LHCb_2014 },
            ValueType{ "B^+->K^+mu^+mu^-::F_H[15.00,22.00]@LHCb-2014", &entries::Bplus_to_Kplus_dimuon_F_H_15_to_22_LHCb_2014 },
            ValueType{ "B^+->K^+l^+l^-::R_K[1.00,6.00]@LHCb-2014", &entries::Bplus_to_Kplus_dilepton_r_k_1_to_6_LHCb_2014C },
            // CMS + LHCb
            ValueType{ "B^0_s->mu^+mu^-::BR@CMS-LHCb-2014", &entries::B_s_to_dimuon_CMS_LHCb_2014A },
            ValueType{ "B^0_d->mu^+mu^-::BR@CMS-LHCb-2014", &entries::Bzero_to_dimuon_CMS_LHCb_2014A },

            /* 2015 */
            // LHCb
            ValueType{ "B^0->K^*0mu^+mu^-::AngularObservables[1.10,2.00]@LHCb-2015A", &entries::Bzero_to_Kstarzero_dimuon_aobs_moments_1dot1_to_2_LHCb_2015A },
            ValueType{ "B^0->K^*0mu^+mu^-::AngularObservables[2.00,3.00]@LHCb-2015A", &entries::Bzero_to_Kstarzero_dimuon_aobs_moments_2_to_3_LHCb_2015A },
            ValueType{ "B^0->K^*0mu^+mu^-::AngularObservables[3.00,4.00]@LHCb-2015A", &entries::Bzero_to_Kstarzero_dimuon_aobs_moments_3_to_4_LHCb_2015A },
            ValueType{ "B^0->K^*0mu^+mu^-::AngularObservables[4.00,5.00]@LHCb-2015A", &entries::Bzero_to_Kstarzero_dimuon_aobs_moments_4_to_5_LHCb_2015A },
            ValueType{ "B^0->K^*0mu^+mu^-::AngularObservables[5.00,6.00]@LHCb-2015A", &entries::Bzero_to_Kstarzero_dimuon_aobs_moments_5_to_6_LHCb_2015A },
            // LHCb
            ValueType{ "Lambda_b->Lambdamu^+mu^-::BR[15.0,20.0]@LHCb-2015B", &entries::Lambdab_to_Lambda_dimuon_br_15_to_20_LHCb_2015B },
            ValueType{ "Lambda_b->Lambdamu^+mu^-::F_0[15.0,20.0]@LHCb-2015B", &entries::Lambdab_to_Lambda_dimuon_f_0_15_to_20_LHCb_2015B },
            ValueType{ "Lambda_b->Lambdamu^+mu^-::A_FB^l[15.0,20.0]@LHCb-2015B", &entries::Lambdab_to_Lambda_dimuon_a_fb_l_15_to_20_LHCb_2015B },
            ValueType{ "Lambda_b->Lambdamu^+mu^-::A_FB^h[15.0,20.0]@LHCb-2015B", &entries::Lambdab_to_Lambda_dimuon_a_fb_h_15_to_20_LHCb_2015B },

            /* Theory Constraints */
            ValueType{ "B->K::f_0+f_++f_T@HPQCD-2013A", &entries::B_to_K_fzero_fplus_ftensor_17_to_23_HPQCD_2013A },
            ValueType{ "B->K^*::V@MILC-2013A", &entries::B_to_Kstar_V_15_to_19dot21_MILC_2013A },
            ValueType{ "B->K^*::A_1@MILC-2013A", &entries::B_to_Kstar_A1_15_to_19dot21_MILC_2013A },
            ValueType{ "B->K^*::A_12@MILC-2013A", &entries::B_to_Kstar_A12_15_to_19dot21_MILC_2013A },

            ValueType{ "B_s->K^*::V@MILC-2013A", &entries::Bs_to_Kstar_V_15_to_19dot21_MILC_2013A },
            ValueType{ "B_s->K^*::A_1@MILC-2013A", &entries::Bs_to_Kstar_A1_15_to_19dot21_MILC_2013A },
            ValueType{ "B_s->K^*::A_12@MILC-2013A", &entries::Bs_to_Kstar_A12_15_to_19dot21_MILC_2013A },

            ValueType{ "B->pi::f_+@IKMvD-2014", &entries::B_to_pi_fp_IKMvD_2014 },

            ValueType{ "Lambda_b->Lambda::f_perp+long^V+A@BFvD2014", &entries::LambdaB_to_Lambda_all_v_and_a_0_BFvD2014 },
            ValueType{ "Lambda_b->Lambda::f_perp^V@BFvD2014", &entries::LambdaB_to_Lambda_fperpV_13dot5_to_20dot3_BFvD2014 },
            ValueType{ "Lambda_b->Lambda::f_perp^A@BFvD2014", &entries::LambdaB_to_Lambda_fperpA_13dot5_to_20dot3_BFvD2014 },
            ValueType{ "Lambda_b->Lambda::f_long^V@BFvD2014", &entries::LambdaB_to_Lambda_flongV_13dot5_to_20dot3_BFvD2014 },
            ValueType{ "Lambda_b->Lambda::f_long^A@BFvD2014", &entries::LambdaB_to_Lambda_flongA_13dot5_to_20dot3_BFvD2014 },

            ValueType{ "B->K^*::V+A_0+A_1+A_2[LCSR]@BSZ2015", &entries::B_to_Kstar_V_A0_A1_A2_0dot1_to_12dot1_BSZ_2015_lcsr_parameters },
            ValueType{ "B->K^*::V+A_0+A_1+A_2[LCSR+Lattice]@BSZ2015", &entries::B_to_Kstar_V_A0_A1_A2_0dot1_to_12dot1_BSZ_2015_lcsrlattice_parameters },

            ValueType{ "B->K^*::V+A_0+A_1+A_12@HLMW2015", &entries::B_to_Kstar_V_A0_A1_A12_11dot9_to_17dot8_HLMW_2015 },

            ValueType{ "B->D::f_++f_0@HPQCD2015A", &entries::B_to_D_f_plus_f_zero_0_to_11dot62_HPQCD_2015A },
            ValueType{ "B->pi::f_++f_0@FNALMILC2015A", &entries::B_to_pi_f_plus_f_zero_18_to_26_FNALMILC_2015A },

            ValueType{ "Lambda_b->Lambda::f_perp+long^V+A+T+T5@DM2016", &entries::Lambda_b_to_Lambda_v_a_t_t5_parameters_DM_2016 },
        };

        return constraint_entries;
    }

    /*
     * Adding a new constraint:
     * 1. Instantiate an existing ConstraintEntry in namespace entries{...}
     * 2. Add an entry to the map in make_constraint_entries
     * 4. Run constraint_TEST and check text output for new constraint
     */
    Constraint
    Constraint::make(const QualifiedName & name, const Options & options)
    {
        auto & entries = make_constraint_entries();

        auto e = entries.find(name);
        if (e == entries.end())
            throw UnknownConstraintError(name);

        return e->second->make(e->first, options);
    }

    template <>
    struct WrappedForwardIteratorTraits<Constraints::ConstraintIteratorTag>
    {
        typedef std::map<QualifiedName, const ConstraintEntry *>::const_iterator UnderlyingIterator;
    };
    template class WrappedForwardIterator<Constraints::ConstraintIteratorTag, const std::pair<const QualifiedName, const ConstraintEntry *>>;

    template<>
    struct Implementation<Constraints>
    {
        const std::map<QualifiedName, const ConstraintEntry *> constraint_entries;

        Implementation() :
            constraint_entries(make_constraint_entries())
        {
        }
    };

    Constraints::Constraints() :
        PrivateImplementationPattern<Constraints>(new Implementation<Constraints>())
    {
    }

    Constraints::~Constraints()
    {
    }

    Constraints::ConstraintIterator
    Constraints::begin() const
    {
        return ConstraintIterator(_imp->constraint_entries.cbegin());
    }

    Constraints::ConstraintIterator
    Constraints::end() const
    {
        return ConstraintIterator(_imp->constraint_entries.cend());
    }
}
