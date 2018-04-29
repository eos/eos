/* vim: set sw=4 sts=4 et foldmethod=marker foldmarker={{{,}}} : */

/*
 * Copyright (c) 2011-2018 Danny van Dyk
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

#include <config.h>

#include <eos/constraint.hh>
#include <eos/statistics/log-likelihood.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/observable_set.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qualified-name.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <yaml-cpp/yaml.h>

#include <cmath>
#include <map>
#include <vector>

namespace fs = boost::filesystem;

namespace eos
{
    UnknownConstraintError::UnknownConstraintError(const QualifiedName & name) :
        Exception("Constraint '" + name.str() + "' is unknown")
    {
    }

    ConstraintDeserializationError::ConstraintDeserializationError(const QualifiedName & name, const std::string & msg) :
        Exception("Could not deserialize entry '" + name.str() + "': " + msg)
    {
    }

    ConstraintInputFileParseError::ConstraintInputFileParseError(const std::string & filename, const std::string & msg) :
        Exception("Could not parse constraint input file '" + filename + "': " + msg)
    {
    }

    template <>
    struct WrappedForwardIteratorTraits<ConstraintEntry::ObservableNameIteratorTag>
    {
        typedef std::vector<QualifiedName>::const_iterator UnderlyingIterator;
    };
    template class WrappedForwardIterator<ConstraintEntry::ObservableNameIteratorTag, const QualifiedName>;

    namespace impl
    {
        static bool less(const std::pair<YAML::Node, YAML::Node> & lhs, const std::pair<YAML::Node, YAML::Node> & rhs)
        {
            return lhs.first.as<std::string>() < rhs.first.as<std::string>();
        }
    }

    /// {{{ ConstraintEntryBase
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

            ConstraintEntryBase(const QualifiedName & name,
                    const std::vector<QualifiedName> & observable_names) :
                _name(name),
                _observable_names(observable_names)
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

            virtual void serialize(YAML::Emitter & out) const
            {
                out << YAML::BeginMap;
                out << YAML::Key << "type";
                out << YAML::Value << "unimplemented";
                out << YAML::EndMap;
            }

            virtual ConstraintEntry * deserialize(const YAML::Node &) const
            {
                return nullptr;
            }
    };
    /// }}}

    /// {{{ GaussianConstraintEntry
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

        virtual void serialize(YAML::Emitter & out) const
        {
            out << YAML::BeginMap;
            out << YAML::Key << "type" << YAML::Value << "Gaussian";
            out << YAML::Key << "observable" << YAML::Value << observable.full();
            out << YAML::Key << "kinematics" << YAML::Value << YAML::Flow << YAML::BeginMap;
            for (const auto & k : kinematics)
            {
                out << YAML::Key << k.name() << YAML::Value << k.evaluate();
            }
            out << YAML::EndMap;
            out << YAML::Key << "options" << YAML::Value << YAML::Flow << YAML::BeginMap;
            for (const auto & o : options)
            {
                out << YAML::Key << o.first << YAML::Value << o.second;
            }
            out << YAML::EndMap;
            out << YAML::Key << "mean" << YAML::Value << central;
            out << YAML::Key << "sigma-stat" << YAML::Flow << YAML::BeginMap;
            out << YAML::Key << "hi" << YAML::Value << sigma_hi_stat;
            out << YAML::Key << "lo" << YAML::Value << sigma_lo_stat;
            out << YAML::EndMap;
            out << YAML::Key << "sigma-sys" << YAML::Flow << YAML::BeginMap;
            out << YAML::Key << "hi" << YAML::Value << sigma_hi_sys;
            out << YAML::Key << "lo" << YAML::Value << sigma_lo_sys;
            out << YAML::EndMap;
            out << YAML::Key << "dof" << YAML::Value << number_of_observations;
            out << YAML::EndMap;
        }

        static ConstraintEntry * deserialize(const QualifiedName & name, const YAML::Node & n)
        {
            static const std::string required_keys[] =
            {
                "observable", "kinematics", "options", "mean", "sigma-stat", "sigma-sys", "dof"
            };

            for (auto && k : required_keys)
            {
                if (! n[k].IsDefined())
                {
                    throw ConstraintDeserializationError(name, "required key '" + k + "' not specified");
                }
            }

            static const std::string scalar_keys[] =
            {
                "observable", "mean", "dof"
            };

            for (auto && k : scalar_keys)
            {
                if (YAML::NodeType::Scalar != n[k].Type())
                {
                    throw ConstraintDeserializationError(name, "required key '" + k + "' not mapped to a scalar value");
                }
            }

            static const std::string map_keys[] =
            {
                "kinematics", "options", "sigma-stat", "sigma-sys"
            };

            for (auto && k : map_keys)
            {
                if (YAML::NodeType::Map != n[k].Type())
                {
                    throw ConstraintDeserializationError(name, "required key '" + k + "' not mapped to a map");
                }
            }

            try
            {
                QualifiedName observable(n["observable"].as<std::string>());
                double mean = n["mean"].as<double>();
                unsigned dof = n["dof"].as<unsigned>();

                Kinematics kinematics;
                std::list<std::pair<YAML::Node, YAML::Node>> kinematics_nodes(n["kinematics"].begin(), n["kinematics"].end());
                // yaml-cpp does not guarantee loading of a map in the order it is written. Circumvent this problem
                // by sorting the entries lexicographically.
                kinematics_nodes.sort(&impl::less);
                std::set<std::string> kinematics_keys;
                for (auto && k : kinematics_nodes)
                {
                    std::string key = k.first.as<std::string>();
                    if (! kinematics_keys.insert(key).second)
                        throw ConstraintDeserializationError(name, "kinematics key '" + key + "' encountered more than once");

                    kinematics.declare(key, k.second.as<double>());
                }

                Options options;
                std::list<std::pair<YAML::Node, YAML::Node>> options_nodes(n["options"].begin(), n["options"].end());
                // yaml-cpp does not guarantee loading of a map in the order it is written. Circumvent this problem
                // by sorting the entries lexicographically.
                options_nodes.sort(&impl::less);
                std::set<std::string> options_keys;
                for (auto && o : options_nodes)
                {
                    std::string key = o.first.as<std::string>();
                    if (! options_keys.insert(key).second)
                        throw ConstraintDeserializationError(name, "options key '" + key + "' encountered more than once");

                    options.set(key, o.second.as<std::string>());
                }

                double sigma_hi_stat = n["sigma-stat"]["hi"].as<double>();
                double sigma_lo_stat = n["sigma-stat"]["lo"].as<double>();
                double sigma_hi_sys  = n["sigma-sys"]["hi"].as<double>();
                double sigma_lo_sys  = n["sigma-sys"]["lo"].as<double>();

                return new GaussianConstraintEntry(name.str(), observable, kinematics, options, mean,
                        sigma_hi_stat, sigma_lo_stat, sigma_hi_sys, sigma_lo_sys, dof);
            }
            catch (QualifiedNameSyntaxError & e)
            {
                throw ConstraintDeserializationError(name, "'" + n["observable"].as<std::string>() + "' is not a valid observable name (" + e.what() + ")");
            }
        }
    };
    /// }}}

    /// {{{  LogGammaConstraintEntry
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

        virtual void serialize(YAML::Emitter & out) const
        {
            out << YAML::BeginMap;
            out << YAML::Key << "type" << YAML::Value << "LogGamma";
            out << YAML::Key << "observable" << YAML::Value << observable.full();
            out << YAML::Key << "kinematics" << YAML::Value << YAML::Flow << YAML::BeginMap;
            for (const auto & k : kinematics)
            {
                out << YAML::Key << k.name() << YAML::Value << k.evaluate();
            }
            out << YAML::EndMap;
            out << YAML::Key << "options" << YAML::Value << YAML::Flow << YAML::BeginMap;
            for (const auto & o : options)
            {
                out << YAML::Key << o.first << YAML::Value << o.second;
            }
            out << YAML::EndMap;
            out << YAML::Key << "mode" << YAML::Value << central;
            out << YAML::Key << "sigma" << YAML::Flow << YAML::BeginMap;
            out << YAML::Key << "hi" << YAML::Value << sigma_hi;
            out << YAML::Key << "lo" << YAML::Value << sigma_lo;
            out << YAML::EndMap;
            out << YAML::Key << "dof" << YAML::Value << number_of_observations;
            out << YAML::EndMap;
        }

        static ConstraintEntry * deserialize(const QualifiedName & name, const YAML::Node & n)
        {
            static const std::string required_keys[] =
            {
                "observable", "kinematics", "options", "mode", "sigma", "dof"
            };

            for (auto && k : required_keys)
            {
                if (! n[k].IsDefined())
                {
                    throw ConstraintDeserializationError(name, "required key '" + k + "' not specified");
                }
            }

            static const std::string scalar_keys[] =
            {
                "observable", "mode", "dof"
            };

            for (auto && k : scalar_keys)
            {
                if (YAML::NodeType::Scalar != n[k].Type())
                {
                    throw ConstraintDeserializationError(name, "required key '" + k + "' not mapped to a scalar value");
                }
            }

            static const std::string map_keys[] =
            {
                "kinematics", "options", "sigma"
            };

            for (auto && k : map_keys)
            {
                if (YAML::NodeType::Map != n[k].Type())
                {
                    throw ConstraintDeserializationError(name, "required key '" + k + "' not mapped to a map");
                }
            }

            try
            {
                QualifiedName observable(n["observable"].as<std::string>());
                double mode = n["mode"].as<double>();
                unsigned dof = n["dof"].as<unsigned>();

                Kinematics kinematics;
                std::list<std::pair<YAML::Node, YAML::Node>> kinematics_nodes(n["kinematics"].begin(), n["kinematics"].end());
                // yaml-cpp does not guarantee loading of a map in the order it is written. Circumvent this problem
                // by sorting the entries lexicographically.
                kinematics_nodes.sort(&impl::less);
                std::set<std::string> kinematics_keys;
                for (auto && k : kinematics_nodes)
                {
                    std::string key = k.first.as<std::string>();
                    if (! kinematics_keys.insert(key).second)
                        throw ConstraintDeserializationError(name, "kinematics key '" + key + "' encountered more than once");

                    kinematics.declare(key, k.second.as<double>());
                }

                Options options;
                std::list<std::pair<YAML::Node, YAML::Node>> options_nodes(n["options"].begin(), n["options"].end());
                // yaml-cpp does not guarantee loading of a map in the order it is written. Circumvent this problem
                // by sorting the entries lexicographically.
                options_nodes.sort(&impl::less);
                std::set<std::string> options_keys;
                for (auto && o : options_nodes)
                {
                    std::string key = o.first.as<std::string>();
                    if (! options_keys.insert(key).second)
                        throw ConstraintDeserializationError(name, "options key '" + key + "' encountered more than once");

                    options.set(key, o.second.as<std::string>());
                }

                double sigma_hi = n["sigma"]["hi"].as<double>();
                double sigma_lo = n["sigma"]["lo"].as<double>();

                return new LogGammaConstraintEntry(name.str(), observable, kinematics, options, mode,
                        sigma_hi, sigma_lo, dof);
            }
            catch (QualifiedNameSyntaxError & e)
            {
                throw ConstraintDeserializationError(name, "'" + n["observable"].as<std::string>() + "' is not a valid observable name (" + e.what() + ")");
            }
        }
    };
    /// }}}

    /// {{{ AmorosoConstraintEntry
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

        virtual void serialize(YAML::Emitter & out) const
        {
            out << YAML::BeginMap;
            out << YAML::Key << "type" << YAML::Value << "Amoroso";
            out << YAML::Key << "observable" << YAML::Value << observable.full();
            out << YAML::Key << "kinematics" << YAML::Value << YAML::Flow << YAML::BeginMap;
            for (const auto & k : kinematics)
            {
                out << YAML::Key << k.name() << YAML::Value << k.evaluate();
            }
            out << YAML::EndMap;
            out << YAML::Key << "options" << YAML::Value << YAML::Flow << YAML::BeginMap;
            for (const auto & o : options)
            {
                out << YAML::Key << o.first << YAML::Value << o.second;
            }
            out << YAML::EndMap;
            out << YAML::Key << "physical-limit" << YAML::Value << physical_limit;
            out << YAML::Key << "theta" << YAML::Value << theta;
            out << YAML::Key << "alpha" << YAML::Value << alpha;
            out << YAML::Key << "beta" << YAML::Value << beta;
            out << YAML::EndMap;
        }

        static ConstraintEntry * deserialize(const QualifiedName & name, const YAML::Node & n)
        {
            static const std::string required_keys[] =
            {
                "observable", "kinematics", "options", "physical-limit", "alpha", "beta", "theta"
            };

            for (auto && k : required_keys)
            {
                if (! n[k].IsDefined())
                {
                    throw ConstraintDeserializationError(name, "required key '" + k + "' not specified");
                }
            }

            static const std::string scalar_keys[] =
            {
                "observable", "physical-limit", "alpha", "beta", "theta"
            };

            for (auto && k : scalar_keys)
            {
                if (YAML::NodeType::Scalar != n[k].Type())
                {
                    throw ConstraintDeserializationError(name, "required key '" + k + "' not mapped to a scalar value");
                }
            }

            static const std::string map_keys[] =
            {
                "kinematics", "options"
            };

            for (auto && k : map_keys)
            {
                if (YAML::NodeType::Map != n[k].Type())
                {
                    throw ConstraintDeserializationError(name, "required key '" + k + "' not mapped to a map");
                }
            }

            try
            {
                QualifiedName observable(n["observable"].as<std::string>());

                Kinematics kinematics;
                std::list<std::pair<YAML::Node, YAML::Node>> kinematics_nodes(n["kinematics"].begin(), n["kinematics"].end());
                // yaml-cpp does not guarantee loading of a map in the order it is written. Circumvent this problem
                // by sorting the entries lexicographically.
                kinematics_nodes.sort(&impl::less);
                std::set<std::string> kinematics_keys;
                for (auto && k : kinematics_nodes)
                {
                    std::string key = k.first.as<std::string>();
                    if (! kinematics_keys.insert(key).second)
                        throw ConstraintDeserializationError(name, "kinematics key '" + key + "' encountered more than once");

                    kinematics.declare(key, k.second.as<double>());
                }

                Options options;
                std::list<std::pair<YAML::Node, YAML::Node>> options_nodes(n["options"].begin(), n["options"].end());
                // yaml-cpp does not guarantee loading of a map in the order it is written. Circumvent this problem
                // by sorting the entries lexicographically.
                options_nodes.sort(&impl::less);
                std::set<std::string> options_keys;
                for (auto && o : options_nodes)
                {
                    std::string key = o.first.as<std::string>();
                    if (! options_keys.insert(key).second)
                        throw ConstraintDeserializationError(name, "options key '" + key + "' encountered more than once");

                    options.set(key, o.second.as<std::string>());
                }

                double physical_limit = n["physical-limit"].as<double>();
                double theta          = n["theta"].as<double>();
                double alpha          = n["alpha"].as<double>();
                double beta           = n["beta"].as<double>();

                return new AmorosoConstraintEntry(name.str(), observable, kinematics, options,
                        physical_limit, theta, alpha, beta);
            }
            catch (QualifiedNameSyntaxError & e)
            {
                throw ConstraintDeserializationError(name, "'" + n["observable"].as<std::string>() + "' is not a valid observable name (" + e.what() + ")");
            }
        }
    };
    /// }}}

    /// {{{ MultivariateGaussianConstraintEntry
    struct MultivariateGaussianConstraintEntry :
        public ConstraintEntryBase
    {
        std::vector<QualifiedName> observable_names;

        std::vector<Kinematics> kinematics;

        std::vector<Options> options;

        std::vector<double> means;

        std::vector<double> sigma_stat_hi;
        std::vector<double> sigma_stat_lo;

        std::vector<double> sigma_sys;

        std::vector<std::vector<double>> correlation;

        unsigned number_of_observations;

        unsigned dim;

        MultivariateGaussianConstraintEntry(const QualifiedName & name,
                const std::vector<QualifiedName> & observable_names,
                const std::vector<Kinematics> & kinematics,
                const std::vector<Options> & options,
                const std::vector<double> & means,
                const std::vector<double> & sigma_stat_hi,
                const std::vector<double> & sigma_stat_lo,
                const std::vector<double> & sigma_sys,
                const std::vector<std::vector<double>> & correlation,
                const unsigned number_of_observations) :
            ConstraintEntryBase(name, observable_names),
            observable_names(observable_names),
            kinematics(kinematics),
            options(options),
            means(means),
            sigma_stat_hi(sigma_stat_hi),
            sigma_stat_lo(sigma_stat_lo),
            sigma_sys(sigma_sys),
            correlation(correlation),
            number_of_observations(number_of_observations),
            dim(observable_names.size())
        {
            if (dim != kinematics.size()) { throw InternalError("MultivariateGaussianConstraintEntry: wrong number of kinematics"); }

            if (dim != options.size()) { throw InternalError("MultivariateGaussianConstraintEntry: wrong number of options"); }

            if (dim != means.size()) { throw InternalError("MultivariateGaussianConstraintEntry: wrong number of means"); }

            if (dim != sigma_stat_hi.size()) { throw InternalError("MultivariateGaussianConstraintEntry: wrong number of sigma-stat(hi)"); }

            if (dim != sigma_stat_lo.size()) { throw InternalError("MultivariateGaussianConstraintEntry: wrong number of sigma-stat(lo)"); }

            if (dim != sigma_sys.size()) { throw InternalError("MultivariateGaussianConstraintEntry: wrong number of sigma-sys"); }

            if (dim < number_of_observations) { throw InternalError("MultivariateGaussianConstraintEntry: wrong number of observations"); }

            if (dim != correlation.size()) { throw InternalError("MultviateGaussianConstraintEntry: wrong number of rows in correlation"); }
            for (auto i = 0u ; i < dim ; ++i)
            {
                if (dim != correlation[i].size()) { throw InternalError("MultviateGaussianConstraintEntry: wrong number of columns in correlation row " + stringify(i)); }
            }
        }

        virtual ~MultivariateGaussianConstraintEntry() = default;

        virtual const std::string & type() const
        {
            static const std::string type("MultivariateGaussian<" + stringify(dim) + "> (using correlation matrix)");

            return type;
        }

        virtual Constraint make(const QualifiedName & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);

            std::vector<ObservablePtr> observables(dim, nullptr);
            for (auto i = 0u ; i < dim ; ++i)
            {
                observables[i] = Observable::make(this->observable_names[i], parameters, this->kinematics[i], this->options[i] + options);
                if (! observables[i].get())
                    throw InternalError("make_multivariate_gaussian_constraint<" + stringify(dim) + ">: " + name.str() + ": '" + this->observable_names[i].str() + "' is not a valid observable name");
            }

            std::vector<double> variances(dim, 0.0);
            if ("symmetric+quadratic" == options.get("uncertainty", "symmetric+quadratic"))
            {
                for (auto i = 0u ; i < dim ; ++i)
                {
                    double combined_lo = power_of<2>(sigma_stat_lo[i]) + power_of<2>(sigma_sys[i]);
                    double combined_hi = power_of<2>(sigma_stat_hi[i]) + power_of<2>(sigma_sys[i]);

                    variances[i] = std::max(combined_lo, combined_hi);
                }
            }

            // create GSL vector for the mean
            gsl_vector * means = gsl_vector_alloc(dim);
            for (auto i = 0u ; i < dim ; ++i)
            {
                gsl_vector_set(means, i, this->means[i]);
            }

            // create GSL matrix for the covariance
            gsl_matrix * covariance = gsl_matrix_alloc(dim, dim);
            for (auto i = 0u ; i < dim ; ++i)
            {
                for (auto j = 0u ; j < dim ; ++j)
                {
                    double value = std::sqrt(variances[i] * variances[j]) * correlation[i][j];
                    gsl_matrix_set(covariance, i, j, value);
                }
            }

            LogLikelihoodBlockPtr block = LogLikelihoodBlock::MultivariateGaussian(cache, observables, means, covariance, number_of_observations);

            return Constraint(name, std::vector<ObservablePtr>(observables.begin(), observables.end()), { block });
        }

        virtual std::ostream & insert(std::ostream & os) const
        {
            os << _name.full() << ":" << std::endl;
            os << "    type: MultivariateGaussian<" << dim << ">" << std::endl;

            return os;
        }

        virtual void serialize(YAML::Emitter & out) const
        {
            out << YAML::BeginMap;
            out << YAML::Key << "type" << YAML::Value << "MultivariateGaussian";
            out << YAML::Key << "dim" << YAML::Value << dim;
            out << YAML::Key << "observables" << YAML::Value << YAML::BeginSeq;
            for (const auto & o : observable_names)
            {
                out << o.full();
            }
            out << YAML::EndSeq;
            out << YAML::Key << "kinematics" << YAML::Value << YAML::BeginSeq;
            for (const auto & k : kinematics)
            {
                out << YAML::Flow << YAML::BeginMap;
                for (const auto & kk : k)
                {
                    out << YAML::Key << kk.name() << YAML::Value << kk.evaluate();
                }
                out << YAML::EndMap;
            }
            out << YAML::EndSeq;
            out << YAML::Key << "options" << YAML::Value << YAML::BeginSeq;
            for (const auto & o : options)
            {
                out << YAML::Flow << YAML::BeginMap;
                for (const auto & oo : o)
                {
                    out << YAML::Key << oo.first << YAML::Value << oo.second;
                }
                out << YAML::EndMap;
            }
            out << YAML::EndSeq;
            out << YAML::Key << "means" << YAML::Value << YAML::Flow << YAML::BeginSeq;
            for (const auto & m : means)
            {
                out << m;
            }
            out << YAML::EndSeq;
            out << YAML::Key << "sigma-stat-hi" << YAML::Value << YAML::Flow << YAML::BeginSeq;
            for (const auto & s : sigma_stat_hi)
            {
                out << s;
            }
            out << YAML::EndSeq;
            out << YAML::Key << "sigma-stat-lo" << YAML::Value << YAML::Flow << YAML::BeginSeq;
            for (const auto & s : sigma_stat_lo)
            {
                out << s;
            }
            out << YAML::EndSeq;
            out << YAML::Key << "sigma-sys" << YAML::Value << YAML::Flow << YAML::BeginSeq;
            for (const auto & s : sigma_sys)
            {
                out << s;
            }
            out << YAML::EndSeq;
            out << YAML::Key << "correlations" << YAML::Value << YAML::BeginSeq;
            for (const auto & row : correlation)
            {
                out << YAML::Flow << YAML::BeginSeq;
                for (const auto & corr : row)
                {
                    out << corr;
                }
                out << YAML::EndSeq;
            }
            out << YAML::EndSeq;
            out << YAML::Key << "dof" << YAML::Value << number_of_observations;
            out << YAML::EndMap;
        }

        static ConstraintEntry * deserialize(const QualifiedName & name, const YAML::Node & n)
        {
            static const std::string required_keys[] =
            {
                "dim", "observables", "kinematics", "options", "means", "sigma-stat-hi", "sigma-stat-lo", "sigma-sys", "correlations", "dof"
            };

            for (auto && k : required_keys)
            {
                if (! n[k].IsDefined())
                {
                    throw ConstraintDeserializationError(name, "required key '" + k + "' not specified");
                }
            }

            static const std::string scalar_keys[] =
            {
                "dim", "dof"
            };

            for (auto && k : scalar_keys)
            {
                if (YAML::NodeType::Scalar != n[k].Type())
                {
                    throw ConstraintDeserializationError(name, "required key '" + k + "' not mapped to a scalar value");
                }
            }

            static const std::string seq_keys[] =
            {
                "observables", "kinematics", "options", "means", "sigma-stat-hi", "sigma-stat-lo", "sigma-sys", "correlations"
            };

            for (auto && k : seq_keys)
            {
                if (YAML::NodeType::Sequence != n[k].Type())
                {
                    throw ConstraintDeserializationError(name, "required key '" + k + "' not mapped to a sequence");
                }
            }

            try
            {
                unsigned dof = n["dof"].as<unsigned>();

                std::vector<QualifiedName> observables;
                for (auto && o : n["observables"])
                {
                    observables.push_back(QualifiedName(o.as<std::string>()));
                }

                std::vector<Kinematics> kinematics;
                for (auto && entry : n["kinematics"])
                {
                    if (! n.IsMap())
                    {
                        throw ConstraintDeserializationError(name, "non-map entry encountered in kinematics sequence");
                    }

                    kinematics.push_back(Kinematics{ });
                    std::list<std::pair<YAML::Node, YAML::Node>> kinematics_nodes(entry.begin(), entry.end());
                    // yaml-cpp does not guarantee loading of a map in the order it is written. Circumvent this problem
                    // by sorting the entries lexicographically.
                    kinematics_nodes.sort(&impl::less);
                    std::set<std::string> kinematics_keys;
                    for (auto && k : kinematics_nodes)
                    {
                        std::string key = k.first.as<std::string>();
                        if (! kinematics_keys.insert(key).second)
                            throw ConstraintDeserializationError(name, "kinematics key '" + key + "' encountered more than once");

                        kinematics.back().declare(key, k.second.as<double>());
                    }
                }

                std::vector<Options> options;
                for (auto && entry : n["options"])
                {
                    if (! n.IsMap())
                    {
                        throw ConstraintDeserializationError(name, "non-map entry encountered in options sequence");
                    }

                    options.push_back(Options{ });
                    std::list<std::pair<YAML::Node, YAML::Node>> options_nodes(entry.begin(), entry.end());
                    // yaml-cpp does not guarantee loading of a map in the order it is written. Circumvent this problem
                    // by sorting the entries lexicographically.
                    options_nodes.sort(&impl::less);
                    std::set<std::string> options_keys;
                    for (auto && o : options_nodes)
                    {
                        std::string key = o.first.as<std::string>();
                        if (! options_keys.insert(key).second)
                            throw ConstraintDeserializationError(name, "options key '" + key + "' encountered more than once");

                        options.back().set(key, o.second.as<std::string>());
                    }
                }

                std::vector<double> means;
                for (auto && v : n["means"])
                {
                    means.push_back(v.as<double>());
                }

                std::vector<double> sigma_stat_hi;
                for (auto && v : n["sigma-stat-hi"])
                {
                    sigma_stat_hi.push_back(v.as<double>());
                }

                std::vector<double> sigma_stat_lo;
                for (auto && v : n["sigma-stat-lo"])
                {
                    sigma_stat_lo.push_back(v.as<double>());
                }

                std::vector<double> sigma_sys;
                for (auto && v : n["sigma-sys"])
                {
                    sigma_sys.push_back(v.as<double>());
                }

                std::vector<std::vector<double>> correlations;
                for (auto && row : n["correlations"])
                {
                    correlations.push_back(std::vector<double>());

                    for (auto && v : row)
                    {
                        correlations.back().push_back(v.as<double>());
                    }
                }

                return new MultivariateGaussianConstraintEntry(name.str(), observables, kinematics, options, means,
                        sigma_stat_hi, sigma_stat_lo, sigma_sys, correlations, dof);
            }
            catch (QualifiedNameSyntaxError & e)
            {
                throw ConstraintDeserializationError(name, "'" + n["observable"].as<std::string>() + "' is not a valid observable name (" + e.what() + ")");
            }
        }
    };
    /// }}}

    /// {{{ MultivariateGaussianCovarianceConstraintEntry
    struct MultivariateGaussianCovarianceConstraintEntry :
        public ConstraintEntryBase
    {
        std::vector<QualifiedName> observables;

        std::vector<Kinematics> kinematics;

        std::vector<Options> options;

        std::vector<double> means;

        std::vector<std::vector<double>> covariance;

        unsigned number_of_observations;

        unsigned dim;

        MultivariateGaussianCovarianceConstraintEntry(const QualifiedName & name,
                const std::vector<QualifiedName> & observables,
                const std::vector<Kinematics> & kinematics,
                const std::vector<Options> & options,
                const std::vector<double> & means,
                const std::vector<std::vector<double>> & covariance,
                const unsigned number_of_observations) :
            ConstraintEntryBase(name, observables),
            observables(observables),
            kinematics(kinematics),
            options(options),
            means(means),
            covariance(covariance),
            number_of_observations(number_of_observations),
            dim(observables.size())
        {
            if (dim != kinematics.size()) { throw InternalError("MultivariateGaussianConstraintEntry: wrong number of kinematics"); }

            if (dim != options.size()) { throw InternalError("MultivariateGaussianConstraintEntry: wrong number of options"); }

            if (dim != means.size()) { throw InternalError("MultivariateGaussianConstraintEntry: wrong number of means"); }

            if (dim < number_of_observations) { throw InternalError("MultivariateGaussianConstraintEntry: wrong number of observations"); }

            if (dim != covariance.size()) { throw InternalError("MultviateGaussianConstraintEntry: wrong number of rows in covariance"); }
            for (auto i = 0u ; i < dim ; ++i)
            {
                if (dim != covariance[i].size()) { throw InternalError("MultviateGaussianConstraintEntry: wrong number of columns in covariance row " + stringify(i)); }
            }
        }

        virtual ~MultivariateGaussianCovarianceConstraintEntry() = default;

        virtual const std::string & type() const
        {
            static const std::string type("MultivariateGaussian<" + stringify(dim) + "> (using covariance matrix)");

            return type;
        }

        virtual Constraint make(const QualifiedName & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);

            std::vector<ObservablePtr> observables(dim, nullptr);
            for (auto i = 0u ; i < dim ; ++i)
            {
                observables[i] = Observable::make(this->observables[i], parameters, this->kinematics[i], this->options[i] + options);
                if (! observables[i].get())
                    throw InternalError("make_multivariate_gaussian_covariance_constraint<" + stringify(dim) + ">: " + name.str() + ": '" + this->observables[i].str() + "' is not a valid observable name");
            }

            // create GSL vector for the mean
            gsl_vector * means = gsl_vector_alloc(dim);
            for (auto i = 0u ; i < dim ; ++i)
            {
                gsl_vector_set(means, i, this->means[i]);
            }

            // create GSL matrix for the covariance
            gsl_matrix * covariance = gsl_matrix_alloc(dim, dim);
            for (auto i = 0u ; i < dim ; ++i)
            {
                for (auto j = 0u ; j < dim ; ++j)
                {
                    gsl_matrix_set(covariance, i, j, this->covariance[i][j]);
                }
            }

            auto block = LogLikelihoodBlock::MultivariateGaussian(cache, observables, means, covariance, number_of_observations);

            return Constraint(name, std::vector<ObservablePtr>(observables.begin(), observables.end()), { block });
        }

        virtual std::ostream & insert(std::ostream & os) const
        {
            os << _name.full() << ":" << std::endl;
            os << "    type: MultivariateGaussianCovariance<" << dim << ">" << std::endl;

            return os;
        }

        virtual void serialize(YAML::Emitter & out) const
        {
            out << YAML::BeginMap;
            out << YAML::Key << "type" << YAML::Value << "MultivariateGaussian(Covariance)";
            out << YAML::Key << "dim" << YAML::Value << dim;
            out << YAML::Key << "observables" << YAML::Value << YAML::BeginSeq;
            for (const auto & o : observables)
            {
                out << o.full();
            }
            out << YAML::EndSeq;
            out << YAML::Key << "kinematics" << YAML::Value << YAML::BeginSeq;
            for (const auto & k : kinematics)
            {
                out << YAML::Flow << YAML::BeginMap;
                for (const auto & kk : k)
                {
                    out << YAML::Key << kk.name() << YAML::Value << kk.evaluate();
                }
                out << YAML::EndMap;
            }
            out << YAML::EndSeq;
            out << YAML::Key << "options" << YAML::Value << YAML::BeginSeq;
            for (const auto & o : options)
            {
                out << YAML::Flow << YAML::BeginMap;
                for (const auto & oo : o)
                {
                    out << YAML::Key << oo.first << YAML::Value << oo.second;
                }
                out << YAML::EndMap;
            }
            out << YAML::EndSeq;
            out << YAML::Key << "means" << YAML::Value << YAML::Flow << YAML::BeginSeq;
            for (const auto & m : means)
            {
                out << m;
            }
            out << YAML::EndSeq;
            out << YAML::Key << "covariance" << YAML::Value << YAML::BeginSeq;
            for (const auto & row : covariance)
            {
                out << YAML::Flow << YAML::BeginSeq;
                for (const auto & cov : row)
                {
                    out << cov;
                }
                out << YAML::EndSeq;
            }
            out << YAML::EndSeq;
            out << YAML::Key << "dof" << YAML::Value << number_of_observations;
            out << YAML::EndMap;
        }

        static ConstraintEntry * deserialize(const QualifiedName & name, const YAML::Node & n)
        {
            static const std::string required_keys[] =
            {
                "dim", "observables", "kinematics", "options", "means", "covariance", "dof"
            };

            for (auto && k : required_keys)
            {
                if (! n[k].IsDefined())
                {
                    throw ConstraintDeserializationError(name, "required key '" + k + "' not specified");
                }
            }

            static const std::string scalar_keys[] =
            {
                "dim", "dof"
            };

            for (auto && k : scalar_keys)
            {
                if (YAML::NodeType::Scalar != n[k].Type())
                {
                    throw ConstraintDeserializationError(name, "required key '" + k + "' not mapped to a scalar value");
                }
            }

            static const std::string seq_keys[] =
            {
                "observables", "kinematics", "options", "means", "covariance"
            };

            for (auto && k : seq_keys)
            {
                if (YAML::NodeType::Sequence != n[k].Type())
                {
                    throw ConstraintDeserializationError(name, "required key '" + k + "' not mapped to a sequence");
                }
            }

            try
            {
                unsigned dof = n["dof"].as<unsigned>();

                std::vector<QualifiedName> observables;
                for (auto && o : n["observables"])
                {
                    observables.push_back(QualifiedName(o.as<std::string>()));
                }

                std::vector<Kinematics> kinematics;
                for (auto && entry : n["kinematics"])
                {
                    if (! n.IsMap())
                    {
                        throw ConstraintDeserializationError(name, "non-map entry encountered in kinematics sequence");
                    }

                    kinematics.push_back(Kinematics{ });
                    std::list<std::pair<YAML::Node, YAML::Node>> kinematics_nodes(entry.begin(), entry.end());
                    // yaml-cpp does not guarantee loading of a map in the order it is written. Circumvent this problem
                    // by sorting the entries lexicographically.
                    kinematics_nodes.sort(&impl::less);
                    std::set<std::string> kinematics_keys;
                    for (auto && k : kinematics_nodes)
                    {
                        std::string key = k.first.as<std::string>();
                        if (! kinematics_keys.insert(key).second)
                            throw ConstraintDeserializationError(name, "kinematics key '" + key + "' encountered more than once");

                        kinematics.back().declare(key, k.second.as<double>());
                    }
                }

                std::vector<Options> options;
                for (auto && entry : n["options"])
                {
                    if (! n.IsMap())
                    {
                        throw ConstraintDeserializationError(name, "non-map entry encountered in options sequence");
                    }

                    options.push_back(Options{ });
                    std::list<std::pair<YAML::Node, YAML::Node>> options_nodes(entry.begin(), entry.end());
                    // yaml-cpp does not guarantee loading of a map in the order it is written. Circumvent this problem
                    // by sorting the entries lexicographically.
                    options_nodes.sort(&impl::less);
                    std::set<std::string> options_keys;
                    for (auto && o : options_nodes)
                    {
                        std::string key = o.first.as<std::string>();
                        if (! options_keys.insert(key).second)
                            throw ConstraintDeserializationError(name, "options key '" + key + "' encountered more than once");

                        options.back().set(key, o.second.as<std::string>());
                    }
                }

                std::vector<double> means;
                for (auto && v : n["means"])
                {
                    means.push_back(v.as<double>());
                }

                std::vector<std::vector<double>> covariance;
                for (auto && row : n["covariance"])
                {
                    covariance.push_back(std::vector<double>());

                    for (auto && v : row)
                    {
                        covariance.back().push_back(v.as<double>());
                    }
                }

                return new MultivariateGaussianCovarianceConstraintEntry(name.str(), observables, kinematics, options, means, covariance, dof);
            }
            catch (QualifiedNameSyntaxError & e)
            {
                throw ConstraintDeserializationError(name, "'" + n["observable"].as<std::string>() + "' is not a valid observable name (" + e.what() + ")");
            }
        }
    };
    /// }}}

    ConstraintEntry::~ConstraintEntry() = default;

    ConstraintEntry *
    ConstraintEntry::FromYAML(const QualifiedName & name, const YAML::Node & n)
    {
        // make sure we deserialize from a map
        if (YAML::NodeType::Map != n.Type())
        {
            throw ConstraintDeserializationError(name, "YAML node is not a map");
        }

        if (! n["type"].IsDefined())
        {
            throw ConstraintDeserializationError(name, "YAML node has not key 'type'");
        }

        static const std::map<std::string, std::function<ConstraintEntry * (const QualifiedName &, const YAML::Node &)>> deserializers
        {
            { "Amoroso",                          &AmorosoConstraintEntry::deserialize                        },
            { "Gaussian",                         &GaussianConstraintEntry::deserialize                       },
            { "LogGamma",                         &LogGammaConstraintEntry::deserialize                       },
            { "MultivariateGaussian",             &MultivariateGaussianConstraintEntry::deserialize           },
            { "MultivariateGaussian(Covariance)", &MultivariateGaussianCovarianceConstraintEntry::deserialize }
        };

        std::string type = n["type"].as<std::string>();
        auto i = deserializers.find(type);
        if (i == deserializers.end())
        {
            throw ConstraintDeserializationError(name, "unsupported type '" + type + "'");
        }

        return i->second(name, n);
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

    std::map<QualifiedName, const ConstraintEntry *>
    load_constraint_entries()
    {
        using ValueType = std::map<QualifiedName, const ConstraintEntry *>::value_type;

        std::map<QualifiedName, const ConstraintEntry *> result;

        fs::path base;
        if (std::getenv("EOS_TESTS_CONSTRAINTS"))
        {
            std::string envvar = std::string(std::getenv("EOS_TESTS_CONSTRAINTS"));
            base = fs::system_complete(envvar);
        }
        else if (std::getenv("EOS_HOME"))
        {
            std::string envvar = std::string(std::getenv("EOS_HOME"));
            base = fs::system_complete(envvar) / "constraints";
        }
        else
        {
            base = fs::system_complete(EOS_DATADIR "/eos/constraints/");
        }

        if (! fs::exists(base))
        {
            throw InternalError("Could not find the constraint input files");
        }

        if (! fs::is_directory(base))
        {
            throw InternalError("Expect '" + base.string() + " to be a directory");
        }

        for (fs::directory_iterator f(base), f_end ; f != f_end ; ++f)
        {
            auto file_path = f->path();

            if (! fs::is_regular_file(status(file_path)))
                continue;

            if (".yaml" != file_path.extension().string())
                continue;

            std::string file = file_path.string();
            try
            {
                YAML::Node node = YAML::LoadFile(file);

                for (auto && p : node)
                {
                    std::string keyname = p.first.Scalar();

                    if ("@metadata@" == keyname)
                        continue;

                    QualifiedName name(keyname);
                    ConstraintEntry * entry = ConstraintEntry::FromYAML(name, p.second);

                    if (! result.insert(ValueType{ name, entry }).second)
                    {
                        throw ConstraintInputFileParseError(file, "encountered duplicate constraint '" + keyname + "'");
                    }
                }
            }
            catch (ConstraintDeserializationError & e)
            {
                throw ConstraintInputFileParseError(file, e.what());
            }
        }

        return result;
    }

    static
    const std::map<QualifiedName, const ConstraintEntry *> &
    make_constraint_entries()
    {
        static const std::map<QualifiedName, const ConstraintEntry *> entries = load_constraint_entries();

        return entries;
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
