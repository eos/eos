/* vim: set sw=4 sts=4 et foldmethod=marker foldmarker={{{,}}} : */

/*
 * Copyright (c) 2011-2021 Danny van Dyk
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
#include <eos/maths/gsl-interface.hh>
#include <eos/maths/power-of.hh>
#include <eos/statistics/log-likelihood.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/log.hh>
#include <eos/utils/observable_set.hh>
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

    ConstraintEntryEncodingError::ConstraintEntryEncodingError(const QualifiedName & name) :
        Exception("Constraint '" + name.str() + "' contains non-ascii characters")
    {
    }

    ConstraintInputFileParseError::ConstraintInputFileParseError(const std::string & filename, const std::string & msg) :
        Exception("Could not parse constraint input file '" + filename + "': " + msg)
    {
    }

    template <>
    struct WrappedForwardIteratorTraits<ConstraintEntry::ObservableNameIteratorTag>
    {
        using UnderlyingIterator = std::vector<QualifiedName>::const_iterator;
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

        GaussianConstraintEntry(const std::string & name,
                const QualifiedName & observable,
                const Kinematics & kinematics, const Options & options,
                const double & central,
                const double & sigma_hi_stat, const double & sigma_lo_stat,
                const double & sigma_hi_sys, const double & sigma_lo_sys) :
            ConstraintEntryBase(name, observable),
            observable(observable),
            kinematics(kinematics),
            options(options),
            central(central),
            sigma_hi_stat(sigma_hi_stat),
            sigma_lo_stat(sigma_lo_stat),
            sigma_hi_sys(sigma_hi_sys),
            sigma_lo_sys(sigma_lo_sys)
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

            for (const auto & [key, value] : this->options)
            {
                if (options.has(key) && (value != options[key]))
                {
                    Log::instance()->message("[GaussianConstraintEntry.make]", ll_debug)
                        << "Constraint '" << name << "' provides option key '" << key << "' with value '" << value << "'; user is overriding this preset with '" << options[key] << "'";
                }
            }

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

        virtual LogPriorPtr make_prior(const Parameters & parameters, const Options & options) const
        {
            throw InternalError("GaussianConstraintEntry::make_prior: not yet implemented");
            return nullptr;
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
            out << YAML::DoublePrecision(9);
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
        }

        static ConstraintEntry * deserialize(const QualifiedName & name, const YAML::Node & n)
        {
            static const std::string required_keys[] =
            {
                "observable", "kinematics", "options", "mean", "sigma-stat", "sigma-sys"
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
                "observable", "mean"
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

                    options.declare(key, o.second.as<std::string>());
                }

                double sigma_hi_stat = n["sigma-stat"]["hi"].as<double>();
                double sigma_lo_stat = n["sigma-stat"]["lo"].as<double>();
                double sigma_hi_sys  = n["sigma-sys"]["hi"].as<double>();
                double sigma_lo_sys  = n["sigma-sys"]["lo"].as<double>();

                return new GaussianConstraintEntry(name.str(), observable, kinematics, options, mean,
                        sigma_hi_stat, sigma_lo_stat, sigma_hi_sys, sigma_lo_sys);
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

        double alpha, lambda;

        LogGammaConstraintEntry(const std::string & name,
                const QualifiedName & observable,
                const Kinematics & kinematics, const Options & options,
                const double & central,
                const double & sigma_hi, const double & sigma_lo,
                const double & alpha, const double & lambda) :
            ConstraintEntryBase(name, observable),
            observable(observable),
            kinematics(kinematics),
            options(options),
            central(central),
            sigma_hi(sigma_hi),
            sigma_lo(sigma_lo),
            alpha(alpha),
            lambda(lambda)
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

            for (const auto & [key, value] : this->options)
            {
                if (options.has(key) && (value != options[key]))
                {
                    Log::instance()->message("[LogGammaConstraintEntry.make]", ll_debug)
                        << "Constraint '" << name << "' provides option key '" << key << "' with value '" << value << "'; user is overriding this preset with '" << options[key] << "'";
                }
            }

            ObservablePtr observable = Observable::make(this->observable, parameters, this->kinematics, this->options + options);
            if (! observable.get())
                throw InternalError("make_LogGamma_constraint: " + name.str() + ": '" + this->observable.str() + "' is not a valid observable name");

            double min = this->central - this->sigma_lo;
            double max = this->central + this->sigma_hi;

            LogLikelihoodBlockPtr block = LogLikelihoodBlock::LogGamma(cache, observable, min, this->central, max, alpha, lambda);

            return Constraint(name, { observable }, { block });
        }

        virtual LogPriorPtr make_prior(const Parameters & parameters, const Options & options) const
        {
            throw InternalError("LogGammaConstraintEntry::make_prior: not yet implemented");
            return nullptr;
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
            out << YAML::DoublePrecision(9);
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
            out << YAML::Key << "alpha" << YAML::Value << alpha;
            out << YAML::Key << "lambda" << YAML::Value << lambda;
            out << YAML::EndMap;
        }

        static ConstraintEntry * deserialize(const QualifiedName & name, const YAML::Node & n)
        {
            static const std::string required_keys[] =
            {
                "observable", "kinematics", "options", "mode", "sigma", "alpha", "lambda"
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
                "observable", "mode", "alpha", "lambda"
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

                    options.declare(key, o.second.as<std::string>());
                }

                double sigma_hi = n["sigma"]["hi"].as<double>();
                double sigma_lo = n["sigma"]["lo"].as<double>();

                double alpha  = n["alpha"].as<double>();
                double lambda = n["lambda"].as<double>();

                return new LogGammaConstraintEntry(name.str(), observable, kinematics, options, mode,
                        sigma_hi, sigma_lo, alpha, lambda);
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
            static const std::string type("Amoroso");

            return type;
        }

        virtual Constraint make(const QualifiedName & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);

            for (const auto & [key, value] : this->options)
            {
                if (options.has(key) && (value != options[key]))
                {
                    Log::instance()->message("[AmorosoConstraintEntry.make]", ll_debug)
                        << "Constraint '" << name << "' provides option key '" << key << "' with value '" << value << "'; user is overriding this preset with '" << options[key] << "'";
                }
            }

            ObservablePtr observable = Observable::make(this->observable, parameters, this->kinematics, this->options + options);
            if (! observable.get())
                throw InternalError("make_amoroso_constraint: " + name.str() + ": '" + this->observable.str() + "' is not a valid observable name");

            LogLikelihoodBlockPtr block = LogLikelihoodBlock::Amoroso(cache, observable, physical_limit, theta, alpha, beta);

            return Constraint(name, { observable }, { block });
        }

        virtual LogPriorPtr make_prior(const Parameters & parameters, const Options & options) const
        {
            throw InternalError("AmorosoConstraintEntry::make_prior: not yet implemented");
            return nullptr;
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
            out << YAML::DoublePrecision(9);
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

                    options.declare(key, o.second.as<std::string>());
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

            if (dim != correlation.size()) { throw InternalError("MultivariateGaussianConstraintEntry: wrong number of rows in correlation"); }
            for (auto i = 0u ; i < dim ; ++i)
            {
                if (dim != correlation[i].size()) { throw InternalError("MultivariateGaussianConstraintEntry: wrong number of columns in correlation row " + stringify(i)); }
            }
        }

        virtual ~MultivariateGaussianConstraintEntry() = default;

        virtual const std::string & type() const
        {
            static const std::string type("MultivariateGaussian");

            return type;
        }

        virtual Constraint make(const QualifiedName & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);

            // If specified, these options allow to specify a subset of the measurements
            const unsigned begin = destringify<unsigned>(options.get("begin", "0"));
            const unsigned end = destringify<unsigned>(options.get("end", stringify(dim)));

            if (end > dim)
                throw InvalidOptionValueError("End of the measurements sub-sample: end", options.get("end", stringify(dim)) , "Cannot use a value of 'end' pointing beyond the number of measurements.");

            if (begin >= end)
                throw InvalidOptionValueError("First measurement of the sub-sample: begin", options.get("begin", "0"), "Cannot use a value for 'begin' equal to or larger than 'end'");

            const unsigned subdim_meas = end - begin;

            std::vector<ObservablePtr> observables(subdim_meas, nullptr);
            for (auto i = begin ; i < end ; ++i)
            {
                for (const auto & [key, value] : this->options[i])
                {
                    if (options.has(key) && (value != options[key]))
                    {
                        Log::instance()->message("[MultivariateGaussianConstraintEntry.make]", ll_debug)
                            << "Constraint '" << name << "' in observable '" << this->observable_names[i] << "' provides option key '" << key << "' with value '" << value << "'; user is overriding this preset with '" << options[key] << "'";
                    }
                }

                observables[i - begin] = Observable::make(this->observable_names[i], parameters, this->kinematics[i], this->options[i] + options);
                if (! observables[i - begin].get())
                    throw InternalError("make_multivariate_gaussian_constraint<" + stringify(dim) + ">: " + name.str() + ": '" + this->observable_names[i].str() + "' is not a valid observable name");
            }

            std::vector<double> variances(subdim_meas, 0.0);
            if ("symmetric+quadratic" == options.get("uncertainty", "symmetric+quadratic"))
            {
                for (auto i = begin ; i < end ; ++i)
                {
                    double combined_lo = power_of<2>(sigma_stat_lo[i]) + power_of<2>(sigma_sys[i]);
                    double combined_hi = power_of<2>(sigma_stat_hi[i]) + power_of<2>(sigma_sys[i]);

                    variances[i - begin] = std::max(combined_lo, combined_hi);
                }
            }

            // create GSL vector for the mean
            gsl_vector * means = gsl_vector_alloc(subdim_meas);
            for (auto i = begin ; i < end ; ++i)
            {
                gsl_vector_set(means, i - begin, this->means[i]);
            }

            // create GSL matrix for the covariance
            gsl_matrix * covariance = gsl_matrix_alloc(subdim_meas, subdim_meas);
            for (auto i = begin ; i < end ; ++i)
            {
                for (auto j = begin ; j < end ; ++j)
                {
                    double value = std::sqrt(variances[i - begin] * variances[j - begin]) * correlation[i][j];
                    gsl_matrix_set(covariance, i - begin, j - begin, value);
                }
            }

            gsl_matrix * response = gsl_matrix_calloc(subdim_meas, subdim_meas);
            gsl_matrix_set_identity(response);

            if (subdim_meas != response->size2)
                throw InternalError("Constraint " + name.full() + ": number of predictions and number of columns in response matrix are not identical");

            // create the block for number_of_observations = subdim_meas
            //   - if we sliced the measurements, we need to adjust the number of observations
            //   - if we did not slice, we made sure that number_of_observations == dim == subdim_meas
            LogLikelihoodBlockPtr block = LogLikelihoodBlock::MultivariateGaussian(cache, observables, means, covariance, response, subdim_meas);

            return Constraint(name, std::vector<ObservablePtr>(observables.begin(), observables.end()), { block });
        }

        virtual LogPriorPtr make_prior(const Parameters & parameters, const Options & options) const
        {
            // If specified, these options allow to specify a subset of the measurements
            unsigned begin = destringify<unsigned>(options.get("begin", "0"));
            unsigned end = destringify<unsigned>(options.get("end", stringify(dim)));

            if (end > dim)
                throw InvalidOptionValueError("End of the measurements sub-sample: end", options.get("end", stringify(dim)) , "Cannot use a value of 'end' pointing beyond the number of measurements.");

            if (begin >= end)
                throw InvalidOptionValueError("First measurement of the sub-sample: begin", options.get("begin", "0"), "Cannot use a value for 'begin' equal to or larger than 'end'");

            unsigned subdim = end - begin;

            // create GSL vector for the mean
            gsl_vector * means = gsl_vector_alloc(subdim);
            for (auto i = begin ; i < end ; ++i)
            {
                gsl_vector_set(means, i - begin, this->means[i]);
            }

            std::vector<double> variances(subdim, 0.0);
            if ("symmetric+quadratic" == options.get("uncertainty", "symmetric+quadratic"))
            {
                for (auto i = begin ; i < end ; ++i)
                {
                    double combined_lo = power_of<2>(sigma_stat_lo[i]) + power_of<2>(sigma_sys[i]);
                    double combined_hi = power_of<2>(sigma_stat_hi[i]) + power_of<2>(sigma_sys[i]);

                    variances[i - begin] = std::max(combined_lo, combined_hi);
                }
            }

            // create GSL matrix for the covariance
            gsl_matrix * covariance = gsl_matrix_alloc(subdim, subdim);
            for (auto i = begin ; i < end ; ++i)
            {
                for (auto j = begin ; j < end ; ++j)
                {
                    double value = std::sqrt(variances[i - begin] * variances[j - begin]) * correlation[i][j];
                    gsl_matrix_set(covariance, i - begin, j - begin, value);
                }
            }

            return LogPrior::MultivariateGaussian(parameters, observable_names, means, covariance);
        }

        virtual std::ostream & insert(std::ostream & os) const
        {
            os << _name.full() << ":" << std::endl;
            os << "    type: MultivariateGaussian<" << dim << ">" << std::endl;

            return os;
        }

        virtual void serialize(YAML::Emitter & out) const
        {
            out << YAML::DoublePrecision(9);
            out << YAML::BeginMap;
            out << YAML::Key << "type" << YAML::Value << "MultivariateGaussian";
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
                "observables", "kinematics", "options", "means", "sigma-stat-hi", "sigma-stat-lo", "sigma-sys", "correlations"
            };

            for (auto && k : required_keys)
            {
                if (! n[k].IsDefined())
                {
                    throw ConstraintDeserializationError(name, "required key '" + k + "' not specified");
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

                        options.back().declare(key, o.second.as<std::string>());
                    }
                }

                std::vector<double> means;
                for (auto && v : n["means"])
                {
                    means.push_back(v.as<double>());
                }

                // Test for the presence of the optional "dof" parameter
                unsigned dof;
                if (n["dof"])
                {
                    if (YAML::NodeType::Scalar != n["dof"].Type())
                    {
                        throw ConstraintDeserializationError(name, "optonal key 'dof' not mapped to a scalar value");
                    }
                    dof = n["dof"].as<unsigned>();
                }
                else
                    dof = means.size();

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

        const GSLVectorPtr means;

        const GSLMatrixPtr covariance;

        const GSLMatrixPtr response;

        unsigned number_of_observations;

        unsigned dim_meas, dim_pred;

        MultivariateGaussianCovarianceConstraintEntry(const QualifiedName & name,
                const std::vector<QualifiedName> & observables,
                const std::vector<Kinematics> & kinematics,
                const std::vector<Options> & options,
                gsl_vector * const means,
                gsl_matrix * const covariance,
                gsl_matrix * const response,
                const unsigned number_of_observations) :
            ConstraintEntryBase(name, observables),
            observables(observables),
            kinematics(kinematics),
            options(options),
            means(means),
            covariance(covariance),
            response(response),
            number_of_observations(number_of_observations),
            dim_meas(means->size),
            dim_pred(observables.size())
        {
            if ((nullptr == response) && (dim_meas != dim_pred)) { throw InternalError("MultivariateGaussianConstraintEntry: number of measurements does not equal number of predictions in absence of a response matrix"); }

            if (dim_meas != covariance->size1) { throw InternalError("MultivariateGaussianConstraintEntry: number of rows in covariance does not equal number of measurements"); }
            if (dim_meas != covariance->size2) { throw InternalError("MultivariateGaussianConstraintEntry: number of columns in covariance does not equal number of measurements"); }

            if ((nullptr != response) && (dim_meas != response->size1)) { throw InternalError("MultivariateGaussianConstraintEntry: number of rows in response does not equal number of measurements");}
            if ((nullptr != response) && (dim_pred != response->size2)) { throw InternalError("MultivariateGaussianConstraintEntry: number of columns in response does not equal number of predictions");}

            if (dim_pred != kinematics.size()) { throw InternalError("MultivariateGaussianConstraintEntry: number of kinematics entries does not equal number of predictions"); }

            if (dim_pred != options.size()) { throw InternalError("MultivariateGaussianConstraintEntry: number of options entries does not equal number of predictions"); }

            if (dim_meas < number_of_observations) { throw InternalError("MultivariateGaussianConstraintEntry: number of observations larger than number of measurements"); }
        }

        virtual ~MultivariateGaussianCovarianceConstraintEntry()
        {
        }

        virtual const std::string & type() const
        {
            static const std::string type("MultivariateGaussian(Covariance)");

            return type;
        }

        virtual Constraint make(const QualifiedName & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);

            std::vector<ObservablePtr> observables(dim_pred, nullptr);
            for (auto i = 0u ; i < dim_pred ; ++i)
            {
                for (const auto & [key, value] : this->options[i])
                {
                    if (options.has(key) && (value != options[key]))
                    {
                        Log::instance()->message("[MultivariateGaussianCovarianceConstraintEntry.make]", ll_debug)
                            << "Constraint '" << name << "' in observable '" << this->observables[i] << "' provides option key '" << key << "' with value '" << value << "'; user is overriding this preset with '" << options[key] << "'";
                    }
                }

                observables[i] = Observable::make(this->observables[i], parameters, this->kinematics[i], this->options[i] + options);
                if (! observables[i].get())
                    throw InternalError("make_multivariate_gaussian_covariance_constraint<measurements=" + stringify(dim_meas) + ",predictions=" + stringify(dim_pred) + ">: " + name.str() + ": '" + this->observables[i].str() + "' is not a valid observable name");
            }

            // If specified, these options allow to specify a subset of the measurements
            unsigned begin = destringify<unsigned>(options.get("begin", "0"));
            unsigned end = destringify<unsigned>(options.get("end", stringify(dim_meas)));

            if (end > dim_meas)
                throw InvalidOptionValueError("End of the measurements sub-sample: end", options.get("end", stringify(dim_meas)), "Cannot use a value of 'end' pointing beyond the number of measurements.");

            if (begin >= end)
                throw InvalidOptionValueError("First measurement of the sub-sample: begin", options.get("begin", "0"), "Cannot use a value for 'begin' equal to or larger than 'end'");

            if ((nullptr != this->response) && (end != dim_meas))
                throw InternalError("Response matrices and begin and end options are mutually incompatible.");

            if ((nullptr != this->response) && (begin != 0))
                throw InternalError("Response matrices and begin and end options are mutually incompatible.");

            unsigned subdim_meas = end - begin;

            // create GSL vector for the mean
            gsl_vector * means = gsl_vector_alloc(subdim_meas);
            gsl_vector_view means_subset = gsl_vector_subvector(this->means.get(), begin, subdim_meas);
            gsl_vector_memcpy(means, &(means_subset.vector));

            // create GSL matrix for the covariance
            gsl_matrix * covariance = gsl_matrix_alloc(subdim_meas, subdim_meas);
            gsl_matrix_view covariance_subset = gsl_matrix_submatrix(this->covariance.get(), begin, begin, subdim_meas, subdim_meas);
            gsl_matrix_memcpy(covariance, &(covariance_subset.matrix));

            if (this->response)
            {
                // create GSL matrix for the response
                gsl_matrix * response = gsl_matrix_calloc(subdim_meas, dim_pred);
                gsl_matrix_view covariance_subset = gsl_matrix_submatrix(this->response.get(), begin, 0, subdim_meas, dim_pred);
                gsl_matrix_memcpy(response, &(covariance_subset.matrix));

                // create the block for number_of_observations = subdim_meas
                //   - if we sliced the measurements, we need to adjust the number of observations
                //   - if we did not slice, we made sure that number_of_observations == dim == subdim_meas
                auto block = LogLikelihoodBlock::MultivariateGaussian(cache, observables, means, covariance, response, subdim_meas);

                return Constraint(name, std::vector<ObservablePtr>(observables.begin(), observables.end()), { block });
            }
            else
            {
                // extract the observables
                std::vector<ObservablePtr>::const_iterator first_obs = observables.begin() + begin;
                std::vector<ObservablePtr>::const_iterator last_obs  = observables.begin() + end;
                std::vector<ObservablePtr> restricted_observables(first_obs, last_obs);

                // create GSL matrix for the response
                gsl_matrix * response = gsl_matrix_calloc(subdim_meas, subdim_meas);
                gsl_matrix_set_identity(response);

                // create the block for number_of_observations = subdim_meas
                //   - if we sliced the measurements, we need to adjust the number of observations
                //   - if we did not slice, we made sure that number_of_observations == dim == subdim_meas
                auto block = LogLikelihoodBlock::MultivariateGaussian(cache, restricted_observables, means, covariance, response, subdim_meas);

                return Constraint(name, std::vector<ObservablePtr>(restricted_observables.begin(), restricted_observables.end()), { block });
            }
        }

        virtual LogPriorPtr make_prior(const Parameters & parameters, const Options & options) const
        {
            // If specified, these options allow to specify a subset of the measurements
            unsigned begin = destringify<unsigned>(options.get("begin", "0"));
            unsigned end = destringify<unsigned>(options.get("end", stringify(dim_meas)));

            if (end > dim_meas)
                throw InvalidOptionValueError("End of the measurements sub-sample: end", options.get("end", stringify(dim_meas)) , "Cannot use a value of 'end' pointing beyond the number of measurements.");

            if (begin >= end)
                throw InvalidOptionValueError("First measurement of the sub-sample: begin", options.get("begin", "0"), "Cannot use a value for 'begin' equal to or larger than 'end'");

            unsigned subdim_meas = end - begin;

            // create GSL vector for the mean
            gsl_vector * means = gsl_vector_alloc(subdim_meas);
            gsl_vector_view means_subset = gsl_vector_subvector(this->means.get(), begin, subdim_meas);
            gsl_vector_memcpy(means, &(means_subset.vector));

            // create GSL matrix for the covariance
            gsl_matrix * covariance = gsl_matrix_alloc(subdim_meas, subdim_meas);
            gsl_matrix_view covariance_subset = gsl_matrix_submatrix(this->covariance.get(), begin, begin, subdim_meas, subdim_meas);
            gsl_matrix_memcpy(covariance, &(covariance_subset.matrix));

            return LogPrior::MultivariateGaussian(parameters, observables, means, covariance);
        }

        virtual std::ostream & insert(std::ostream & os) const
        {
            os << _name.full() << ":" << std::endl;
            os << "    type: MultivariateGaussianCovariance<measurements=" << dim_meas << ",predictions=" << dim_pred << ">" << std::endl;

            return os;
        }

        virtual void serialize(YAML::Emitter & out) const
        {
            out << YAML::DoublePrecision(9);
            out << YAML::BeginMap;
            out << YAML::Key << "type" << YAML::Value << "MultivariateGaussian(Covariance)";
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
            for (auto i = 0u ; i < dim_meas ; ++i)
            {
                out << gsl_vector_get(means.get(), i);
            }
            out << YAML::EndSeq;
            out << YAML::Key << "covariance" << YAML::Value << YAML::BeginSeq;
            for (auto i = 0u ; i < dim_meas ; ++i)
            {
                out << YAML::Flow << YAML::BeginSeq;
                for (auto j = 0u ; j < dim_meas ; ++j)
                {
                    out << gsl_matrix_get(covariance.get(), i, j);
                }
                out << YAML::EndSeq;
            }
            out << YAML::EndSeq;
            if (response)
            {
                out << YAML::Key << "response" << YAML::Value << YAML::BeginSeq;
                for (auto i = 0u ; i < dim_meas ; ++i)
                {
                    out << YAML::Flow << YAML::BeginSeq;
                    for (auto j = 0u ; j < dim_pred ; ++j)
                    {
                        out << gsl_matrix_get(response.get(), i, j);
                    }
                    out << YAML::EndSeq;
                }
                out << YAML::EndSeq;
            }
            out << YAML::Key << "dof" << YAML::Value << number_of_observations;
            out << YAML::EndMap;
        }

        static ConstraintEntry * deserialize(const QualifiedName & name, const YAML::Node & n)
        {
            static const std::string required_keys[] =
            {
                "observables", "kinematics", "options", "means", "covariance"
            };

            for (auto && k : required_keys)
            {
                if (! n[k].IsDefined())
                {
                    throw ConstraintDeserializationError(name, "required key '" + k + "' not specified");
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

                        options.back().declare(key, o.second.as<std::string>());
                    }
                }

                std::vector<double> _means;
                for (auto && v : n["means"])
                {
                    _means.push_back(v.as<double>());
                }
                gsl_vector * means = gsl_vector_alloc(_means.size());
                for (auto i = 0u ; i < _means.size() ; ++i)
                {
                    gsl_vector_set(means, i, _means[i]);
                }

                // Test for the presence of the optional "dof" parameter
                unsigned dof;
                if (n["dof"])
                {
                    if (YAML::NodeType::Scalar != n["dof"].Type())
                    {
                        throw ConstraintDeserializationError(name, "optonal key 'dof' not mapped to a scalar value");
                    }
                    dof = n["dof"].as<unsigned>();
                }
                else
                    dof = _means.size();

                std::vector<std::vector<double>> _covariance;
                for (auto && row : n["covariance"])
                {
                    _covariance.push_back(std::vector<double>());

                    for (auto && v : row)
                    {
                        _covariance.back().push_back(v.as<double>());
                    }
                }
                gsl_matrix * covariance = gsl_matrix_alloc(_covariance.size(), _covariance.size());
                for (auto i = 0u ; i < _covariance.size() ; ++i)
                {
                    const auto & _row = _covariance[i];

                    if (_row.size() != _covariance.size())
                    {
                        throw ConstraintDeserializationError(name, "covariance matrix is not square; row " + stringify(i) + " has " + stringify(_row.size()) + " columns; expected " + stringify(_covariance.size()));
                    }

                    for (auto j = 0u ; j < _row.size() ; ++j)
                    {
                        gsl_matrix_set(covariance, i, j, _row[j]);
                    }
                }

                gsl_matrix * response = nullptr;
                if (n["response"].IsDefined())
                {
                    std::vector<std::vector<double>> _response;
                    for (auto && row : n["response"])
                    {
                        _response.push_back(std::vector<double>());

                        for (auto && v : row)
                        {
                            _response.back().push_back(v.as<double>());
                        }
                    }

                    unsigned _response_rows = _response.size();
                    unsigned _response_cols = _response[0].size();

                    response = gsl_matrix_alloc(_response_rows, _response_cols);
                    for (auto i = 0u ; i < _response_rows ; ++i)
                    {
                        const auto & _row = _response[i];

                        if (_row.size() != _response_cols)
                        {
                            throw ConstraintDeserializationError(name, "response matrix is invalid; row " + stringify(i) + " has " + stringify(_row.size()) + " columns; expected " + stringify(_response_cols));
                        }

                        for (auto j = 0u ; j < _row.size() ; ++j)
                        {
                            gsl_matrix_set(response, i, j, _row[j]);
                        }
                    }
                }

                return new MultivariateGaussianCovarianceConstraintEntry(name.str(), observables, kinematics, options, means, covariance, response, dof);
            }
            catch (QualifiedNameSyntaxError & e)
            {
                throw ConstraintDeserializationError(name, "'" + n["observable"].as<std::string>() + "' is not a valid observable name (" + e.what() + ")");
            }
        }
    };
    /// }}}

    /// {{{ UniformBoundConstraintEntry
    struct UniformBoundConstraintEntry :
        public ConstraintEntryBase
    {
        std::vector<QualifiedName> observable_names;

        std::vector<Kinematics> kinematics;

        std::vector<Options> options;

        const double bound;
        const double uncertainty;

        unsigned number_of_observables;

        UniformBoundConstraintEntry(const std::string & name,
                const std::vector<QualifiedName> & observable_names,
                const std::vector<Kinematics> & kinematics,
                const std::vector<Options> & options,
                const double & bound,
                const double & uncertainty) :
            ConstraintEntryBase(name, observable_names),
            observable_names(observable_names),
            kinematics(kinematics),
            options(options),
            bound(bound),
            uncertainty(uncertainty),
            number_of_observables(observable_names.size())
        {
            if (number_of_observables != kinematics.size()) { throw InternalError("UniformBoundConstraintEntry: wrong number of kinematics"); }

            if (number_of_observables != options.size()) { throw InternalError("UniformBoundConstraintEntry: wrong number of options"); }
        }

        virtual ~UniformBoundConstraintEntry() = default;

        virtual const std::string & type() const
        {
            static const std::string type("UniformBound");

            return type;
        }

        virtual Constraint make(const QualifiedName & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);


            std::vector<ObservablePtr> observables(number_of_observables, nullptr);
            for (auto i = 0u ; i < number_of_observables ; ++i)
            {
                for (const auto & [key, value] : this->options[i])
                {
                    if (options.has(key) && (value != options[key]))
                    {
                        Log::instance()->message("[UniformBoundConstraintEntry.make]", ll_debug)
                            << "Constraint '" << name << "' in observable '" << this->observable_names[i] << "' provides option key '" << key << "' with value '" << value << "'; user is overriding this preset with '" << options[key] << "'";
                    }
                }

                observables[i] = Observable::make(this->observable_names[i], parameters, this->kinematics[i], this->options[i] + options);
                if (! observables[i].get())
                    throw InternalError("make_uniform_bound_constraint<" + stringify(number_of_observables) + ">: " + name.str() + ": '" + this->observable_names[i].str() + "' is not a valid observable name");
            }

            LogLikelihoodBlockPtr block = LogLikelihoodBlock::UniformBound(cache, observables, this->bound, this->uncertainty);

            return Constraint(name, std::vector<ObservablePtr>(observables.begin(), observables.end()), { block });
        }

        virtual LogPriorPtr make_prior(const Parameters & parameters, const Options & options) const
        {
            throw InternalError("UniformBoundConstraintEntry::make_prior: not yet implemented");
            return nullptr;
        }

        virtual std::ostream & insert(std::ostream & os) const
        {
            os << _name.full() << ":" << std::endl;
            os << "    type: UniformBound<" << number_of_observables << ">" << std::endl;

            return os;
        }

        virtual void serialize(YAML::Emitter & out) const
        {
            out << YAML::DoublePrecision(9);
            out << YAML::BeginMap;
            out << YAML::Key << "type" << YAML::Value << "UniformBound";
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
            out << YAML::Key << "bound" << YAML::Value << bound;
            out << YAML::Key << "uncertainty" << YAML::Value << uncertainty;
            out << YAML::EndMap;
        }

        static ConstraintEntry * deserialize(const QualifiedName & name, const YAML::Node & n)
        {
            static const std::string required_keys[] =
            {
                "observables", "kinematics", "options", "bound", "uncertainty"
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
                "bound", "uncertainty"
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
                "observables", "kinematics", "options"
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

                        options.back().declare(key, o.second.as<std::string>());
                    }
                }

                double bound       = n["bound"].as<double>();
                double uncertainty = n["uncertainty"].as<double>();

                return new UniformBoundConstraintEntry(name.str(), observables, kinematics, options, bound, uncertainty);
            }
            catch (QualifiedNameSyntaxError & e)
            {
                throw ConstraintDeserializationError(name, "'" + n["observable"].as<std::string>() + "' is not a valid observable name (" + e.what() + ")");
            }
        }
    };
    /// }}}

    /// {{{ MixtureConstraintEntry
    struct MixtureConstraintEntry :
        public ConstraintEntryBase
    {
        std::vector<QualifiedName> observables;

        std::vector<Kinematics> kinematics;

        std::vector<Options> options;

        std::vector<GSLVectorPtr> means;

        std::vector<GSLMatrixPtr> covariances;

        std::vector<double> weights;

        std::vector<std::array<double, 2>> test_stat;

        unsigned number_of_observations;

        unsigned dim_meas, dim_pred;

        explicit
        MixtureConstraintEntry(const QualifiedName & name,
                const std::vector<QualifiedName> & observables,
                const std::vector<Kinematics> & kinematics,
                const std::vector<Options> & options,
                std::vector<GSLVectorPtr> && _means,
                std::vector<GSLMatrixPtr> && _covariances,
                const std::vector<double> & weights,
                const std::vector<std::array<double, 2>> & test_stat,
                const unsigned number_of_observations) :
            ConstraintEntryBase(name, observables),
            observables(observables),
            kinematics(kinematics),
            options(options),
            means(std::move(_means)),
            covariances(std::move(_covariances)),
            weights(weights),
            test_stat(test_stat),
            number_of_observations(number_of_observations),
            dim_pred(observables.size())
        {
            if (means.size() != covariances.size()) { throw InternalError("MixtureConstraintEntry: number of components does not agree between means and covariances"); }
            if (means.size() != weights.size()) { throw InternalError("MixtureConstraintEntry: number of components does not agree between means and weights"); }

            if (means.empty()) { throw InternalError("MixtureConstraintEntry: need at least one component"); }

            dim_meas = means.front()->size;
            for (unsigned int i = 0 ; i < means.size() ; ++i)
            {
                if (dim_meas != means[i]->size) { throw InternalError("MixtureConstraintEntry: mean vectors are not all equal in size"); }
            }

            if (dim_meas != dim_pred) { throw InternalError("MixtureConstraintEntry: number of measurements does not equal number of predictions"); }

            for (unsigned int i = 0 ; i < means.size() ; ++i)
            {
                if (dim_meas != covariances[i]->size1) { throw InternalError("MixtureConstraintEntry: number of rows in at least one covariance does not equal number of measurements"); }
                if (dim_meas != covariances[i]->size2) { throw InternalError("MixtureConstraintEntry: number of columns in at least one covariance does not equal number of measurements"); }
            }

            if (dim_pred != kinematics.size()) { throw InternalError("MixtureConstraintEntry: number of kinematics entries does not equal number of predictions"); }

            if (dim_pred != options.size()) { throw InternalError("MixtureConstraintEntry: number of options entries does not equal number of predictions"); }

            if (dim_meas < number_of_observations) { throw InternalError("MixtureConstraintEntry: number of observations larger than number of measurements"); }
        }

        MixtureConstraintEntry(const MixtureConstraintEntry &) = delete;
        MixtureConstraintEntry & operator= (const MixtureConstraintEntry &) = delete;

        virtual ~MixtureConstraintEntry()
        {
        }

        virtual const std::string & type() const
        {
            static const std::string type("Mixture");

            return type;
        }

        virtual Constraint make(const QualifiedName & name, const Options & options) const
        {
            Parameters parameters(Parameters::Defaults());
            ObservableCache cache(parameters);

            std::vector<ObservablePtr> observables(dim_pred, nullptr);
            for (auto i = 0u ; i < dim_pred ; ++i)
            {
                for (const auto & [key, value] : this->options[i])
                {
                    if (options.has(key) && (value != options[key]))
                    {
                        Log::instance()->message("[MixtureConstraintEntry.make]", ll_debug)
                            << "Constraint '" << name << "' in observable '" << this->observables[i] << "' provides option key '" << key << "' with value '" << value << "'; user is overriding this preset with '" << options[key] << "'";
                    }
                }

                observables[i] = Observable::make(this->observables[i], parameters, this->kinematics[i], this->options[i] + options);
                if (! observables[i].get())
                    throw InternalError("make_multivariate_gaussian_covariance_constraint<measurements=" + stringify(dim_meas) + ",predictions=" + stringify(dim_pred) + ">: " + name.str() + ": '" + this->observables[i].str() + "' is not a valid observable name");
            }

            std::vector<LogLikelihoodBlockPtr> components;

            for (unsigned i = 0 ; i < this->means.size() ; ++i)
            {
                // create GSL vector for the mean of each component
                gsl_vector * mean = gsl_vector_alloc(dim_meas);
                gsl_vector_memcpy(mean, this->means[i].get());

                // create GSL matrix for the covariance of each component
                gsl_matrix * covariance = gsl_matrix_alloc(dim_meas, dim_meas);
                gsl_matrix_memcpy(covariance, this->covariances[i].get());

                gsl_matrix * response = gsl_matrix_alloc(dim_meas, dim_pred);
                gsl_matrix_set_identity(response);

                components.push_back(LogLikelihoodBlock::MultivariateGaussian(cache, observables, mean, covariance, response, dim_meas));
            }

            auto block = LogLikelihoodBlock::Mixture(components, weights, test_stat);

            return Constraint(name, std::vector<ObservablePtr>(observables.begin(), observables.end()), { block });
        }

        virtual LogPriorPtr make_prior(const Parameters & parameters, const Options & options) const
        {
            throw InternalError("MixtureConstraintEntry::make_prior: not yet implemented");
            return nullptr;
        }

        virtual std::ostream & insert(std::ostream & os) const
        {
            os << _name.full() << ":" << std::endl;
            os << "    type: Mixture<components=" << means.size() << ",measurements=" << dim_meas << ",predictions=" << dim_pred << ">" << std::endl;

            return os;
        }

        virtual void serialize(YAML::Emitter & out) const
        {
            out << YAML::DoublePrecision(9);
            out << YAML::BeginMap;
            out << YAML::Key << "type" << YAML::Value << "Mixture";
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
            out << YAML::Key << "components" << YAML::Value << YAML::BeginSeq;
            for (unsigned n = 0 ; n < means.size() ; ++n)
            {
                out << YAML::BeginMap;
                out << YAML::Key << "means" << YAML::Value << YAML::Flow << YAML::BeginSeq;
                for (auto i = 0u ; i < dim_meas ; ++i)
                {
                    out << gsl_vector_get(means[n].get(), i);
                }
                out << YAML::EndSeq;
                out << YAML::Key << "covariance" << YAML::Value << YAML::BeginSeq;
                for (auto i = 0u ; i < dim_meas ; ++i)
                {
                    out << YAML::Flow << YAML::BeginSeq;
                    for (auto j = 0u ; j < dim_meas ; ++j)
                    {
                        out << gsl_matrix_get(covariances[n].get(), i, j);
                    }
                    out << YAML::EndSeq;
                }
                out << YAML::EndSeq;
                out << YAML::EndMap;
            }
            out << YAML::EndSeq;
            out << YAML::Key << "weights" << YAML::Value << YAML::Flow << YAML::BeginSeq;
            for (const auto & w : weights)
            {
                out << w;
            }
            out << YAML::EndSeq;
            out << YAML::Key << "test statistics" << YAML::Value << YAML::BeginMap;
            out << YAML::Key << "sigma" << YAML::Value << YAML::Flow << YAML::BeginSeq;
            if (test_stat.size() != 0)
            {
                for (auto i = 0u ; i < test_stat[0].size() ; ++i)
                {
                    out << test_stat[0][i];
                }
            }
            out << YAML::EndSeq;
            out << YAML::Key << "densities" << YAML::Value << YAML::Flow << YAML::BeginSeq;
            if (test_stat.size() != 0)
            {
                for (auto i = 0u ; i < test_stat[0].size() ; ++i)
                {
                    out << test_stat[1][i];
                }
            }
            out << YAML::EndSeq;
            out << YAML::EndMap;
            out << YAML::Key << "dof" << YAML::Value << number_of_observations;
            out << YAML::EndMap;
        }

        static ConstraintEntry * deserialize(const QualifiedName & name, const YAML::Node & n)
        {
            static const std::string required_keys[] =
            {
                "observables", "kinematics", "options", "components", "weights", "test statistics"
            };

            for (auto && k : required_keys)
            {
                if (! n[k].IsDefined())
                {
                    throw ConstraintDeserializationError(name, "required key '" + k + "' not specified");
                }
            }

            static const std::string seq_keys[] =
            {
                "observables", "kinematics", "options", "components", "weights"
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

                        options.back().declare(key, o.second.as<std::string>());
                    }
                }

                std::vector<GSLVectorPtr> means;
                std::vector<GSLMatrixPtr> covariances;
                for (auto && c : n["components"])
                {
                    if (YAML::NodeType::Sequence != c["means"].Type())
                    {
                        throw ConstraintDeserializationError(name, "required key 'means' not mapped to a sequence");
                    }

                    std::vector<double> _mean;
                    for (auto && v : c["means"])
                    {
                        _mean.push_back(v.as<double>());
                    }
                    GSLVectorPtr mean = make_gsl_vector(_mean.size());
                    for (auto i = 0u ; i < _mean.size() ; ++i)
                    {
                        gsl_vector_set(mean.get(), i, _mean[i]);
                    }
                    means.emplace_back(std::move(mean));

                    if (YAML::NodeType::Sequence != c["covariance"].Type())
                    {
                        throw ConstraintDeserializationError(name, "required key 'covariance' not mapped to a sequence");
                    }

                    std::vector<std::vector<double>> _covariance;
                    for (auto && row : c["covariance"])
                    {
                        _covariance.push_back(std::vector<double>());

                        for (auto && v : row)
                        {
                            _covariance.back().push_back(v.as<double>());
                        }
                    }
                    GSLMatrixPtr covariance = make_gsl_matrix(_covariance.size(), _covariance.size());
                    for (auto i = 0u ; i < _covariance.size() ; ++i)
                    {
                        const auto & _row = _covariance[i];

                        if (_row.size() != _covariance.size())
                        {
                            throw ConstraintDeserializationError(name, "covariance matrix is not square; row " + stringify(i) + " has " + stringify(_row.size()) + " columns; expected " + stringify(_covariance.size()));
                        }

                        for (auto j = 0u ; j < _row.size() ; ++j)
                        {
                            gsl_matrix_set(covariance.get(), i, j, _row[j]);
                        }
                    }
                    covariances.emplace_back(std::move(covariance));
                }

                // Test for the presence of the optional "dof" parameter or infer it from the first component
                unsigned dof;
                if (n["dof"])
                {
                    if (YAML::NodeType::Scalar != n["dof"].Type())
                    {
                        throw ConstraintDeserializationError(name, "optonal key 'dof' not mapped to a scalar value");
                    }
                    dof = n["dof"].as<unsigned>();
                }
                else
                    dof = means[0]->size;

                std::vector<double> weights;
                for (auto && v : n["weights"])
                {
                    weights.push_back(v.as<double>());
                }

                std::vector<double> sigma;
                for (auto && v : n["test statistics"]["sigma"])
                {
                    sigma.push_back(v.as<double>());
                }

                std::vector<double> densities;
                for (auto && v : n["test statistics"]["densities"])
                {
                    densities.push_back(v.as<double>());
                }

                if (sigma.size() != densities.size())
                {
                    throw ConstraintDeserializationError(name, "'sigma' and 'densities' have different size in 'test statistics'");
                }

                std::vector<std::array<double, 2>> test_stat;
                auto sigma_it = sigma.begin();
                for (auto value = densities.cbegin(); value != densities.cend() ; ++value, ++sigma_it)
                {
                    test_stat.push_back(std::array<double, 2>{*sigma_it, *value});
                }

                return new MixtureConstraintEntry(name.str(), observables, kinematics, options, std::move(means), std::move(covariances), std::move(weights), std::move(test_stat), dof);
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
    ConstraintEntry::FromYAML(const QualifiedName & name, const std::string & s)
    {
        // valid ASCII characters are limited to 0 <= c <= 0x7f (127)
        if (s.cend() != std::find_if(s.cbegin(), s.cend(), [&] (unsigned char c) -> bool { return (c > 0x7f); }))
            throw ConstraintEntryEncodingError(name);

        YAML::Node node = YAML::Load(s);

        return ConstraintEntry::FromYAML(name, node);
    }

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
            { "MultivariateGaussian(Covariance)", &MultivariateGaussianCovarianceConstraintEntry::deserialize },
            { "UniformBound",                     &UniformBoundConstraintEntry::deserialize                   },
            { "Mixture",                          &MixtureConstraintEntry::deserialize,                       },
        };

        std::string type = n["type"].as<std::string>();
        auto i = deserializers.find(type);
        if (i == deserializers.end())
        {
            throw ConstraintDeserializationError(name, "unsupported type '" + type + "'");
        }

        return i->second(name, n);
    }

    std::string
    ConstraintEntry::serialize() const
    {
        YAML::Emitter out;
        this->serialize(out);

        return { out.c_str() };
    }

    /* Constraint */
    template <>
    struct WrappedForwardIteratorTraits<Constraint::BlockIteratorTag>
    {
        using UnderlyingIterator = std::vector<LogLikelihoodBlockPtr>::iterator;
    };
    template class WrappedForwardIterator<Constraint::BlockIteratorTag, LogLikelihoodBlockPtr>;

    template <>
    struct WrappedForwardIteratorTraits<Constraint::ObservableIteratorTag>
    {
        using UnderlyingIterator = ObservableSet::Iterator;
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

    using ConstraintFactory = std::function<Constraint (const QualifiedName &, const Options & options)>;

    template <typename Factory_>
    ConstraintFactory make_factory(const Factory_ & f)
    {
        return std::bind(&Factory_::make, f, std::placeholders::_1, std::placeholders::_2);
    }

    std::map<QualifiedName, std::shared_ptr<const ConstraintEntry>>
    load_constraint_entries()
    {
        Context context("When loading constraint entries:");

        using ValueType = std::map<QualifiedName, std::shared_ptr<const ConstraintEntry>>::value_type;

        std::map<QualifiedName, std::shared_ptr<const ConstraintEntry>> result;

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
                Context context("When parsing file '" + file_path.string() + "':");

                YAML::Node node;
                try
                {
                    node = YAML::LoadFile(file);
                }
                catch (YAML::ParserException & e)
                {
                    throw ConstraintInputFileParseError(file, e.what());
                }

                for (auto && p : node)
                {
                    std::string keyname = p.first.Scalar();

                    if ("@metadata@" == keyname)
                        continue;

                    Context context("When parsing constraint '" + keyname + "':");

                    QualifiedName name(keyname);
                    std::shared_ptr<const ConstraintEntry> entry{ ConstraintEntry::FromYAML(name, p.second) };

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

    class ConstraintEntries :
        public InstantiationPolicy<ConstraintEntries, Singleton>
    {
        private:
            std::map<QualifiedName, std::shared_ptr<const ConstraintEntry>> _entries;

            ConstraintEntries() :
                _entries(load_constraint_entries())
            {
            }

            ~ConstraintEntries() = default;

        public:
            friend class InstantiationPolicy<ConstraintEntries, Singleton>;

            const std::map<QualifiedName, std::shared_ptr<const ConstraintEntry>> & entries() const
            {
                return _entries;
            }

            void insert(const QualifiedName & key, const std::shared_ptr<const ConstraintEntry> & value)
            {
                _entries[key] = value;
            }
    };

    Constraint
    Constraint::make(const QualifiedName & name, const Options & options)
    {
        auto & entries = ConstraintEntries::instance()->entries();

        auto e = entries.find(name);
        if (e == entries.end())
            throw UnknownConstraintError(name);

        return e->second->make(e->first, name.options() + options); // options supersede name.options
    }

    template <>
    struct WrappedForwardIteratorTraits<Constraints::ConstraintIteratorTag>
    {
        using UnderlyingIterator = std::map<QualifiedName, std::shared_ptr<const ConstraintEntry>>::const_iterator;
    };
    template class WrappedForwardIterator<Constraints::ConstraintIteratorTag, const std::pair<const QualifiedName, std::shared_ptr<const ConstraintEntry>>>;

    template<>
    struct Implementation<Constraints>
    {
        const std::map<QualifiedName, std::shared_ptr<const ConstraintEntry>> constraint_entries;

        Implementation() :
            constraint_entries(ConstraintEntries::instance()->entries())
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

    std::shared_ptr<const ConstraintEntry>
    Constraints::operator[] (const QualifiedName & name) const
    {
        auto i = _imp->constraint_entries.find(name);

        if (_imp->constraint_entries.end() == i)
        {
            return {};
        }

        return i->second;
    }

    std::shared_ptr<const ConstraintEntry>
    Constraints::insert(const QualifiedName & name, const std::string & entry) const
    {
        std::shared_ptr<const ConstraintEntry> _entry(ConstraintEntry::FromYAML(name, entry));

        if (nullptr == _entry.get())
            throw InternalError("ConstraintEntry::FromYAML should never return nullptr");

        ConstraintEntries::instance()->insert(name, _entry);

        return _entry;
    }
}
