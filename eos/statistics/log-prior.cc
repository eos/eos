/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Frederik Beaujean
 * Copyright (c) 2017-2024 Danny van Dyk
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

#include <eos/statistics/log-prior.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/log.hh>
#include <eos/maths/power-of.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <cmath>
#include <limits>
#include <numeric>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_spline.h>

#include <config.h>

#ifdef EOS_USE_GSL_LINALG_CHOLESKY_DECOMP
#  if (EOS_USE_GSL_LINALG_CHOLESKY_DECOMP == 1)
#    define GSL_LINALG_CHOLESKY_DECOMP gsl_linalg_cholesky_decomp
#  else
#    define GSL_LINALG_CHOLESKY_DECOMP gsl_linalg_cholesky_decomp1
#  endif
#else
#  error EOS_USE_GSL_LINALG_CHOLESKY_DECOMP not defined.
#endif

namespace eos
{
    template <>
    struct WrappedForwardIteratorTraits<LogPrior::IteratorTag>
    {
        using UnderlyingIterator = std::vector<Parameter>::iterator;
    };

    namespace priors
    {
        struct RangeError :
            public Exception
        {
            RangeError(const std::string & message) throw () :
                Exception("Range Error: " + message)
            {
            }
        };

        struct UnknownPriorError :
            public Exception
        {
            UnknownPriorError(const std::string & message) throw () :
                Exception("Unknown prior error: " + message)
            {
            }
        };

        /*!
         * Flat or uniform prior
         */
        class Flat :
            public LogPrior
        {
            private:
                Parameter _parameter;

                std::string _name;

                const double _min, _max;

                // The flat prior always returns this value
                double _value;

            public:
                Flat(const Parameters & parameters, const std::string & name, const double & min, const double & max) :
                    LogPrior(parameters),
                    _parameter(parameters[name]),
                    _name(name),
                    _min(min),
                    _max(max),
                    _value(std::log(1.0 / (_max - _min)))
                {
                    if (_min >= _max)
                    {
                        throw RangeError("LogPrior::Flat(" + _name +"): minimum (" + stringify(_min)
                                          + ") must be smaller than maximum (" + stringify(_max) + ")");
                    }
                    _varied_parameters.push_back(_parameter);
                }

                virtual ~Flat()
                {
                }

                virtual std::string as_string() const
                {
                    std::string result = "Parameter: " + _name + ", prior type: flat, range: [" + stringify(_min) + "," + stringify(_max) + "]";
                    return result;
                }

                virtual double operator()() const
                {
                    return _value;
                }

                virtual LogPriorPtr clone(const Parameters & parameters) const
                {
                    return LogPriorPtr(new priors::Flat(parameters, _name, _min, _max));
                }

                virtual void sample()
                {
                    _parameter.set(_parameter.evaluate_generator() * (_max - _min) + _min);
                }

                virtual void compute_cdf()
                {
                    _parameter.set_generator((_parameter.evaluate() - _min) / (_max - _min));
                }

                virtual bool informative() const
                {
                    return false;
                }
        };

        /*!
         * [asymmetric] Gaussian or Normal prior distribution with finite support
         */
        class CurtailedGauss :
            public LogPrior
        {
            private:
                Parameter _parameter;

                const std::string _name;

                const double _min, _max;

                const double _lower, _central, _upper;

                const double _sigma_lower, _sigma_upper;

                // coefficients needed for sampling from asymmetric Gaussian on finite support
                // the pdf is a piecewise function of y given x^{+a}_{-b}
                // P(y|x,a,b) = \theta(y-x) c_a \mathcal{N}(y|x,a)  + \theta(x-y) c_b \mathcal{N}(y|x,b)
                // Fix c_a, c_b by requiring that the PDF
                // a) be continuous at x
                // b) integrate to one over the range
                const double _c_a, _c_b;

                // The probability covered to the left of the central value. It is just
                // c_b (1/2 - \Phi(y_-|x,b)
                // but for the sampling it is faster to precompute the call to Gaussian cumulative.
                const double _prob_lower;

                // PDF normalization factor precomputed for operator()
                const double _norm_lower, _norm_upper;

            public:
                CurtailedGauss(const Parameters & parameters, const std::string & name, const double & min, const double & max,
                        const double & lower, const double & central, const double & upper) :
                    LogPrior(parameters),
                    _parameter(parameters[name]),
                    _name(name),
                    _min(min),
                    _max(max),
                    _lower(lower),
                    _central(central),
                    _upper(upper),
                    _sigma_lower(central - lower),
                    _sigma_upper(upper - central),
                    _c_a(1.0 / ((_sigma_lower / _sigma_upper) * (0.5 - gsl_cdf_gaussian_P(min - central, _sigma_lower))
                                 + gsl_cdf_gaussian_P(max - central, _sigma_upper) - 0.5)),
                    _c_b(_sigma_lower/ _sigma_upper * _c_a),
                    _prob_lower(_c_b * (0.5 - gsl_cdf_gaussian_P(min - central, _sigma_lower))),
                    _norm_lower(std::log(_c_b / std::sqrt(2 * M_PI) / _sigma_lower)),
                    _norm_upper(std::log(_c_a / std::sqrt(2 * M_PI) / _sigma_upper))
                {
                    if (min >= max)
                    {
                        throw RangeError("LogPrior::Gauss(" + _name +"): minimum (" + stringify(_min)
                                          + ") must be smaller than maximum (" + stringify(_max) + ")");
                    }
                    _varied_parameters.push_back(_parameter);
                }

                virtual ~CurtailedGauss()
                {
                }

                virtual std::string as_string() const
                {
                    std::string result = "Parameter: " + _name + ", prior type: Gaussian, range: [" + stringify(_min) + "," + stringify(_max) + "]";

                    result += ", x = " + stringify(_central);
                    if (std::abs(_sigma_upper - _sigma_lower) < 1e-15)
                    {
                        result += " +- " + stringify(_sigma_upper);
                    }
                    else
                    {
                        result += " + " + stringify(_sigma_upper) + " - " + stringify(_sigma_lower);
                    }
                    return result;
                }

                virtual double operator()() const
                {
                    double sigma = 0.0, norm = 0.0;

                    // read parameter's current value
                    double x = _parameter.evaluate();

                    if (x < _central)
                    {
                        sigma = _sigma_lower;
                        norm = _norm_lower;
                    }
                    else
                    {
                        sigma = _sigma_upper;
                        norm = _norm_upper;
                    }
                    return norm - 0.5 * power_of<2>((x - _central) / sigma);
                }

                virtual LogPriorPtr clone(const Parameters & parameters) const
                {
                    return LogPriorPtr(new priors::CurtailedGauss(parameters, _name, _min, _max, _lower, _central, _upper));
                }

                virtual void sample()
                {
                    // CDF = c \Phi(x - x_{central} / \sigma) + b

                    // find out if sample in upper or lower part
                    const auto p = _parameter.evaluate_generator();

                    if (p < _prob_lower)
                       _parameter.set(gsl_cdf_gaussian_Pinv((p - _prob_lower) / _c_b + 0.5,  _sigma_lower) + _central);
                    else
                       _parameter.set(gsl_cdf_gaussian_Pinv((p - _prob_lower) / _c_a + 0.5,  _sigma_upper) + _central);
                }

                virtual void compute_cdf()
                {
                    // CDF = c \Phi(x - x_{central} / \sigma) + b

                    // find out if sample in upper or lower part
                    const auto x = _parameter.evaluate();

                    if (x < _central)
                        _parameter.set_generator(_c_b * gsl_cdf_gaussian_P((x - _central) / _sigma_lower, 1.0) + _prob_lower);
                    else
                        _parameter.set_generator(_c_a * gsl_cdf_gaussian_P((x - _central) / _sigma_upper, 1.0) + _prob_lower);
                }

                virtual bool informative() const
                {
                    return true;
                }
        };

        /*!
         * Prior distribution for renormalization scales
         */
        class Scale :
            public LogPrior
        {
            private:
                Parameter _parameter;

                const std::string _name;

                const double _mu_0, _lambda;

                const double _min, _max;

                const double _ln_lambda;

            public:
                Scale(const Parameters & parameters, const std::string & name, const double & min, const double & max, const double & mu_0, const double & lambda) :
                    LogPrior(parameters),
                    _parameter(parameters[name]),
                    _name(name),
                    _mu_0(mu_0),
                    _lambda(lambda),
                    _min(std::max(min, mu_0 / lambda)),
                    _max(std::min(max, mu_0 * lambda)),
                    _ln_lambda(std::log(lambda))
                {
                    _varied_parameters.push_back(_parameter);
                }

                virtual ~Scale()
                {
                }

                virtual std::string as_string() const
                {
                    std::string result = "Parameter: " + _name + ", prior type: Scale, range: [" + stringify(_mu_0 / _lambda) + "," + stringify(_mu_0 * _lambda) + "]";
                    result += ", mu_0 = " + stringify(_mu_0) + ", lambda = " + stringify(_lambda);

                    return result;
                }

                virtual double operator()() const
                {
                    // read parameter's current value
                    double x = _parameter.evaluate();

                    if ((x < _min) || (_max < x))
                        return -std::numeric_limits<double>::infinity();

                    return 1.0 / (2.0 * _ln_lambda * x);
                }

                virtual LogPriorPtr clone(const Parameters & parameters) const
                {
                    return LogPriorPtr(new priors::Scale(parameters, _name, _min, _max, _mu_0, _lambda));
                }

                virtual void sample()
                {
                    // CDF: p = [\ln x - \ln \mu_0 + \ln \lambda] / (2.0 \ln \lambda)
                    // inverse CDF: x = \mu_0 * \lambda^(2 p - 1)

                    _parameter.set(_mu_0 * std::pow(_lambda, 2.0 * _parameter.evaluate_generator() - 1.0));
                }

                virtual void compute_cdf()
                {
                    // CDF: p = [\ln x - \ln \mu_0 + \ln \lambda] / (2.0 \ln \lambda)
                    // inverse CDF: x = \mu_0 * \lambda^(2 p - 1)

                    _parameter.set_generator((std::log(_parameter.evaluate() / _mu_0) + _ln_lambda) / (2.0 * _ln_lambda));
                }

                virtual bool informative() const
                {
                    return true;
                }
        };

        /*!
         * Gaussian prior distribution
         */
        class Gaussian :
            public LogPrior
        {
            private:
                Parameter _parameter;

                QualifiedName _name;

                const double _mu;

                const double _sigma;

                const double _ln_norm;

            public:
                Gaussian(const Parameters & parameters, const QualifiedName & name, const double & mu, const double & sigma) :
                    LogPrior(parameters),
                    _parameter(parameters[name]),
                    _name(name),
                    _mu(mu),
                    _sigma(sigma),
                    _ln_norm(-0.5 * std::log(2 * M_PI) - std::log(sigma))
                {
                }

                ~Gaussian() = default;

                virtual std::string as_string() const
                {
                    std::string result = "Parameter: " + _name.full() + ", prior type: gaussian";
                    return result;
                }

                virtual double operator()() const
                {
                    double x = _parameter.evaluate();

                    return _ln_norm - 0.5 * power_of<2>((x - _mu) / _sigma);
                }

                virtual LogPriorPtr clone(const Parameters & parameters) const
                {
                    return LogPriorPtr(new priors::Gaussian(parameters, _name, _mu, _sigma));
                }

                virtual void sample()
                {
                    const double u = _parameter.evaluate_generator();
                    const double x = gsl_cdf_gaussian_Pinv(u, _sigma) + _mu;
                    _parameter.set(x);
                }

                virtual void compute_cdf()
                {
                    const double x = _parameter.evaluate();
                    const double u = gsl_cdf_gaussian_P(x - _mu, _sigma);
                    _parameter.set_generator(u);
                }

                virtual bool informative() const
                {
                    return true;
                }
        };

        /*!
         * Multivariate Gaussian prior distribution
         */
        class MultivariateGaussian :
            public LogPrior
        {
            private:
                std::vector<Parameter> _parameters;

                std::vector<QualifiedName> _names;

                const unsigned _dim;

                // inputs
                gsl_vector * const _mean;
                gsl_matrix * const _covariance;

                // the normalization constant of the density
                const double _norm;

                // cholesky matrix of covariance, inverse of the cholesky matrix, and inverse of covariance
                gsl_matrix * _chol;
                gsl_matrix * _chol_inv;
                gsl_matrix * _covariance_inv;

                // temporary storage for evaluation
                gsl_vector * _observables;
                gsl_vector * _measurements;
                gsl_vector * _measurements_2;

            public:
                MultivariateGaussian(const Parameters & parameters, const std::vector<QualifiedName> & names, gsl_vector * mean, gsl_matrix * covariance) :
                    LogPrior(parameters),
                    _names(names),
                    _dim(names.size()),
                    _mean(mean),
                    _covariance(covariance),
                    _norm(compute_norm()),
                    _chol(gsl_matrix_alloc(covariance->size1, covariance->size2)),
                    _chol_inv(gsl_matrix_alloc(covariance->size1, covariance->size2)),
                    _covariance_inv(gsl_matrix_alloc(covariance->size1, covariance->size2)),
                    _observables(gsl_vector_alloc(_dim)),
                    _measurements(gsl_vector_alloc(_dim)),
                    _measurements_2(gsl_vector_alloc(_dim))
                {
                    if (_covariance->size1 != _covariance->size2)
                        throw InternalError("priors::MultivariateGaussian: covariance matrix is not a square matrix");

                    if (_covariance->size1 != _mean->size)
                        throw InternalError("priors::MultivariateGaussian: number of parameters and dimension of covariance matrix are not identical");

                    if (_dim != _mean->size)
                        throw InternalError("priors::MultivariateGaussian: number of parameters and dimension of mean vector are not identical");

                    for (auto & n : names)
                    {
                        const auto & param = parameters[n];
                        _parameters.push_back(param);
                        _varied_parameters.push_back(param);
                    }

                    // cholesky decomposition (informally: the sqrt of the covariance matrix)
                    // the GSL matrix contains both the cholesky and its transpose, see GSL reference, ch. 14.5
                    cholesky();
                    invert_covariance();

                    // keep only the lower and diagonal parts, set upper parts to zero
                    for (unsigned i = 0; i < _dim ; ++i)
                    {
                        for (unsigned j = i + 1 ; j < _dim ; ++j)
                        {
                            gsl_matrix_set(_chol, i, j, 0.0);
                        }
                    }
                }

                virtual ~MultivariateGaussian()
                {
                    gsl_matrix_free(_covariance_inv);
                    gsl_matrix_free(_chol_inv);
                    gsl_matrix_free(_chol);
                    gsl_matrix_free(_covariance);

                    gsl_vector_free(_measurements_2);
                    gsl_vector_free(_measurements);
                    gsl_vector_free(_observables);
                    gsl_vector_free(_mean);
                }

                virtual std::string as_string() const
                {
                    throw InternalError("priors::MultivariateGaussian::as_string() not implemented");

                    return "";
                }

                // compute the normalization constant on log scale
                // -k/2 * log 2 Pi - 1/2 log(abs(det(V^{-1})))
                double compute_norm()
                {
                    // copy covariance matrix
                    gsl_matrix * M = gsl_matrix_alloc(_dim, _dim);
                    gsl_matrix_memcpy(M, _covariance);

                    // find LU decomposition
                    int signum = 0;
                    gsl_permutation * p = gsl_permutation_alloc(_dim);
                    gsl_linalg_LU_decomp(M, p, &signum);

                    // calculate determinant
                    const double log_det = gsl_linalg_LU_lndet(M);

                    gsl_permutation_free(p);
                    gsl_matrix_free(M);

                    return -0.5 * _dim * std::log(2 * M_PI) - 0.5 * log_det;
                }

                // compute cholesky decomposition of covariance matrix
                void cholesky()
                {
                    // copy covariance matrix
                    gsl_matrix_memcpy(_chol, _covariance);
                    if (GSL_SUCCESS != GSL_LINALG_CHOLESKY_DECOMP(_chol))
                    {
                        throw InternalError("priors::MultivariateGaussian: Cholesky decomposition failed");
                    }
                }

                // invert covariance matrix based on previously obtained Cholesky decomposition
                void invert_covariance()
                {
                    // copy cholesky matrix
                    gsl_matrix_memcpy(_covariance_inv, _chol);

                    // compute inverse matrix from cholesky
                    if (GSL_SUCCESS != gsl_linalg_cholesky_invert(_covariance_inv))
                    {
                        throw InternalError("priors::MultivariateGaussian: Cholesky inversion failed");
                    }

                    // delete upper part of _chol, which contains the elements of the original covariance matrix
                    for (unsigned i = 0u; i < _dim ; ++i)
                    {
                        for (unsigned j = i + 1u ; j < _dim ; ++j)
                        {
                            gsl_matrix_set(_chol, i, j, 0.0);
                        }
                    }

                    // copy cholesky matrix to _chol_inv and invert in place
                    // inverse will be stored in upper triangle of _chol_inv
                    gsl_matrix_memcpy(_chol_inv, _chol);
                    gsl_linalg_tri_invert(CblasUpper, CblasNonUnit, _chol_inv);

                    // delete lower part of _chol_inv, which contains the elements of the original cholesky matrix
                    for (unsigned i = 0u; i < _dim ; ++i)
                    {
                        for (unsigned j = 0u ; j < i ; ++j)
                        {
                            gsl_matrix_set(_chol_inv, i, j, 0.0);
                        }
                    }
                }

                virtual double operator()() const
                {
                    // read parameters
                    for (auto i = 0u ; i < _dim ; ++i)
                    {
                        gsl_vector_set(_observables, i, _parameters[i].evaluate());
                    }

                    // prepare for centering
                    //   measurements <- mean
                    gsl_vector_memcpy(_measurements, _mean);

                    // center the gaussian:
                    //   measurements <- measurements - observables
                    gsl_blas_daxpy(-1.0, _observables, _measurements);

                    // measurements_2 <- inv(covariance) * measurements
                    gsl_blas_dgemv(CblasNoTrans, 1.0, _covariance_inv, _measurements, 0.0, _measurements_2);

                    double chi_square = 0.0;
                    gsl_blas_ddot(_measurements, _measurements_2, &chi_square);

                    return _norm - 0.5 * chi_square;
                }

                virtual LogPriorPtr clone(const Parameters & parameters) const
                {
                    gsl_vector * mean = gsl_vector_alloc(_dim);
                    gsl_vector_memcpy(mean, _mean);

                    gsl_matrix * covariance = gsl_matrix_alloc(_dim, _dim);
                    gsl_matrix_memcpy(covariance, _covariance);

                    return LogPriorPtr(new priors::MultivariateGaussian(parameters, _names, mean, covariance));
                }

                virtual void sample()
                {
                    // generate standard normals in _measurements
                    for (auto i = 0u ; i < _dim ; ++i)
                    {
                        const double u = _parameters[i].evaluate_generator();
                        const double z = gsl_cdf_ugaussian_Pinv(u);
                        gsl_vector_set(_measurements, i, z);
                    }

                    // transform: measurements_2 <- _chol * measurements
                    gsl_blas_dgemv(CblasNoTrans, 1.0, _chol, _measurements, 0.0, _measurements_2);

                    // transform: measurements_2 <- measurements_2 + _mean
                    gsl_vector_add(_measurements_2, _mean);

                    // set parameters
                    for (auto i = 0u ; i < _dim ; ++i)
                    {
                        _parameters[i].set(gsl_vector_get(_measurements_2, i));
                    }
                }

                virtual void compute_cdf()
                {
                    // get parameters
                    for (auto i = 0u ; i < _dim ; ++i)
                    {
                        gsl_vector_set(_measurements_2, i, _parameters[i].evaluate());
                    }

                    // transform: measurements_2 <- measurements_2 - _mean
                    gsl_vector_sub(_measurements_2, _mean);

                    // transform: measurements <- _chol_inv * measurements_2
                    gsl_blas_dgemv(CblasNoTrans, 1.0, _chol_inv, _measurements_2, 0.0, _measurements);

                    // compute cdfs for standard normals in _measurements
                    for (auto i = 0u ; i < _dim ; ++i)
                    {
                        const double z = gsl_vector_get(_measurements, i);
                        const double u = gsl_cdf_ugaussian_P(z);
                        _parameters[i].set_generator(u);
                    }
                }

                virtual bool informative() const
                {
                    return true;
                }
        };

        /*!
         * Poisson prior distribution
         */
        class Poisson :
            public LogPrior
        {
            private:
                Parameter _parameter;

                std::string _name;

                const double _k;

                const double _ln_norm;

            public:
                Poisson(const Parameters & parameters, const std::string & name, const double & k) :
                    LogPrior(parameters),
                    _parameter(parameters[name]),
                    _name(name),
                    _k(k),
                    _ln_norm(-1.0 * std::log(gsl_sf_gamma(_k + 1.0)) + std::log(_k))
                    {
                        _varied_parameters.push_back(_parameter);
                    }

                virtual ~Poisson()
                {
                }

                virtual std::string as_string() const
                {
                    std::string result = "Parameter: " + _name + ", prior type: poisson";
                    return result;
                }

                virtual double operator()() const
                {
                    double lambda = _parameter.evaluate() * _k;

                    return _ln_norm - lambda + _k * std::log(lambda);
                }

                virtual LogPriorPtr clone(const Parameters & parameters) const
                {
                    return LogPriorPtr(new priors::Poisson(parameters, _name, _k));
                }

                virtual void sample()
                {
                    const double u = _parameter.evaluate_generator();
                    const double lambda = gsl_cdf_gamma_Pinv(u, _k + 1.0, 1.0);
                    _parameter.set(lambda / _k);
                }

                virtual void compute_cdf()
                {
                    const double lambda = _parameter.evaluate() * _k;
                    const double u = gsl_cdf_gamma_P(lambda, _k + 1.0, 1.0);
                    _parameter.set_generator(u);
                }

                virtual bool informative() const
                {
                    return true;
                }
        };
    }

    LogPrior::LogPrior(const Parameters & parameters) :
        _parameters(parameters)
    {
    }

    LogPrior::Iterator
    LogPrior::begin()
    {
        return LogPrior::Iterator(_varied_parameters.begin());
    }

    LogPrior::Iterator
    LogPrior::end()
    {
        return LogPrior::Iterator(_varied_parameters.end());
    }

    LogPriorPtr
    LogPrior::Flat(const Parameters & parameters, const std::string & name, const double & min, const double & max)
    {
        LogPriorPtr prior = std::make_shared<eos::priors::Flat>(parameters, name, min, max);

        return prior;
    }

    LogPriorPtr
    LogPrior::CurtailedGauss(const Parameters & parameters, const std::string & name, const double & min, const double & max,
            const double & lower, const double & central, const double & upper)
    {
        // check input
        if (lower >= central)
            throw InternalError("LogPrior::Gauss: lower value (" + stringify(lower) + ") >= central value (" + stringify(central) + ")");

        if (upper <= central)
            throw InternalError("LogPrior::Gauss: upper value (" + stringify(upper) + ") <= central value (" + stringify(central) + ")");

        LogPriorPtr prior = std::make_shared<eos::priors::CurtailedGauss>(parameters, name, min, max, lower, central, upper);

        return prior;
    }

    LogPriorPtr
    LogPrior::Scale(const Parameters & parameters, const std::string & name, const double & min, const double & max,
            const double & mu_0, const double & lambda)
    {
        // check input
        if (mu_0 <= 0.0)
            throw InternalError("LogPrior::Scale: default value mu_0 must be strictly positive");

        if (lambda <= 1.0)
            throw InternalError("LogPrior::Scale: scale factor lambda must be strictly larger than 1");

        LogPriorPtr prior = std::make_shared<eos::priors::Scale>(parameters, name, min, max, mu_0, lambda);

        return prior;
    }

    LogPriorPtr
    LogPrior::Gaussian(const Parameters & parameters, const QualifiedName & name, const double & mu, const double & sigma)
    {
        LogPriorPtr prior = std::make_shared<eos::priors::Gaussian>(parameters, name, mu, sigma);

        return prior;
    }

    LogPriorPtr
    LogPrior::MultivariateGaussian(const Parameters & parameters, const std::vector<QualifiedName> & names,
            gsl_vector * mean, gsl_matrix * covariance)
    {
        LogPriorPtr prior = std::make_shared<eos::priors::MultivariateGaussian>(parameters, names, mean, covariance);

        return prior;
    }

    LogPriorPtr
    LogPrior::Poisson(const Parameters & parameters, const std::string & name, const double & k)
    {
        LogPriorPtr prior = std::make_shared<eos::priors::Poisson>(parameters, name, k);

        return prior;
    }

    LogPriorPtr
    LogPrior::Make(const Parameters & parameters, const std::string & s)
    {
        // extract parameter name
        auto loc1 = s.find_first_of(':');
        auto loc2 = s.find_first_of(',');

        std::string par_name = s.substr(loc1 + 2, loc2 - loc1 - 2);

        // extract type
        loc1 = s.find_first_of(':', loc2 + 1);
        loc2 = s.find_first_of(',', loc2 + 1);

        std::string prior_type = s.substr(loc1 + 2, loc2 - loc1 - 2);

        // extract range
        loc1 = s.find_first_of('[', loc2 + 1);
        loc2 = s.find_first_of(',', loc2 + 1);
        std::string range_min = s.substr(loc1 + 1, loc2 - loc1 - 1);

        loc1 = s.find_first_of(',', loc1 + 1);
        loc2 = s.find_first_of(']', loc1 + 1);
        std::string range_max = s.substr(loc1 + 1, loc2 - loc1 - 1);

        const auto & [min, max] = std::tuple(destringify<double>(range_min), destringify<double>(range_max));

        if (prior_type == "flat")
        {
            return Flat(parameters, par_name, min, max);
        }
        if (prior_type == "Gaussian")
        {
            // extract central value
            loc1 = s.find_first_of('=', loc2 + 1);
            loc2 = s.find_first_of('+', loc2 + 1);

            double central = destringify<double>(s.substr(loc1 + 2, loc2 - loc1 - 3));

            double sigma_upper, sigma_lower;

            // extract sigma_upper, lower
            if ( s[loc2 + 1] == '-')
            {
                sigma_upper = destringify<double>(s.substr(loc2 + 3));
                sigma_lower = sigma_upper;
            }
            else
            {
                loc1 = loc2;
                loc2 = s.find_first_of('-', loc2 + 1);

                sigma_upper = destringify<double>(s.substr(loc1 + 2, loc2 - loc1 - 3));

                loc1 = loc2;
                loc2 = s.find_first_of(',', loc2 + 1);

                // Gaussian
                if (loc2 == std::string::npos)
                {
                    sigma_lower = destringify<double>(s.substr(loc1 + 1));
                }
                // LogGamma
                // don't parse until end of string, but until next comma
                else
                {
                    sigma_lower = destringify<double>(s.substr(loc1 + 2, loc2 - loc1 - 2));
                }
            }

            if (prior_type == "Gaussian")
                return CurtailedGauss(parameters, par_name, min, max, central - sigma_lower, central, central + sigma_upper);
        }

        throw priors::UnknownPriorError("Cannot construct prior from '" + s + "'");
    }

    template class WrappedForwardIterator<LogPrior::IteratorTag, Parameter>;
}
