/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2014 Stephan Jahn
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

#ifndef EOS_GUARD_SRC_UTILS_ANALYSIS_C_HH
#define EOS_GUARD_SRC_UTILS_ANALYSIS_C_HH 1

#include <eos/utils/analysis.hh>

/// C like interface to handle an Analysic object
extern "C" {
    using namespace eos;

    /// constructor
    Analysis * EOS_Analysis_new(LogLikelihood * log_likelihood);

    /// destructor
    void EOS_Analysis_delete(Analysis * ana);

    /*!
     * Update parameter values in the underlying Parameters objects
     * and evaluate the Analysis.log_posterior.
     */
    double EOS_Analysis_log_posterior(Analysis * ana, const double * par_vals);

    /*!
     * Get info about analysis object.
     *
     * @note The returned pointer must explicitly be deleted by using \c free.
     */
    char * EOS_Analysis_info(Analysis * ana);

    /*!
     * Analysis.add Flat Prior
     *
     * @note The returned pointer must explicitly be deleted by using \c free.
     */
    char * EOS_Analysis_add_Flat(Analysis * ana, const char * par_name,
                                 const double range_min, const double range_max,
                                 bool nuisance);

    /*!
     * Analysis.add Gaussian Prior
     *
     * @note The returned pointer must explicitly be deleted by using \c free.
     */
    char * EOS_Analysis_add_Gauss(Analysis * ana, const char * par_name,
                                  const double range_min, const double range_max,
                                  const double lower, const double central, const double upper,
                                  bool nuisance);
    /*!
     * Analysis.add LogGamma Prior
     *
     * @note The returned pointer must explicitly be deleted by using \c free.
     */
    char * EOS_Analysis_add_LogGamma(Analysis * ana, const char * par_name,
                                     const double range_min, const double range_max,
                                     const double lower, const double central, const double upper,
                                     bool nuisance);
}

#endif
