/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_UTILS_SIGNAL_PDF_HH
#define EOS_GUARD_EOS_UTILS_SIGNAL_PDF_HH 1

#include <eos/utils/density.hh>
#include <eos/utils/exception.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>

#include <string>
#include <memory>

namespace eos
{
    class KinematicRange
    {
        public:
            const char * name;

            const double min, max;

            operator const char * () const
            {
                return name;
            }
    };

    class SignalPDF;

    typedef std::shared_ptr<SignalPDF> SignalPDFPtr;

    class SignalPDF :
        public Density
    {
        public:
            virtual const std::string & name() const = 0;

            virtual double evaluate() const = 0;

            virtual Kinematics kinematics() = 0;

            virtual Parameters parameters() = 0;

            virtual Options options() = 0;

            virtual DensityPtr clone() const = 0;

            virtual DensityPtr clone(const Parameters & parameters) const = 0;

            static SignalPDFPtr make(const std::string & name, const Parameters & parameters, const Kinematics & kinematics, const Options & options);
    };

    class SignalPDFFactory
    {
        public:
            SignalPDFFactory();

            virtual ~SignalPDFFactory();

            virtual SignalPDFPtr make(const Parameters &, const Kinematics &, const Options &) const = 0;
    };

    /*!
     * SignalPDFNameError is thrown when SignalPDF::make encounters a malformed observable name.
     */
    struct SignalPDFNameError :
        public Exception
    {
        ///@name Basic Functions
        ///@{
        /*!
         * Constructor.
         *
         * @param name The offending malformed observable name.
         */
        SignalPDFNameError(const std::string & name);
        ///@}
    };
}

#endif
