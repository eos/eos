/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2025 Marta Burgos
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

#ifndef EOS_GUARD_EOS_NONLEPTONIC_AMPLITUDES_QCDF_COEFFICIENTS_HH
#define EOS_GUARD_EOS_NONLEPTONIC_AMPLITUDES_QCDF_COEFFICIENTS_HH 1

#include <eos/maths/complex.hh>
#include <eos/maths/power-of.hh>
#include <eos/models/model.hh>
//#include <eos/models/wilson-coefficients.hh>
#include <eos/nonleptonic-amplitudes/nonleptonic-amplitudes-fwd.hh>
#include <eos/nonleptonic-amplitudes/nonleptonic-amplitudes.hh>
#include <eos/utils/parameters.hh>
#include <eos/form-factors/psd-lcdas.hh>


#include <array>
#include <cmath>

namespace eos
{
    template <typename Transition_> class QCDFCoefficients;

    template <> class QCDFCoefficients<PToPP> : public NonleptonicAmplitudes<PToPP>
    {
        private:
            std::shared_ptr<Model> model;
            QuarkFlavorOption opt_q;
            LightMesonOption opt_p1;
            LightMesonOption opt_p2;


            UsedParameter mB;
            UsedParameter mb;
            UsedParameter mu;
            UsedParameter mus;
            UsedParameter mB_q_0;
            UsedParameter mP1;
            UsedParameter mP2;
            UsedParameter FP1;
            UsedParameter FP2;
            UsedParameter fB;
            UsedParameter fP1;
            UsedParameter fP2;

            std::shared_ptr<PseudoscalarLCDAs> lcdasP1;
            std::shared_ptr<PseudoscalarLCDAs> lcdasP2;

            eos::WilsonCoefficients<BToS> wc2;

            complex<double> lcda;

            complex<double> a1();
            complex<double> a2();
            complex<double> a3();
            complex<double> a4();
            complex<double> a5();
            complex<double> a6();
            complex<double> a7();
            complex<double> a8();
            complex<double> a9();
            complex<double> a10();


            complex<double> b1();
            complex<double> b2();
            complex<double> b3u();
            complex<double> b3c();
            complex<double> b4u();
            complex<double> b4c();
            complex<double> b3EWu();
            complex<double> b3EWc();
            complex<double> b4EWu();
            complex<double> b4EWc();

            static const std::vector<OptionSpecification> options;

        public:

            QCDFCoefficients(const Parameters & p, const Options & o);

            ~QCDFCoefficients() {};

            inline void
            update() const
            {

            }

            static NonleptonicAmplitudes<PToPP> * make(const Parameters &, const Options &);

            // Helper functions

            //Veretex functions
            complex<double> gvertex(const double & x) const;
            std::array<complex<double>, 11> Vertex(const double & x) const;

            // Hard scattering functions
            std::array<complex<double>, 11> HardSpec(const double & x, const double & y) const;

            // Penguin functions
            complex<double> GM2(const double & s) const;
            complex<double> GM2hat(const double & s) const;
            complex<double> Gsx(const double & s,const double & x) const;
            std::array<complex<double>, 11> Penguin(const double & x) const;

            //wilson coefficients: erase later
            //std::array<complex<double>, 11> WilsonCoefficients(const double & x) const;

    };

}

#endif
