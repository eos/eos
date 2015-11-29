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

#include <eos/signal-pdf.hh>
#include <eos/b-decays/bs-to-kstar-l-nu.hh>
#include <eos/rare-b-decays/exclusive-b-to-s-dilepton-large-recoil.hh>
#include <eos/rare-b-decays/exclusive-b-to-s-dilepton-low-recoil.hh>
#include <eos/utils/concrete-signal-pdf.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <map>

namespace eos
{
    namespace test
    {
        // PDF = (1/2 L_0 + 1/3 L_1 + 1/4 L_2) / 2
        class Legendre1DPDF
        {
            public:
                Legendre1DPDF(const Parameters &, const Options &)
                {
                }

                double result(const double & z) const
                {
                    return (9.0 + 8.0 * z + 9.0 * z * z) / 24.0;
                }
        };
    }

    SignalPDFNameError::SignalPDFNameError(const std::string & name) :
        Exception("SignalPDF name '" + name + "' is malformed")
    {
    }

    SignalPDFFactory::SignalPDFFactory()
    {
    }

    SignalPDFFactory::~SignalPDFFactory()
    {
    }

    template <typename Decay_, typename Tuple_, typename ... Args_>
    std::pair<std::string, SignalPDFFactory *> make_signal_pdf(const char * name,
            double (Decay_::* function)(const Args_ & ...) const,
            const Tuple_ & kinematic_ranges)
    {
        std::string sname(name);

        return std::make_pair(sname, make_concrete_observable_factory(sname, function, kinematic_ranges));
    }

    SignalPDFPtr
    SignalPDF::make(const std::string & _name, const Parameters & parameters, const Kinematics & kinematics, const Options & _options)
    {
        static const std::map<std::string, SignalPDFFactory *> signal_pdfs
        {
            /* Internal Tests */

            make_signal_pdf("Test::Legendre1D",
                    &test::Legendre1DPDF::result,
                    std::make_tuple( KinematicRange{ "z", -1.0, +1.0 } )),

            /* Exclusive Decays */

            /* Exclusive B Decays */

            make_signal_pdf("B_s->K^*lnu::d^4Gamma",
                &BsToKstarLeptonNeutrino::four_differential_decay_width,
                std::make_tuple( KinematicRange{ "s", 0.02, 19.71 },
                    KinematicRange{ "cos(theta_l)", -1.0, +1.0 },
                    KinematicRange{ "cos(theta_k)", -1.0, +1.0 },
                    KinematicRange{ "phi", 0.0, 2.0 * M_PI })),

            /* Exclusive Rare B Decays */

            make_signal_pdf("B->Kll::d^2Gamma@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::two_differential_decay_width,
                    std::make_tuple( KinematicRange{ "s", 0.04, 22.87 },
                        KinematicRange{ "cos(theta_l)", -1.0, +1.0 })),

            make_signal_pdf("B->Kll::d^2Gamma@LowRecoil",
                    &BToKDilepton<LowRecoil>::two_differential_decay_width,
                    std::make_tuple( KinematicRange{ "s", 0.00, 22.87 },
                        KinematicRange{ "cos(theta_l)", -1.0, +1.0 })),

            make_signal_pdf("B->K^*ll::d^4Gamma@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::four_differential_decay_width,
                    std::make_tuple( KinematicRange{ "s", 1.00, 6.00 },
                        KinematicRange{ "cos(theta_l)", -1.0, +1.0 },
                        KinematicRange{ "cos(theta_k)", -1.0, +1.0 },
                        KinematicRange{ "phi", 0.0, 2.0 * M_PI })),

            make_signal_pdf("B->K^*ll::d^4Gamma@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::four_differential_decay_width,
                    std::make_tuple( KinematicRange{ "s", 0.00, 19.21 },
                        KinematicRange{ "cos(theta_l)", -1.0, +1.0 },
                        KinematicRange{ "cos(theta_k)", -1.0, +1.0 },
                        KinematicRange{ "phi", 0.0, 2.0 * M_PI })),
        };

        Options options;
        std::string name(_name);

        std::string::size_type pos;
        while (std::string::npos != (pos = name.rfind(',')))
        {
            std::string::size_type sep(name.find('=', pos + 1));
            if (std::string::npos == sep)
                throw SignalPDFNameError(_name);

            std::string key(name.substr(pos + 1, sep - pos - 1));
            std::string value(name.substr(sep + 1));

            options.set(key, value);
            name.erase(pos);
        }

        auto i = signal_pdfs.find(name);
        if (signal_pdfs.end() == i)
            return SignalPDFPtr();

        return i->second->make(parameters, kinematics, options + _options);
    }
}
