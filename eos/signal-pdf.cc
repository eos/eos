/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015, 2016, 2017 Danny van Dyk
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
#include <eos/b-decays/b-to-d-l-nu.hh>
#include <eos/b-decays/b-to-d-l-x-nu.hh>
#include <eos/b-decays/b-to-pi-l-x-nu.hh>
#include <eos/b-decays/b-to-pi-pi-l-nu.hh>
#include <eos/b-decays/bs-to-kstar-l-nu.hh>
#include <eos/b-decays/lambdab-to-lambdac2625-l-nu.hh>
#include <eos/rare-b-decays/exclusive-b-to-s-dilepton-large-recoil.hh>
#include <eos/rare-b-decays/exclusive-b-to-s-dilepton-low-recoil.hh>
#include <eos/utils/concrete-signal-pdf.hh>
#include <eos/utils/density.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <map>
#include <ostream>

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

                static const std::string description;
        };

        const std::string Legendre1DPDF::description = "1D PDF up to 2nd order in z; used for unit tests only.";
    }

    SignalPDFNameError::SignalPDFNameError(const std::string & name) :
        Exception("SignalPDF name '" + name + "' is malformed")
    {
    }

    SignalPDFEntry::SignalPDFEntry()
    {
    }

    SignalPDFEntry::~SignalPDFEntry() = default;

    template <>
    struct WrappedForwardIteratorTraits<SignalPDFEntry::KinematicRangeIteratorTag>
    {
        typedef const KinematicRange * UnderlyingIterator;
    };
    template class WrappedForwardIterator<SignalPDFEntry::KinematicRangeIteratorTag, const KinematicRange>;

    std::ostream &
    SignalPDFEntry::insert(std::ostream & os) const
    {
        os << "<empty SignalPDF description>" << std::endl;
        return os;
    }

    template <typename Decay_, typename ... FunctionArgs_, typename ... KinematicRanges_>
    std::pair<QualifiedName, SignalPDFEntry *> make_signal_pdf(const char * name,
            double (Decay_::* function)(const FunctionArgs_ & ...) const,
            const Options & default_options,
            const KinematicRanges_ & ... kinematic_ranges)
    {
        static_assert(sizeof...(FunctionArgs_) == sizeof...(KinematicRanges_), "Need as many function arguments as kinematics ranges!");

        QualifiedName qn(name);

        return std::make_pair(qn, make_concrete_signal_pdf_entry(qn, function, default_options, kinematic_ranges...));
    }

    const std::map<QualifiedName, const SignalPDFEntry *> &
    make_signal_pdf_entries()
    {
        static const std::map<QualifiedName, const SignalPDFEntry *> signal_pdf_entries
        {
            /* Internal Tests */

            make_signal_pdf("Test::Legendre1D",
                    &test::Legendre1DPDF::result,
                    Options{ },
                    KinematicRange{ "z", -1.0, +1.0, "" }),

            /* Exclusive Decays */

            /* Exclusive B Decays */

            make_signal_pdf("B->pipimunu::d^3Gamma@QCDF",
                    &BToPiPiLeptonNeutrino::triple_differential_branching_ratio,
                    Options{ },
                    KinematicRange{ "q2", 0.01, 0.93859, BToPiPiLeptonNeutrino::kinematics_description_q2 },
                    KinematicRange{ "k2", 18.582, 27.872, BToPiPiLeptonNeutrino::kinematics_description_k2 },
                    KinematicRange{ "cos(theta)", -1.0, +1.0, BToPiPiLeptonNeutrino::kinematics_description_z }),

            make_signal_pdf("B->Dmunu::dGamma",
                    &BToDLeptonNeutrino::differential_branching_ratio,
                    Options{ },
                    KinematicRange{ "s", 0.0, 11.62, BToDLeptonNeutrino::kinematics_description_s }),

            make_signal_pdf("B->Dmunu::d^2Gamma",
                    &BToDLeptonInclusiveNeutrinos::normalized_differential_decay_width_1nu,
                    Options{ },
                    KinematicRange{ "s", 0.0, 19.71, BToDLeptonInclusiveNeutrinos::kinematics_description_s },
                    KinematicRange{ "cos(theta)", -1.0, +1.0, BToDLeptonInclusiveNeutrinos::kinematics_description_c_theta}),

            make_signal_pdf("B->Dmu3nu::d^5Gamma",
                    &BToDLeptonInclusiveNeutrinos::normalized_differential_decay_width_3nu,
                    Options{ },
                    KinematicRange{ "s", 3.16, 19.71, BToDLeptonInclusiveNeutrinos::kinematics_description_s },
                    KinematicRange{ "snunubar", 0.0, 3.16, BToDLeptonInclusiveNeutrinos::kinematics_description_snunubar },
                    KinematicRange{ "cos(theta_tau)", -1.0, +1.0, BToDLeptonInclusiveNeutrinos::kinematics_description_c_theta_tau },
                    KinematicRange{ "phi", 0.0, 2.0 * M_PI, BToDLeptonInclusiveNeutrinos::kinematics_description_phi },
                    KinematicRange{ "cos(theta_mu^*)", -1.0, +1.0, BToDLeptonInclusiveNeutrinos::kinematics_description_c_theta_mu_star }),

            make_signal_pdf("B->pimunu::d^2Gamma",
                    &BToPiLeptonInclusiveNeutrinos::normalized_differential_decay_width_1nu,
                    Options{ },
                    KinematicRange{ "s", 0.0, 26.41, BToPiLeptonInclusiveNeutrinos::kinematics_description_s},
                    KinematicRange{ "cos(theta)", -1.0, +1.0, BToPiLeptonInclusiveNeutrinos::kinematics_description_c_theta}),

            make_signal_pdf("B->pimu3nu::d^5Gamma",
                    &BToPiLeptonInclusiveNeutrinos::normalized_differential_decay_width_3nu,
                    Options{ },
                    KinematicRange{ "s", 3.16, 26.41, BToPiLeptonInclusiveNeutrinos::kinematics_description_s },
                    KinematicRange{ "snunubar", 0.0, 3.16, BToPiLeptonInclusiveNeutrinos::kinematics_description_snunubar },
                    KinematicRange{ "cos(theta_tau)", -1.0, +1.0, BToPiLeptonInclusiveNeutrinos::kinematics_description_c_theta_tau },
                    KinematicRange{ "phi", 0.0, 2.0 * M_PI, BToPiLeptonInclusiveNeutrinos::kinematics_description_phi },
                    KinematicRange{ "cos(theta_mu^*)", -1.0, +1.0, BToPiLeptonInclusiveNeutrinos::kinematics_description_c_theta_mu_star }),

            make_signal_pdf("B_s->K^*lnu::d^4Gamma",
                    &BsToKstarLeptonNeutrino::four_differential_decay_width,
                    Options{ },
                    KinematicRange{ "s", 0.02, 19.71, BsToKstarLeptonNeutrino::kinematics_description_s },
                    KinematicRange{ "cos(theta_l)", -1.0, +1.0, BsToKstarLeptonNeutrino::kinematics_description_c_theta_l },
                    KinematicRange{ "cos(theta_k)", -1.0, +1.0, BsToKstarLeptonNeutrino::kinematics_description_c_theta_k },
                    KinematicRange{ "phi", 0.0, 2.0 * M_PI, BsToKstarLeptonNeutrino::kinematics_description_phi }),

            make_signal_pdf("Lambda_b->Lambda_c(2625)lnu::dGamma",
                    &LambdaBToLambdaC2625LeptonNeutrino::differential_branching_ratio,
                    Options{ },
                    KinematicRange{ "s", 0.011, 8.9478, LambdaBToLambdaC2625LeptonNeutrino::kinematics_description_s }),

            make_signal_pdf("Lambda_b->Lambda_c(2625)lnu::d^2Gamma",
                    &LambdaBToLambdaC2625LeptonNeutrino::double_differential_branching_ratio,
                    Options{ },
                    KinematicRange{ "s", 0.011, 8.9478, LambdaBToLambdaC2625LeptonNeutrino::kinematics_description_s },
                    KinematicRange{ "cos(theta_l)", -1.0, +1.0, LambdaBToLambdaC2625LeptonNeutrino::kinematics_description_c_theta_l }),

            /* Exclusive Rare B Decays */

            make_signal_pdf("B->Kll::d^2Gamma@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::two_differential_decay_width,
                    Options{ },
                    KinematicRange{ "s", 1.00, 6.00, BToKDilepton<LargeRecoil>::kinematics_description_s },
                    KinematicRange{ "cos(theta_l)", -1.0, +1.0, BToKDilepton<LargeRecoil>::kinematics_description_c_theta_l }),

            make_signal_pdf("B->Kll::d^2Gamma@LowRecoil",
                    &BToKDilepton<LowRecoil>::two_differential_decay_width,
                    Options{ },
                    KinematicRange{ "s", 15.00, 22.87, BToKDilepton<LowRecoil>::kinematics_description_s },
                    KinematicRange{ "cos(theta_l)", -1.0, +1.0, BToKDilepton<LowRecoil>::kinematics_description_c_theta_l }),

            make_signal_pdf("B->K^*ll::d^4Gamma@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::four_differential_decay_width,
                    Options{ },
                    KinematicRange{ "s", 1.00, 6.00, BToKstarDilepton<LargeRecoil>::kinematics_description_s },
                    KinematicRange{ "cos(theta_l)", -1.0, +1.0, BToKstarDilepton<LargeRecoil>::kinematics_description_c_theta_l },
                    KinematicRange{ "cos(theta_k)", -1.0, +1.0, BToKstarDilepton<LargeRecoil>::kinematics_description_c_theta_k },
                    KinematicRange{ "phi", 0.0, 2.0 * M_PI, BToKstarDilepton<LargeRecoil>::kinematics_description_phi }),

            make_signal_pdf("B->K^*ll::d^4Gamma@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::four_differential_decay_width,
                    Options{ },
                    KinematicRange{ "s", 15.00, 19.21, BToKstarDilepton<LowRecoil>::kinematics_description_s },
                    KinematicRange{ "cos(theta_l)", -1.0, +1.0, BToKstarDilepton<LowRecoil>::kinematics_description_c_theta_l },
                    KinematicRange{ "cos(theta_k)", -1.0, +1.0, BToKstarDilepton<LowRecoil>::kinematics_description_c_theta_k },
                    KinematicRange{ "phi", 0.0, 2.0 * M_PI, BToKstarDilepton<LowRecoil>::kinematics_description_phi }),
        };

        return signal_pdf_entries;
    }

    SignalPDFPtr
    SignalPDF::make(const QualifiedName & name, const Parameters & parameters, const Kinematics & kinematics, const Options & _options)
    {
        static const std::map<QualifiedName, const SignalPDFEntry *> signal_pdf_entries = make_signal_pdf_entries();

        // check if 'name' matches any of the implemented Signal PDFs
        {
            auto i = signal_pdf_entries.find(name);
            if (signal_pdf_entries.end() != i)
                return i->second->make(parameters, kinematics, name.options() + _options);
        }

        return SignalPDFPtr();
    }

    template <>
    struct WrappedForwardIteratorTraits<SignalPDFs::SignalPDFIteratorTag>
    {
        typedef std::map<QualifiedName, const SignalPDFEntry *>::iterator UnderlyingIterator;
    };
    template class WrappedForwardIterator<SignalPDFs::SignalPDFIteratorTag, std::pair<const QualifiedName, const SignalPDFEntry *>>;

    template<>
    struct Implementation<SignalPDFs>
    {
        std::map<QualifiedName, const SignalPDFEntry *> signal_pdf_entries;

        Implementation() :
            signal_pdf_entries(make_signal_pdf_entries())
        {
        }
    };

    SignalPDFs::SignalPDFs() :
        PrivateImplementationPattern<SignalPDFs>(new Implementation<SignalPDFs>())
    {
    }

    SignalPDFs::~SignalPDFs()
    {
    }

    SignalPDFs::SignalPDFIterator
    SignalPDFs::begin() const
    {
        return SignalPDFIterator(_imp->signal_pdf_entries.begin());
    }

    SignalPDFs::SignalPDFIterator
    SignalPDFs::end() const
    {
        return SignalPDFIterator(_imp->signal_pdf_entries.end());
    }
}
