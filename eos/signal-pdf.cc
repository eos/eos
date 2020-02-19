/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015, 2016, 2017 Danny van Dyk
 * Copyright (c) 2019 Ahmet Kokulu
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
#include <eos/b-decays/b-to-d-l-x-nu.hh>
#include <eos/b-decays/b-to-dstar-l-nu.hh>
#include <eos/b-decays/b-to-psd-l-nu.hh>
#include <eos/b-decays/b-to-pi-l-x-nu.hh>
#include <eos/b-decays/b-to-pi-pi-l-nu.hh>
#include <eos/b-decays/bs-to-kstar-l-nu.hh>
#include <eos/b-decays/lambdab-to-lambdac-l-nu.hh>
#include <eos/b-decays/lambdab-to-lambdac2625-l-nu.hh>
#include <eos/rare-b-decays/exclusive-b-to-s-dilepton-large-recoil.hh>
#include <eos/rare-b-decays/exclusive-b-to-s-dilepton-low-recoil.hh>
#include <eos/utils/density.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <functional>
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

                double pdf(const double & z) const
                {
                    return (9.0 + 8.0 * z + 9.0 * z * z);
                }

                double norm(const double & z_min, const double & z_max) const
                {
                    using std::pow;

                    return (9.0 * (z_max - z_min) + 4.0 * (pow(z_max, 2.0) - pow(z_min, 2.0)) + 3.0 * (pow(z_max, 3.0) - pow(z_min, 3.0)));
                }

                static const std::string description;
        };

        const std::string Legendre1DPDF::description = "1D PDF up to 2nd order in z; used for unit tests only.";
    }
}

#include <eos/utils/concrete-signal-pdf.hh>

namespace eos
{
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

    template <typename Decay_, typename ... PDFArgs_, typename ... PDFKinematicRanges_, typename ... NormArgs_,  typename ... NormKinematicNames_>
    std::pair<QualifiedName, SignalPDFEntry *> make_signal_pdf(const char * name,
            const Options & default_options,
            double (Decay_::* pdf)(const PDFArgs_ & ...) const,
            const std::tuple<PDFKinematicRanges_ ...> & pdf_kinematic_ranges,
            double (Decay_::* norm)(const NormArgs_ & ...) const,
            const std::tuple<NormKinematicNames_ ...> & norm_kinematic_names)
    {
        static_assert(sizeof...(PDFArgs_) == sizeof...(PDFKinematicRanges_), "Need as many function arguments for the PDF as kinematics ranges!");
        static_assert(sizeof...(NormArgs_) == sizeof...(NormKinematicNames_), "Need as many function arguments for the normalization as kinematics names!");

        QualifiedName qn(name);

        return std::make_pair(qn, make_concrete_signal_pdf_entry(qn, default_options, pdf, pdf_kinematic_ranges, norm, norm_kinematic_names));
    }

    const std::map<QualifiedName, const SignalPDFEntry *> &
    make_signal_pdf_entries()
    {
        static const std::map<QualifiedName, const SignalPDFEntry *> signal_pdf_entries
        {
            /* Internal Tests */

            make_signal_pdf("Test::Legendre1D",
                    Options{ },
                    &test::Legendre1DPDF::pdf,
                    std::make_tuple(
                        KinematicRange{ "z", -1.0, +1.0, "" }
                    ),
                    &test::Legendre1DPDF::norm,
                    std::make_tuple(
                        "z_min",
                        "z_max"
                    )
                ),

            /* Exclusive Decays */

            /* Exclusive B Decays */

            make_signal_pdf("B->pipimunu::d^3Gamma@QCDF",
                    Options{ },
                    &BToPiPiLeptonNeutrino::triple_differential_branching_ratio,
                    std::make_tuple(
                        KinematicRange{ "q2", 0.01, 0.93859, BToPiPiLeptonNeutrino::kinematics_description_q2 },
                        KinematicRange{ "k2", 18.582, 27.872, BToPiPiLeptonNeutrino::kinematics_description_k2 },
                        KinematicRange{ "cos(theta)", -1.0, +1.0, BToPiPiLeptonNeutrino::kinematics_description_z }
                    ),
                    &BToPiPiLeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "q2_min",
                        "q2_max",
                        "k2_min",
                        "k2_max",
                        "cos(theta)_min",
                        "cos(theta)_max"
                    )
                ),

            make_signal_pdf("B->pilnu::dGamma/dq2",
                    Options{ { "U", "u" } },
                    &BToPseudoscalarLeptonNeutrino::differential_branching_ratio,
                    std::make_tuple(
                        KinematicRange{ "q2", 0.0, 26.41, BToPseudoscalarLeptonNeutrino::kinematics_description_q2 }
                    ),
                    &BToPseudoscalarLeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "q2_min",
                        "q2_max"
                    )
                ),

            make_signal_pdf("B->pilnu::d^2Gamma/dq2/dcos(theta_l)",
                    Options{ {"U", "u"} },
                    &BToPseudoscalarLeptonNeutrino::normalized_two_differential_decay_width,
                    std::make_tuple(
                        KinematicRange{ "q2", 0.0, 26.41, BToPseudoscalarLeptonNeutrino::kinematics_description_q2 },
                        KinematicRange{ "cos(theta_l)", -1.0, +1.0, BToPseudoscalarLeptonNeutrino::kinematics_description_c_theta_l}
                    ),
                    &BToPseudoscalarLeptonNeutrino::normalized_integrated_decay_width,
                    std::make_tuple(
                        "q2_min",
                        "q2_max"
                    )
                ),

            make_signal_pdf("B->Dlnu::dGamma/dq2",
                    Options{ { "U", "c" } },
                    &BToPseudoscalarLeptonNeutrino::differential_branching_ratio,
                    std::make_tuple(
                        KinematicRange{ "q2", 0.0, 11.62, BToPseudoscalarLeptonNeutrino::kinematics_description_q2 }
                    ),
                    &BToPseudoscalarLeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "q2_min",
                        "q2_max"
                    )
                ),

            make_signal_pdf("B->Dlnu::d^2Gamma/dq2/dcos(theta_l)",
                    Options{ {"U", "c"} },
                    &BToPseudoscalarLeptonNeutrino::normalized_two_differential_decay_width,
                    std::make_tuple(
                        KinematicRange{ "q2", 0.0, 11.62, BToPseudoscalarLeptonNeutrino::kinematics_description_q2 },
                        KinematicRange{ "cos(theta_l)", -1.0, +1.0, BToPseudoscalarLeptonNeutrino::kinematics_description_c_theta_l}
                    ),
                    &BToPseudoscalarLeptonNeutrino::normalized_integrated_decay_width,
                    std::make_tuple(
                        "q2_min",
                        "q2_max"
                    )
                ),

            make_signal_pdf("B->D^*munu::dBR",
                    Options{ {"l", "mu"} },
                    &BToDstarLeptonNeutrino::differential_branching_ratio,
                    std::make_tuple(
                        KinematicRange{ "q2", 0.0, 10.68, BToDstarLeptonNeutrino::kinematics_description_q2 }
                    ),
                    &BToDstarLeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "q2_min",
                        "q2_max"
                    )
                ),

            make_signal_pdf("B->D^*munu::d^4Gamma",
                    Options{ {"l", "mu"} },
                    &BToDstarLeptonNeutrino::normalized_four_differential_decay_width,
                    std::make_tuple(
                        KinematicRange{ "q2", 0.0, 10.68, BToDstarLeptonNeutrino::kinematics_description_q2 },
                        KinematicRange{ "cos(theta_l)", -1.0, +1.0, BToDstarLeptonNeutrino::kinematics_description_c_theta_l },
                        KinematicRange{ "cos(theta_d)", -1.0, +1.0, BToDstarLeptonNeutrino::kinematics_description_c_theta_d },
                        KinematicRange{ "phi", 0.0, 2.0 * M_PI, BToDstarLeptonNeutrino::kinematics_description_phi }
                    ),
                    &BToDstarLeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "q2_min",
                        "q2_max"
                    )
                ),

            make_signal_pdf("B->Dmu1nu::d^2Gamma",
                    Options{ },
                    &BToDLeptonInclusiveNeutrinos::differential_decay_width_1nu,
                    std::make_tuple(
                        KinematicRange{ "s", 0.0, 19.71, BToDLeptonInclusiveNeutrinos::kinematics_description_s },
                        KinematicRange{ "cos(theta)", -1.0, +1.0, BToDLeptonInclusiveNeutrinos::kinematics_description_c_theta}
                    ),
                    &BToDLeptonInclusiveNeutrinos::integrated_decay_width_1nu,
                    std::make_tuple(
                        "s_min",
                        "s_max"
                    )
                ),

            make_signal_pdf("B->Dmu3nu::d^5Gamma",
                    Options{ },
                    &BToDLeptonInclusiveNeutrinos::differential_decay_width_3nu,
                    std::make_tuple(
                        KinematicRange{ "s", 3.16, 19.71, BToDLeptonInclusiveNeutrinos::kinematics_description_s },
                        KinematicRange{ "snunubar", 0.0, 3.16, BToDLeptonInclusiveNeutrinos::kinematics_description_snunubar },
                        KinematicRange{ "cos(theta_tau)", -1.0, +1.0, BToDLeptonInclusiveNeutrinos::kinematics_description_c_theta_tau },
                        KinematicRange{ "phi", 0.0, 2.0 * M_PI, BToDLeptonInclusiveNeutrinos::kinematics_description_phi },
                        KinematicRange{ "cos(theta_mu^*)", -1.0, +1.0, BToDLeptonInclusiveNeutrinos::kinematics_description_c_theta_mu_star }
                    ),
                    &BToDLeptonInclusiveNeutrinos::integrated_decay_width_3nu,
                    std::make_tuple(
                        "s_min",
                        "s_max"
                    )
                ),

            make_signal_pdf("B->pimu1nu::d^2Gamma",
                    Options{ },
                    &BToPiLeptonInclusiveNeutrinos::differential_decay_width_1nu,
                    std::make_tuple(
                        KinematicRange{ "s", 0.0, 26.41, BToPiLeptonInclusiveNeutrinos::kinematics_description_s},
                        KinematicRange{ "cos(theta)", -1.0, +1.0, BToPiLeptonInclusiveNeutrinos::kinematics_description_c_theta}
                    ),
                    &BToPiLeptonInclusiveNeutrinos::integrated_decay_width_1nu,
                    std::make_tuple(
                        "s_min",
                        "s_max"
                    )
                ),

            make_signal_pdf("B->pimu3nu::d^5Gamma",
                    Options{ },
                    &BToPiLeptonInclusiveNeutrinos::differential_decay_width_3nu,
                    std::make_tuple(
                        KinematicRange{ "s", 3.16, 26.41, BToPiLeptonInclusiveNeutrinos::kinematics_description_s },
                        KinematicRange{ "snunubar", 0.0, 3.16, BToPiLeptonInclusiveNeutrinos::kinematics_description_snunubar },
                        KinematicRange{ "cos(theta_tau)", -1.0, +1.0, BToPiLeptonInclusiveNeutrinos::kinematics_description_c_theta_tau },
                        KinematicRange{ "phi", 0.0, 2.0 * M_PI, BToPiLeptonInclusiveNeutrinos::kinematics_description_phi },
                        KinematicRange{ "cos(theta_mu^*)", -1.0, +1.0, BToPiLeptonInclusiveNeutrinos::kinematics_description_c_theta_mu_star }
                    ),
                    &BToPiLeptonInclusiveNeutrinos::integrated_decay_width_3nu,
                    std::make_tuple(
                        "s_min",
                        "s_max"
                    )
                ),

            make_signal_pdf("B_s->K^*lnu::d^4Gamma",
                    Options{ },
                    &BsToKstarLeptonNeutrino::four_differential_decay_width,
                    std::make_tuple(
                        KinematicRange{ "s", 0.02, 19.71, BsToKstarLeptonNeutrino::kinematics_description_s },
                        KinematicRange{ "cos(theta_l)", -1.0, +1.0, BsToKstarLeptonNeutrino::kinematics_description_c_theta_l },
                        KinematicRange{ "cos(theta_k)", -1.0, +1.0, BsToKstarLeptonNeutrino::kinematics_description_c_theta_k },
                        KinematicRange{ "phi", 0.0, 2.0 * M_PI, BsToKstarLeptonNeutrino::kinematics_description_phi }
                    ),
                    &BsToKstarLeptonNeutrino::integrated_decay_width,
                    std::make_tuple(
                        "s_min",
                        "s_max"
                    )
                ),

            make_signal_pdf("Lambda_b->Lambda_clnu::dGamma",
                    Options{ },
                    &LambdaBToLambdaCLeptonNeutrino::differential_branching_ratio,
                    std::make_tuple(
                        KinematicRange{ "s", 0.011, 11.1, LambdaBToLambdaCLeptonNeutrino::kinematics_description_q2 }
                    ),
                    &LambdaBToLambdaCLeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "s_min",
                        "s_max"
                    )
                ),

            make_signal_pdf("Lambda_b->Lambda_c(2625)lnu::dGamma",
                    Options{ },
                    &LambdaBToLambdaC2625LeptonNeutrino::differential_branching_ratio,
                    std::make_tuple(
                        KinematicRange{ "s", 0.011, 8.9478, LambdaBToLambdaC2625LeptonNeutrino::kinematics_description_s }
                    ),
                    &LambdaBToLambdaC2625LeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "s_min",
                        "s_max"
                    )
                ),

            make_signal_pdf("Lambda_b->Lambda_c(2625)lnu::d^2Gamma",
                    Options{ },
                    &LambdaBToLambdaC2625LeptonNeutrino::double_differential_branching_ratio,
                    std::make_tuple(
                        KinematicRange{ "s", 0.011, 8.9478, LambdaBToLambdaC2625LeptonNeutrino::kinematics_description_s },
                        KinematicRange{ "cos(theta_l)", -1.0, +1.0, LambdaBToLambdaC2625LeptonNeutrino::kinematics_description_c_theta_l }
                    ),
                    &LambdaBToLambdaC2625LeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple(
                        "s_min",
                        "s_max"
                    )
                ),

            /* Exclusive Rare B Decays */

            make_signal_pdf("B->Kll::d^2Gamma@LargeRecoil",
                    Options{ },
                    &BToKDilepton<LargeRecoil>::two_differential_decay_width,
                    std::make_tuple(
                        KinematicRange{ "s",                  1.00,  6.00, BToKDilepton<LargeRecoil>::kinematics_description_s },
                        KinematicRange{ "cos(theta_l)^LHCb", -1.0,  +1.0,  BToKDilepton<LargeRecoil>::kinematics_description_c_theta_l }
                    ),
                    &BToKDilepton<LargeRecoil>::integrated_decay_width,
                    std::make_tuple(
                        "s_min",
                        "s_max"
                    )
                ),

            make_signal_pdf("B->Kll::d^2Gamma@LowRecoil",
                    Options{ },
                    &BToKDilepton<LowRecoil>::two_differential_decay_width,
                    std::make_tuple(
                        KinematicRange{ "s",                 15.00, 22.87, BToKDilepton<LowRecoil>::kinematics_description_s },
                        KinematicRange{ "cos(theta_l)^LHCb", -1.0,  +1.0,  BToKDilepton<LowRecoil>::kinematics_description_c_theta_l }
                    ),
                    &BToKDilepton<LowRecoil>::integrated_decay_width,
                    std::make_tuple(
                        "s_min",
                        "s_max"
                    )
                ),

            make_signal_pdf("B->K^*ll::d^4Gamma@LargeRecoil",
                    Options{ },
                    &BToKstarDilepton<LargeRecoil>::four_differential_decay_width,
                    std::make_tuple(
                        KinematicRange{ "s",                  1.00,  6.00,       BToKstarDilepton<LargeRecoil>::kinematics_description_s         },
                        KinematicRange{ "cos(theta_l)^LHCb", -1.0,  +1.0,        BToKstarDilepton<LargeRecoil>::kinematics_description_c_theta_l },
                        KinematicRange{ "cos(theta_k)^LHCb", -1.0,  +1.0,        BToKstarDilepton<LargeRecoil>::kinematics_description_c_theta_k },
                        KinematicRange{ "phi^LHCb",           0.0,   2.0 * M_PI, BToKstarDilepton<LargeRecoil>::kinematics_description_phi       }
                    ),
                    &BToKstarDilepton<LargeRecoil>::integrated_decay_width,
                    std::make_tuple(
                        "s_min",
                        "s_max"
                    )
                ),

            make_signal_pdf("B->K^*ll::d^4Gamma@LowRecoil",
                    Options{ },
                    &BToKstarDilepton<LowRecoil>::four_differential_decay_width,
                    std::make_tuple(
                        KinematicRange{ "s",                  15.00,  19.21,      BToKstarDilepton<LowRecoil>::kinematics_description_s },
                        KinematicRange{ "cos(theta_l)^LHCb", -1.0,   +1.0,        BToKstarDilepton<LowRecoil>::kinematics_description_c_theta_l },
                        KinematicRange{ "cos(theta_k)^LHCb", -1.0,   +1.0,        BToKstarDilepton<LowRecoil>::kinematics_description_c_theta_k },
                        KinematicRange{ "phi^LHCb",           0.0,    2.0 * M_PI, BToKstarDilepton<LowRecoil>::kinematics_description_phi }
                    ),
                    &BToKstarDilepton<LowRecoil>::integrated_decay_width,
                    std::make_tuple(
                        "s_min",
                        "s_max"
                    )
                ),

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
