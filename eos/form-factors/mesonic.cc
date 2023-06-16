/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013, 2014, 2015, 2016, 2018 Danny van Dyk
 * Copyright (c) 2015 Christoph Bobeth
 * Copyright (c) 2018 Ahmet Kokulu
 * Copyright (c) 2019 Nico Gubernari
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

#include <eos/form-factors/analytic-b-to-gamma-qcdf.hh>
#include <eos/form-factors/analytic-b-to-psd-dkmmo2008.hh>
#include <eos/form-factors/analytic-b-to-pi-pi.hh>
#include <eos/form-factors/analytic-b-to-p-lcsr.hh>
#include <eos/form-factors/analytic-b-to-v-lcsr.hh>
#include <eos/form-factors/form-factors.hh>
#include <eos/form-factors/parametric-bcl2008.hh>
#include <eos/form-factors/parametric-bfw2010.hh>
#include <eos/form-factors/parametric-bgl1997.hh>
#include <eos/form-factors/parametric-bgjvd2019.hh>
#include <eos/form-factors/parametric-bsz2015.hh>
#include <eos/form-factors/parametric-fvdv2018.hh>
#include <eos/form-factors/parametric-kkvdz2022.hh>
#include <eos/form-factors/parametric-kmpw2010.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/qualified-name.hh>

#include <cmath>
#include <limits>
#include <map>

namespace eos
{
    /* P -> V Processes */

    FormFactors<PToV>::~FormFactors()
    {
    }

    const std::map<FormFactorFactory<PToV>::KeyType, FormFactorFactory<PToV>::ValueType>
    FormFactorFactory<PToV>::form_factors
    {
        { "B->omega::BSZ2015",    &BSZ2015FormFactors<BToOmega,   PToV>::make         },
        { "B->rho::BSZ2015",      &BSZ2015FormFactors<BToRho,     PToV>::make         },
        { "B->K^*::KMPW2010",     &KMPW2010FormFactors<PToV>::make                    },
        { "B->K^*::BSZ2015",      &BSZ2015FormFactors<BToKstar,   PToV>::make         },
        { "B->K^*::BFW2010",      &BFW2010FormFactors<BToKstar,   PToV>::make         },
        { "B->D^*::BSZ2015",      &BSZ2015FormFactors<BToDstar,   PToV>::make         },
        { "B->D^*::BGJvD2019",    &HQETFormFactors<BToDstar,      PToV>::make         },
        { "B->D^*::BGL1997",      &BGL1997FormFactors<BToDstar>::make                 },
        { "B_s->K^*::BSZ2015",    &BSZ2015FormFactors<BsToKstar,  PToV>::make         },
        { "B_s->D_s^*::BSZ2015",  &BSZ2015FormFactors<BsToDsstar, PToV>::make         },
        { "B_s->D_s^*::BGJvD2019",&HQETFormFactors<BsToDsstar,    PToV>::make         },
        { "B_s->phi::BSZ2015",    &BSZ2015FormFactors<BsToPhi,    PToV>::make         },
        { "B_s->phi::BFW2010",    &BFW2010FormFactors<BsToPhi,    PToV>::make         },
        // analytic computations
        { "B->K^*::B-LCSR",       &AnalyticFormFactorBToVLCSR<lcsr::BToKstar>::make   },
        { "B->D^*::B-LCSR",       &AnalyticFormFactorBToVLCSR<lcsr::BToDstar>::make   },
        { "B->rho::B-LCSR",       &AnalyticFormFactorBToVLCSR<lcsr::BToRho>::make     },
        { "B_s->K^*::B-LCSR",     &AnalyticFormFactorBToVLCSR<lcsr::BsToKstar>::make  },
        { "B_s->phi::B-LCSR",     &AnalyticFormFactorBToVLCSR<lcsr::BsToPhi>::make    },
        { "B_s->D_s^*::B-LCSR",   &AnalyticFormFactorBToVLCSR<lcsr::BsToDsstar>::make }
    };

    complex<double>
    FormFactors<PToV>::v(const complex<double> &) const
    {
        throw InternalError("P->V form factor V for complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    complex<double>
    FormFactors<PToV>::a_0(const complex<double> &) const
    {
        throw InternalError("P->V form factor A_0 for complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    complex<double>
    FormFactors<PToV>::a_1(const complex<double> &) const
    {
        throw InternalError("P->V form factors A_1 complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    complex<double>
    FormFactors<PToV>::a_12(const complex<double> &) const
    {
        throw InternalError("P->V form factor A_12 for complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    complex<double>
    FormFactors<PToV>::a_2(const complex<double> &) const
    {
        throw InternalError("P->V form factor A_2 for complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    complex<double>
    FormFactors<PToV>::t_1(const complex<double> &) const
    {
        throw InternalError("P->V form factor T_1 for complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    complex<double>
    FormFactors<PToV>::t_2(const complex<double> &) const
    {
        throw InternalError("P->V form factor T_2 for complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    complex<double>
    FormFactors<PToV>::t_23(const complex<double> &) const
    {
        throw InternalError("P->V form factor T_23 for complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    std::shared_ptr<FormFactors<PToV>>
    FormFactorFactory<PToV>::create(const QualifiedName & name, const Parameters & parameters, const Options & options)
    {
        Context ctx("When creating a P->V form factor");

        std::shared_ptr<FormFactors<PToV>> result;

        auto i = form_factors.find(name);
        if (form_factors.end() != i)
        {
            result.reset(i->second(parameters, name.options() + options));
            return result;
        }

        throw NoSuchFormFactorError(name.prefix_part().str(), name.name_part().str());
        return result;
    }

    OptionSpecification
    FormFactorFactory<PToV>::option_specification(const qnp::Prefix & process)
    {
        OptionSpecification result { "form-factors", {}, "" };
        for (const auto & ff : FormFactorFactory<PToV>::form_factors)
        {
            if (process == std::get<0>(ff).prefix_part())
                result.allowed_values.push_back(std::get<0>(ff).name_part().str());
        }

        return result;
    }

    OptionSpecification
    FormFactorFactory<PToV>::option_specification()
    {
        std::set<std::string> allowed_values;
        for (const auto & ff : FormFactorFactory<PToV>::form_factors)
        {
            allowed_values.insert(std::get<0>(ff).name_part().str());
        }

        OptionSpecification result { "form-factors", { allowed_values.cbegin(), allowed_values.cend() }, "" };
        return result;
    }

    /* B_{u,d} -> omega */

    /* B_{u,d} -> K^* */

    /* B_s -> phi */

    /* P -> gamma Processes */

    FormFactors<PToGamma>::~FormFactors()
    {
    }

    const std::map<FormFactorFactory<PToGamma>::KeyType, FormFactorFactory<PToGamma>::ValueType>
    FormFactorFactory<PToGamma>::form_factors
    {
        { KeyType("B->gamma::FLvD2022QCDF"), &AnalyticFormFactorBToGammaQCDF::make }
    };

    std::shared_ptr<FormFactors<PToGamma>>
    FormFactorFactory<PToGamma>::create(const QualifiedName & name, const Parameters & parameters, const Options & options)
    {
        Context ctx("When creating a P->gamma form factor");

        std::shared_ptr<FormFactors<PToGamma>> result;

        auto & form_factors = FormFactorFactory<PToGamma>::form_factors;
        auto i = form_factors.find(name);
        if (form_factors.end() != i)
        {
            result.reset(i->second(parameters, name.options() + options));
            return result;
        }

        throw NoSuchFormFactorError(name.prefix_part().str(), name.name_part().str());
        return result;
    }

    OptionSpecification
    FormFactorFactory<PToGamma>::option_specification(const qnp::Prefix & process)
    {
        OptionSpecification result { "form-factors", {}, "" };
        for (const auto & ff : FormFactorFactory<PToGamma>::form_factors)
        {
            if (process == std::get<0>(ff).prefix_part())
                result.allowed_values.push_back(std::get<0>(ff).name_part().str());
        }

        return result;
    }

    /* P -> gamma^* Processes */

    FormFactors<PToGammaOffShell>::~FormFactors() = default;

    const std::map<FormFactorFactory<PToGammaOffShell>::KeyType, FormFactorFactory<PToGammaOffShell>::ValueType>
    FormFactorFactory<PToGammaOffShell>::form_factors
    {
        { KeyType("B->gamma^*::KKvDZ2022"), &KKvDZ2022FormFactors::make }
    };

    std::shared_ptr<FormFactors<PToGammaOffShell>>
    FormFactorFactory<PToGammaOffShell>::create(const QualifiedName & name, const Parameters & parameters, const Options & options)
    {
        Context ctx("When creating a P->gamma^* form factor");

        std::shared_ptr<FormFactors<PToGammaOffShell>> result;

        auto & form_factors = FormFactorFactory<PToGammaOffShell>::form_factors;
        auto i = form_factors.find(name);
        if (form_factors.end() != i)
        {
            result.reset(i->second(parameters, name.options() + options));
            return result;
        }

        throw NoSuchFormFactorError(name.prefix_part().str(), name.name_part().str());
        return result;
    }

    OptionSpecification
    FormFactorFactory<PToGammaOffShell>::option_specification(const qnp::Prefix & process)
    {
        OptionSpecification result { "form-factors", {}, "" };
        for (const auto & ff : FormFactorFactory<PToGammaOffShell>::form_factors)
        {
            if (process == std::get<0>(ff).prefix_part())
                result.allowed_values.push_back(std::get<0>(ff).name_part().str());
        }

        return result;
    }

    /* P -> P Processes */

    /* B_{u,d} -> pi */

    /* B_{u,d} -> D */

    FormFactors<PToP>::~FormFactors() = default;

    double FormFactors<PToP>::f_m(const double & /*s*/) const
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double FormFactors<PToP>::f_p_d1(const double & s) const
    {
        using namespace std::placeholders;

        std::function<double (const double &)> f = [&, this](const double & q2) -> double { return this->f_p(q2); };

        return derivative<1u, deriv::TwoSided>(f, s);
    }

    double FormFactors<PToP>::f_p_d2(const double & s) const
    {
        using namespace std::placeholders;

        std::function<double (const double &)> f = [&, this](const double & q2) -> double { return this->f_p(q2); };

        return derivative<2u, deriv::TwoSided>(f, s);
    }

    const std::map<FormFactorFactory<PToP>::KeyType, FormFactorFactory<PToP>::ValueType>
    FormFactorFactory<PToP>::form_factors
    {
        // parametrizations
        // b -> s
        { "B->K::BCL2008",       &BCL2008FormFactors<BToK, 3u>::make                                                                                          },
        { "B->K::KMPW2010",      &KMPW2010FormFactors<PToP>::make                                                                                             },
        { "B->K::BSZ2015",       &BSZ2015FormFactors<BToK,   PToP>::make                                                                                      },
        { "B->K::BFW2010",       &BFW2010FormFactors<BToK,   PToP>::make                                                                                      },
        // b -> u
        { "B->pi::BCL2008",      &BCL2008FormFactors<BToPi, 3u>::make                                                                                         },
        { "B->pi::BCL2008-4",    &BCL2008FormFactors<BToPi, 4u>::make                                                                                         },
        { "B->pi::BCL2008-5",    &BCL2008FormFactors<BToPi, 5u>::make                                                                                         },
        { "B->pi::BSZ2015",      &BSZ2015FormFactors<BToPi,  PToP>::make                                                                                      },
        { "B_s->K::BFW2010",     &BFW2010FormFactors<BsToK,  PToP>::make                                                                                      },
        { "B_s->K::BSZ2015",     &BSZ2015FormFactors<BsToK,  PToP>::make                                                                                      },
        // b -> c
        { "B->D::BCL2008",       &BCL2008FormFactors<BToD, 3u>::make                                                                                          },
        { "B->D::BSZ2015",       &BSZ2015FormFactors<BToD,   PToP>::make                                                                                      },
        { "B->D::BGJvD2019",     &HQETFormFactors<BToD,      PToP>::make                                                                                      },
        { "B->D::BGL1997",       &BGL1997FormFactors<BToD>::make                                                                                              },
        { "B_s->D_s::BSZ2015",   &BSZ2015FormFactors<BsToDs, PToP>::make                                                                                      },
        { "B_s->D_s::BGJvD2019", &HQETFormFactors<BsToDs,    PToP>::make                                                                                      },
        // c -> d
        { "D->pi::BSZ2015",      &BSZ2015FormFactors<DToPi,  PToP>::make                                                                                      },
        { "D_s->K::BSZ2015",     &BSZ2015FormFactors<DsToK,  PToP>::make                                                                                      },
        // c -> s
        { "D->K::BSZ2015",       &BSZ2015FormFactors<DToK,   PToP>::make                                                                                      },
        // analytic computations
        { "B->pi::DKMMO2008",    &AnalyticFormFactorBToPseudoscalarDKMMO2008<QuarkFlavor::bottom, QuarkFlavor::up, QuarkFlavor::down>::make                   },
        { "B_s->K::DKMMO2008",   &AnalyticFormFactorBToPseudoscalarDKMMO2008<QuarkFlavor::bottom, QuarkFlavor::up, QuarkFlavor::strange>::make                },
        { "B->pi::B-LCSR",       &AnalyticFormFactorBToPLCSR<lcsr::BToPi>::make                                                                               },
        { "B->K::B-LCSR",        &AnalyticFormFactorBToPLCSR<lcsr::BToK>::make                                                                                },
        { "B->D::B-LCSR",        &AnalyticFormFactorBToPLCSR<lcsr::BToD>::make                                                                                },
        { "B_s->K::B-LCSR",      &AnalyticFormFactorBToPLCSR<lcsr::BsToK>::make                                                                               },
        { "B_s->D_s::B-LCSR",    &AnalyticFormFactorBToPLCSR<lcsr::BsToDs>::make                                                                              }
    };

    complex<double>
    FormFactors<PToP>::f_p(const complex<double> &) const
    {
        throw InternalError("P->P form factor f_+ for complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    complex<double>
    FormFactors<PToP>::f_0(const complex<double> &) const
    {
        throw InternalError("P->P form factor f_0 for complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    complex<double>
    FormFactors<PToP>::f_t(const complex<double> &) const
    {
        throw InternalError("P->P form factor f_t for complex q2 is not implemented for this parametrisation");
        return complex<double>(std::numeric_limits<double>::signaling_NaN());
    }

    std::shared_ptr<FormFactors<PToP>>
    FormFactorFactory<PToP>::create(const QualifiedName & name, const Parameters & parameters, const Options & options)
    {
        Context ctx("When creating a P->P form factor");

        std::shared_ptr<FormFactors<PToP>> result;

        auto i = FormFactorFactory<PToP>::form_factors.find(name);
        if (FormFactorFactory<PToP>::form_factors.end() != i)
        {
            result.reset(i->second(parameters, name.options() + options));
            return result;
        }

        throw NoSuchFormFactorError(name.prefix_part().str(), name.name_part().str());
        return result;
    }

    OptionSpecification
    FormFactorFactory<PToP>::option_specification(const qnp::Prefix & process)
    {
        OptionSpecification result { "form-factors", {}, "" };
        for (const auto & ff : FormFactorFactory<PToP>::form_factors)
        {
            if (process == std::get<0>(ff).prefix_part())
                result.allowed_values.push_back(std::get<0>(ff).name_part().str());
        }

        return result;
    }

    OptionSpecification
    FormFactorFactory<PToP>::option_specification()
    {
        std::set<std::string> allowed_values;
        for (const auto & ff : FormFactorFactory<PToP>::form_factors)
        {
            allowed_values.insert(std::get<0>(ff).name_part().str());
        }

        OptionSpecification result { "form-factors", { allowed_values.cbegin(), allowed_values.cend() }, "" };
        return result;
    }

    /* B_{u,d} -> pi */

    /* P -> PP Processes */

    FormFactors<PToPP>::~FormFactors()
    {
    }

    const std::map<FormFactorFactory<PToPP>::KeyType, FormFactorFactory<PToPP>::ValueType>
    FormFactorFactory<PToPP>::form_factors
    {
        // analytic computations
        { "B->pipi::BFvD2016",            &AnalyticFormFactorBToPiPiBFvD2016::make   },
        { "B->pipi::FvDV2018-Dispersive", &AnalyticFormFactorBToPiPiFvDV2018::make   },
        { "B->pipi::FvDV2018",            &FvDV2018FormFactors<BToPiPi>::make        },
    };

    std::shared_ptr<FormFactors<PToPP>>
    FormFactorFactory<PToPP>::create(const QualifiedName & name, const Parameters & parameters, const Options & options)
    {
        Context ctx("When creating a P->PP form factor");

        std::shared_ptr<FormFactors<PToPP>> result;

        auto i = FormFactorFactory<PToPP>::form_factors.find(name);
        if (FormFactorFactory<PToPP>::form_factors.end() != i)
        {
            result.reset(i->second(parameters, name.options() + options));
            return result;
        }

        throw NoSuchFormFactorError(name.prefix_part().str(), name.name_part().str());
        return result;
    }

    OptionSpecification
    FormFactorFactory<PToPP>::option_specification(const qnp::Prefix & process)
    {
        OptionSpecification result { "form-factors", {}, "" };
        for (const auto & ff : FormFactorFactory<PToPP>::form_factors)
        {
            if (process == std::get<0>(ff).prefix_part())
                result.allowed_values.push_back(std::get<0>(ff).name_part().str());
        }

        return result;
    }

    /* V -> P Processes */

    FormFactors<VToP>::~FormFactors()
    {
    }

    const std::map<FormFactorFactory<VToP>::KeyType, FormFactorFactory<VToP>::ValueType>
    FormFactorFactory<VToP>::form_factors
    {
        // parametrizations
        // b -> c
        { "B^*->D::BGJvD2019",        &HQETFormFactors<BstarToD, VToP>::make          },
    };

    std::shared_ptr<FormFactors<VToP>>
    FormFactorFactory<VToP>::create(const QualifiedName & name, const Parameters & parameters, const Options & options)
    {
        Context ctx("When creating a V->P form factor");

        std::shared_ptr<FormFactors<VToP>> result;

        auto i = FormFactorFactory<VToP>::form_factors.find(name);
        if (FormFactorFactory<VToP>::form_factors.end() != i)
        {
            result.reset(i->second(parameters, name.options() + options));
            return result;
        }

        throw NoSuchFormFactorError(name.prefix_part().str(), name.name_part().str());
        return result;
    }

    OptionSpecification
    FormFactorFactory<VToP>::option_specification(const qnp::Prefix & process)
    {
        OptionSpecification result { "form-factors", {}, "" };
        for (const auto & ff : FormFactorFactory<VToP>::form_factors)
        {
            if (process == std::get<0>(ff).prefix_part())
                result.allowed_values.push_back(std::get<0>(ff).name_part().str());
        }

        return { result };
    }

    /* V -> V Processes */

    FormFactors<VToV>::~FormFactors()
    {
    }

    const std::map<FormFactorFactory<VToV>::KeyType, FormFactorFactory<VToV>::ValueType>
    FormFactorFactory<VToV>::form_factors
    {
        // parametrizations
        // b -> c
        // not yet supported
        { "B^*->D^*::BGJvD2019",      &HQETFormFactors<BstarToDstar, VToV>::make      },
    };

    std::shared_ptr<FormFactors<VToV>>
    FormFactorFactory<VToV>::create(const QualifiedName & name, const Parameters & parameters, const Options & options)
    {
        Context ctx("When creating a V->V form factor");

        std::shared_ptr<FormFactors<VToV>> result;

        auto i = FormFactorFactory<VToV>::form_factors.find(name);
        if (FormFactorFactory<VToV>::form_factors.end() != i)
        {
            result.reset(i->second(parameters, name.options() + options));
            return result;
        }

        throw NoSuchFormFactorError(name.prefix_part().str(), name.name_part().str());
        return result;
    }

    OptionSpecification
    FormFactorFactory<VToV>::option_specification(const qnp::Prefix & process)
    {
        OptionSpecification result { "form-factors", {}, "" };
        for (const auto & ff : FormFactorFactory<VToV>::form_factors)
        {
            if (process == std::get<0>(ff).prefix_part())
                result.allowed_values.push_back(std::get<0>(ff).name_part().str());
        }

        return result;
    }
}
