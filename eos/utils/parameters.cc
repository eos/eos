/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2014, 2015 Danny van Dyk
 * Copyright (c) 2010 Christian Wacker
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

#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/stringify.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

#include <cmath>
#include <map>
#include <random>
#include <vector>

namespace eos
{
    struct Parameter::Template
    {
        std::string name;

        double min, central, max;
    };

    struct Parameter::Data :
        Parameter::Template
    {
        double value;

        Parameter::Id id;

        Data(const Parameter::Template & t, const Parameter::Id & i) :
            Parameter::Template(t),
            value(t.central),
            id(i)
        {
        }
    };

    struct Parameters::Data
    {
        std::vector<Parameter::Data> data;
    };

    template class WrappedForwardIterator<Parameters::IteratorTag, Parameter>;

    template <>
    struct Implementation<Parameters>
    {
        std::shared_ptr<Parameters::Data> parameters_data;

        std::map<std::string, unsigned> parameters_map;

        std::vector<Parameter> parameters;

        Implementation(const std::initializer_list<Parameter::Template> & list) :
            parameters_data(new Parameters::Data)
        {
            unsigned idx(0);
            for (auto i(list.begin()), i_end(list.end()) ; i != i_end ; ++i, ++idx)
            {
                parameters_data->data.push_back(Parameter::Data(*i, idx));
                parameters_map[i->name] = idx;
                parameters.push_back(Parameter(parameters_data, idx));
            }
        }

        Implementation(const Implementation & other) :
            parameters_data(new Parameters::Data(*other.parameters_data)),
            parameters_map(other.parameters_map)
        {
            parameters.reserve(other.parameters.size());
            for (unsigned i = 0 ; i != parameters.size() ; ++i)
            {
                parameters.push_back(Parameter(parameters_data, i));
            }
        }
    };

    Parameters::Parameters(Implementation<Parameters> * imp) :
        PrivateImplementationPattern<Parameters>(imp)
    {
    }

    Parameters::~Parameters()
    {
    }

    Parameters
    Parameters::clone() const
    {
        return Parameters(new Implementation<Parameters>(*_imp));
    }

    Parameter
    Parameters::operator[] (const std::string & name) const
    {
        auto i(_imp->parameters_map.find(name));

        if (_imp->parameters_map.end() == i)
            throw UnknownParameterError(name);

        return Parameter(_imp->parameters_data, i->second);
    }

    Parameter
    Parameters::operator[] (const Parameter::Id & id) const
    {
        if (id >= _imp->parameters.size())
            throw InternalError("Parameters::operator[] (Parameter::Id): invalid id '" + stringify(id) + "'");

        return _imp->parameters[id];
    }

    void
    Parameters::set(const std::string & name, const double & value)
    {
        auto i(_imp->parameters_map.find(name));

        if (_imp->parameters_map.end() == i)
            throw UnknownParameterError(name);

        _imp->parameters_data->data[i->second].value = value;
    }

    Parameters::Iterator
    Parameters::begin() const
    {
        return Parameters::Iterator(_imp->parameters.begin());
    }

    Parameters::Iterator
    Parameters::end() const
    {
        return Parameters::Iterator(_imp->parameters.end());
    }

    bool
    Parameters::operator!= (const Parameters & rhs) const
    {
        return rhs._imp.get() != this->_imp.get();
    }

    Parameters
    Parameters::Defaults()
    {
        return Parameters(new Implementation<Parameters>({
            Parameter::Template{"hbar",                                          +6.58211913e-25, +6.58211928e-25, +6.58211943e-25}, // GeV s, cf. [PDG2012]:  p. 4

            // Wilson coefficients C1 - C6 at mu = 4.2 GeV to NNLL accuary as calculated in the SM by EOS's StandardModel class.
            // For the calculations, cf. [BMU1999].
            Parameter::Template{"c1",                                            -0.29063621,     -0.29063621,     -0.29063621    },
            Parameter::Template{"c2",                                            +1.01029623,     +1.01029623,     +1.01029623    },
            Parameter::Template{"c3",                                            -0.00616220,     -0.00616220,     -0.00616220    },
            Parameter::Template{"c4",                                            -0.08730376,     -0.08730376,     -0.08730376    },
            Parameter::Template{"c5",                                            +0.00042854,     +0.00042854,     +0.00042854    },
            Parameter::Template{"c6",                                            +0.00115807,     +0.00115807,     +0.00115807    },
            Parameter::Template{"Abs{c7}",                                       +0.0,            +0.33726473,     +1.0           },
            Parameter::Template{"Arg{c7}",                                       +0.0,            +M_PI,           +2.0 * M_PI    },
            Parameter::Template{"Re{c7}",                                        -1.0,            -0.33726473,     +1.0           },
            Parameter::Template{"Im{c7}",                                        -1.0,            +0.0,            +1.0           },
            Parameter::Template{"c8",                                            -0.18288898,     -0.18288898,     -0.18288898    },
            Parameter::Template{"Abs{c9}",                                       +0.0,            +4.27342842,    +15.0           },
            Parameter::Template{"Arg{c9}",                                       +0.0,            +0.0,            +2.0 * M_PI    },
            Parameter::Template{"Re{c9}",                                       -15.0,            +4.27342842,    +15.0           },
            Parameter::Template{"Im{c9}",                                       -15.0,            +0.0,           +15.0           },
            Parameter::Template{"Abs{c10}",                                      +0.0,            +4.16611761,    +15.0           },
            Parameter::Template{"Arg{c10}",                                      +0.0,            +M_PI,           +2.0 * M_PI    },
            Parameter::Template{"Re{c10}",                                      -15.0,            -4.16611761,    +15.0           },
            Parameter::Template{"Im{c10}",                                      -15.0,            +0.0,           +15.0           },
            // Primed Wilson coefficients are negligible in the SM
            Parameter::Template{"Abs{c7'}",                                      +0.0,            +0.0,            +1.0           },
            Parameter::Template{"Arg{c7'}",                                      +0.0,            +0.0,            +2.0 * M_PI    },
            Parameter::Template{"Re{c7'}",                                       -1.0,            +0.0,            +1.0           },
            Parameter::Template{"Im{c7'}",                                       -1.0,            +0.0,            +1.0           },
            Parameter::Template{"c8'",                                           +0.0,            +0.0,            +0.0           },
            Parameter::Template{"Abs{c9'}",                                      +0.0,            +0.0,           +15.0           },
            Parameter::Template{"Arg{c9'}",                                      +0.0,            +0.0,            +2.0 * M_PI    },
            Parameter::Template{"Re{c9'}",                                      -15.0,            +0.0,           +15.0           },
            Parameter::Template{"Im{c9'}",                                      -15.0,            +0.0,           +15.0           },
            Parameter::Template{"Abs{c10'}",                                     +0.0,            +0.0,           +15.0           },
            Parameter::Template{"Arg{c10'}",                                     +0.0,            +0.0,            +2.0 * M_PI    },
            Parameter::Template{"Re{c10'}",                                     -15.0,            +0.0,           +15.0           },
            Parameter::Template{"Im{c10'}",                                     -15.0,            +0.0,           +15.0           },
            // scalar, pseudoscalar, and tensor b -> s ll coefficients
            Parameter::Template{"Abs{cS}",                                       +0.0,            +0.0,            +1.0           },
            Parameter::Template{"Arg{cS}",                                       +0.0,            +0.0,            +2.0 * M_PI    },
            Parameter::Template{"Re{cS}",                                        -1.0,            +0.0,            +1.0           },
            Parameter::Template{"Im{cS}",                                        -1.0,            +0.0,            +1.0           },
            Parameter::Template{"Abs{cS'}",                                      +0.0,            +0.0,            +1.0           },
            Parameter::Template{"Arg{cS'}",                                      +0.0,            +0.0,            +2.0 * M_PI    },
            Parameter::Template{"Re{cS'}",                                       -1.0,            +0.0,            +1.0           },
            Parameter::Template{"Im{cS'}",                                       -1.0,            +0.0,            +1.0           },
            Parameter::Template{"Abs{cP}",                                       +0.0,            +0.0,            +1.0           },
            Parameter::Template{"Arg{cP}",                                       +0.0,            +0.0,            +2.0 * M_PI    },
            Parameter::Template{"Re{cP}",                                        -1.0,            +0.0,            +1.0           },
            Parameter::Template{"Im{cP}",                                        -1.0,            +0.0,            +1.0           },
            Parameter::Template{"Abs{cP'}",                                      +0.0,            +0.0,            +1.0           },
            Parameter::Template{"Arg{cP'}",                                      +0.0,            +0.0,            +2.0 * M_PI    },
            Parameter::Template{"Re{cP'}",                                       -1.0,            +0.0,            +1.0           },
            Parameter::Template{"Im{cP'}",                                       -1.0,            +0.0,            +1.0           },
            Parameter::Template{"Abs{cT}",                                       +0.0,            +0.0,            +1.0           },
            Parameter::Template{"Arg{cT}",                                       +0.0,            +0.0,            +2.0 * M_PI    },
            Parameter::Template{"Re{cT}",                                        -1.0,            +0.0,            +1.0           },
            Parameter::Template{"Im{cT}",                                        -1.0,            +0.0,            +1.0           },
            Parameter::Template{"Abs{cT5}",                                      +0.0,            +0.0,            +1.0           },
            Parameter::Template{"Arg{cT5}",                                      +0.0,            +0.0,            +2.0 * M_PI    },
            Parameter::Template{"Re{cT5}",                                       -1.0,            +0.0,            +1.0           },
            Parameter::Template{"Im{cT5}",                                       -1.0,            +0.0,            +1.0           },
            // b->u Wilson coefficients
            Parameter::Template{"Re{cVL}",                                        1.0,             1.0,             1.0           },
            Parameter::Template{"Im{cVL}",                                        0.0,             0.0,             0.0           },
            Parameter::Template{"Re{cVR}",                                        0.0,             0.0,             0.0           },
            Parameter::Template{"Im{cVR}",                                        0.0,             0.0,             0.0           },
            Parameter::Template{"Re{cSL}",                                        0.0,             0.0,             0.0           },
            Parameter::Template{"Im{cSL}",                                        0.0,             0.0,             0.0           },
            Parameter::Template{"Re{cSR}",                                        0.0,             0.0,             0.0           },
            Parameter::Template{"Im{cSR}",                                        0.0,             0.0,             0.0           },
            Parameter::Template{"Re{cT}",                                         0.0,             0.0,             0.0           },
            Parameter::Template{"Im{cT}",                                         0.0,             0.0,             0.0           },
            // Factorization scale
            Parameter::Template{"mu",                                            +2.4,            +4.2,            +9.6           },
            // GSW Parameter
            Parameter::Template{"GSW::sin^2(theta)",                             +0.23104,        +0.23116,        +0.23128       }, // cf. [PDG2012]:  p. 4
            // Wolfenstein parameters of CKM, cf. [UTFIT2013]
            Parameter::Template{"CKM::A",                                        +0.814,          +0.827,          +0.840         },
            Parameter::Template{"CKM::lambda",                                   +0.22470,        +0.22535,        +0.22600       },
            Parameter::Template{"CKM::rhobar",                                   +0.111,          +0.132,          +0.153         },
            Parameter::Template{"CKM::etabar",                                   +0.336,          +0.350,          +0.364         },
            // CKM matrix elements for use by CKMScanModel
            Parameter::Template{"CKM::abs(V_ud)",                                +1.00,           +1.00,           +1.00          },
            Parameter::Template{"CKM::arg(V_ud)",                                +1.00,           +1.00,           +1.00          },
            Parameter::Template{"CKM::abs(V_us)",                                +0.22,           +0.22,           +0.22          },
            Parameter::Template{"CKM::arg(V_us)",                                +0.22,           +0.22,           +0.22          },
            Parameter::Template{"CKM::abs(V_ub)",                                +0.00,           +0.00,           +0.00          },
            Parameter::Template{"CKM::arg(V_ub)",                                +0.00,           +0.00,           +0.00          },
            Parameter::Template{"CKM::abs(V_cd)",                                +1.00,           +1.00,           +1.00          },
            Parameter::Template{"CKM::arg(V_cd)",                                +1.00,           +1.00,           +1.00          },
            Parameter::Template{"CKM::abs(V_cs)",                                +0.22,           +0.22,           +0.22          },
            Parameter::Template{"CKM::arg(V_cs)",                                +0.22,           +0.22,           +0.22          },
            Parameter::Template{"CKM::abs(V_cb)",                                +0.00,           +0.00,           +0.00          },
            Parameter::Template{"CKM::arg(V_cb)",                                +0.00,           +0.00,           +0.00          },
            Parameter::Template{"CKM::abs(V_td)",                                +0.00,           +0.00,           +0.00          },
            Parameter::Template{"CKM::arg(V_td)",                                +0.00,           +0.00,           +0.00          },
            Parameter::Template{"CKM::abs(V_ts)",                                +0.04,           +0.04,           +0.04          },
            Parameter::Template{"CKM::arg(V_ts)",                                +0.04,           +0.04,           +0.04          },
            Parameter::Template{"CKM::abs(V_tb)",                                +1.00,           +1.00,           +1.00          },
            Parameter::Template{"CKM::arg(V_tb)",                                +1.00,           +1.00,           +1.00          },
            // QED inputs
            Parameter::Template{"QED::alpha_e(m_b)",                             +1.0/133.0,      +1.0/133.0,      +1.0/128.0     }, // alpha_e(m_b) .. alpha_e(m_W)
            // QCD inputs
            Parameter::Template{"QCD::alpha_s(MZ)",                              +0.1191,         +0.1184,         +0.1177        }, // cf. [PDG2012], p. 1
            Parameter::Template{"QCD::mu_t",                                     +170.0,          +170.0,          +170.0         },
            Parameter::Template{"QCD::mu_b",                                     +4.2,            +4.2,            +4.2           },
            Parameter::Template{"QCD::mu_c",                                     +1.0,            +1.0,            +1.0           },
            Parameter::Template{"QCD::Lambda",                                   +0.5,            +0.5,            +0.5           },
            Parameter::Template{"QCD::cond_GG",                                  +0.000,          +0.012,          +0.018         },
            Parameter::Template{"QCD::m_0",                                      +0.7746,         +0.89443,        +1.0           }, // GeV, cf. [DKMMO2008]
            Parameter::Template{"QCD::r_vac",                                    +0.1,            +0.55,           +1.0           }, // cf. [DKMMO2008]
            // G_Fermi
            Parameter::Template{"G_Fermi",                                       +1.1663781e-5,   +1.1663787e-5,   +1.1663793e-5  }, // cf. [PDG2012], p. 5

            /* Masses in GeV */
            // Lepton masses
            Parameter::Template{"mass::e",                                       +5.10999e-4,     +5.10999e-4,     +5.10999e-4    }, // cf. [PDG2012], p. 13
            Parameter::Template{"mass::mu",                                      +1.05658e-1,     +1.05658e-1,     +1.05658e-1    }, // cf. [PDG2012], p. 13
            Parameter::Template{"mass::tau",                                     +1.77666,        +1.77682,        +1.77698       }, // cf. [PDG2012], p. 14
            // Quark masses
            Parameter::Template{"mass::d(2GeV)",                                 +0.0045,         +0.0048,         +0.0055        }, // cf. [PDG2012], p. 21
            Parameter::Template{"mass::ud(2GeV)",                                +0.008,          +0.008,          +0.008         },
            Parameter::Template{"mass::s(2GeV)",                                 +0.090,          +0.095,          +0.010         }, // cf. [PDG2012], p. 21
            Parameter::Template{"mass::c",                                       +1.250,          +1.275,          +1.300         }, // cf. [PDG2012], p. 21
            Parameter::Template{"mass::b(MSbar)",                                +4.15,           +4.18,           +4.21          }, // cf. [PDG2012], p. 21
            Parameter::Template{"mass::t(pole)",                                 +172.5,          +173.5,          +174.5         }, // cf. [PDG2012], p. 22

            // Pi meson masses
            Parameter::Template{"mass::pi^0",                                    +0.1349760,      +0.1349766,      +0.1349772     }, // cf. [PDG-Live], May 2014
            Parameter::Template{"mass::pi^+",                                    +0.13956983,     +0.13957018,     +0.13957053    },
            // K meson masses
            Parameter::Template{"mass::K_d",                                     +0.497590,       +0.497614,       +0.497638      }, // cf. [PDG2012], p. 41
            Parameter::Template{"mass::K_u",                                     +0.493661,       +0.493677,       +0.493693      }, // cf. [PDG2012], p. 39
            Parameter::Template{"mass::K^*_d",                                   +0.89572,        +0.89594,        +0.89616       }, // cf. [PDG2012], p. 45
            Parameter::Template{"mass::K^*_u",                                   +0.89140,        +0.89166,        +0.89192       }, // cf. [PDG2012], p. 45
            // B meson masses
            Parameter::Template{"mass::B_d",                                     +5.27941,        +5.27958,        +5.27975       }, // cf. [PDG2012], p. 84
            Parameter::Template{"mass::B_u",                                     +5.27908,        +5.27925,        +5.27942       }, // cf. [PDG2012], p. 70
            Parameter::Template{"mass::B_s",                                     +5.36653,        +5.36677,        +5.36701       }, // cf. [PDG2012], p. 105
            // s baryon masses
            Parameter::Template{"mass::Lambda",                                  +1.115677,       +1.115683,       +1.115689      }, // cf. [PDG2012], p. 149
            // b baryon masses
            Parameter::Template{"mass::Lambda_b",                                +5.6187,         +5.6194,         +5.6201        }, // cf. [PDG2012], p. 165
            // Gauge boson masses
            Parameter::Template{"mass::W",                                       +80.370,         +80.385,         +80.400        }, // cf. [PDG2012], p. 8
            Parameter::Template{"mass::Z",                                       +91.1855,        +91.1876,        +91.1897       }, // cf. [PDG2012], p. 9

            /* Decay constants */
            Parameter::Template{"decay-constant::B_d",                           +0.1859,         +0.1906,         +0.1953        }, // GeV, cf. [LATAVG2011]
            Parameter::Template{"decay-constant::B_u",                           +0.1859,         +0.1906,         +0.1953        }, // GeV, cf. [LATAVG2011]
            Parameter::Template{"decay-constant::B_s",                           +0.2226,         +0.2276,         +0.2326        }, // GeV, cf. [LATAVG2011]
            Parameter::Template{"decay-constant::K_d",                           +0.155,          +0.1561,         +0.1572        }, // GeV, cf. [LATAVG2011]
            Parameter::Template{"decay-constant::K_u",                           +0.155,          +0.1561,         +0.1572        }, // GeV, cf. [LATAVG2011]
            Parameter::Template{"decay-constant::pi",                            +0.1288,         +0.1302,         +0.1326        }, // GeV, cf. [FLAG2013], eq. (43)

            /* Pion LCDA parameters */
            Parameter::Template{"pi::a2@1GeV",                                   +0.16,           +0.16,           +0.16          }, // cf. [DKMMO2008]
            Parameter::Template{"pi::a4@1GeV",                                   +0.04,           +0.04,           +0.04          }, // cf. [DKMMO2008]
            Parameter::Template{"pi::f3@1GeV",                                   +0.0030,         +0.0045,         +0.0060        }, // GeV^2, cf. [DKMMO2008]
            Parameter::Template{"pi::omega3@1GeV",                               -2.2,            -1.5,            -0.8           }, // cf. [DKMMO2008]
            Parameter::Template{"pi::omega4@1GeV",                               +0.2,            +0.2,            +0.2           }, // cf. [DKMMO2008]
            Parameter::Template{"pi::delta^2@1GeV",                              +0.12,           +0.18,           +0.24          }, // GeV^2, cf. [DKMMO2008]

            /* b->s matching parameters */
            Parameter::Template{"b->s::mu_0c",                                   +80.0,           +80.0,           +80.0          },
            Parameter::Template{"b->s::mu_0t",                                   +120.0,          +120.0,          +120.0         },

            /* Mean life times */
            Parameter::Template{"life_time::B_d",                                +1.512e-12,      +1.519e-12,      +1.526e-12     }, // cf. [PDG2012], p. 84
            Parameter::Template{"life_time::B_u",                                +1.633e-12,      +1.641e-12,      +1.649e-12     }, // cf. [PDG2012], p. 70
            Parameter::Template{"life_time::B_s",                                +1.505e-12,      +1.516e-12,      +1.527e-12     }, // cf. [PDG-Live], September 2013
            Parameter::Template{"life_time::B@Y(4S)",                            +1.573e-12,      +1.580e-12,      +1.589e-12     }, // based on 50:50 admixture of B^+B^- and B^0Bbar^0 @ Y(4S)
            Parameter::Template{"life_time::Lambda_b",                           +1.438e-12,      +1.451e-12,      +1.464e-12     }, // cf. PDG2014 update by [HFAG]

            /* Decay width differences in neutral meson systems */
            Parameter::Template{"life_time::Delta_B_d",                          +0.0,            +0.0,            +0.0           }, // smaller than 1eV, i.e. zero for all intents and purposes, cf. [PDG2012]
            Parameter::Template{"life_time::Delta_B_s",                          +0.070e+12,      +0.081e+12,      +0.092e+12     }, // cf. [PDG2012+]

            /* Decay parameter of the Lambda hyperon */
            Parameter::Template{"Lambda::alpha",                                 +0.629,          +0.642,          +0.655         }, // cf. [PDG2012+]

            // Form factor uncertainties
            Parameter::Template{"formfactors::xi_perp_uncertainty",              +0.89,           +1.0,            +1.11          },
            Parameter::Template{"formfactors::xi_par_uncertainty",               +0.86,           +1.0,            +1.14          },

            // form factor parameters for LCSR calculation according to [DKMMO2008]
            Parameter::Template{"B->pi::zeta(NNLO)@DKMMO2008",                   -1.0,             0.0,            +1.0           },
            Parameter::Template{"B->pi::M^2@DKMMO2008",                         +15.0,           +18.0,           +20.0           }, // GeV^2, cf. [DKMMO2008]
            Parameter::Template{"B->pi::Mp^2@DKMMO2008",                         +5.0,            +5.0,            +5.0           }, // GeV^2, cf. [DKMMO2008]
            Parameter::Template{"B->pi::mu@DKMMO2008",                           +3.0,            +3.0,            +3.0           }, // GeV, cf. [DKMMO2008]
            Parameter::Template{"B->pi::s_0^B@DKMMO2008",                       +35.75,          +35.75,          +35.75          }, // GeV^2, cf. [DKMMO2008]
            Parameter::Template{"B->pi::sp_0^B@DKMMO2008",                      +35.5,           +35.6,           +36.0           }, // GeV^2, cf. [DKMMO2008]

            // form factor parameters for B->pi, B->K in [BCL2008] parametrization
            Parameter::Template{"B->pi::f_+(0)@BCL2008",                         +0.22,           +0.26,           +0.30          },
            Parameter::Template{"B->pi::b_+^1@BCL2008",                          +0.00,           +0.00,           +0.00          },
            Parameter::Template{"B->pi::b_+^2@BCL2008",                          +0.00,           +0.00,           +0.00          },
            Parameter::Template{"B->pi::b_0^1@BCL2008",                          +0.00,           +0.00,           +0.00          },
            Parameter::Template{"B->pi::b_0^2@BCL2008",                          +0.00,           +0.00,           +0.00          },
            Parameter::Template{"B->pi::f_T(0)@BCL2008",                         +0.00,           +0.00,           +0.00          },
            Parameter::Template{"B->pi::b_T^1@BCL2008",                          +0.00,           +0.00,           +0.00          },
            Parameter::Template{"B->pi::b_T^2@BCL2008",                          +0.00,           +0.00,           +0.00          },

            // form factor parameters for B->K^* according to [BZ2004] (approximate)
            Parameter::Template{"B->K^*::a0_uncertainty@BZ2004",                 +0.85,           +1.0,            +1.15          },
            Parameter::Template{"B->K^*::a1_uncertainty@BZ2004",                 +0.85,           +1.0,            +1.15          },
            Parameter::Template{"B->K^*::a2_uncertainty@BZ2004",                 +0.85,           +1.0,            +1.15          },
            Parameter::Template{"B->K^*::v_uncertainty@BZ2004",                  +0.85,           +1.0,            +1.15          },

            // form factor parameters for B->K^* according to [KMPW2010], Table 4, p. 31
            Parameter::Template{"B->K^*::F^V(0)@KMPW2010",                       +0.24,           +0.36,           +0.59          },
            Parameter::Template{"B->K^*::F^A0(0)@KMPW2010",                      +0.22,           +0.29,           +0.39          },
            Parameter::Template{"B->K^*::F^A1(0)@KMPW2010",                      +0.15,           +0.25,           +0.41          },
            Parameter::Template{"B->K^*::F^A2(0)@KMPW2010",                      +0.13,           +0.23,           +0.42          },
            Parameter::Template{"B->K^*::b^V_1@KMPW2010",                        -5.2,            -4.8,            -4.0           },
            Parameter::Template{"B->K^*::b^A0_1@KMPW2010",                      -21.2,            -18.2,          -16.9           },
            Parameter::Template{"B->K^*::b^A1_1@KMPW2010",                       -0.46,           +0.34,           +1.2           },
            Parameter::Template{"B->K^*::b^A2_1@KMPW2010",                       -2.2,            -0.85,           +2.03          },

            // form factor parameter for B_s->K^* according to [FMvD2015].
            Parameter::Template{"B_s->K^*::F_perp(0)@FMvD2015",                   0.349,           0.349,           0.349         },
            Parameter::Template{"B_s->K^*::F_para(0)@FMvD2015",                   0.379,           0.379,           0.379         },
            Parameter::Template{"B_s->K^*::F_long(0)@FMvD2015",                   0.314,           0.314,           0.314         },
            Parameter::Template{"B_s->K^*::F_time(0)@FMvD2015",                   0.314,           0.314,           0.314         },
            Parameter::Template{"B_s->K^*::beta_perp^1@FMvD2015",                -4.980,          -4.980,          -4.980         },
            Parameter::Template{"B_s->K^*::beta_para^1@FMvD2015",                +0.065,          +0.065,          +0.065         },
            Parameter::Template{"B_s->K^*::beta_time^1@FMvD2015",                +0.065,          +0.065,          +0.065         },

            // form factor parameters for B->K^* simple series expansion (SSE) based on
            // LCSR according to [BFW2010], table 10, p. 33, corrected values from Aoife Bharucha
            Parameter::Template{"B->K^*::beta^V0_0@BFW2010",                     +0.31,           +0.31,           +0.31,         },
            Parameter::Template{"B->K^*::beta^V0_1@BFW2010",                     +0.74,           +0.74,           +0.74,         },
            Parameter::Template{"B->K^*::beta^V1_0@BFW2010",                     +0.62,           +0.62,           +0.62,         },
            Parameter::Template{"B->K^*::beta^V1_1@BFW2010",                     -1.4,            -1.4,            -1.4,          },
            Parameter::Template{"B->K^*::beta^V2_0@BFW2010",                     +0.45,           +0.45,           +0.45,         },
            Parameter::Template{"B->K^*::beta^V2_1@BFW2010",                     +0.35,           +0.35,           +0.35,         },
            Parameter::Template{"B->K^*::beta^Vt_0@BFW2010",                     +0.49,           +0.49,           +0.49,         },
            Parameter::Template{"B->K^*::beta^Vt_1@BFW2010",                     -1.4,            -1.4,            -1.4,          },

            // form factor parameters for B->K according to [BZ2004v2] (approximate)
            Parameter::Template{"B->K::fp_uncertainty@BZ2004v2",                 +0.85,           +1.0,            +1.15          },
            Parameter::Template{"B->K::f0_uncertainty@BZ2004v2",                 +0.85,           +1.0,            +1.15          },
            Parameter::Template{"B->K::ft_uncertainty@BZ2004v2",                 +0.85,           +1.0,            +1.15          },

            // form factor parameters for B->K according to [KMPW2010], Table 4, p. 31
            Parameter::Template{"B->K::F^p(0)@KMPW2010",                         +0.32,           +0.34,           +0.39          },
            Parameter::Template{"B->K::F^t(0)@KMPW2010",                         +0.36,           +0.39,           +0.44          },
            Parameter::Template{"B->K::b^p_1@KMPW2010",                          -3.7,            -2.1,            -1.2           },
            Parameter::Template{"B->K::b^0_1@KMPW2010",                          -5.2,            -4.3,            -3.5           },
            Parameter::Template{"B->K::b^t_1@KMPW2010",                          -4.2,            -2.2,            -1.2           },

            // form factor parameters for B->K simple series expansion (SSE) based on
            // LCSR according to [BFW2010], table 6, p. 22, corrected values from Aoife Bharucha
            Parameter::Template{"B->K::alpha^V0_0@BFW2010",                      +0.48,           +0.48,           +0.48,          },
            Parameter::Template{"B->K::alpha^V0_1@BFW2010",                      -1.05,           -1.05,           -1.05,          },
            Parameter::Template{"B->K::alpha^Vt_0np@BFW2010",                    +0.52,           +0.52,           +0.52,          },
            Parameter::Template{"B->K::alpha^Vt_1np@BFW2010",                    -1.4,            -1.4,            -1.4,           },
            Parameter::Template{"B->K::alpha^T0_0@BFW2010",                      +0.48,           +0.48,           +0.48,          },
            Parameter::Template{"B->K::alpha^T0_1@BFW2010",                      -1.09,           -1.09,           -1.09,          },

            // form factor parameters for Lambda_b -> Lambda decays according to
            // [BFvD2014], table ?, p. ??
            Parameter::Template{"Lambda_b->Lambda::f_0^V(0)@BFvD2014",           +0.22,           +0.28,           +0.34,         },
            Parameter::Template{"Lambda_b->Lambda::b_1_0^V@BFvD2014",           -19.00,          -11.00,           -6.00,         },
            Parameter::Template{"Lambda_b->Lambda::f_0^A(0)@BFvD2014",           +0.21,           +0.27,           +0.33,         },
            Parameter::Template{"Lambda_b->Lambda::b_1_0^A@BFvD2014",           -16.00,           -8.00,           -5.00,         },
            Parameter::Template{"Lambda_b->Lambda::f_perp^V(0)@BFvD2014",        +0.22,           +0.28,           +0.34,         },
            Parameter::Template{"Lambda_b->Lambda::b_1_perp^V@BFvD2014",        -19.00,          -11.00,           -6.00,         },
            Parameter::Template{"Lambda_b->Lambda::f_perp^A(0)@BFvD2014",        +0.24,           +0.30,           +0.36,         },
            Parameter::Template{"Lambda_b->Lambda::b_1_perp^A@BFvD2014",        -12.00,           -6.00,           -3.00,         },
            Parameter::Template{"Lambda_b->Lambda::sigma_perp^V@BFvD2014",       -0.20,            0.00,           +0.20,         },
            Parameter::Template{"Lambda_b->Lambda::sigma_perp^A@BFvD2014",       -0.30,            0.00,           +0.30,         },
            Parameter::Template{"Lambda_b->Lambda::sigma_long^V@BFvD2014",       -0.21,            0.00,           +0.21,         },
            Parameter::Template{"Lambda_b->Lambda::sigma_long^A@BFvD2014",       -0.30,            0.00,           +0.30,         },

            // B LCDA parameters
            Parameter::Template{"lambda_B_p",                                    +0.370,          +0.485,          +0.600         }, // GeV, cf. [BHvD2010], Table I
            // B->K LCDA Parameter
            Parameter::Template{"B->K::a_1@1GeV",                                +0.03,           +0.06,           +0.09          }, // cf. [BBL2006], Table 3
            Parameter::Template{"B->K::a_2@1GeV",                                +0.10,           +0.25,           +0.4           }, // cf. [BBL2006], Table 3
            Parameter::Template{"B->K::a_4@1GeV",                                -0.115,          -0.015,          +0.085         }, // cf. [BZ2004v3], Eq. (24)
            Parameter::Template{"B->K::a_1@2.2GeV",                              +0.024,          +0.048,          +0.071         }, // cf. [BBL2006], Table 3 and cf. [BHP2007] App. A, pp. 24-25
            Parameter::Template{"B->K::a_2@2.2GeV",                              +0.070,          +0.174,          +0.278         }, // cf. [BBL2006], Table 3 and cf. [BHP2007] App. A, pp. 24-25
            Parameter::Template{"B->K::a_4@2.2GeV",                              -0.0679,         -0.0089,         +0.0502        }, // cf. [BZ2004v3], Eq. (24) and cf. [BHP2007] App. A, pp. 24-25
            // B->K^*, K^* LCDA parameters
            Parameter::Template{"B->K^*::a_1_par",                               +0.03,           +0.1,            +0.17          },
            Parameter::Template{"B->K^*::a_2_par",                               +0.0,            +0.1,            +0.2           },
            Parameter::Template{"B->K^*::a_1_perp",                              +0.03,           +0.1,            +0.17          },
            Parameter::Template{"B->K^*::a_2_perp",                              +0.0,            +0.1,            +0.2           },
            Parameter::Template{"B->K^*::f_Kstar_par",                           +0.212,          +0.217,          +0.222         }, // GeV, cf. [BHvD2010], Table I
            Parameter::Template{"B->K^*::f_Kstar_perp@2GeV",                     +0.168,          +0.173,          +0.178         }, // GeV @2 Gev, 0.185 +/-0.005 GeV, cf. [BHvD2010], Table I
            // B->K^*ll uncertainties from subleading terms for Large Recoil
            Parameter::Template{"B->K^*ll::sl_uncertainty@LargeRecoil",          +0.95,           +1.0,            +1.05          },
            Parameter::Template{"B->K^*ll::A_0^L_uncertainty@LargeRecoil",       +0.95,           +1.0,            +1.05          },
            Parameter::Template{"B->K^*ll::A_0^R_uncertainty@LargeRecoil",       +0.95,           +1.0,            +1.05          },
            Parameter::Template{"B->K^*ll::A_par^L_uncertainty@LargeRecoil",     +0.95,           +1.0,            +1.05          },
            Parameter::Template{"B->K^*ll::A_par^R_uncertainty@LargeRecoil",     +0.95,           +1.0,            +1.05          },
            Parameter::Template{"B->K^*ll::A_perp^L_uncertainty@LargeRecoil",    +0.95,           +1.0,            +1.05          },
            Parameter::Template{"B->K^*ll::A_perp^R_uncertainty@LargeRecoil",    +0.95,           +1.0,            +1.05          },
            // B->Pll uncertainties at subleading order at Large Recoil
            Parameter::Template{"B->Pll::Lambda_pseudo@LargeRecoil",             -0.5,            +0.0,            +0.5           },
            Parameter::Template{"B->Pll::sl_phase_pseudo@LargeRecoil",           -M_PI/2.0,       +0.0,            +M_PI/2.0      },
            // B->Pll uncertainties at subleading order at Low Recoil
            Parameter::Template{"B->Pll::Lambda_pseudo@LowRecoil",               -0.5,            +0.0,            +0.5           },
            Parameter::Template{"B->Pll::sl_phase_pseudo@LowRecoil",             -M_PI/2.0,       +0.0,            +M_PI/2.0      },
            // B->Vll uncertainties at subleading order at Low Recoil
            Parameter::Template{"B->Vll::Lambda@LowRecoil",                      -0.5,            +0.0,            +0.5           },
            Parameter::Template{"B->Vll::Lambda_0@LowRecoil",                    -0.5,            +0.0,            +0.5           },
            Parameter::Template{"B->Vll::Lambda_pa@LowRecoil",                   -0.5,            +0.0,            +0.5           },
            Parameter::Template{"B->Vll::Lambda_pp@LowRecoil",                   -0.5,            +0.0,            +0.5           },
            Parameter::Template{"B->Vll::sl_phase@LowRecoil",                    -M_PI/2.0,       +0.0,            +M_PI/2.0      },
            Parameter::Template{"B->Vll::sl_phase_0@LowRecoil",                  -M_PI/2.0,       +0.0,            +M_PI/2.0      },
            Parameter::Template{"B->Vll::sl_phase_pa@LowRecoil",                 -M_PI/2.0,       +0.0,            +M_PI/2.0      },
            Parameter::Template{"B->Vll::sl_phase_pp@LowRecoil",                 -M_PI/2.0,       +0.0,            +M_PI/2.0      },
            // HQET parameters
            Parameter::Template{"HQET::lambda_1",                                -0.20,           -0.20,           -0.20          }, // cf. [ALGH2001], Table 2, p. 13
            Parameter::Template{"HQET::lambda_2",                                +0.12,           +0.12,           +0.12          }, // cf. [ALGH2001], Table 2, p. 13
            // Heavy Quark Expansion parameters for hadronic matrix elements ~ <B|O|B>
            Parameter::Template{"B->B::mu_pi^2@1GeV",                            +0.35,           +0.45,           +0.55          }, // cf. [BBMU2003], Eq. (19), p. 9
            Parameter::Template{"B->B::mu_G^2@1GeV",                             +0.33,           +0.35,           +0.38          }, // cf. [BBMU2003], Eq. (17), p. 9
            Parameter::Template{"B->B::rho_D^3@1GeV",                            +0.10,           +0.20,           +0.30          }, // cf. [BBMU2003], between Eqs. (19),(20), p. 9
            Parameter::Template{"B->B::rho_LS^3@1GeV",                           -0.30,           -0.15,           -0.00          }, // cf. [BBMU2003], Eq. (20), p. 9
            // B->X_s gamma SM theory uncertainty
            Parameter::Template{"B->X_sgamma::uncertainty",                      -1.0,            +0.0,            +1.0           },

            // todo all exp::* could be in constraint.cc as well
            // Experimental Input
            Parameter::Template{"exp::BR(B->X_clnu)",                            +0.1005,         +0.1033,         +0.1061        }, // cf. [PDG2012], p. 88
            // todo which numbers to take?
            Parameter::Template{"exp::C(B->X_clnu, B->X_ulnu)",                  +0.57,           +0.58,           +0.59          },
            Parameter::Template{"exp::CKM(B->X_sll, B->X_clnu)",                 +0.975218,       +0.98549,        +0.995277      },
            // Parameterize unknown admixture of l=e, l=mu in B->X_sll
            Parameter::Template{"exp::Admixture-BR(B->X_sll)",                   +0.95,            +1.0,            +1.05         }, // BR varies by up to +/-5% between l=mu and l=e
        }));
    }

    Parameter::Parameter(const std::shared_ptr<Parameters::Data> & parameters_data, unsigned index) :
        _parameters_data(parameters_data),
        _index(index)
    {
    }

    Parameter::Parameter(const Parameter & other) :
        _parameters_data(other._parameters_data),
        _index(other._index)
    {
    }

    Parameter::~Parameter()
    {
    }

    MutablePtr
    Parameter::clone() const
    {
        return MutablePtr(new Parameter(_parameters_data, _index));
    }

    Parameter::operator double () const
    {
        return _parameters_data->data[_index].value;
    }

    double
    Parameter::operator() () const
    {
        return _parameters_data->data[_index].value;
    }

    const Parameter &
    Parameter::operator= (const double & value)
    {
        _parameters_data->data[_index].value = value;

        return *this;
    }

    const double &
    Parameter::central() const
    {
        return _parameters_data->data[_index].central;
    }

    const double &
    Parameter::max() const
    {
        return _parameters_data->data[_index].max;
    }

    const double &
    Parameter::min() const
    {
        return _parameters_data->data[_index].min;
    }

    const std::string &
    Parameter::name() const
    {
        return _parameters_data->data[_index].name;
    }

    Parameter::Id
    Parameter::id() const
    {
        return _parameters_data->data[_index].id;
    }

    template class WrappedForwardIterator<ParameterUser::ConstIteratorTag, const Parameter::Id>;

    ParameterUser::ConstIterator
    ParameterUser::begin() const
    {
        return ConstIterator(_ids.cbegin());
    }

    ParameterUser::ConstIterator
    ParameterUser::end() const
    {
        return ConstIterator(_ids.cend());
    }

    void
    ParameterUser::uses(const Parameter::Id & id)
    {
        _ids.insert(id);
    }

    void
    ParameterUser::uses(const ParameterUser & other)
    {
        _ids.insert(other._ids.cbegin(), other._ids.cend());
    }

    UsedParameter::UsedParameter(const Parameter & parameter, ParameterUser & user) :
        Parameter(parameter)
    {
        user.uses(parameter.id());
    }

    UnknownParameterError::UnknownParameterError(const std::string & name) throw () :
        Exception("Unknown parameter: '" + name + "'")
    {
    }
}
