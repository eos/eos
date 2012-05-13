/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013 Danny van Dyk
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
            Parameter::Template{"hbar",                                          +6.58211883e-25, +6.58211899e-25, +6.58211915e-25}, // GeV s, cf. [PDG2010]:  p. 4

            // Wilson coefficients C1 - C6 at mu = 4.2 GeV to NNLL accuary as calculated in the SM by EOS's StandardModel class.
            // For the calculations, cf. [BMU1999].
            Parameter::Template{"c1",                                            -0.28768333,     -0.28768333,     -0.28768333    },
            Parameter::Template{"c2",                                            +1.01013250,     +1.01013250,     +1.01013250    },
            Parameter::Template{"c3",                                            -0.00600697,     -0.00600697,     -0.00600697    },
            Parameter::Template{"c4",                                            -0.08597076,     -0.08597076,     -0.08597076    },
            Parameter::Template{"c5",                                            +0.00041824,     +0.00041824,     +0.00041824    },
            Parameter::Template{"c6",                                            +0.00112410,     +0.00112410,     +0.00112410    },
            Parameter::Template{"Abs{c7}",                                       +0.0,            +0.33613067,     +1.0           },
            Parameter::Template{"Arg{c7}",                                       +0.0,            +M_PI,           +2.0 * M_PI    },
            Parameter::Template{"Re{c7}",                                        -1.0,            -0.33613067,     +1.0           },
            Parameter::Template{"Im{c7}",                                        -1.0,            +0.0,            +1.0           },
            Parameter::Template{"c8",                                            -0.18205267,     -0.18205267,     -0.18205267    },
            Parameter::Template{"Abs{c9}",                                       +0.0,            +4.27450580,    +15.0           },
            Parameter::Template{"Arg{c9}",                                       +0.0,            +0.0,            +2.0 * M_PI    },
            Parameter::Template{"Re{c9}",                                       -15.0,            +4.27450580,    +15.0           },
            Parameter::Template{"Im{c9}",                                       -15.0,            +0.0,           +15.0           },
            Parameter::Template{"Abs{c10}",                                      +0.0,            +4.16020202,    +15.0           },
            Parameter::Template{"Arg{c10}",                                      +0.0,            +M_PI,           +2.0 * M_PI    },
            Parameter::Template{"Re{c10}",                                      -15.0,            -4.16020202,    +15.0           },
            Parameter::Template{"Im{c10}",                                      -15.0,            +0.0,           +15.0           },
            // Primed Wilson coefficients are negligible in the SM
            Parameter::Template{"Abs{c7'}",                                      +0.0,            +0.0,            +1.0           },
            Parameter::Template{"Arg{c7'}",                                      +0.0,            +0.0,            +2.0 * M_PI    },
            Parameter::Template{"Re{c7'}",                                       -1.0,            +0.0,            +1.0           },
            Parameter::Template{"Im{c7'}",                                       -1.0,            +0.0,            +1.0           },
            Parameter::Template{"Abs{c9'}",                                      +0.0,            +0.0,           +15.0           },
            Parameter::Template{"Arg{c9'}",                                      +0.0,            +0.0,            +2.0 * M_PI    },
            Parameter::Template{"Re{c9'}",                                      -15.0,            +0.0,           +15.0           },
            Parameter::Template{"Im{c9'}",                                      -15.0,            +0.0,           +15.0           },
            Parameter::Template{"Abs{c10'}",                                     +0.0,            +0.0,           +15.0           },
            Parameter::Template{"Arg{c10'}",                                     +0.0,            +0.0,            +2.0 * M_PI    },
            Parameter::Template{"Re{c10'}",                                     -15.0,            +0.0,           +15.0           },
            Parameter::Template{"Im{c10'}",                                     -15.0,            +0.0,           +15.0           },
            // Factorization scale
            Parameter::Template{"mu",                                            +2.4,            +4.2,            +9.6           },
            // GSW Parameter
            Parameter::Template{"GSW::sin^2(theta)",                             +0.23103,        +0.23116,        +0.23129       },
            // Wolfenstein parameters of CKM, cf. [CKMfitter04] Table 2, p. 48 and ICHEP10 result s +/- 1 sigma
            Parameter::Template{"CKM::A",                                        +0.785,          +0.812,          +0.825         },
            Parameter::Template{"CKM::lambda",                                   +0.22466,        +0.22543,        +0.22620       },
            Parameter::Template{"CKM::rhobar",                                   +0.119,          +0.144,          +0.169         },
            Parameter::Template{"CKM::etabar",                                   +0.327,          +0.342,          +0.358         },
            // QED inputs
            Parameter::Template{"QED::alpha_e(m_b)",                             +1.0/133.0,      +1.0/133.0,      +1.0/128.0     }, // alpha_e(m_b) .. alpha_e(m_W)
            // QCD inputs
            Parameter::Template{"QCD::alpha_s(MZ)",                              +0.117620,       +0.117620,       +0.117620      },
            Parameter::Template{"QCD::mu_t",                                     +170.0,          +170.0,          +170.0         },
            Parameter::Template{"QCD::mu_b",                                     +4.2,            +4.2,            +4.2           },
            Parameter::Template{"QCD::mu_c",                                     +1.0,            +1.0,            +1.0           },
            Parameter::Template{"QCD::Lambda",                                   +0.5,            +0.5,            +0.5           },
            // G_Fermi
            Parameter::Template{"G_Fermi",                                       +1.16636e-5,     +1.16637e-5,     +1.16638e-5    }, // cf. [PDG2010], p. 5

            /* Masses in GeV */
            // Lepton masses
            Parameter::Template{"mass::e",                                       +5.10999e-4,     +5.10999e-4,     +5.10999e-4    }, // cf. [PDG2010], p. 13
            Parameter::Template{"mass::mu",                                      +1.05658e-1,     +1.05658e-1,     +1.05658e-1    }, // cf. [PDG2010], p. 13
            Parameter::Template{"mass::tau",                                     +1.77666,        +1.77682,        +1.77698       }, // cf. [PDG2010], p. 14
            // Quark masses
            Parameter::Template{"mass::s(2GeV)",                                 +0.080,          +0.101,          +0.130         }, // cf. [PDG2010], p. 21
            Parameter::Template{"mass::c",                                       +1.18,           +1.27,           +1.34          }, // cf. [PDG2010], p. 21
            Parameter::Template{"mass::b(MSbar)",                                +4.13,           +4.19,           +4.37          }, // cf. [PDG2010], p. 21
            Parameter::Template{"mass::t(pole)",                                 +172.0,          +173.3,          +174.6         }, // cf. [TEWWG2010], p. 21
            // K meson masses
            Parameter::Template{"mass::K0",                                      +0.497590,       +0.497614,       +0.497638      }, // cf. [PDG2010], p. 41
            Parameter::Template{"mass::K^*0",                                    +0.89572,        +0.89594,        +0.89616       }, // cf. [PDG2010], p. 45
            // B meson masses
            Parameter::Template{"mass::B_d",                                     +5.2792,         +5.2795,         +5.2798        }, // cf. [PDG2010], p. 82
            Parameter::Template{"mass::B_u",                                     +5.27888,        +5.27917,        +5.27946       }, // cf. [PDG2010], p. 69
            Parameter::Template{"mass::B_s",                                     +5.3657,         +5.3663,         +5.3669        }, // cf. [PDG2010], p. 103
            // Gauge boson masses
            Parameter::Template{"mass::W",                                       +80.376,         +80.399,         +80.422        }, // cf. [PDG2010], p. 8
            Parameter::Template{"mass::Z",                                       +91.1855,        +91.1876,        +91.1897       }, // cf. [PDG2010], p. 9

            /* Decay constants */
            Parameter::Template{"decay-constant::B_d",                           +0.200,          +0.212,          +0.224         }, // GeV, cf. [S2010], p. 3
            Parameter::Template{"decay-constant::B_u",                           +0.200,          +0.212,          +0.224         }, // GeV, cf. [S2010], p. 3
            Parameter::Template{"decay-constant::B_s",                           +0.244,          +0.256,          +0.268         }, // GeV, cf. [S2010], p. 3
            Parameter::Template{"decay-constant::K_d",                           +0.1549,         +0.1561,         +0.1573        }, // GeV, cf. [PDGBOOK2010], p. 864, Eq. (7)
            Parameter::Template{"decay-constant::K_u",                           +0.1549,         +0.1561,         +0.1573        }, // GeV, cf. [PDGBOOK2010], p. 864, Eq. (7)
            /* b->s matching parameters */
            Parameter::Template{"b->s::mu_0c",                                   +80.0,           +80.0,           +80.0          },
            Parameter::Template{"b->s::mu_0t",                                   +120.0,          +120.0,          +120.0         },

            /* Mean life times */
            Parameter::Template{"life_time::B_d",                                +1.516e-12,      +1.525e-12,      +1.534e-12     }, // cf. [PDG2010], p. 82
            Parameter::Template{"life_time::B_u",                                +1.627e-12,      +1.638e-12,      +1.649e-12     }, // cf. [PDG2010], p. 69
            Parameter::Template{"life_time::B_s",                                +1.446e-12,      +1.472e-12,      +1.496e-12     }, // cf. [PDG2010], p. 103

            // Form factor uncertainties
            Parameter::Template{"formfactors::xi_perp_uncertainty",              +0.89,           +1.0,            +1.11          },
            Parameter::Template{"formfactors::xi_par_uncertainty",               +0.86,           +1.0,            +1.14          },

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

            // form factor parameters for B->K according to [BZ2004v2] (approximate)
            Parameter::Template{"B->K::fp_uncertainty@BZ2004v2",                 +0.85,           +1.0,            +1.15          },
            Parameter::Template{"B->K::f0_uncertainty@BZ2004v2",                 +0.85,           +1.0,            +1.15          },
            Parameter::Template{"B->K::ft_uncertainty@BZ2004v2",                 +0.85,           +1.0,            +1.15          },

            // form factor parameters for B->K according to [KMPW2010], Table 4, p. 31
            Parameter::Template{"B->K::F^p(0)@KMPW2010",                         +0.32,           +0.34,           +0.39          },
            Parameter::Template{"B->K::F^0(0)@KMPW2010",                         +0.32,           +0.34,           +0.39          },
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
            // B->Pll uncertainties at subleasing order at Large Recoil
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
            // B->X_s HQET parameters
            Parameter::Template{"B->X_s::lambda_1",                              -0.20,           -0.20,           -0.20          }, // cf. [ALGH2001], Table 2, p. 13
            Parameter::Template{"B->X_s::lambda_2",                              +0.12,           +0.12,           +0.12          }, // cf. [ALGH2001], Table 2, p. 13
            // Heavy Quark Expansion parameters for hadronic matrix elements ~ <B|O|B>
            Parameter::Template{"B->B::mu_pi^2@1GeV",                            +0.35,           +0.45,           +0.55          }, // cf. [BBMU2003], Eq. (19), p. 9
            Parameter::Template{"B->B::mu_G^2@1GeV",                             +0.33,           +0.35,           +0.38          }, // cf. [BBMU2003], Eq. (17), p. 9
            Parameter::Template{"B->B::rho_D^3@1GeV",                            +0.10,           +0.20,           +0.30          }, // cf. [BBMU2003], between Eqs. (19),(20), p. 9
            Parameter::Template{"B->B::rho_LS^3@1GeV",                           -0.30,           -0.15,           -0.00          }, // cf. [BBMU2003], Eq. (20), p. 9
            // B->X_s gamma SM theory uncertainty
            Parameter::Template{"B->X_sgamma::uncertainty",                      -1.0,            +0.0,            +1.0           },
            // Experimental Input
            Parameter::Template{"exp::BR(B->X_clnu)",                            +0.1042,         +0.1057,         +0.1072        }, // cf. [PDG2008], p. 82
            Parameter::Template{"exp::C(B->X_clnu, B->X_ulnu)",                  +0.57,           +0.58,           +0.59          },
            Parameter::Template{"exp::CKM(B->X_sll, B->X_clnu)",                 +0.975218,       +0.98549,        +0.995277      },
            // Parametrise unknown admixture of l=e, l=mu in B->X_sll
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
