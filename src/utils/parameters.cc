/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/parameters.hh>
#include <src/utils/private_implementation_pattern-impl.hh>

#include <cmath>
#include <map>
#include <vector>

namespace eos
{
    template <>
    struct Implementation<Parameters>
    {
        std::vector<Parameters::NameValuePair> parameters;

        std::map<std::string, unsigned> parameter_map;

        Implementation(const std::initializer_list<Parameters::NameValuePair> & list) :
            parameters(list)
        {
            unsigned idx(0);
            for (auto i(list.begin()), i_end(list.end()) ; i != i_end ; ++i, ++idx)
            {
                parameter_map[i->name] = idx;
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
        auto i(_imp->parameter_map.find(name));

        if (_imp->parameter_map.end() == i)
            throw UnknownParameterError(name);

        return Parameter(_imp, i->second);
    }

    void
    Parameters::set(const std::string & name, const double & value)
    {
        auto i(_imp->parameter_map.find(name));

        if (_imp->parameter_map.end() == i)
            throw UnknownParameterError(name);

        _imp->parameters[i->second].value = value;
    }

    Parameters
    Parameters::FromList(const std::initializer_list<Parameters::NameValuePair> & list)
    {
        return Parameters(new Implementation<Parameters>(list));
    }

    Parameters
    Parameters::Defaults()
    {
        return Parameters::FromList({
            // Wilson coefficients C1 - C6 at mu = 4.2 GeV to NLL accuary, based on [CMM1997]
            Parameters::NameValuePair{"c1", -0.323},
            Parameters::NameValuePair{"c2", +1.00931},
            Parameters::NameValuePair{"c3", -0.00522869},
            Parameters::NameValuePair{"c4", -0.0879473},
            Parameters::NameValuePair{"c5", +0.000374755},
            Parameters::NameValuePair{"c6", +0.00105859},
            // Wilson coefficients C7 - c10 at mu = 4.2 GeV to NNLL, based on ?
            Parameters::NameValuePair{"Abs{c7}",  0.331}, // c7eff = -0.306
            Parameters::NameValuePair{"Arg{c7}",  M_PI},
            Parameters::NameValuePair{"c8", -0.181}, // c8eff = -0.168
            Parameters::NameValuePair{"Abs{c9}",  4.27},
            Parameters::NameValuePair{"Arg{c9}",  0.00},
            Parameters::NameValuePair{"Abs{c10}", 4.17},
            Parameters::NameValuePair{"Arg{c10}", M_PI},
            // Primed Wilson coefficients are negligible in the SM
            Parameters::NameValuePair{"c7prime", 0.0},
            Parameters::NameValuePair{"c9prime", 0.0},
            Parameters::NameValuePair{"c10prime", 0.0},
            // Factorization scale
            Parameters::NameValuePair{"mu", 2.4, 4.2, 9.6},
            // CKM matrix elements, cf. [PDG2008], Eqs. (11.4), (11.5), p. 169 and Eq. (11.26), p. 174
            Parameters::NameValuePair{"CKM::A", 0.793, 0.814, 0.835},
            Parameters::NameValuePair{"CKM::lambda", 0.2247, 0.2257, 0.2266},
            Parameters::NameValuePair{"CKM::rhobar", 0.119, 0.135, 0.166},
            Parameters::NameValuePair{"CKM::etabar", 0.332, 0.349, 0.364},
            Parameters::NameValuePair{"CKM::|V_cb|", 0.0404, 0.0417, 0.0430},
            // QCD inputs
            Parameters::NameValuePair{"QCD::alpha_s(MZ)", 0.117620},
            Parameters::NameValuePair{"QCD::mu_t", 170.0},
            Parameters::NameValuePair{"QCD::mu_b",   4.2},
            Parameters::NameValuePair{"QCD::mu_c",   1.0},
            // Masses in GeV
            Parameters::NameValuePair{"mass::b(MSbar)", 4.13, 4.20, 4.37}, // cf. [PDG2008], p. 21
            Parameters::NameValuePair{"mass::c", 1.16, 1.27, 1.34}, // cf. [PDG2008], p. 21
            Parameters::NameValuePair{"mass::s", 0.00}, // we neglect m_s, cf. [BHvD2010], Eq. (??)
            Parameters::NameValuePair{"mass::t", 169.1, 171.2, 173.3}, // cf. [PDG2008], p. 21
            Parameters::NameValuePair{"mass::tau", 1.77684}, // cf. [PDG2008], p. 14
            Parameters::NameValuePair{"mass::B0", 5.27920, 5.27953, 5.27986}, // cf. [PDG2008], p. 79
            Parameters::NameValuePair{"mass::K^*0", 0.896}, // cf. [PDG2008], p. 44
            Parameters::NameValuePair{"mass::W", 80.373, 80.398, 80.423}, // cf. [PDG2008], p. 8
            Parameters::NameValuePair{"mass::Z", 91.1855, 91.1876, 91.1897}, // cf. [PDG2008], p. 9
            // Form factor uncertainties
            Parameters::NameValuePair{"formfactors::a0_uncertainty", 0.85, 1.0, 1.15},
            Parameters::NameValuePair{"formfactors::a1_uncertainty", 0.85, 1.0, 1.15},
            Parameters::NameValuePair{"formfactors::a2_uncertainty", 0.85, 1.0, 1.15},
            Parameters::NameValuePair{"formfactors::v_uncertainty",  0.85, 1.0, 1.15},
            Parameters::NameValuePair{"formfactors::xi_perp_uncertainty", 0.89, 1.0, 1.11},
            Parameters::NameValuePair{"formfactors::xi_par_uncertainty",  0.86, 1.0, 1.14},
            // B LCDA parameters
            Parameters::NameValuePair{"f_B",        0.17,  0.20,  0.23}, // GeV, cf. [BHvD2010], Table I
            Parameters::NameValuePair{"lambda_B_p", 0.370, 0.485, 0.600}, // GeV, cf. [BHvD2010], Table I
            // B->K^*, K^* LCDA parameters
            Parameters::NameValuePair{"B->K^*::a_1_par",   0.03,  0.1,   0.17},
            Parameters::NameValuePair{"B->K^*::a_2_par",   0.0,   0.1,   0.2},
            Parameters::NameValuePair{"B->K^*::a_1_perp",  0.03,  0.1,   0.17},
            Parameters::NameValuePair{"B->K^*::a_2_perp",  0.0,   0.1,   0.2},
            Parameters::NameValuePair{"B->K^*::f_Kstar_par",       0.212, 0.217, 0.222}, // GeV, cf. [BHvD2010], Table I
            Parameters::NameValuePair{"B->K^*::f_Kstar_perp@2GeV", 0.168, 0.173, 0.178}, // GeV @2 Gev, 0.185 +/-0.005 GeV, cf. [BHvD2010], Table I
            // B->K^*ll uncertainties from subleading terms for Large Recoil
            Parameters::NameValuePair{"B->K^*ll::A_0^L_uncertainty@LargeRecoil",    0.95, 1.0, 1.05},
            Parameters::NameValuePair{"B->K^*ll::A_0^R_uncertainty@LargeRecoil",    0.95, 1.0, 1.05},
            Parameters::NameValuePair{"B->K^*ll::A_par^L_uncertainty@LargeRecoil",  0.95, 1.0, 1.05},
            Parameters::NameValuePair{"B->K^*ll::A_par^R_uncertainty@LargeRecoil",  0.95, 1.0, 1.05},
            Parameters::NameValuePair{"B->K^*ll::A_perp^L_uncertainty@LargeRecoil", 0.95, 1.0, 1.05},
            Parameters::NameValuePair{"B->K^*ll::A_perp^R_uncertainty@LargeRecoil", 0.95, 1.0, 1.05},
             // B->K^*ll uncertainties from subleading terms for Low Recoil
            Parameters::NameValuePair{"B->K^*ll::A_0^L_uncertainty@LowRecoil",    0.95, 1.0, 1.05},
            Parameters::NameValuePair{"B->K^*ll::A_0^R_uncertainty@LowRecoil",    0.95, 1.0, 1.05},
            Parameters::NameValuePair{"B->K^*ll::A_par^L_uncertainty@LowRecoil",  0.95, 1.0, 1.05},
            Parameters::NameValuePair{"B->K^*ll::A_par^R_uncertainty@LowRecoil",  0.95, 1.0, 1.05},
            Parameters::NameValuePair{"B->K^*ll::A_perp^L_uncertainty@LowRecoil", 0.95, 1.0, 1.05},
            Parameters::NameValuePair{"B->K^*ll::A_perp^R_uncertainty@LowRecoil", 0.95, 1.0, 1.05},
            // B->K^*ll uncertainties from improved Isgur-Wise relations at Low Recoil
            Parameters::NameValuePair{"B->K^*ll::IW_long_uncertainty", 0.8, 1.0, 1.2},
            Parameters::NameValuePair{"B->K^*ll::IW_par_uncertainty",  0.8, 1.0, 1.2},
            Parameters::NameValuePair{"B->K^*ll::IW_perp_uncertainty", 0.8, 1.0, 1.2},
            // B->X_s HQET parameters
            Parameters::NameValuePair{"B->X_s::lambda_1", -0.20}, // cf. [ALGH2001], Table 2, p. 13
            Parameters::NameValuePair{"B->X_s::lambda_2", +0.12}, // cf. [ALGH2001], Table 2, p. 13
            // B->X_s gamma SM theory uncertainty
            Parameters::NameValuePair{"B->X_sgamma::uncertainty", -1.0, 0.0, +1.0},
            // Experimental Input
            Parameters::NameValuePair{"exp::BR(B->X_clnu)", 0.1042, 0.1057, 0.1072}, // cf. [PDG2008], p. 82
            Parameters::NameValuePair{"exp::C(B->X_clnu, B->X_ulnu)", 0.57, 0.58, 0.59},
            Parameters::NameValuePair{"exp::CKM(B->X_sll, B->X_clnu)", 0.975218, 0.98549, 0.995277},
        });
    }

    Parameter::Parameter(const std::tr1::shared_ptr<Implementation<Parameters>> & imp, unsigned index) :
        _imp(imp),
        _index(index)
    {
    }

    Parameter::~Parameter()
    {
    }

    Parameter::operator double () const
    {
        return _imp->parameters[_index].value;
    }

    double
    Parameter::operator() () const
    {
        return _imp->parameters[_index].value;
    }

    const Parameter &
    Parameter::operator= (const double & value)
    {
        _imp->parameters[_index].value = value;

        return *this;
    }

    const double &
    Parameter::max() const
    {
        return _imp->parameters[_index].max;
    }

    const double &
    Parameter::min() const
    {
        return _imp->parameters[_index].min;
    }

    const std::string &
    Parameter::name() const
    {
        return _imp->parameters[_index].name;
    }

    UnknownParameterError::UnknownParameterError(const std::string & name) throw () :
        Exception("Unknown parameter: '" + name + "'")
    {
    }
}
