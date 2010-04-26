/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/parameters.hh>
#include <src/utils/private_implementation_pattern-impl.hh>

#include <map>
#include <vector>

namespace wf
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
            Parameters::NameValuePair{"c7", -0.33},
            Parameters::NameValuePair{"c8", -0.166},
            Parameters::NameValuePair{"c9", +4.23},
            Parameters::NameValuePair{"c10", -4.17},
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
            // Masses in GeV
            Parameters::NameValuePair{"mass::b(MSbar)", 4.13, 4.20, 4.37}, // cf. [PDG2008], p. 21
            Parameters::NameValuePair{"mass::c", 1.16, 1.27, 1.34}, // cf. [PDG2008], p. 21
            Parameters::NameValuePair{"mass::s", 0.00}, // we neglect m_s, cf. [BHvD2010], Eq. (??)
            Parameters::NameValuePair{"mass::t", 169.1, 171.2, 173.3}, // cf. [PDG2008], p. 21
            Parameters::NameValuePair{"mass::B0", 5.27920, 5.27953, 5.27986}, // cf. [PDG2008], p. 79
            Parameters::NameValuePair{"mass::K^*0", 0.896}, // cf. [PDG2008], p. 44
            Parameters::NameValuePair{"mass::W", 80.373, 80.398, 80.423}, // cf. [PDG2008], p. 8
            Parameters::NameValuePair{"mass::Z", 91.1855, 91.1876, 91.1897}, // cf. [PDG2008], p. 9
            // Form factor uncertainties
            Parameters::NameValuePair{"formfactors::a0_uncertainty", 0.85, 1.0, 1.15},
            Parameters::NameValuePair{"formfactors::a1_uncertainty", 0.85, 1.0, 1.15},
            Parameters::NameValuePair{"formfactors::a2_uncertainty", 0.85, 1.0, 1.15},
            Parameters::NameValuePair{"formfactors::v_uncertainty", 0.85, 1.0, 1.15},
            // B->K^*, K^* LCDA parameters
            Parameters::NameValuePair{"B->K^*::a_1_par", 0.03, 0.1, 0.17},
            Parameters::NameValuePair{"B->K^*::a_2_par", 0.0, 0.1, 0.2},
            Parameters::NameValuePair{"B->K^*::a_1_perp", 0.03, 0.1, 0.17},
            Parameters::NameValuePair{"B->K^*::a_2_perp", 0.0, 0.1, 0.2},
            // B->K^*ll uncertainties from subleading terms
            Parameters::NameValuePair{"B->K^*ll::A_0^L_uncertainty", 0.9, 1.0, 1.1},
            Parameters::NameValuePair{"B->K^*ll::A_0^R_uncertainty", 0.9, 1.0, 1.1},
            Parameters::NameValuePair{"B->K^*ll::A_par^L_uncertainty", 0.9, 1.0, 1.1},
            Parameters::NameValuePair{"B->K^*ll::A_par^R_uncertainty", 0.9, 1.0, 1.1},
            Parameters::NameValuePair{"B->K^*ll::A_perp^L_uncertainty", 0.9, 1.0, 1.1},
            Parameters::NameValuePair{"B->K^*ll::A_perp^R_uncertainty", 0.9, 1.0, 1.1},
            // Experimental Input
            Parameters::NameValuePair{"exp::BR(B->X_clnu)", 0.099, 0.101, 0.105}, // cf. [PDG2008], p. 82
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
