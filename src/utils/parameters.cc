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
            // Wilson coefficients at mu = 4.8 GeV to NNLL accuary, cf. ABBBSW2008, p. 6, Table 2.
            // We use C_{7,8} = C_{7,8}^eff - C'_{7,8}^eff.
            Parameters::NameValuePair{"c1", -0.257},
            Parameters::NameValuePair{"c2", +1.009},
            Parameters::NameValuePair{"c3", +0.005},
            Parameters::NameValuePair{"c4", -0.078},
            Parameters::NameValuePair{"c5", +0.000},
            Parameters::NameValuePair{"c6", +0.001},
            Parameters::NameValuePair{"c7", -0.298},
            Parameters::NameValuePair{"c8", -0.164},
            Parameters::NameValuePair{"c9", +4.211},
            Parameters::NameValuePair{"c10", -4.103},
            // Primed Wilson coefficients are negligible in the SM
            Parameters::NameValuePair{"c7prime", 0.0},
            Parameters::NameValuePair{"c9prime", 0.0},
            Parameters::NameValuePair{"c10prime", 0.0},
            // Factorization scale
            Parameters::NameValuePair{"mu", 4.8},
            // CKM matrix elements, cf. [PDG2008], Eqs. (11.4), (11.5), p. 169 and Eq. (11.26), p. 174
            Parameters::NameValuePair{"CKM::A", 0.793, 0.814, 0.835},
            Parameters::NameValuePair{"CKM::lambda", 0.2247, 0.2257, 0.2266},
            Parameters::NameValuePair{"CKM::rhobar", 0.119, 0.135, 0.166},
            Parameters::NameValuePair{"CKM::etabar", 0.332, 0.349, 0.364},
            // Masses in GeV
            Parameters::NameValuePair{"mass::b(MSbar)", 4.13, 4.20, 4.37}, // cf. [PDG2008], p. 21
            Parameters::NameValuePair{"mass::c", 1.16, 1.27, 1.34}, // cf. [PDG2008], p. 21
            Parameters::NameValuePair{"mass::B0", 5.27920, 5.27953, 5.27986}, // cf. [PDG2008], p. 79
            Parameters::NameValuePair{"mass::K^*0", 0.896}, // cf. [PDG2008], p. 44
            // Form factor uncertainties
            Parameters::NameValuePair{"formfactors::a0_uncertainty", 0.85, 1.0, 1.15},
            Parameters::NameValuePair{"formfactors::a1_uncertainty", 0.85, 1.0, 1.15},
            Parameters::NameValuePair{"formfactors::a2_uncertainty", 0.85, 1.0, 1.15},
            Parameters::NameValuePair{"formfactors::v_uncertainty", 0.85, 1.0, 1.15},
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

    UnknownParameterError::UnknownParameterError(const std::string & name) throw () :
        Exception("Unknown parameter: '" + name + "'")
    {
    }
}
