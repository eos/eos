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
        std::vector<double> parameters;

        std::map<std::string, unsigned> parameter_map;

        Implementation(const std::initializer_list<Parameters::NameValuePair> & list)
        {
            for (auto i(list.begin()), i_end(list.end()) ; i != i_end ; ++i)
            {
                parameter_map[i->first] = parameters.size();
                parameters.push_back(i->second);
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

    Parameter
    Parameters::declare(const std::string & name, const double & value)
    {
        auto i(_imp->parameter_map.find(name));

        if (_imp->parameter_map.end() != i) // TODO: Should this throw an exception instead?
        {
            _imp->parameters[i->second] = value;
            return Parameter(_imp, i->second);
        }
        else
        {
            _imp->parameter_map[name] = _imp->parameters.size();
            Parameter result(_imp, _imp->parameters.size());
            _imp->parameters.push_back(value);
            return result;
        }
    }

    void
    Parameters::set(const std::string & name, const double & value)
    {
        auto i(_imp->parameter_map.find(name));

        if (_imp->parameter_map.end() == i)
            throw UnknownParameterError(name);

        _imp->parameters[i->second] = value;
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
            Parameters::NameValuePair{"CKM::A", 0.814},
            Parameters::NameValuePair{"CKM::lambda", 0.2257},
            Parameters::NameValuePair{"CKM::rhobar", 0.135},
            Parameters::NameValuePair{"CKM::etabar", 0.349},
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
        return _imp->parameters[_index];
    }

    double
    Parameter::operator() () const
    {
        return _imp->parameters[_index];
    }

    const Parameter &
    Parameter::operator= (const double & value)
    {
        _imp->parameters[_index] = value;
    }

    UnknownParameterError::UnknownParameterError(const std::string & name) throw () :
        Exception("Unknown parameter: '" + name + "'")
    {
    }
}
