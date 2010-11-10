/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/kinematic.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/stringify.hh>

#include <map>

namespace eos
{
    template <>
    struct Implementation<Kinematics>
    {
        std::map<std::string, double> variables;
    };

    Kinematics::Kinematics() :
        PrivateImplementationPattern<Kinematics>(new Implementation<Kinematics>)
    {
    }

    Kinematics::~Kinematics()
    {
    }

    double
    Kinematics::operator[] (const std::string & variable) const
    {
        auto i(_imp->variables.find(variable));

        if (_imp->variables.end() == i)
            throw UnknownKinematicVariableError(variable);

        return i->second;
    }

    std::string
    Kinematics::as_string () const
    {
        std::string result;

        auto i(_imp->variables.cbegin()), i_end(_imp->variables.cend());
        if (i != i_end)
        {
            result = i->first + '=' + stringify(i->second);
            ++i;
        }


        for ( ; i != i_end ; ++i)
        {
            result += ", " + i->first + '=' + stringify(i->second);
        }

        return result;
    }

    void
    Kinematics::declare(const std::string & variable)
    {
        auto i(_imp->variables.find(variable));

        if (_imp->variables.end() == i)
            _imp->variables[variable] = 0.0;
    }

    void
    Kinematics::set(const std::string & variable, const double & value)
    {
        auto i(_imp->variables.find(variable));

        if (_imp->variables.end() == i)
            throw UnknownKinematicVariableError(variable);

        i->second = value;
    }

    UnknownKinematicVariableError::UnknownKinematicVariableError(const std::string & variable) throw () :
        Exception("Unknown kinematic variable: '" + variable + "'")
    {
    }
}
