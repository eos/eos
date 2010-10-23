/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/kinematic.hh>
#include <src/utils/private_implementation_pattern-impl.hh>

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
