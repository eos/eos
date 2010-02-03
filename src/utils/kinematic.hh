/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_UTILS_KINEMATIC_HH
#define WFITTER_GUARD_SRC_UTILS_KINEMATIC_HH 1

#include <src/utils/exception.hh>
#include <src/utils/private_implementation_pattern.hh>

namespace wf
{
    struct UnknownKinematicVariableError :
        public Exception
    {
        UnknownKinematicVariableError(const std::string & variable) throw ();
    };

    class Kinematics :
        public PrivateImplementationPattern<Kinematics>
    {
        public:
            Kinematics();

            ~Kinematics();

            double operator[] (const std::string & variable) const;

            void declare(const std::string & variable);

            void set(const std::string & variable, const double & value);
    };

    inline double lambda(const double & a, const double & b, const double & c)
    {
        return a * a + b * b + c * c - 2.0 * (a * b + a * c + b * c);
    }
}

#endif
