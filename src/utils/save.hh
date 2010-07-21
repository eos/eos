/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_UTILS_SAVE_HH
#define WFITTER_GUARD_SRC_UTILS_SAVE_HH 1

namespace wf
{
    template <typename T_> class Save
    {
        private:
            T_ & _variable;

            T_ _original;

        public:
            Save(T_ & variable, const T_ & replacement) :
                _variable(variable),
                _original(variable)
            {
                _variable = replacement;
            }

            ~Save()
            {
                _variable = _original;
            }
    };
}

#endif
