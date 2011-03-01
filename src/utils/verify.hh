/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_SRC_UTILS_VERIFY_HH
#define EOS_GUARD_SRC_UTILS_VERIFY_HH 1

#include <src/utils/exception.hh>

namespace eos
{
    class VerifiedRangeError :
        public Exception
    {
        protected:
            VerifiedRangeError(const std::string & message);
    };

    struct VerifiedRangeOverflow :
        public VerifiedRangeError
    {
        VerifiedRangeOverflow(const std::string & value, const std::string & maximum);
    };

    struct VerifiedRangeUnderflow :
        public VerifiedRangeError
    {
        VerifiedRangeUnderflow(const std::string & value, const std::string & minimum);
    };

    template <typename T_>
    class VerifiedRange
    {
        private:
            T_ _min;

            T_ _max;

            T_ _value;

            const T_ & _verify(const T_ & t)
            {
                if (t < _min)
                    throw VerifiedRangeUnderflow(stringify(t), stringify(_min));

                if (t > _max)
                    throw VerifiedRangeOverflow(stringify(t), stringify(_max));

                return t;
            }

        public:
            VerifiedRange(const T_ & min, const T_ & max, const T_ & value) :
                _min(min),
                _max(max),
                _value(_verify(value))
            {
            }

            operator T_ () const
            {
                return _value;
            }

            const T_ & operator= (const T_ & value)
            {
                _value = _verify(value);
            }
    };
}

#endif
