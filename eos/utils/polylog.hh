/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef EOS_GUARD_EOS_UTILS_POLYLOG_HH
#define EOS_GUARD_EOS_UTILS_POLYLOG_HH 1

#include <eos/utils/complex.hh>

namespace eos
{
    complex<double> dilog(const complex<double> & z) __attribute__ ((pure));

    complex<double> trilog(const complex<double> & z) __attribute__ ((pure));
}

#endif
