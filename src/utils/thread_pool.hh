/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef WFITTER_GUARD_SRC_UTILS_THREAD_POOL_HH
#define WFITTER_GUARD_SRC_UTILS_THREAD_POOL_HH 1

#include <src/utils/instantiation_policy.hh>
#include <src/utils/private_implementation_pattern.hh>
#include <src/utils/ticket.hh>

#include <tr1/functional>

namespace wf
{
    class ThreadPool :
        public InstantiationPolicy<ThreadPool, Singleton>,
        public PrivateImplementationPattern<ThreadPool>
    {
        public:
            ThreadPool();

            ~ThreadPool();

            Ticket enqueue(const std::tr1::function<void (void)> & work);

            static ThreadPool * instance();
    };
}

#endif
