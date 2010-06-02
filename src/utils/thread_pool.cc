/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/condition_variable.hh>
#include <src/utils/instantiation_policy-impl.hh>
#include <src/utils/lock.hh>
#include <src/utils/mutex.hh>
#include <src/utils/private_implementation_pattern-impl.hh>
#include <src/utils/thread.hh>
#include <src/utils/thread_pool.hh>

#include <list>

#include <unistd.h>

#include <iostream>

namespace wf
{
    template <>
    struct Implementation<ThreadPool>
    {
        unsigned number_of_threads;

        // Thread termination
        Mutex * const terminate_mutex;

        bool terminate;

        // Job handling
        Mutex * const job_mutex;

        ConditionVariable * const job_arrival;

        unsigned waiting_jobs;
        unsigned done_jobs;

        std::list<std::pair<Ticket, std::tr1::function<void (void)> *>> queue;

        std::list<Thread *> threads;

        void thread_function()
        {
            std::tr1::function<void (void)> * job;
            Ticket ticket;

            do
            {
                {
                    Lock l(*terminate_mutex);
                    if (terminate)
                        break;
                }

                job = 0;

                {
                    Lock l(*job_mutex);

                    if (queue.empty())
                    {
                        waiting_jobs += 1;
                        job_arrival->wait(*job_mutex);
                        waiting_jobs -= 1;
                    }

                    if (queue.empty())
                    {
                        Lock l(*terminate_mutex);
                        if (terminate)
                            break;
                        else
                            continue;
                    }

                    ticket = queue.front().first;
                    job = queue.front().second;
                    queue.pop_front();
                }

                (*job)();
                ticket.mark();
            }
            while (true);
        }

        Implementation() :
            terminate_mutex(new Mutex),
            terminate(false),
            job_mutex(new Mutex),
            job_arrival(new ConditionVariable),
            waiting_jobs(0),
            done_jobs(0),
            number_of_threads(sysconf(_SC_NPROCESSORS_CONF))
        {
            for (int i(0) ; i < number_of_threads ; ++i)
            {
                threads.push_back(new Thread(std::tr1::bind(&Implementation<ThreadPool>::thread_function, this)));
            }
        }

        ~Implementation()
        {
            {
                Lock l(*terminate_mutex);
                terminate = true;
            }

            {
                Lock l(*job_mutex);
                job_arrival->broadcast();
            }

            for (auto t(threads.begin()), t_end(threads.end()) ; t != t_end ; ++t)
            {
                delete *t;
            }
        }
    };

    ThreadPool::ThreadPool() :
        InstantiationPolicy<ThreadPool, Singleton>(),
        PrivateImplementationPattern<ThreadPool>(new Implementation<ThreadPool>)
    {
    }

    ThreadPool::~ThreadPool()
    {
    }

    Ticket
    ThreadPool::enqueue(const std::tr1::function<void (void)> & job)
    {
        std::pair<Ticket, std::tr1::function<void (void)> *> item(Ticket(), new std::tr1::function<void (void)>(job));

        {
            Lock l(*_imp->job_mutex);
            _imp->queue.push_back(item);

            if (_imp->waiting_jobs > 0)
                _imp->job_arrival->signal();
        }

        return item.first;
    }

    ThreadPool *
    ThreadPool::instance()
    {
        return InstantiationPolicy<ThreadPool, Singleton>::instance();
    }
}

