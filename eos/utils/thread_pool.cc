/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2021 Danny van Dyk
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <eos/utils/condition_variable.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/instantiation_policy-impl.hh>
#include <eos/utils/lock.hh>
#include <eos/utils/mutex.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/thread.hh>
#include <eos/utils/thread_pool.hh>

#include <list>
#include <unistd.h>

namespace eos
{
    template <> struct Implementation<ThreadPool>
    {
            unsigned      number_of_threads;
            unsigned long nominal_capacity;
            unsigned long stop_capacity;

            // Thread termination
            Mutex * const terminate_mutex;

            bool terminate;

            // Job handling
            Mutex * const job_mutex;

            ConditionVariable * const job_arrival;
            ConditionVariable * const job_capacity;

            unsigned long waiting_for_jobs;
            unsigned long pending_jobs;

            std::list<std::pair<Ticket, std::function<void(void)> *>> queue;

            std::list<Thread *> threads;

            void
            thread_function()
            {
                std::function<void(void)> * job;
                Ticket                      ticket;

                do
                {
                    {
                        Lock l(*terminate_mutex);
                        if (terminate)
                        {
                            break;
                        }
                    }

                    job = 0;

                    {
                        Lock l(*job_mutex);

                        if (queue.empty())
                        {
                            waiting_for_jobs += 1;
                            job_arrival->wait(*job_mutex);
                            waiting_for_jobs -= 1;
                        }

                        if (queue.empty())
                        {
                            Lock l(*terminate_mutex);
                            if (terminate)
                            {
                                break;
                            }
                            else
                            {
                                continue;
                            }
                        }

                        ticket = queue.front().first;
                        job    = queue.front().second;
                        queue.pop_front();
                    }

                    (*job)();
                    {
                        Lock l(*job_mutex);
                        pending_jobs -= 1;

                        if (pending_jobs == nominal_capacity)
                        {
                            job_capacity->signal();
                        }
                    }
                    ticket.mark();
                }
                while (true);
            }

            static unsigned
            _number_of_threads()
            {
                // by default, use as many threads as configured processors available to the system
                unsigned result = sysconf(_SC_NPROCESSORS_CONF);

                // limit number of threads to user configured value
                const char * env_max_threads = std::getenv("EOS_MAX_THREADS");
                if (env_max_threads)
                {
                    unsigned max_threads = destringify<unsigned>(env_max_threads);
                    result               = std::min(result, max_threads);
                }

                return result;
            }

            Implementation() :
                number_of_threads(_number_of_threads()),
                nominal_capacity(number_of_threads * 10),
                stop_capacity(nominal_capacity * 2),
                terminate_mutex(new Mutex),
                terminate(false),
                job_mutex(new Mutex),
                job_arrival(new ConditionVariable),
                job_capacity(new ConditionVariable),
                waiting_for_jobs(0),
                pending_jobs(0)
            {
                for (unsigned i(0); i < number_of_threads; ++i)
                {
                    threads.push_back(new Thread(std::bind(&Implementation<ThreadPool>::thread_function, this)));
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

                for (auto & thread : threads)
                {
                    delete thread;
                }
            }
    };

    ThreadPool::ThreadPool() :
        InstantiationPolicy<ThreadPool, Singleton>(),
        PrivateImplementationPattern<ThreadPool>(new Implementation<ThreadPool>)
    {
    }

    ThreadPool::~ThreadPool() {}

    Ticket
    ThreadPool::enqueue(const std::function<void(void)> & job)
    {
        std::pair<Ticket, std::function<void(void)> *> item(Ticket(), new std::function<void(void)>(job));

        {
            Lock l(*_imp->job_mutex);
            _imp->queue.push_back(item);
            _imp->pending_jobs += 1;

            if (_imp->waiting_for_jobs > 0)
            {
                _imp->job_arrival->signal();
            }
        }

        return item.first;
    }

    ThreadPool *
    ThreadPool::instance()
    {
        return InstantiationPolicy<ThreadPool, Singleton>::instance();
    }

    void
    ThreadPool::wait_for_free_capacity()
    {
        Lock l(*_imp->job_mutex);

        if (_imp->pending_jobs < _imp->stop_capacity)
        {
            return;
        }

        _imp->job_capacity->wait(*_imp->job_mutex);
    }

    unsigned
    ThreadPool::number_of_threads() const
    {
        return _imp->number_of_threads;
    }
} // namespace eos
