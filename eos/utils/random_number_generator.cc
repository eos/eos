/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011 Danny van Dyk
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

#include <eos/utils/exception.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/random_number_generator.hh>

#include <gsl/gsl_rng.h>

namespace eos
{
    template <>
    struct Implementation<RandomNumberGenerator>
    {
        gsl_rng * rng;

        Implementation(const unsigned long & seed) :
            rng(gsl_rng_alloc(gsl_rng_mt19937))
        {
            if (! rng)
                throw InternalError("RandomNumberGenerator: Could not allocate new gsl_rng");

            gsl_rng_set(rng, seed);
        }

        ~Implementation()
        {
            gsl_rng_free(rng);
        }
    };

    RandomNumberGenerator::RandomNumberGenerator(const unsigned long & seed) :
        PrivateImplementationPattern<RandomNumberGenerator>(new Implementation<RandomNumberGenerator>(seed))
    {
    }

    RandomNumberGenerator::~RandomNumberGenerator()
    {
    }

    double
    RandomNumberGenerator::operator() (void)
    {
        return gsl_rng_uniform(_imp->rng);
    }
}
