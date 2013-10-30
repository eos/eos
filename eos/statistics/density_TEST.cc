/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2013 Frederik Beaujean
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

#include <eos/statistics/density_TEST.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>
#include <test/test.hh>

#include <cmath>
#include <map>

using namespace test;
using namespace eos;

namespace
{
    /*!
     * A multivariate normal distribution with
     * zero mean and unit covariance on the log scale
     */
    double multivariate_unit_normal_pdf(const std::vector<double> parameters)
    {
        double result = 0;
        for (const auto & p : parameters)
        {
            result += -log(std::sqrt(2 * M_PI)) - power_of<2>(p) / 2.0;
        }
        return result;
    }
}

namespace eos
{
    TestParameter::TestParameter(const std::string & name, double value) :
        _name(name),
        _value(value)
    {
    }

    TestParameter::~TestParameter()
    {
    }

    MutablePtr
    TestParameter::clone() const
    {
        return MutablePtr(new TestParameter(_name, _value));
    }

    TestParameter::operator double () const
    {
        return _value;
    }

    double
    TestParameter::operator() () const
    {
        return _value;
    }

    double
    TestParameter::evaluate() const
    {
        return _value;
    }

    const Mutable &
    TestParameter::operator= (const double & value)
    {
        _value = value;

        return *this;
    }

    void
    TestParameter::set(const double & value)
    {
        _value = value;
    }

    const std::string &
    TestParameter::name() const
    {
        return _name;
    }

    template class WrappedForwardIterator<Parameters::IteratorTag, ParameterDescription>;

    template <>
    struct Implementation<SimpleParameters>
    {
        // forbid parameters with same name
        std::map<std::string, SimpleParameter::Index> parameters_map;
        std::shared_ptr<std::vector<double>> values;
        std::vector<ParameterDescription> defs;

        Implementation() :
            values(new std::vector<double>)
        {

        }

        SimpleParameter & declare(const std::string & name, const double & min,
                                const double & max, bool nuisance=false)
        {
            auto i(parameters_map.find(name));

            if (parameters_map.end() == i)
            {
                SimpleParameter::Index id = defs.size();
                parameters_map[name] = id;
                values->push_back(0);
                SimpleParameter * p = new SimpleParameter(name, id, values);
                defs.push_back(ParameterDescription{ MutablePtr(p), min, max, nuisance });
                return *p;
            }
            else
            {
                return *(static_cast<SimpleParameter *>(defs[i->second].parameter.get()));
            }
        }
    };

    SimpleParameter::SimpleParameter(const std::string & name, const Index & index,
                                     const std::shared_ptr<std::vector<double>> & parameters) :
        _name(name),
        _index(index),
        _parameters(parameters)
    {
    }

    SimpleParameter::~SimpleParameter()
    {
    }

    MutablePtr
    SimpleParameter::clone() const
    {
        return MutablePtr(new SimpleParameter(_name, _index, _parameters));
    }

    SimpleParameter::operator double () const
    {
        return (*_parameters)[_index];
    }

    double
    SimpleParameter::operator() () const
    {
        return (*_parameters)[_index];
    }

    double
    SimpleParameter::evaluate() const
    {
        return (*_parameters)[_index];
    }

    const Mutable &
    SimpleParameter::operator= (const double & value)
    {
        (*_parameters)[_index] = value;

        return *this;
    }

    void
    SimpleParameter::set(const double & value)
    {
        (*_parameters)[_index] = value;
    }

    const std::string &
    SimpleParameter::name() const
    {
        return _name;
    }

    SimpleParameters::SimpleParameters() :
        PrivateImplementationPattern<SimpleParameters>(new Implementation<SimpleParameters>())
    {
    }

    SimpleParameters::~SimpleParameters()
    {
    }

    SimpleParameters
    SimpleParameters::clone() const
    {
        SimpleParameters result;

        // copy parameters
        for (auto & d : _imp->defs)
        {
            result.declare(d.parameter->name(), d.min, d.max, d.nuisance);
        }

        // copy values
        std::copy(_imp->values->begin(), _imp->values->end(), result._imp->values->begin());

        return result;
    }

    SimpleParameters::Iterator
    SimpleParameters::begin() const
    {
        return _imp->defs.begin();
    }

    SimpleParameters::Iterator
    SimpleParameters::end() const
    {
        return _imp->defs.end();
    }

    SimpleParameter &
    SimpleParameters::declare(const std::string & name, const double & min,
                              const double & max, bool nuisance)
    {
        return _imp->declare(name, min, max, nuisance);
    }

    SimpleParameter &
    SimpleParameters::operator[] (const std::string & name) const
    {
        auto i(_imp->parameters_map.find(name));

        if (_imp->parameters_map.end() == i)
            throw UnknownParameterError(name);

        return *(static_cast<SimpleParameter *>(_imp->defs[i->second].parameter.get()));
    }

    SimpleParameter &
    SimpleParameters::operator[] (const SimpleParameter::Index & id) const
    {
        if (id >= _imp->defs.size())
            throw InternalError("Parameters::operator[] (Parameter::Id): invalid id '" + stringify(id) + "'");

        return *(static_cast<SimpleParameter *>(_imp->defs[id].parameter.get()));
    }

    bool
    SimpleParameters::operator!= (const SimpleParameters & rhs) const
    {
        return _imp->values.get() != rhs._imp->values.get();
    }

    const std::vector<double> &
    SimpleParameters::values() const
    {
        return *_imp->values;
    }

    TestDensity::TestDensity(const WrappedDensity & density) :
        _density(density)
    {
    }

    TestDensity::~TestDensity()
    {
    }

    void
    TestDensity::add(const ParameterDescription & def)
    {
        _defs.push_back(def);
    }

     double
     TestDensity::evaluate() const
     {
         // copy values
         _parameter_values.resize(_defs.size(), 0.0);
         unsigned i = 0;
         for (auto & d : _defs)
         {
             _parameter_values[i] = *(d.parameter);
             ++i;
         }

         return _density(_parameter_values);
     }

     DensityPtr
     TestDensity::clone() const
     {
         TestDensity * density = new TestDensity(_density);
         for (auto & d : _defs)
         {
             density->add(ParameterDescription{ d.parameter->clone(), d.min, d.max, d.nuisance });
         }
         return DensityPtr(density);
     }

     Density::Iterator
     TestDensity::begin() const
     {
         return _defs.begin();
     }

     Density::Iterator
     TestDensity::end() const
     {
         return _defs.end();
     }

     TestDensity
     make_multivariate_unit_normal(const unsigned & ndim)
     {
         TestDensity::WrappedDensity wrapped_density(::multivariate_unit_normal_pdf);
         TestDensity density(wrapped_density);

         for (unsigned i = 0 ; i < ndim ; ++i)
         {
             TestParameter p(std::string("par") + stringify(i));
             density.add(ParameterDescription{ p.clone(), -5, 5, false });
         }

         return density;
     }
}

class DensityTest :
    public TestCase
{
    public:
        DensityTest() :
            TestCase("density_test")
        {
        }

        void simple() const
        {
            // create, access, and modify
            {
                SimpleParameters p;
                auto mH = p.declare("mH", 120, 130);
                TEST_CHECK_EQUAL(mH.name(), "mH");

                mH.set(125);
                TEST_CHECK_EQUAL(double(mH), 125);

                p["mH"] = 129;
                TEST_CHECK_EQUAL(mH, 129);

                p[0] = 128;
                TEST_CHECK_EQUAL(mH, 128);

                TEST_CHECK_EQUAL(p.values().size(), 1);
                TEST_CHECK_EQUAL(p.values().front(), 128);

                TEST_CHECK( (p != p) == false);

                TEST_CHECK_EQUAL(p.begin()->min, 120);
                TEST_CHECK_EQUAL(p.begin()->max, 130);
                TEST_CHECK_EQUAL(p.begin()->nuisance, false);
            }

            // cloning
            {
                SimpleParameters p1;
                p1.declare("mH", 120, 130);
                p1.declare("mt", 170, 180);

                p1[0] = 125;
                p1[1] = 174;

                SimpleParameters p2 = p1.clone();

                TEST_CHECK(p1 != p2);
                TEST_CHECK_EQUAL(p1[0], p2[0]);

                // now modify p1, does p2 change?
                p1[0] = 126;
                TEST_CHECK_EQUAL(p2[0], 125);

                p2[1] = 173;
                TEST_CHECK_EQUAL(p1[1], 174);
            }
        }

        void test() const
        {
            static const double eps = 1e-13;

            // create
            {
                TestDensity::WrappedDensity wrapped_density(::multivariate_unit_normal_pdf);
                TestDensity density(wrapped_density);
                Parameters p = Parameters::Defaults();

                Parameter p1 = p.declare("x", 1.5);
                density.add(ParameterDescription{ p1.clone(), -5, 5, false });
                Parameter p2 = p.declare("y", -0.3);
                density.add(ParameterDescription{ p2.clone(), -5, 5, false });

                static const double result = -3.0078770664093453;
                TEST_CHECK_RELATIVE_ERROR(density.evaluate(), result, eps);

                // copy
                TestDensity density_copy = density;
                TEST_CHECK_RELATIVE_ERROR(density.evaluate(), result, eps);

                // clone
                DensityPtr density_clone = density.clone();
                TEST_CHECK_RELATIVE_ERROR(density_clone->evaluate(), result, eps);
            }

            // modify
            {
                TestDensity density(make_multivariate_unit_normal(2));
                density.begin()->parameter->set(1.5);
                (++density.begin())->parameter->set(-0.3);

                TEST_CHECK_RELATIVE_ERROR(density.evaluate(), -3.0078770664093453, eps);
            }
        }

        virtual void run() const
        {
            test();
            simple();
        }
} density_test;
