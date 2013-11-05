#include <eos/statistics/density-wrapper.hh>
#include <eos/utils/wrapped_forward_iterator-impl.hh>

namespace eos
{
    DensityWrapper::DensityWrapper(const WrappedDensity & density) :
        _density(density)
    {
    }

    DensityWrapper::DensityWrapper(RawDensity func) :
        _density(WrappedDensity(func))
    {
    }

    DensityWrapper::~DensityWrapper()
    {
    }

    void
    DensityWrapper::add(const ParameterDescription & def)
    {
        _defs.push_back(def);
    }

     double
     DensityWrapper::evaluate() const
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
     DensityWrapper::clone() const
     {
         DensityWrapper * density = new DensityWrapper(_density);
         for (auto & d : _defs)
         {
             density->add(ParameterDescription{ d.parameter->clone(), d.min, d.max, d.nuisance });
         }
         return DensityPtr(density);
     }

     Density::Iterator
     DensityWrapper::begin() const
     {
         return _defs.begin();
     }

     Density::Iterator
     DensityWrapper::end() const
     {
         return _defs.end();
     }
}
