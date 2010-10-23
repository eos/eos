/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/observable.hh>
#include <src/utils/private_implementation_pattern-impl.hh>

#include <map>

namespace eos
{
    Observable::Observable()
    {
    }

    Observable::~Observable()
    {
    }

    template <>
    struct Implementation<ObservableOptions>
    {
        std::map<std::string, std::string> options;
    };

    ObservableOptions::ObservableOptions() :
        PrivateImplementationPattern<ObservableOptions>(new Implementation<ObservableOptions>)
    {
    }

    ObservableOptions::~ObservableOptions()
    {
    }

    const std::string &
    ObservableOptions::operator[] (const std::string & key) const
    {
        auto i(_imp->options.find(key));
        if (_imp->options.end() == i)
            throw UnknownOptionError(key);

        return i->second;
    }

    bool
    ObservableOptions::has(const std::string & key) const
    {
        return _imp->options.end() != _imp->options.find(key);
    }

    void
    ObservableOptions::set(const std::string & key, const std::string & value)
    {
        auto i(_imp->options.find(key));
        if (_imp->options.end() != i)
        {
            i->second = value;
        }
        else
        {
            _imp->options[key] = value;
        }
    }

    ObservableFactory::ObservableFactory()
    {
    }

    ObservableFactory::~ObservableFactory()
    {
    }

    UnknownOptionError::UnknownOptionError(const std::string & key) throw () :
        Exception("Unknown option: '" + key + "'")
    {
    }
}
