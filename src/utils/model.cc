/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <src/utils/model.hh>
#include <src/utils/standard-model.hh>

#include <map>

namespace eos
{
    Model::~Model()
    {
    }

    std::shared_ptr<Model>
    Model::make(const std::string & name, const Parameters & parameters)
    {
        typedef std::function<std::shared_ptr<Model> (const Parameters &)> ModelMaker;
        static const std::map<std::string, ModelMaker> model_makers
        {
            std::make_pair("SM", &StandardModel::make)
        };

        auto i = model_makers.find(name);

        if (model_makers.cend() == i)
            throw NoSuchModelError(name);

        return i->second(parameters);
    }

    NoSuchModelError::NoSuchModelError(const std::string & name) :
        Exception("No such model: '" + name + "'")
    {
    }
}
