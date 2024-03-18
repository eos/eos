import yaml as _yaml

class Deserializable:
    def __init__(self):
        pass

    @classmethod
    def from_yaml(cls, yaml_data:str):
        kwargs = _yaml.safe_load(yaml_data)
        return cls.from_dict(**kwargs)

    @staticmethod
    def make(cls, **kwargs):
        if not isinstance(cls, type):
            raise ValueError(f'Argument cls=\'{cls}\' is not a type')

        if not issubclass(cls, Deserializable):
            raise ValueError(f'Class {cls.__name__} is not a subclass of Deserializable')

        try:
            result = cls(**kwargs)
        except Exception as e:
            raise ValueError(f'When creating {cls.__name__} from {kwargs}: {e}')

        return result

    @classmethod
    def from_dict(cls, **kwargs):
        return Deserializable.make(cls, **kwargs)
