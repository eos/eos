import eos
from jinja_util import print_template, qn_to_link_map

def get_units():
   unit_names = {m if m[0].isupper() else 'Undefined' for m in dir(eos.Unit)}
   return [getattr(eos.Unit, n)() for n in unit_names]

if __name__ == '__main__':

    print_template(__file__,
        units = [(u.__str__(), r':math:`'+u.latex()+'`') for u in get_units()]
    )
