import eos
import re
import yaml
from jinja_util import print_template, qn_to_link_map

def latex_to_rst(s):
    return(re.sub(r'\$([^\$]*)\$', r':math:`\1`', s))

def make_constraints():
    # All observables (pseudo or otherwise), to be used to see if a constraint is an observable or a parameter
    obs = [qn for qn, _ in eos.Observables()]

    result = []
    for qn, entry in eos.Constraints():
        data = yaml.safe_load(entry.serialize())

        if 'observable' in data:
            unique_observables = { str(data['observable']) }
        elif 'observables' in data:
            unique_observables = { str(o) for o in data['observables'] }

        observable_entries = []
        for oqn in unique_observables:
            anchor = oqn.translate(qn_to_link_map).lower()
            if oqn in obs:
                # We are looking at a (pseudo-)observable so link to the observables page
                observable_entries.append(f'`{ oqn } <observables.html#{ anchor }>`_')
            else:
                # We are looking at a parameter, so link there
                observable_entries.append(f'`{ oqn } <parameters.html#{ anchor }>`_')

        constraint = {
            'observables': ', '.join(observable_entries)
        }
        result.append((qn, constraint))

    return result

if __name__ == '__main__':

    print_template(__file__,
        version = eos.__version__,
        constraints = make_constraints(),
        len = len,
    )
