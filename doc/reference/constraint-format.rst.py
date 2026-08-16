import eos
import yaml
from jinja_util import print_template


# Constraint names pinned for documentation purposes. The named entry is serialized live from
# eos.Constraints() at doc-build time, so the shown example always matches the exact current on-disk
# format (key names, kinematics/options shape, etc.) rather than a hand-transcribed copy that could
# silently go stale relative to the schema table above it. Pinning the identity still avoids drift
# caused by the constraint database gaining, losing, or reordering entries: if a pinned name ever
# disappears or changes type, find_examples() below fails loudly instead of silently substituting a
# different example.
EXAMPLE_NAMES = {
    'Gaussian':                        'B^0_s->mu^+mu^-::BR@CMS:2016A',
    'Amoroso':                         'B^0_s->mu^+mu^-::BR@CMS:2013B',
    'MultivariateGaussian':            'B^0->K^*0gamma::S_K+C_K@BaBar:2008A',
    'MultivariateGaussian(Covariance)': 'B->pi::f_0+f_+@FLAG:2024A',
    'UniformBound':                    'b->u::MesonicBound[1^+_A]',
    'Mixture':                         'ublnu::P(WET)@LMNRvD:2023A',
}

# Types with no shipped instance to pin at all: nothing to serialize, so this is necessarily static
# text, hand-crafted purely to illustrate the format.
STATIC_EXAMPLES = {
    'LogGamma': '''\
# No constraint of type LogGamma currently ships with EOS; the following illustrates the format.
Example::Constraint@LogGamma:2026A:
    type: LogGamma
    observable: mass::b(MSbar)
    kinematics: {}
    options: {}
    mode: 0.53
    sigma: { hi: 0.1, lo: 0.19 }
    alpha: 0.383056
    lambda: 0.0687907
''',
}


def format_example(qn, body):
    """Render a constraint entry's serialized body as a full, named database entry."""
    indented_body = '\n'.join(('  ' + line if line else line) for line in body.splitlines())

    return f'{qn}:\n{indented_body}'


def find_examples():
    constraints = eos.Constraints()

    examples = {}
    for type_name, qn in EXAMPLE_NAMES.items():
        try:
            entry = constraints[qn]
        except RuntimeError as e:
            raise RuntimeError(
                f"constraint-format.rst.py: the example pinned for '{type_name}' ('{qn}') no longer "
                f"exists in the constraint database ({e}); update EXAMPLE_NAMES with a new example"
            )

        if entry.type() != type_name:
            raise RuntimeError(
                f"constraint-format.rst.py: the example pinned for '{type_name}' ('{qn}') is now of "
                f"type '{entry.type()}'; update EXAMPLE_NAMES with a new example"
            )

        examples[type_name] = format_example(qn, entry.serialize())

    examples.update(STATIC_EXAMPLES)

    return examples


def make_constraint_types():
    examples = find_examples()

    rows = []
    for type_name in sorted(eos.ConstraintEntry.known_types()):
        fields  = yaml.safe_load(eos.ConstraintEntry.schema_as_yaml(type_name))
        example = examples.get(type_name)
        rows.append((type_name, fields, example))

    return rows


if __name__ == '__main__':

    print_template(__file__,
        version = eos.__version__,
        constraint_types = make_constraint_types(),
    )
