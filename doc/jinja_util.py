"""
A helper to avoid (boilerplate) code duplication.

Usage: Create *.rst.py and *.rst.jinja file in the same directory.
The python file renders the jinja template.

In *.rst.py:
```
[...other imports...]
from jinja_util import print_template

[...processing for template variables...]

print_template(__file__,
    [...kwargs...]
)

```
"""
import jinja2, os

def get_template(abs_path):
    """Get the jinja2 template for the *.rst.py file.

    The template must be in the same directory with the same
    file name prefix.
    """
    dir_path  = os.path.dirname(abs_path)
    fname = os.path.basename(abs_path)
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(searchpath=dir_path),
        undefined=jinja2.StrictUndefined, # error for undefined variables in
                                          # template
        trim_blocks=True,
        lstrip_blocks=True)
    template_fname = fname.replace('.py', '.jinja')
    template = env.get_template(template_fname)

    return template


def print_template(rst_py__file__, **kwargs):

    template = get_template(os.path.abspath(rst_py__file__))
    print(template.render(**kwargs))


# Define the string mapping to convert an EOS QualifiedName to something suitable
# for reStructuredText and HTML
qn_to_link_map = {
    ord(':'): 'co', ord('@'): 'at', ord('/'): 'sl', ord('_'): 'un',
    ord('('): 'po', ord(')'): 'pc', ord('+'): 'pp', ord('-'): 'mm',
    ord('>'): 'to', ord('^'): 'ca', ord('*'): 'as', ord('['): 'so',
    ord(']'): 'sc', ord('{'): 'bo', ord('}'): 'bc', ord('\''): 'pr',
    ord(','): 'cm'
}


def make_parameter_link_keys():
    """Map the qualified name of every parameter to its link target on the parameters page.

    RST target names are matched without regard to case, so parameters whose names differ
    only in case need disambiguation.
    """
    import eos

    keys, used = {}, set()
    for section in eos.Parameters().sections():
        for group in section:
            for parameter in group:
                qn  = str(parameter.name())
                key = qn.translate(qn_to_link_map).lower()

                candidate, index = key, 0
                while candidate in used:
                    index    += 1
                    candidate = f'{key}-{index}'

                used.add(candidate)
                keys[qn] = candidate

    return keys
