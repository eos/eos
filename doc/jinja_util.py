"""
A helper to avoid (boilerplate) code duplication.

Usage: Create *.rst.py and *.rst.jinja file in the same directory.
The python file renderes the jinja template.

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
        trim_blocks=True,
        lstrip_blocks=True)
    template_fname = fname.replace('.py', '.jinja')
    template = env.get_template(template_fname)

    return template


def print_template(rst_py__file__, **kwargs):

    template = get_template(os.path.abspath(rst_py__file__))
    print(template.render(**kwargs))

