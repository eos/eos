import eos
from jinja_util import print_template

# get docstrings from classes
plot_types = {}
for key, PlotterClass in eos.plot.Plotter.plot_types.items():
    # legacy: plotter class has not yet implemented Specs?

    # temporary: skip non-converted classes
    has_specs = hasattr(PlotterClass, 'specs_prototype')
    if not has_specs:
        continue

    oneline = PlotterClass.__doc__.split('\n')[0]
    content = PlotterClass._api_doc if '_api_doc' in dir(PlotterClass) else ''

    # extract arguments
    s = PlotterClass.specs_prototype
    args = []
    for name in s.names:
        opt = ''
        if s[name].optional:
            opt = ' (optional) '
        args.append(f"{name}{opt}, {s[name].validate_set.__doc__}")
    plot_types.update({
        key: {
            'oneline': oneline,
            'content': content,
            'arguments': args
        }
    })


print_template(__file__,
    plot_types = plot_types
)

