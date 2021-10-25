import eos
from jinja_util import print_template

plot_types_attrs = eos.plot.Plotter.plot_types

# get docstrings from classes
plot_type_descs = dict()
for key in plot_types_attrs:
    PlotterClass = plot_types_attrs[key]
    desc = PlotterClass.__doc__.split('\n')[0]
    plot_type_descs.update({key: desc})

print_template(__file__,
    plot_type_descs = plot_type_descs
)

