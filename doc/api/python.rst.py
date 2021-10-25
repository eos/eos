import eos
from jinja_util import print_template

plot_types_attribute = eos.plot.Plotter.plot_types

# get docstrings from classes
plot_type_descrs = dict()
for key in plot_types_attribute:
    PlotterClass = plot_types_attribute[key]
    descr = PlotterClass.__doc__.split('\n')[0]
    plot_type_descrs.update({key: descr})

print_template(__file__,
    plot_type_descrs = plot_type_descrs
)

