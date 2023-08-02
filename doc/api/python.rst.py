import eos
import eos.figure
from jinja_util import print_template


# get docstrings from classes
plot_types = {}
for key, PlotterClass in eos.plot.Plotter.plot_types.items():
    oneline = PlotterClass.__doc__.split('\n')[0]
    content = PlotterClass._api_doc if '_api_doc' in dir(PlotterClass) else ''
    plot_types.update({
        key: {
            'oneline': oneline,
            'content': content
        }
    })


# Get item types
figure_item_types = [] # tuple of (key, class, description)
reg = eos.figure.item.ItemFactory.registry
for item_key, item_class in reg.items():
    description = item_class.__doc__.splitlines()[0] # First line of docstring
    figure_item_types.append((item_key, f"eos.figure.item.{item_class.__name__}", description))


# Get plot types
figure_plot_types = [] # tuple of (key, class, description)
reg = eos.figure.plot.PlotFactory.registry
for plot_key, plot_class in reg.items():
    doc = plot_class.__doc__
    if doc is not None:
        description = doc.splitlines()[0] # First line of docstring
        figure_plot_types.append((plot_key, f"eos.figure.plot.{plot_class.__name__}", description))
    else:
        # Do not show in docs when no doc string is present
        pass


print_template(__file__,
    plot_types = plot_types,
    figure_item_types = figure_item_types,
    figure_plot_types = figure_plot_types
)
