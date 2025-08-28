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


# Get figure types
figure_figure_types = [] # tuple of (key, class, description)
reg = eos.figure.figure.FigureFactory.registry
for figure_key, figure_class in reg.items():
    description = figure_class.__doc__.splitlines()[0] # First line of docstring
    figure_figure_types.append((figure_key, f"{figure_class.__name__}", description))


# Get plot types
figure_plot_types = [] # tuple of (key, class, description)
reg = eos.figure.plot.PlotFactory.registry
for plot_key, plot_class in reg.items():
    description = plot_class.__doc__.splitlines()[0] # First line of docstring
    figure_plot_types.append((plot_key, f"{plot_class.__name__}", description))


# Get item types
figure_item_types = [] # tuple of (key, class, description)
reg = eos.figure.item.ItemFactory.registry
for item_key, item_class in reg.items():
    description = item_class.__doc__.splitlines()[0] # First line of docstring
    figure_item_types.append((item_key, f"{item_class.__name__}", description))


# Document eos.tasks automatically
excluded_tasks = []
task_names = [task.__name__ for task in eos.tasks._tasks.values() if task.__name__ not in excluded_tasks]
task_names = sorted(task_names)

print_template(__file__,
    plot_types = plot_types,
    figure_figure_types = figure_figure_types,
    figure_item_types = figure_item_types,
    figure_plot_types = figure_plot_types,
    task_names = task_names
)
