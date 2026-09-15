import eos
import eos.figure
from eos._api import API_BASIC_CLASSES, API_COMMON_CLASSES
from jinja_util import print_template


# Get figure types
figure_types = [] # tuple of (key, class, description)
reg = eos.figure.figure.FigureFactory.registry
for figure_key, figure_class in reg.items():
    description = figure_class.__doc__.splitlines()[0] # First line of docstring
    figure_types.append((figure_key, f"{figure_class.__name__}", description))


# Get plot types
plot_types = [] # tuple of (key, class, description)
reg = eos.figure.plot.PlotFactory.registry
for plot_key, plot_class in reg.items():
    description = plot_class.__doc__.splitlines()[0] # First line of docstring
    plot_types.append((plot_key, f"{plot_class.__name__}", description))


# Get item types
item_types = [] # tuple of (key, class, description)
reg = eos.figure.item.ItemFactory.registry
for item_key, item_class in reg.items():
    description = item_class.__doc__.splitlines()[0] # First line of docstring
    item_types.append((item_key, f"{item_class.__name__}", description))


# Document eos.tasks automatically
excluded_tasks = []
task_names = [task.__name__ for task in eos.tasks._tasks.values() if task.__name__ not in excluded_tasks]
task_names = sorted(task_names)

print_template(__file__,
    api_basic_classes = API_BASIC_CLASSES,
    api_common_classes = API_COMMON_CLASSES,
    figure_types = figure_types,
    item_types = item_types,
    plot_types = plot_types,
    task_names = task_names
)
