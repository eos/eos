import eos
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

print_template(__file__,
    plot_types = plot_types
)
