import eos.figure
from jinja_util import print_template


def _types(registry):
    """Turn a factory's registry into (selector, class name, one-line description) tuples."""
    rows = []
    for selector, cls in registry.items():
        doc = (cls.__doc__ or '').strip()
        oneline = doc.splitlines()[0].strip() if doc else ''
        rows.append((selector, cls.__name__, oneline))
    return rows


# Get the recognized figure, plot, and item types from their factories' registries.
figure_types = _types(eos.figure.figure.FigureFactory.registry)
plot_types   = _types(eos.figure.plot.PlotFactory.registry)
item_types   = _types(eos.figure.item.ItemFactory.registry)

print_template(__file__,
    figure_types = figure_types,
    plot_types   = plot_types,
    item_types   = item_types,
)
