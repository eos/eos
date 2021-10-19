import eos
import jinja2, os
path = os.path.dirname(os.path.abspath(__file__))
env      = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath=path))
template = env.get_template('python.rst.jinja')

plot_types_attribute = eos.plot.Plotter.plot_types

# get docstrings from classes
plot_type_descrs = dict()
for key in plot_types_attribute:
    PlotterClass = plot_types_attribute[key]
    descr = PlotterClass.__doc__.split('\n')[0]
    plot_type_descrs.update({key: descr})

rendered_str = template.render(
    plot_type_descrs = plot_type_descrs
)

print(rendered_str)

