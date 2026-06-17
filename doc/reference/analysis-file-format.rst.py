import eos
import eos.analysis_file_description as afd
from jinja_util import print_template


def _types(registry):
    """Turn a dispatcher's registry into (selector, class name, one-line description) tuples."""
    rows = []
    for selector, cls in registry.items():
        doc = (cls.__doc__ or '').strip()
        oneline = doc.splitlines()[0].strip() if doc else ''
        rows.append((selector, cls.__name__, oneline))
    return rows


# Get the recognized prior and mask description types from their dispatchers' registries.
# Drop the 'flat', 'gauss', and 'gauss (with min/max)' entries from the table: 'flat' and 'gauss'
# are merely aliases for 'uniform' and 'gaussian', while a curtailed Gaussian is selected by adding
# 'min'/'max' to a 'gaussian' prior rather than by a distinct type.
prior_types = [t for t in _types(afd.PriorDescription.registry) if t[0] not in ('flat', 'gauss', 'gauss (with min/max)')]
mask_types  = _types(afd.MaskDescription.registry)

print_template(__file__,
    prior_types = prior_types,
    mask_types  = mask_types,
)
