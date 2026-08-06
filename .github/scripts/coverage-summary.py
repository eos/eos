#!/usr/bin/env python3
# vim: set sw=4 sts=4 et tw=120 :

"""Compute the overall C++, Python, and combined code coverage, append the
result to a CSV file, generate an SVG badge from the most recent entry, and
optionally write a per-file summary.

This script does *not* create a full coverage report; it only processes the
intermediate summary files produced by ``gcovr`` and ``coverage.py``.
"""

import argparse
import csv
import json
import os


def relative_to_root(path, root):
    """Return ``path`` relative to ``root`` if it lies below it, else unchanged.

    ``gcovr`` reports paths relative to its ``--root``, while ``coverage.py``
    reports absolute paths for sources outside the build directory. Making both
    relative to the source checkout yields comparable, host-independent names in
    the per-file summary.
    """
    if root is None:
        return path
    absolute = os.path.abspath(path)
    prefix   = os.path.abspath(root)
    if absolute == prefix or absolute.startswith(prefix + os.sep):
        return os.path.relpath(absolute, prefix)
    return path


def cpp_line_counts(report_path, root=None):
    """Return per-file and overall physical C++ line counts from a gcovr report.

    The result is a pair ``(per_file, (covered, total))``, where ``per_file``
    maps a file name to its own ``(covered, total)`` pair.

    A template header included in many translation units is instantiated once
    per unit, and ``gcov`` emits a separate coverage record for each
    instantiation *at the same physical source line*. ``gcovr``'s
    ``--json-summary`` sums those records, so a 67-line header instantiated in
    40 units is reported as ~2700 "lines" -- vastly over-weighting template-
    heavy headers in the overall percentage.

    To avoid this, we consume the detailed ``gcovr --json`` report (which lists
    the individual per-line records) and deduplicate by ``(file, line)``: a
    physical line counts once and is considered covered if it is executed in
    *any* instantiation.

    For robustness we fall back to the top-level ``line_total``/``line_covered``
    fields when the report contains no per-line data (e.g. a ``--json-summary``
    file); those numbers carry the inflation described above, and no per-file
    data is available in that case.
    """
    with open(report_path) as f:
        data = json.load(f)

    files = data.get('files')
    if not files or 'lines' not in files[0]:
        total   = int(data.get('line_total', 0))
        covered = int(data.get('line_covered', 0))
        return {}, (covered, total)

    # Maximum hit count seen for each physical line of each file.
    best = {}
    for entry in files:
        fname    = relative_to_root(entry['file'], root)
        per_line = best.setdefault(fname, {})
        for line in entry.get('lines', []):
            number = line['line_number']
            count  = line.get('count', 0)
            if count > per_line.get(number, -1):
                per_line[number] = count

    per_file = {
        fname: (sum(1 for count in lines.values() if count > 0), len(lines))
        for fname, lines in best.items()
    }
    covered = sum(entry[0] for entry in per_file.values())
    total   = sum(entry[1] for entry in per_file.values())
    return per_file, (covered, total)


def python_line_counts(json_path, root=None):
    """Return per-file and overall Python statement counts from a coverage.py file.

    The result is a pair ``(per_file, (covered, total))``, where ``per_file``
    maps a file name to its own ``(covered, total)`` pair. The overall counts
    are taken from the report's ``totals`` block where available, so that they
    are unaffected by any normalisation of the file names.
    """
    with open(json_path) as f:
        data = json.load(f)

    per_file = {}
    for fname, entry in (data.get('files') or {}).items():
        summary = entry.get('summary', {})
        per_file[relative_to_root(fname, root)] = (
            int(summary.get('covered_lines', 0)),
            int(summary.get('num_statements', 0)),
        )

    totals = data.get('totals') or {}
    if totals:
        covered = int(totals.get('covered_lines', 0))
        total   = int(totals.get('num_statements', 0))
    else:
        covered = sum(entry[0] for entry in per_file.values())
        total   = sum(entry[1] for entry in per_file.values())
    return per_file, (covered, total)


def percentage(covered, total):
    """Return the coverage percentage, guarding against division by zero."""
    return 100.0 * covered / total if total > 0 else 0.0


def append_csv(csv_path, date, cpp, python, combined):
    """Append a row to the coverage CSV, creating the header if necessary."""
    is_new = not os.path.exists(csv_path) or os.path.getsize(csv_path) == 0
    with open(csv_path, 'a', newline='') as f:
        writer = csv.writer(f)
        if is_new:
            writer.writerow(['date', 'cpp', 'python', 'combined'])
        writer.writerow([date, f'{cpp:.2f}', f'{python:.2f}', f'{combined:.2f}'])


def write_files_csv(csv_path, cpp_per_file, python_per_file):
    """Write the per-file coverage summary, replacing any previous version.

    In contrast to the historical ``coverage.csv``, this file is overwritten on
    every run: the time series lives in the git history of the coverage
    repository, and the file itself stays small (a few tens of kiB) and yields
    readable diffs. Rows are sorted so that the diffs remain stable.
    """
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['language', 'file', 'covered', 'total', 'percentage'])
        for language, per_file in (('cpp', cpp_per_file), ('python', python_per_file)):
            for fname in sorted(per_file):
                covered, total = per_file[fname]
                writer.writerow([language, fname, covered, total, f'{percentage(covered, total):.2f}'])


def most_recent_entry(csv_path):
    """Return the most recent (cpp, python, combined) tuple from the CSV."""
    with open(csv_path, newline='') as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f'no entries found in {csv_path}')
    last = rows[-1]
    return float(last['cpp']), float(last['python']), float(last['combined'])


def colour_for(value):
    """Return a shields.io-style colour for the given coverage percentage."""
    if value >= 90.0:
        return '#4c1'    # brightgreen
    if value >= 75.0:
        return '#97ca00' # green
    if value >= 60.0:
        return '#a4a61d' # yellowgreen
    if value >= 45.0:
        return '#dfb317' # yellow
    if value >= 30.0:
        return '#fe7d37' # orange
    return '#e05d44'     # red


def text_width(text):
    """A crude monospace-ish width estimate (in px) for the badge text."""
    return int(len(text) * 6.5) + 10


def make_badge(badge_path, cpp, python, combined):
    """Write an SVG badge displaying the three coverage percentages."""
    # The badge consists of a grey label segment followed by one coloured
    # segment per coverage figure.
    segments = [('coverage', '#555')]
    for label, value in (('C++', cpp), ('py', python), ('all', combined)):
        segments.append((f'{label} {value:.0f}%', colour_for(value)))

    widths  = [text_width(text) for text, _ in segments]
    total   = sum(widths)
    height  = 20

    parts = []
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'xmlns:xlink="http://www.w3.org/1999/xlink" width="{total}" height="{height}" '
        f'role="img" aria-label="code coverage">'
    )
    parts.append(
        f'<title>coverage: C++ {cpp:.0f}%, py {python:.0f}%, all {combined:.0f}%</title>'
    )
    # Rounded-corner clip and subtle gradient for the classic badge look.
    parts.append(
        f'<linearGradient id="s" x2="0" y2="100%">'
        f'<stop offset="0" stop-color="#bbb" stop-opacity=".1"/>'
        f'<stop offset="1" stop-opacity=".1"/></linearGradient>'
    )
    parts.append(
        f'<clipPath id="r"><rect width="{total}" height="{height}" rx="3" fill="#fff"/></clipPath>'
    )
    parts.append('<g clip-path="url(#r)">')

    x = 0
    for (text, colour), width in zip(segments, widths):
        parts.append(f'<rect x="{x}" width="{width}" height="{height}" fill="{colour}"/>')
        x += width
    parts.append(f'<rect width="{total}" height="{height}" fill="url(#s)"/>')
    parts.append('</g>')

    parts.append(
        '<g fill="#fff" text-anchor="middle" '
        'font-family="Verdana,Geneva,DejaVu Sans,sans-serif" '
        'text-rendering="geometricPrecision" font-size="11">'
    )
    x = 0
    for (text, _), width in zip(segments, widths):
        cx = x + width / 2
        # Drop shadow for legibility, then the actual text.
        parts.append(
            f'<text x="{cx:.1f}" y="15" fill="#010101" fill-opacity=".3">{text}</text>'
        )
        parts.append(f'<text x="{cx:.1f}" y="14">{text}</text>')
        x += width
    parts.append('</g>')
    parts.append('</svg>')

    with open(badge_path, 'w') as f:
        f.write('\n'.join(parts) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--date',        required=True, help='date of this coverage run (YYYY-MM-DD)')
    parser.add_argument('--cpp-summary', required=True,
                        help='path to the detailed gcovr JSON report (gcovr --json)')
    parser.add_argument('--python-json', required=True, help='path to the coverage.py JSON file')
    parser.add_argument('--csv',         required=True, help='path to the coverage CSV file to append to')
    parser.add_argument('--badge',       required=True, help='path to the SVG badge file to (re)generate')
    parser.add_argument('--files-csv',   required=False,
                        help='path to the per-file coverage CSV file to (re)generate')
    parser.add_argument('--source-root', required=False,
                        help='path to the source checkout; file names are reported relative to it')
    args = parser.parse_args()

    cpp_per_file,    (cpp_covered, cpp_total)       = cpp_line_counts(args.cpp_summary, args.source_root)
    python_per_file, (python_covered, python_total) = python_line_counts(args.python_json, args.source_root)

    cpp_pct      = percentage(cpp_covered, cpp_total)
    python_pct   = percentage(python_covered, python_total)
    combined_pct = percentage(cpp_covered + python_covered, cpp_total + python_total)

    print(f'C++      coverage: {cpp_pct:6.2f}%  ({cpp_covered}/{cpp_total} lines)')
    print(f'Python   coverage: {python_pct:6.2f}%  ({python_covered}/{python_total} statements)')
    print(f'Combined coverage: {combined_pct:6.2f}%  '
          f'({cpp_covered + python_covered}/{cpp_total + python_total})')

    append_csv(args.csv, args.date, cpp_pct, python_pct, combined_pct)

    if args.files_csv:
        write_files_csv(args.files_csv, cpp_per_file, python_per_file)
        print(f'per-file summary: {len(cpp_per_file)} C++ and {len(python_per_file)} Python files '
              f'written to {args.files_csv}')

    # The badge always reflects the most recent entry in the CSV file.
    cpp_latest, python_latest, combined_latest = most_recent_entry(args.csv)
    make_badge(args.badge, cpp_latest, python_latest, combined_latest)


if __name__ == '__main__':
    main()
