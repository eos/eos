#!/usr/bin/env python3
# vim: set sw=4 sts=4 et tw=120 :

"""Collect the compiler diagnostics emitted by ``g++ -fdiagnostics-format=sarif-file``.

During the build, GCC writes one ``<source>.sarif`` file per translation unit
into the directory in which the compiler was invoked. In an out-of-tree (VPATH)
build that directory lives under the build tree, and the source ``uri`` recorded
inside each SARIF file is relative to it -- neither of which GitHub can map back
to files in the repository.

This script serves two purposes, controlled by the arguments:

* ``--sarif-out``/``--list-out``: walk the build tree, rewrite every diagnostic
  location to a repository-root-relative path, merge all runs into a single
  SARIF document (suitable for ``github/codeql-action/upload-sarif``), and emit a
  normalised, sorted, de-duplicated plain-text list of warnings.

* ``--baseline``: compare a previously generated list against a committed
  baseline, print any *new* warning as a GitHub ``::warning`` workflow command,
  and exit non-zero if the set of warnings has grown.

A warning is keyed on ``path<TAB>ruleId<TAB>message`` -- deliberately excluding
the line number, so that unrelated edits which merely shift line numbers do not
register as new warnings.
"""

import argparse
import json
import os
import sys


SEP = '\t'


def _iter_results(sarif):
    """Yield ``(run, result)`` pairs for every result in a SARIF document."""
    for run in sarif.get('runs', []):
        for result in run.get('results', []):
            yield run, result


def _primary_location(result):
    """Return the first physicalLocation of a result, or ``None``."""
    for location in result.get('locations', []):
        physical = location.get('physicalLocation')
        if physical is not None:
            return physical
    return None


def _resolve_uri(sarif_path, uri):
    """Resolve a SARIF ``uri`` to an absolute filesystem path.

    GCC records the ``uri`` relative to the directory in which it ran the
    compiler, which is precisely the directory that holds the ``.sarif`` file.
    """
    if os.path.isabs(uri):
        return os.path.normpath(uri)
    return os.path.normpath(os.path.join(os.path.dirname(sarif_path), uri))


def collect(build_dir, src_dir):
    """Return ``(merged_sarif, warnings)`` gathered from ``build_dir``.

    ``merged_sarif`` is a single SARIF document whose locations are rewritten to
    be relative to ``src_dir`` (the repository root). ``warnings`` is the sorted,
    de-duplicated list of ``path<TAB>rule<TAB>message`` keys. Diagnostics whose
    source lies outside ``src_dir`` (e.g. system headers) are dropped, since they
    cannot be annotated in the repository and are not the project's warnings.
    """
    src_root = os.path.abspath(src_dir)
    merged_runs = []
    warnings = set()

    sarif_paths = []
    for root, _, files in os.walk(build_dir):
        for name in files:
            if name.endswith('.sarif'):
                sarif_paths.append(os.path.join(root, name))
    sarif_paths.sort()

    for sarif_path in sarif_paths:
        with open(sarif_path, encoding='utf-8') as f:
            sarif = json.load(f)

        for run in sarif.get('runs', []):
            kept_results = []
            for result in run.get('results', []):
                # Track compiler *warnings* only. GCC SARIF can also carry
                # ``error`` and ``note`` results; including them would pollute
                # the baseline and let non-warning diagnostics fail the check.
                # A missing ``level`` defaults to "warning" per SARIF 2.1.0, so
                # keep it and drop only explicitly non-warning results.
                if result.get('level') not in (None, 'warning'):
                    continue
                physical = _primary_location(result)
                if physical is None:
                    continue
                artifact = physical.get('artifactLocation', {})
                uri = artifact.get('uri')
                if uri is None:
                    continue

                abs_path = _resolve_uri(sarif_path, uri)
                rel_path = os.path.relpath(abs_path, src_root)
                # Drop diagnostics outside the repository (system headers, etc.).
                if rel_path.startswith(os.pardir):
                    continue

                # Rewrite every artifactLocation in this result to be
                # repository-root-relative and strip the build-tree base id.
                for location in result.get('locations', []):
                    phys = location.get('physicalLocation')
                    if phys is None:
                        continue
                    loc_artifact = phys.get('artifactLocation')
                    if loc_artifact is None:
                        continue
                    loc_uri = loc_artifact.get('uri')
                    if loc_uri is None:
                        continue
                    loc_abs = _resolve_uri(sarif_path, loc_uri)
                    loc_artifact['uri'] = os.path.relpath(loc_abs, src_root)
                    loc_artifact.pop('uriBaseId', None)

                kept_results.append(result)

                rule = result.get('ruleId', '')
                message = result.get('message', {}).get('text', '')
                warnings.add(SEP.join((rel_path, rule, message)))

            if kept_results:
                run['results'] = kept_results
                # The rewritten uris are relative to the repository root, so the
                # build-tree base directory no longer applies.
                run.pop('originalUriBaseIds', None)
                merged_runs.append(run)

    merged_sarif = {
        '$schema': 'https://docs.oasis-open.org/sarif/sarif/v2.1.0/errata01/os/schemas/sarif-schema-2.1.0.json',
        'version': '2.1.0',
        'runs': merged_runs,
    }
    return merged_sarif, sorted(warnings)


def read_list(path):
    """Read a warnings list file into a sorted list, ignoring blank lines."""
    if not os.path.exists(path):
        return []
    with open(path, encoding='utf-8') as f:
        return sorted(line.rstrip('\n') for line in f
                      if line.strip() and not line.startswith('#'))


BASELINE_HEADER = (
    '# Baseline of accepted g++ compiler warnings, one "path<TAB>rule<TAB>message"\n'
    '# entry per line. Regenerate from a build tree with:\n'
    '#   python3 .github/scripts/collect-warnings.py \\\n'
    '#     --build-dir <build> --src-dir <src> --write-baseline .github/warnings-baseline.txt\n'
)


def write_list(path, warnings, header=None):
    with open(path, 'w', encoding='utf-8') as f:
        if header:
            f.write(header)
        for warning in warnings:
            f.write(warning + '\n')


def _escape_data(value):
    """Escape the message portion of a GitHub workflow command."""
    return value.replace('%', '%25').replace('\r', '%0D').replace('\n', '%0A')


def _escape_property(value):
    """Escape a property value (e.g. ``file=``) of a GitHub workflow command.

    In addition to the data escapes, ``:`` and ``,`` are reserved because they
    separate and terminate the property list.
    """
    return _escape_data(value).replace(':', '%3A').replace(',', '%2C')


def check_baseline(warnings, baseline_path):
    """Emit ``::warning`` commands for new warnings and return the exit code.

    New warnings (present now, absent from the baseline) fail the build.
    Warnings present in the baseline but no longer emitted are reported as
    notices so that the baseline can be pruned, but do not fail the build.
    """
    baseline = set(read_list(baseline_path))
    current = set(warnings)

    new_warnings = sorted(current - baseline)
    resolved_warnings = sorted(baseline - current)

    for warning in new_warnings:
        path, _, rest = warning.partition(SEP)
        rule, _, message = rest.partition(SEP)
        # GitHub renders this as an inline annotation even without SARIF upload.
        # Values must be escaped so reserved characters (``%``, newlines, and --
        # for properties -- ``:`` and ``,``) cannot malform the command.
        print(f'::warning file={_escape_property(path)}::'
              f'[{_escape_data(rule)}] {_escape_data(message)}')

    if resolved_warnings:
        print(f'::notice::{len(resolved_warnings)} baselined warning(s) are no longer '
              'emitted; consider regenerating .github/warnings-baseline.txt', file=sys.stderr)

    if new_warnings:
        print(f'\n{len(new_warnings)} new compiler warning(s) introduced:', file=sys.stderr)
        for warning in new_warnings:
            print('  ' + warning.replace(SEP, '  |  '), file=sys.stderr)
        return 1

    print(f'No new compiler warnings ({len(current)} known, baseline has {len(baseline)}).',
          file=sys.stderr)
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--build-dir', help='build tree to scan for *.sarif files')
    parser.add_argument('--src-dir', help='repository root the source uris are made relative to')
    parser.add_argument('--sarif-out', help='write the merged SARIF document here')
    parser.add_argument('--list-out', help='write the normalised warnings list here')
    parser.add_argument('--list-in', help='read a previously generated warnings list from here')
    parser.add_argument('--baseline', help='baseline list to compare against; fail on new warnings')
    parser.add_argument('--write-baseline', help='(re)generate the baseline list at this path')
    args = parser.parse_args()

    # Collection mode: scan the build tree and emit the merged SARIF and list.
    if args.build_dir or args.sarif_out or args.list_out:
        if not (args.build_dir and args.src_dir):
            parser.error('--build-dir and --src-dir are required to collect warnings')
        merged_sarif, warnings = collect(args.build_dir, args.src_dir)
        if args.sarif_out:
            with open(args.sarif_out, 'w', encoding='utf-8') as f:
                json.dump(merged_sarif, f, ensure_ascii=False, indent=1)
        if args.list_out:
            write_list(args.list_out, warnings)
        if args.write_baseline:
            write_list(args.write_baseline, warnings, header=BASELINE_HEADER)
        print(f'Collected {len(warnings)} distinct warning(s) from {args.build_dir}.',
              file=sys.stderr)
    elif args.list_in:
        warnings = read_list(args.list_in)
        if args.write_baseline:
            write_list(args.write_baseline, warnings, header=BASELINE_HEADER)
    else:
        parser.error('nothing to do: pass --build-dir (collect) or --list-in (check)')

    if args.baseline:
        return check_baseline(warnings, args.baseline)
    return 0


if __name__ == '__main__':
    sys.exit(main())
