# Copyright (c) 2026 Danny van Dyk
#
# This file is part of the EOS project. EOS is free software;
# you can redistribute it and/or modify it under the terms of the GNU General
# Public License version 2, as published by the Free Software Foundation.
#
# EOS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 59 Temple
# Place, Suite 330, Boston, MA  02111-1307  USA

from __future__ import annotations

import argparse
from collections.abc import Sequence
from contextlib import contextmanager, redirect_stderr, redirect_stdout
import logging
from pathlib import Path
import sys
from typing import TextIO

from .data_checks import (
    CheckContext,
    CheckFactory,
    CheckScope,
    PlainTextRenderer,
    register_basic_checks,
    run_checks,
)
from .data_checks_dataset import register_dataset_checks
from .data_github import GitHubClient, GitHubError
from .data_release import (
    LocalGitClient, ReleasePlanningError, SourceCheckError,
    plan_create, plan_publish, render_create_plan, render_publish_plan,
)
from .data_release_execution import (
    ReleaseExecutionError, execute_create, execute_publish,
)
from .data_register import (
    RegistrationError, execute_register, plan_register, render_register_plan,
)


def create_check_factory() -> CheckFactory:
    """Return the built-in check registry."""
    factory = CheckFactory()
    register_basic_checks(factory)
    register_dataset_checks(factory)
    return factory


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Download and manage public EOS datasets')
    common_subparser = argparse.ArgumentParser(add_help=False)
    common_subparser.add_argument(
        '-v', '--verbose',
        help='Increase the verbosity of the script',
        dest='verbose', action='count', default=0,
    )
    subparsers = parser.add_subparsers(title='commands', dest='command')

    parser_download = subparsers.add_parser(
        'download',
        parents=[common_subparser],
        description='Download a single public EOS dataset from Github or Zenodo.',
        help='Download a public EOS dataset.',
    )
    parser_download.add_argument('id', metavar='ID', help='The unique id of the public EOS dataset.')
    parser_download.set_defaults(cmd=cmd_download)

    parser_list = subparsers.add_parser(
        'list',
        parents=[common_subparser],
        description='List the ids of all public EOS datasets.',
        help='List the ids of all public EOS datasets.',
    )
    parser_list.set_defaults(cmd=cmd_list)

    parser_update = subparsers.add_parser(
        'update',
        parents=[common_subparser],
        description='Update the list of public EOS datasets from Github.',
        help='Update the list of public EOS datasets from Github.',
    )
    parser_update.set_defaults(cmd=cmd_update)

    parser_check = subparsers.add_parser(
        'check',
        parents=[common_subparser],
        description='Check an EOS dataset candidate.',
        help='Check an EOS dataset candidate.',
    )
    parser_check.add_argument(
        '--analysis-file',
        action='append',
        default=None,
        metavar='PATH',
        help='Select an analysis file relative to DIRECTORY; repeat for multiple files.',
    )
    parser_check.add_argument(
        '--main-analysis-file',
        metavar='PATH',
        help='Select the main analysis file from the analysis files.',
    )
    parser_check.add_argument(
        '-i', '--interactive',
        action='store_true',
        help='Use interactive checking (not yet available).',
    )
    parser_check.add_argument(
        'directory',
        nargs='?',
        default='.',
        metavar='DIRECTORY',
        help='Dataset directory to check (default: current directory).',
    )
    parser_check.set_defaults(cmd=cmd_check)

    parser_create = subparsers.add_parser(
        'create',
        parents=[common_subparser],
        description='Create an EOS dataset commit and transfer its annotated tag.',
        help='Create and transfer an EOS dataset tag.',
    )
    parser_create.add_argument('base_id', metavar='BASE_ID')
    parser_create.add_argument('source', metavar='URL_TO_ANALYSIS_REPO')
    parser_create.add_argument('--analysis-file', action='append', default=None, metavar='PATH')
    parser_create.add_argument('--main-analysis-file', metavar='PATH')
    choice = parser_create.add_mutually_exclusive_group()
    choice.add_argument('--revision', action='store_true')
    choice.add_argument('--replace', action='store_true')
    parser_create.add_argument('-i', '--interactive', action='store_true')
    parser_create.add_argument('--dry-run', action='store_true')
    parser_create.set_defaults(cmd=cmd_create)

    parser_publish = subparsers.add_parser(
        'publish',
        parents=[common_subparser],
        description='Publish an EOS dataset tag and then update its release branch.',
        help='Publish an EOS dataset tag.',
    )
    parser_publish.add_argument('data_id', metavar='DATA_ID')
    parser_publish.add_argument('--dry-run', action='store_true')
    parser_publish.set_defaults(cmd=cmd_publish)

    parser_register = subparsers.add_parser(
        'register',
        parents=[common_subparser],
        description='Register a completed Zenodo release in datasets.yaml.',
        help='Register a completed Zenodo release.',
    )
    parser_register.add_argument('data_id', metavar='ID')
    parser_register.add_argument('--doi', required=True, metavar='DOI')
    parser_register.add_argument('--keyword', required=True, action='append', metavar='KEYWORD')
    parser_register.add_argument('--eos-version', required=True, metavar='VERSION')
    parser_register.add_argument('--likelihood', action='append', default=None, metavar='NAME:FILENAME:FILETYPE')
    parser_register.add_argument('--repo', default='eos/data', metavar='OWNER/REPO')
    parser_register.add_argument('--dry-run', action='store_true')
    parser_register.set_defaults(cmd=cmd_register)

    return parser


@contextmanager
def _configured_logging(verbosity: int, stderr: TextIO):
    import eos

    levels = {
        0: logging.ERROR,
        1: logging.WARNING,
        2: logging.INFO,
        3: logging.DEBUG,
    }
    handler = eos.stderr_handler
    previous_level = handler.level
    previous_stream = handler.stream
    handler.setLevel(levels[min(verbosity, 3)])
    handler.setStream(stderr)
    try:
        yield
    finally:
        handler.setStream(previous_stream)
        handler.setLevel(previous_level)


def _datasets():
    from ..datasets import DataSets

    return DataSets()


def cmd_download(args: argparse.Namespace, **_kwargs) -> int:
    _datasets().download(args.id)
    return 0


def cmd_list(args: argparse.Namespace, *, stdout: TextIO | None = None, **_kwargs) -> int:
    stdout = sys.stdout if stdout is None else stdout
    for dataset_id, dataset in _datasets().datasets():
        print(f'{dataset_id:<9}  -  {dataset.title}', file=stdout)
        print(f'{"":<9}     {", ".join(dataset.authors)}', file=stdout)
    return 0


def cmd_update(args: argparse.Namespace, **_kwargs) -> int:
    _datasets().update()
    return 0


def cmd_check(
    args: argparse.Namespace,
    *,
    stdout: TextIO | None = None,
    stderr: TextIO | None = None,
    check_factory: CheckFactory | None = None,
    **_kwargs,
) -> int:
    stdout = sys.stdout if stdout is None else stdout
    stderr = sys.stderr if stderr is None else stderr
    if args.interactive:
        print('eos-data check: interactive operation is not yet available', file=stderr)
        return 2

    context = CheckContext(
        dataset_root=Path(args.directory),
        analysis_paths=tuple(Path(path) for path in (args.analysis_file or ())),
        main_analysis_path=(
            Path(args.main_analysis_file) if args.main_analysis_file is not None else None
        ),
    )
    factory = check_factory if check_factory is not None else create_check_factory()
    result = run_checks(factory, context, scope=CheckScope.COMPLETE)
    PlainTextRenderer().write(result, stdout)
    return result.exit_status


def cmd_create(
    args: argparse.Namespace,
    *,
    stdout: TextIO | None = None,
    stderr: TextIO | None = None,
    check_factory: CheckFactory | None = None,
    git_client=None,
    github_client=None,
) -> int:
    stdout = sys.stdout if stdout is None else stdout
    stderr = sys.stderr if stderr is None else stderr
    if args.interactive:
        print('eos-data create: interactive operation is not yet available', file=stderr)
        return 2
    git = LocalGitClient() if git_client is None else git_client
    github = GitHubClient() if github_client is None else github_client
    factory = create_check_factory() if check_factory is None else check_factory
    try:
        plan = plan_create(
            args.base_id, args.source, Path.cwd(),
            analysis_files=args.analysis_file or (),
            main_analysis_file=args.main_analysis_file,
            revision=args.revision,
            replace=args.replace,
            git=git,
            github=github,
            check_factory=factory,
        )
        if args.dry_run:
            stdout.write(render_create_plan(plan))
            return 0
        result = execute_create(plan, git=git, github=github, check_factory=factory)
        print(f'Created {result.data_id} at commit {result.commit_id}', file=stdout)
        print(f'Local/remote annotated tag object: {result.local_tag_object_id}', file=stdout)
        if result.old_tag_object_id is not None:
            print(
                f'Replaced tag object {result.old_tag_object_id} with {result.remote_tag_object_id}; '
                f'commit {result.old_commit_id} with {result.commit_id}',
                file=stdout,
            )
        return 0
    except SourceCheckError as error:
        print(error, file=stderr)
        return int(error.result.exit_status)
    except (ReleasePlanningError, ReleaseExecutionError, GitHubError) as error:
        print(f'eos-data create: {error}', file=stderr)
        return 2


def cmd_publish(
    args: argparse.Namespace,
    *,
    stdout: TextIO | None = None,
    stderr: TextIO | None = None,
    check_factory: CheckFactory | None = None,
    git_client=None,
    github_client=None,
) -> int:
    stdout = sys.stdout if stdout is None else stdout
    stderr = sys.stderr if stderr is None else stderr
    git = LocalGitClient() if git_client is None else git_client
    github = GitHubClient() if github_client is None else github_client
    factory = create_check_factory() if check_factory is None else check_factory
    try:
        plan = plan_publish(
            args.data_id, Path.cwd(), git=git, github=github, check_factory=factory,
        )
        if args.dry_run:
            stdout.write(render_publish_plan(plan))
            return 0
        result = execute_publish(plan, git=git, github=github, check_factory=factory)
        print(f'Published {result.data_id}: {result.release_url}', file=stdout)
        print(f'Remote branch {result.branch_name} -> {result.branch_commit_id}', file=stdout)
        return 0
    except SourceCheckError as error:
        print(error, file=stderr)
        return int(error.result.exit_status)
    except (ReleasePlanningError, ReleaseExecutionError, GitHubError) as error:
        print(f'eos-data publish: {error}', file=stderr)
        return 2


def cmd_register(
    args: argparse.Namespace,
    *,
    stdout: TextIO | None = None,
    stderr: TextIO | None = None,
    check_factory: CheckFactory | None = None,
    git_client=None,
    github_client=None,
) -> int:
    stdout = sys.stdout if stdout is None else stdout
    stderr = sys.stderr if stderr is None else stderr
    git = LocalGitClient() if git_client is None else git_client
    github = GitHubClient(repository=args.repo) if github_client is None else github_client
    factory = create_check_factory() if check_factory is None else check_factory
    try:
        plan = plan_register(
            args.data_id, doi=args.doi, keywords=args.keyword,
            eos_version=args.eos_version, likelihoods=args.likelihood or (),
            repository=args.repo, github=github, check_factory=factory,
        )
        if args.dry_run:
            stdout.write(render_register_plan(plan))
            return 0
        result = execute_register(plan, git=git, github=github, check_factory=factory)
        print(
            f'Registered {result.data_id} on {result.branch_name} at {result.commit_id}',
            file=stdout,
        )
        print(
            f"Create a pull request from branch '{result.branch_name}' to complete registration.",
            file=stdout,
        )
        return 0
    except SourceCheckError as error:
        print(error, file=stderr)
        return int(error.result.exit_status)
    except (RegistrationError, ReleasePlanningError, GitHubError) as error:
        print(f'eos-data register: {error}', file=stderr)
        return 2


def main(
    argv: Sequence[str] | None = None,
    *,
    stdout: TextIO | None = None,
    stderr: TextIO | None = None,
    check_factory: CheckFactory | None = None,
    git_client=None,
    github_client=None,
) -> int:
    stdout = sys.stdout if stdout is None else stdout
    stderr = sys.stderr if stderr is None else stderr
    parser = _parser()
    try:
        with redirect_stdout(stdout), redirect_stderr(stderr):
            args = parser.parse_args(argv)
    except SystemExit as error:
        return int(error.code)

    if not hasattr(args, 'cmd') or not callable(args.cmd):
        parser.print_help(file=stdout)
        return 0

    with _configured_logging(args.verbose, stderr):
        return args.cmd(
            args,
            stdout=stdout,
            stderr=stderr,
            check_factory=check_factory,
            git_client=git_client,
            github_client=github_client,
        )
