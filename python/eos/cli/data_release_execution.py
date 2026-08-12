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

"""Execution of validated EOS dataset create and publish plans."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from .data_github import GitHubRef, GitHubRelease, GitHubTag
from .data_release import (
    CreatePlan, LocalGitClient, PublishPlan, plan_create, plan_publish,
)


class ReleaseExecutionError(RuntimeError):
    def __init__(self, message: str, completed: tuple[str, ...] = ()) -> None:
        self.message = message
        self.completed = completed
        suffix = f"\nCompleted operations: {', '.join(completed)}" if completed else ''
        super().__init__(message + suffix)


class PartialPublishError(ReleaseExecutionError):
    pass


@dataclass(frozen=True)
class CreateExecutionResult:
    data_id: str
    commit_id: str
    local_tag_object_id: str
    remote_tag_object_id: str
    old_tag_object_id: str | None
    old_commit_id: str | None
    completed: tuple[str, ...]


@dataclass(frozen=True)
class PublishExecutionResult:
    data_id: str
    release_url: str
    branch_name: str
    branch_commit_id: str
    completed: tuple[str, ...]


class GitHubWrites(Protocol):
    def list_tag_refs(self) -> tuple[GitHubRef, ...]: ...
    def get_tag(self, ref: GitHubRef) -> GitHubTag: ...
    def get_release(self, tag_name: str) -> GitHubRelease | None: ...
    def get_branch(self, branch_name: str) -> str | None: ...
    def commit_parents(self, commit_id: str) -> tuple[str, ...]: ...
    def is_ancestor(self, ancestor: str, descendant: str) -> bool: ...
    def transfer_local_tag(self, root, remote_name: str, tag_name: str, expected_old: str | None) -> None: ...
    def create_release(self, tag_name: str, title: str, notes: str) -> GitHubRelease: ...
    def create_branch(self, branch_name: str, commit_id: str) -> None: ...
    def update_branch(self, branch_name: str, commit_id: str, expected_old: str) -> None: ...


def _create_signature(plan: CreatePlan) -> tuple:
    return (
        plan.source_description, plan.source_commit_id, plan.target, plan.data_id,
        plan.base_branch, plan.parent_commit_id, plan.analysis_files,
        plan.main_analysis_file, plan.commit_message, plan.annotation,
        plan.old_tag_object_id, plan.old_commit_id, plan.operations,
    )


def _publish_signature(plan: PublishPlan) -> tuple:
    return (
        plan.target, plan.data_id, plan.base_branch, plan.tag_object_id,
        plan.tag_commit_id, plan.annotation, plan.analysis_files,
        plan.main_analysis_file, plan.release_title, plan.release_notes,
        plan.branch_commit_id, plan.operations,
    )


def execute_create(
    plan: CreatePlan,
    *,
    git: LocalGitClient,
    github: GitHubWrites,
    check_factory,
) -> CreateExecutionResult:
    """Revalidate and execute one create plan, recording every completed write."""
    revision = plan.data_id != plan.base_branch
    replace = plan.old_tag_object_id is not None
    fresh = plan_create(
        plan.base_branch,
        plan.source_description,
        plan.target.root,
        analysis_files=plan.analysis_files,
        main_analysis_file=plan.main_analysis_file,
        revision=revision,
        replace=replace,
        git=git,
        github=github,
        check_factory=check_factory,
    )
    if _create_signature(fresh) != _create_signature(plan):
        raise ReleaseExecutionError('create preconditions changed after planning; no writes performed')

    completed: list[str] = []
    try:
        current_branch = git.local_branch(plan.target, plan.base_branch)
        git.prepare_local_branch(
            plan.target, plan.base_branch, current_branch, plan.parent_commit_id,
        )
        completed.append(f'prepare-local-branch:{plan.base_branch}')
        with git.checkout_commit(plan.source_description, plan.source_commit_id) as source_root:
            commit_id = git.create_dataset_commit(
                plan.target, source_root, plan.commit_message, plan.parent_commit_id,
            )
        parents = git.commit_parents(plan.target, commit_id)
        expected_parents = () if plan.parent_commit_id is None else (plan.parent_commit_id,)
        if parents != expected_parents:
            raise ReleaseExecutionError(
                f'new commit {commit_id} has parents {parents}, expected {expected_parents}',
                tuple(completed),
            )
        completed.append(f'create-commit:{commit_id}')

        git.update_local_branch(
            plan.target, plan.base_branch, commit_id, current_branch,
        )
        if git.local_branch(plan.target, plan.base_branch) != commit_id:
            raise ReleaseExecutionError('local base branch failed postcondition verification', tuple(completed))
        completed.append(f'move-local-branch:{plan.base_branch}:{commit_id}')

        local_object, peeled = git.create_local_tag(
            plan.target, plan.data_id, commit_id, plan.annotation.serialize(),
            plan.old_tag_object_id,
        )
        if peeled != commit_id:
            raise ReleaseExecutionError('local tag peeled commit failed verification', tuple(completed))
        completed.append(f'create-annotated-tag:{plan.data_id}:{local_object}')

        github.transfer_local_tag(
            plan.target.root, plan.target.remote_name, plan.data_id, plan.old_tag_object_id,
        )
        remote_ref = next(
            (ref for ref in github.list_tag_refs() if ref.name == plan.data_id), None,
        )
        if remote_ref is None or remote_ref.object_type != 'tag' or remote_ref.object_id != local_object:
            raise ReleaseExecutionError('remote tag ref failed postcondition verification', tuple(completed))
        remote_tag = github.get_tag(remote_ref)
        if (
            remote_tag.commit_id != commit_id
            or remote_tag.message.rstrip('\n') != plan.annotation.serialize()
        ):
            raise ReleaseExecutionError('remote annotated tag failed postcondition verification', tuple(completed))
        if tuple(github.commit_parents(commit_id)) != expected_parents:
            raise ReleaseExecutionError('remote tag lineage failed postcondition verification', tuple(completed))
        completed.append(f'transfer-tag:{plan.data_id}:{remote_ref.object_id}')
        return CreateExecutionResult(
            plan.data_id, commit_id, local_object, remote_ref.object_id,
            plan.old_tag_object_id, plan.old_commit_id, tuple(completed),
        )
    except Exception as error:
        done = error.completed if isinstance(error, ReleaseExecutionError) else tuple(completed)
        cause = error.message if isinstance(error, ReleaseExecutionError) else str(error)
        try:
            local_branch = git.local_branch(plan.target, plan.base_branch) or 'absent'
            local_tag = git.local_tag_refs(plan.target).get(plan.data_id)
        except Exception as observation_error:
            local_branch = f'unavailable ({observation_error})'
            local_tag = 'unavailable'
        try:
            remote_ref = next(
                (ref for ref in github.list_tag_refs() if ref.name == plan.data_id), None,
            )
            remote_tag = remote_ref.object_id if remote_ref is not None else 'absent'
        except Exception as observation_error:
            remote_tag = f'unavailable ({observation_error})'
        state = (
            f'Observed local branch={local_branch}, local tag={local_tag or "absent"}, '
            f'remote tag={remote_tag}. Safe recovery: inspect these refs and the completed '
            'operations, reconcile any race, then run create again; do not delete a published object.'
        )
        raise ReleaseExecutionError(f'create execution failed: {cause}. {state}', done) from error


def execute_publish(
    plan: PublishPlan,
    *,
    git: LocalGitClient,
    github: GitHubWrites,
    check_factory,
) -> PublishExecutionResult:
    """Revalidate and publish a tag before creating or advancing its branch."""
    fresh = plan_publish(
        plan.data_id, plan.target.root, git=git, github=github, check_factory=check_factory,
    )
    if _publish_signature(fresh) != _publish_signature(plan):
        raise ReleaseExecutionError('publish preconditions changed after planning; no writes performed')

    completed: list[str] = []
    try:
        release = github.create_release(plan.data_id, plan.release_title, plan.release_notes)
        completed.append(f'create-release:{plan.data_id}')
    except Exception as error:
        raise ReleaseExecutionError(
            f'release creation failed; remote branch was not touched: {error}', tuple(completed),
        ) from error

    try:
        verified_release = github.get_release(plan.data_id)
        if verified_release is None or verified_release.tag_name != plan.data_id:
            raise RuntimeError('GitHub release failed postcondition verification')
        tag_ref = next(
            (ref for ref in github.list_tag_refs() if ref.name == plan.data_id), None,
        )
        if (
            tag_ref is None
            or tag_ref.object_type != 'tag'
            or tag_ref.object_id != plan.tag_object_id
        ):
            observed = 'absent' if tag_ref is None else (
                f'{tag_ref.object_type} {tag_ref.object_id}'
            )
            raise RuntimeError(
                f'post-release tag ref mismatch: expected tag {plan.tag_object_id}, observed {observed}'
            )
        verified_tag = github.get_tag(tag_ref)
        if (
            verified_tag.commit_id != plan.tag_commit_id
            or verified_tag.message.rstrip('\n') != plan.annotation.serialize()
        ):
            raise RuntimeError(
                'post-release annotated tag object, peeled commit, or metadata changed'
            )
        completed.append(f'verify-release:{plan.data_id}')
    except Exception as error:
        raise PartialPublishError(
            f"release '{plan.data_id}' exists, but its exact tag failed post-release "
            f'verification and the remote branch was not touched. Inspect release and tag '
            f'{plan.data_id} before any ref update. Verification failure: {error}',
            tuple(completed),
        ) from error

    try:
        if plan.branch_commit_id is None:
            github.create_branch(plan.base_branch, plan.tag_commit_id)
            completed.append(f'create-remote-branch:{plan.base_branch}:{plan.tag_commit_id}')
        elif plan.branch_commit_id != plan.tag_commit_id:
            if not github.is_ancestor(plan.branch_commit_id, plan.tag_commit_id):
                raise ReleaseExecutionError('remote branch is no longer an ancestor of the tag commit')
            github.update_branch(plan.base_branch, plan.tag_commit_id, plan.branch_commit_id)
            completed.append(f'advance-remote-branch:{plan.base_branch}:{plan.tag_commit_id}')
        observed = github.get_branch(plan.base_branch)
        if observed != plan.tag_commit_id:
            raise ReleaseExecutionError(
                f'remote branch postcondition failed: observed {observed or "absent"}', tuple(completed),
            )
    except Exception as error:
        try:
            observed = github.get_branch(plan.base_branch)
        except Exception as observation_error:
            observed = f'unavailable ({observation_error})'
        recovery = (
            f"release '{plan.data_id}' exists, but branch 'refs/heads/{plan.base_branch}' "
            f"must target {plan.tag_commit_id}; observed {observed or 'absent'}. "
            f"After verifying ancestry, set that ref explicitly to {plan.tag_commit_id}."
        )
        raise PartialPublishError(f'{recovery} Branch failure: {error}', tuple(completed)) from error

    return PublishExecutionResult(
        plan.data_id, release.url, plan.base_branch, plan.tag_commit_id, tuple(completed),
    )
