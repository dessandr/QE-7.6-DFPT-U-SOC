#!/usr/bin/env python3
"""Collect QE PH dynamical matrices and responses into a qe2pert ``phdir``.

The collector does not transform data.  It maps the irreducible-q index from
a PH ``dyn0`` file onto the image-local scratch layout and emits the filenames
already consumed by qe2pert::read_dynmat/load_dvscf/load_dns.
PH displacement-pattern metadata is preserved for qe2pert::read_patterns.
"""

from __future__ import annotations

import argparse
import filecmp
import os
import re
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


ARTIFACT_KINDS = (
    "dvscf",
    "dnsbare",
    "dnsbare_pattern",
    "dnsscf",
)
TEXT_KINDS = ARTIFACT_KINDS[1:]
NON_HUBBARD_KINDS = ("dvscf",)
HUBBARD_MODES = ("required", "none")


class CollectionError(ValueError):
    """Raised when the PH response set is ambiguous or incomplete."""


@dataclass(frozen=True)
class QGrid:
    """Irreducible-q order read from a PH dyn0 file."""

    mesh: tuple[int, int, int]
    points: tuple[tuple[float, float, float], ...]


@dataclass(frozen=True)
class Artifact:
    """One source response and its flat qe2pert destination name."""

    q_index: int
    kind: str
    source: Path
    destination: str
    size: int


@dataclass(frozen=True)
class DynmatArtifact:
    """One native QE dynmat file and its qe2pert destination name."""

    file_index: int
    source: Path
    destination: str
    size: int


@dataclass(frozen=True)
class PatternArtifact:
    """One PH displacement-pattern file and its replicated source copies."""

    q_index: int
    source: Path
    copies: tuple[Path, ...]
    destination: str
    size: int


def fortran_float(value: str) -> float:
    """Parse a Fortran real token."""

    return float(value.replace("D", "E").replace("d", "e"))


def meaningful_lines(path: Path) -> list[str]:
    """Return nonempty, comment-free records from a small text file."""

    lines = []
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.split("!", 1)[0].strip()
        if line:
            lines.append(line)
    return lines


def parse_dyn0(path: Path) -> QGrid:
    """Read the mesh and irreducible-q ordering from a QE dyn0 file."""

    try:
        lines = meaningful_lines(path)
    except FileNotFoundError as error:
        raise CollectionError(f"dyn0 file not found: {path}") from error
    if len(lines) < 2:
        raise CollectionError(f"{path}: truncated dyn0 header")

    try:
        mesh_values = tuple(int(token) for token in lines[0].split())
        q_count = int(lines[1].split()[0])
    except ValueError as error:
        raise CollectionError(f"{path}: invalid dyn0 header") from error
    if len(mesh_values) != 3 or any(value < 1 for value in mesh_values):
        raise CollectionError(f"{path}: invalid q mesh {mesh_values}")
    if q_count < 1 or len(lines) < q_count + 2:
        raise CollectionError(f"{path}: invalid irreducible-q count {q_count}")

    points = []
    for index, line in enumerate(lines[2 : q_count + 2], 1):
        tokens = line.split()
        if len(tokens) < 3:
            raise CollectionError(f"{path}: truncated q point {index}")
        try:
            points.append(tuple(fortran_float(token) for token in tokens[:3]))
        except ValueError as error:
            raise CollectionError(f"{path}: invalid q point {index}") from error
    return QGrid(mesh=mesh_values, points=tuple(points))


def response_kinds(hubbard: str) -> tuple[str, ...]:
    """Return the response-file contract selected for qe2pert."""

    if hubbard == "required":
        return ARTIFACT_KINDS
    if hubbard == "none":
        return NON_HUBBARD_KINDS
    raise CollectionError(f"unsupported Hubbard mode: {hubbard}")


def discover_dynmat_artifacts(
    dyn0: Path, prefix: str, q_count: int
) -> tuple[DynmatArtifact, ...]:
    """Discover ``dyn0`` and the numbered XML files qe2pert actually reads."""

    if not dyn0.is_file():
        raise CollectionError(f"dyn0 file not found: {dyn0}")
    if not dyn0.name.endswith("0"):
        raise CollectionError(
            f"cannot derive numbered dynmat files from dyn0 name: {dyn0.name}"
        )

    source_base = dyn0.with_name(dyn0.name[:-1])
    artifacts = []
    for file_index in range(q_count + 1):
        if file_index == 0:
            candidates = (dyn0,)
        else:
            base_candidate = Path(f"{source_base}{file_index}")
            xml_candidate = Path(f"{base_candidate}.xml")
            if base_candidate.is_file() and not xml_candidate.is_file():
                raise CollectionError(
                    f"legacy text dynmat {base_candidate} cannot be read by "
                    "qe2pert; export numbered .xml files first"
                )
            candidates = (xml_candidate,) if xml_candidate.is_file() else ()
        if not candidates:
            raise CollectionError(
                f"missing dynmat file {source_base}{file_index}.xml"
            )

        path = candidates[0].resolve()
        artifacts.append(
            DynmatArtifact(
                file_index=file_index,
                source=path,
                destination=(f"{prefix}.dyn{file_index}.xml" if file_index
                             else f"{prefix}.dyn0"),
                size=path.stat().st_size,
            )
        )
    return tuple(artifacts)


def source_name(prefix: str, kind: str) -> str:
    """Return the native PH filename for one response kind."""

    return f"{prefix}.dvscf1" if kind == "dvscf" else f"{prefix}.{kind}"


def destination_name(prefix: str, kind: str, q_index: int) -> str:
    """Return the established flat qe2pert filename."""

    return f"{prefix}.{kind}_q{q_index}"


def image_directories(source: Path) -> tuple[Path, ...]:
    """Return deterministic PH image directories."""

    images = tuple(
        image for image in sorted(source.glob("_ph*")) if image.is_dir()
    )
    if not images:
        raise CollectionError(f"{source}: no _ph* image directories found")
    return images


def q_directories(source: Path, prefix: str) -> tuple[tuple[int, Path], ...]:
    """Enumerate q owners in all image directories without recursion."""

    q_pattern = re.compile(rf"^{re.escape(prefix)}\.q_(\d+)$")
    records = []
    for image in image_directories(source):
        records.append((1, image))
        for child in sorted(image.iterdir()):
            if not child.is_dir():
                continue
            match = q_pattern.match(child.name)
            if match is not None:
                records.append((int(match.group(1)), child))
    return tuple(records)


def discover_pattern_artifacts(
    source: Path, prefix: str, q_count: int
) -> tuple[PatternArtifact, ...]:
    """Find one complete byte-consistent replicated PH pattern set."""

    pattern_name = re.compile(r"^patterns\.(\d+)\.xml$")
    found: dict[int, list[Path]] = {}
    for image in image_directories(source):
        phsave = image / f"{prefix}.phsave"
        if not phsave.is_dir():
            continue
        for path in sorted(phsave.iterdir()):
            if not path.is_file():
                continue
            match = pattern_name.match(path.name)
            if match is not None:
                found.setdefault(int(match.group(1)), []).append(path.resolve())

    unexpected = sorted(index for index in found if index < 1 or index > q_count)
    if unexpected:
        raise CollectionError(
            f"pattern files outside dyn0 q range: {unexpected}"
        )

    artifacts = []
    for q_index in range(1, q_count + 1):
        copies = tuple(found.get(q_index, ()))
        if not copies:
            raise CollectionError(
                f"missing patterns.{q_index}.xml under {source}/_ph*"
            )
        reference = copies[0]
        mismatches = tuple(
            path
            for path in copies[1:]
            if not filecmp.cmp(reference, path, shallow=False)
        )
        if mismatches:
            joined = ", ".join(str(path) for path in copies)
            raise CollectionError(
                f"mismatched replicated pattern q{q_index}: {joined}"
            )
        artifacts.append(
            PatternArtifact(
                q_index=q_index,
                source=reference,
                copies=copies,
                destination=(
                    f"{prefix}.phsave/patterns.{q_index}.xml"
                ),
                size=reference.stat().st_size,
            )
        )
    return tuple(artifacts)


def discover_artifacts(
    source: Path,
    prefix: str,
    q_count: int,
    kinds: Sequence[str] = ARTIFACT_KINDS,
) -> tuple[Artifact, ...]:
    """Find exactly one complete response set for every q index."""

    found: dict[tuple[int, str], list[Path]] = {}
    for q_index, directory in q_directories(source, prefix):
        for kind in kinds:
            path = directory / source_name(prefix, kind)
            if path.is_file():
                found.setdefault((q_index, kind), []).append(path.resolve())

    artifacts = []
    for q_index in range(1, q_count + 1):
        for kind in kinds:
            paths = found.get((q_index, kind), [])
            if not paths:
                raise CollectionError(
                    f"missing q{q_index} {kind} under {source}/_ph*"
                )
            if len(paths) != 1:
                joined = ", ".join(str(path) for path in paths)
                raise CollectionError(f"ambiguous q{q_index} {kind}: {joined}")
            path = paths[0]
            artifacts.append(
                Artifact(
                    q_index=q_index,
                    kind=kind,
                    source=path,
                    destination=destination_name(prefix, kind, q_index),
                    size=path.stat().st_size,
                )
            )

    unexpected = sorted(
        (q_index, kind)
        for q_index, kind in found
        if q_index < 1 or q_index > q_count
    )
    if unexpected:
        raise CollectionError(f"artifacts outside dyn0 q range: {unexpected}")
    return tuple(artifacts)


def count_text_records(path: Path) -> int:
    """Count list-directed complex records in a formatted Hubbard file."""

    return sum(
        len(line.split())
        for line in path.read_text(encoding="utf-8", errors="strict").splitlines()
    )


def validate_artifacts(
    artifacts: Sequence[Artifact],
    q_count: int,
    kinds: Sequence[str] = ARTIFACT_KINDS,
) -> tuple[int, int, int]:
    """Validate file counts, direct-access sizes, and text record counts."""

    expected_files = q_count * len(kinds)
    if len(artifacts) != expected_files:
        raise CollectionError(
            f"expected {expected_files} artifacts, found {len(artifacts)}"
        )
    if any(artifact.size <= 0 for artifact in artifacts):
        raise CollectionError("one or more response files are empty")

    dvscf_sizes = {
        artifact.size for artifact in artifacts if artifact.kind == "dvscf"
    }
    if len(dvscf_sizes) != 1:
        raise CollectionError(f"nonuniform dvscf sizes: {sorted(dvscf_sizes)}")

    text_counts = {
        count_text_records(artifact.source)
        for artifact in artifacts
        if artifact.kind in TEXT_KINDS
    }
    if text_counts and len(text_counts) != 1:
        raise CollectionError(
            f"nonuniform Hubbard text record counts: {sorted(text_counts)}"
        )
    return next(iter(dvscf_sizes)), next(iter(text_counts), 0), sum(
        artifact.size for artifact in artifacts
    )


def validate_dynmat_artifacts(
    artifacts: Sequence[DynmatArtifact], q_count: int
) -> int:
    """Validate complete, nonempty, uniquely named dynmat input files."""

    expected_indices = set(range(q_count + 1))
    actual_indices = {artifact.file_index for artifact in artifacts}
    if actual_indices != expected_indices or len(artifacts) != q_count + 1:
        raise CollectionError(
            "incomplete dynmat indices: "
            f"expected {sorted(expected_indices)}, found {sorted(actual_indices)}"
        )
    if any(artifact.size <= 0 for artifact in artifacts):
        raise CollectionError("one or more dynmat files are empty")
    destinations = {artifact.destination for artifact in artifacts}
    if len(destinations) != len(artifacts):
        raise CollectionError("duplicate dynmat destination names")
    return sum(artifact.size for artifact in artifacts)


def validate_pattern_artifacts(
    artifacts: Sequence[PatternArtifact], q_count: int
) -> tuple[int, int]:
    """Validate complete, nonempty, uniquely named PH pattern metadata."""

    expected_indices = set(range(1, q_count + 1))
    actual_indices = {artifact.q_index for artifact in artifacts}
    if actual_indices != expected_indices or len(artifacts) != q_count:
        raise CollectionError(
            "incomplete pattern indices: "
            f"expected {sorted(expected_indices)}, found {sorted(actual_indices)}"
        )
    if any(artifact.size <= 0 for artifact in artifacts):
        raise CollectionError("one or more pattern files are empty")
    destinations = {artifact.destination for artifact in artifacts}
    if len(destinations) != len(artifacts):
        raise CollectionError("duplicate pattern destination names")
    return (
        sum(artifact.size for artifact in artifacts),
        sum(len(artifact.copies) for artifact in artifacts),
    )


def manifest_text(
    source: Path,
    dyn0: Path,
    prefix: str,
    grid: QGrid,
    artifacts: Sequence[Artifact],
    dynmat_artifacts: Sequence[DynmatArtifact],
    pattern_artifacts: Sequence[PatternArtifact],
    hubbard: str,
) -> str:
    """Build a stable source-to-destination manifest."""

    by_key = {(item.q_index, item.kind): item for item in artifacts}
    lines = [
        f"# source\t{source}",
        f"# dyn0\t{dyn0}",
        f"# prefix\t{prefix}",
        f"# hubbard\t{hubbard}",
        "q_index\tqx\tqy\tqz\tkind\tbytes\tdestination\tsource",
    ]
    dynmat_by_index = {
        item.file_index: item for item in dynmat_artifacts
    }
    pattern_by_index = {
        item.q_index: item for item in pattern_artifacts
    }
    dyn0_item = dynmat_by_index[0]
    lines.append(
        f"0\t-\t-\t-\tdynmat\t{dyn0_item.size}\t"
        f"{dyn0_item.destination}\t{dyn0_item.source}"
    )
    for q_index, point in enumerate(grid.points, 1):
        dynmat_item = dynmat_by_index[q_index]
        lines.append(
            f"{q_index}\t{point[0]:.16e}\t{point[1]:.16e}\t"
            f"{point[2]:.16e}\tdynmat\t{dynmat_item.size}\t"
            f"{dynmat_item.destination}\t{dynmat_item.source}"
        )
        pattern_item = pattern_by_index[q_index]
        try:
            pattern_source = pattern_item.source.relative_to(source)
        except ValueError:
            pattern_source = pattern_item.source
        lines.append(
            f"{q_index}\t{point[0]:.16e}\t{point[1]:.16e}\t"
            f"{point[2]:.16e}\tpattern\t{pattern_item.size}\t"
            f"{pattern_item.destination}\t{pattern_source}"
        )
        for kind in response_kinds(hubbard):
            item = by_key[(q_index, kind)]
            try:
                relative_source = item.source.relative_to(source)
            except ValueError:
                relative_source = item.source
            lines.append(
                f"{q_index}\t{point[0]:.16e}\t{point[1]:.16e}\t"
                f"{point[2]:.16e}\t{kind}\t{item.size}\t"
                f"{item.destination}\t{relative_source}"
            )
    return "\n".join(lines) + "\n"


def materialize_file(source: Path, destination: Path, mode: str) -> None:
    """Create one collected artifact without altering its source."""

    if mode == "copy":
        shutil.copy2(source, destination)
    elif mode == "hardlink":
        os.link(source, destination)
    elif mode == "symlink":
        destination.symlink_to(os.path.relpath(source, destination.parent))
    else:  # pragma: no cover - protected by argparse and callers
        raise CollectionError(f"unsupported collection mode: {mode}")


def collect(
    output: Path,
    mode: str,
    source: Path,
    dyn0: Path,
    prefix: str,
    grid: QGrid,
    artifacts: Sequence[Artifact],
    dynmat_artifacts: Sequence[DynmatArtifact],
    pattern_artifacts: Sequence[PatternArtifact],
    hubbard: str,
) -> None:
    """Atomically create a non-overwriting flat response directory."""

    if output.exists():
        raise CollectionError(f"refusing to overwrite existing output: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{output.name}.tmp.", dir=output.parent)
    )
    try:
        for artifact in (
            *dynmat_artifacts,
            *pattern_artifacts,
            *artifacts,
        ):
            destination = staging / artifact.destination
            destination.parent.mkdir(parents=True, exist_ok=True)
            materialize_file(artifact.source, destination, mode)
            if destination.stat().st_size != artifact.size:
                raise CollectionError(
                    f"size mismatch after collecting {artifact.destination}"
                )
        (staging / "response_manifest.tsv").write_text(
            manifest_text(
                source,
                dyn0,
                prefix,
                grid,
                artifacts,
                dynmat_artifacts,
                pattern_artifacts,
                hubbard,
            ),
            encoding="utf-8",
        )
        if hubbard == "required":
            contract_text = (
                "Hubbard mode: required.\n"
                "dynmat, patterns, dvscf, dnsbare_pattern, and dnsscf "
                "are inputs.\n"
                "dnsbare is retained for diagnostics only.\n"
            )
        else:
            contract_text = (
                "Hubbard mode: none.\n"
                "dynmat, patterns, and dvscf are qe2pert phonon inputs.\n"
            )
        (staging / "README.txt").write_text(
            "Flat QE PH files for qe2pert.\n"
            f"{contract_text}"
            f"Collection mode: {mode}\n",
            encoding="utf-8",
        )
        staging.rename(output)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def build_argument_parser() -> argparse.ArgumentParser:
    """Construct the command-line interface."""

    parser = argparse.ArgumentParser(
        description=(
            "Collect QE dynmat/pattern/response files into the names consumed "
            "by qe2pert, with separate Hubbard and non-Hubbard contracts."
        )
    )
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--dyn0", required=True, type=Path)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--hubbard",
        choices=HUBBARD_MODES,
        default="required",
        help=(
            "response contract: 'required' keeps strict SOC+U checks; "
            "'none' requires dynmat, patterns, and dvscf only "
            "(default: required)"
        ),
    )
    parser.add_argument(
        "--mode",
        choices=("hardlink", "copy", "symlink"),
        default="hardlink",
        help="materialization mode (default: hardlink; requires one filesystem)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="validate and print the mapping without creating output",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point."""

    args = build_argument_parser().parse_args(argv)
    try:
        source = args.source.expanduser().resolve()
        dyn0 = args.dyn0.expanduser().resolve()
        output = args.output.expanduser().resolve()
        grid = parse_dyn0(dyn0)
        kinds = response_kinds(args.hubbard)
        dynmat_artifacts = discover_dynmat_artifacts(
            dyn0, args.prefix, len(grid.points)
        )
        pattern_artifacts = discover_pattern_artifacts(
            source, args.prefix, len(grid.points)
        )
        artifacts = discover_artifacts(
            source, args.prefix, len(grid.points), kinds
        )
        dvscf_size, text_records, total_size = validate_artifacts(
            artifacts, len(grid.points), kinds
        )
        dynmat_size = validate_dynmat_artifacts(
            dynmat_artifacts, len(grid.points)
        )
        pattern_size, pattern_copies = validate_pattern_artifacts(
            pattern_artifacts, len(grid.points)
        )
        marker = (
            "SOC_U_PH_COLLECT"
            if args.hubbard == "required"
            else "QE2PERT_PH_COLLECT"
        )

        print(
            f"{marker}_CHECK "
            f"{grid.mesh[0]} {grid.mesh[1]} {grid.mesh[2]} "
            f"{len(grid.points)} {len(artifacts)} {total_size}"
        )
        print(f"{marker}_CONTRACT {args.hubbard}")
        print(
            f"{marker}_DYNMAT_CHECK {len(dynmat_artifacts)} {dynmat_size}"
        )
        print(
            f"{marker}_PATTERN_CHECK {len(pattern_artifacts)} "
            f"{pattern_copies} {pattern_size}"
        )
        print(f"{marker}_DVSCF_CHECK {dvscf_size}")
        if text_records:
            print(f"{marker}_TEXT_CHECK {text_records}")
        for q_index, point in enumerate(grid.points, 1):
            owner = next(
                item.source.parent
                for item in artifacts
                if item.q_index == q_index and item.kind == "dvscf"
            )
            print(
                f"{marker}_Q "
                f"{q_index} {point[0]:.12e} {point[1]:.12e} "
                f"{point[2]:.12e} {owner}"
            )

        if args.dry_run:
            print(f"{marker}_RESULT PASS DRY_RUN")
            return 0
        collect(
            output,
            args.mode,
            source,
            dyn0,
            args.prefix,
            grid,
            artifacts,
            dynmat_artifacts,
            pattern_artifacts,
            args.hubbard,
        )
        print(f"{marker}_RESULT PASS {args.mode} {output}")
        return 0
    except (CollectionError, OSError, UnicodeError) as error:
        print(f"collect_qe2pert_ph.py: error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
