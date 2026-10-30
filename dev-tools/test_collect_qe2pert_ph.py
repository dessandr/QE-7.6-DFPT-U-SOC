#!/usr/bin/env python3
"""Regression tests for collect_qe2pert_ph.py."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from collect_qe2pert_ph import (
    ARTIFACT_KINDS,
    CollectionError,
    collect,
    discover_artifacts,
    discover_dynmat_artifacts,
    discover_pattern_artifacts,
    parse_dyn0,
    response_kinds,
    source_name,
    validate_artifacts,
    validate_dynmat_artifacts,
    validate_pattern_artifacts,
)


DYN0 = """\
2 1 1
2
0.0 0.0 0.0
0.5 0.0 0.0
"""


class CollectQe2pertPhTest(unittest.TestCase):
    """Check q ownership, flat naming, validation, and hard-link output."""

    def make_fixture(
        self, root: Path, *, hubbard: bool = True
    ) -> tuple[Path, Path]:
        source = root / "scratch"
        q1 = source / "_ph0"
        q2 = source / "_ph1" / "toy.q_2"
        q1.mkdir(parents=True)
        q2.mkdir(parents=True)
        for directory, marker in ((q1, "q1"), (q2, "q2")):
            kinds = ARTIFACT_KINDS if hubbard else response_kinds("none")
            for kind in kinds:
                path = directory / source_name("toy", kind)
                if kind == "dvscf":
                    path.write_bytes((marker + "-dvscf").encode().ljust(32, b"0"))
                else:
                    path.write_text(
                        "\n".join(f"({index}.0,0.0)" for index in range(8))
                        + "\n",
                        encoding="utf-8",
                    )
        dyn0 = root / "toy.dyn0"
        dyn0.write_text(DYN0, encoding="utf-8")
        (root / "toy.dyn1.xml").write_text("q1 dynmat\n", encoding="utf-8")
        (root / "toy.dyn2.xml").write_text("q2 dynmat\n", encoding="utf-8")
        for image in (source / "_ph0", source / "_ph1"):
            phsave = image / "toy.phsave"
            phsave.mkdir()
            for q_index in (1, 2):
                (phsave / f"patterns.{q_index}.xml").write_text(
                    f"<Root><q_index>{q_index}</q_index></Root>\n",
                    encoding="utf-8",
                )
        return source, dyn0

    def test_complete_hardlink_collection(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, dyn0 = self.make_fixture(root)
            grid = parse_dyn0(dyn0)
            artifacts = discover_artifacts(source, "toy", len(grid.points))
            self.assertEqual(validate_artifacts(artifacts, 2), (32, 8, 544))
            dynmats = discover_dynmat_artifacts(dyn0, "toy", 2)
            self.assertEqual(
                validate_dynmat_artifacts(dynmats, 2),
                sum(item.size for item in dynmats),
            )
            patterns = discover_pattern_artifacts(source, "toy", 2)
            self.assertEqual(
                validate_pattern_artifacts(patterns, 2),
                (sum(item.size for item in patterns), 4),
            )

            output = root / "flat"
            collect(
                output,
                "hardlink",
                source,
                dyn0,
                "toy",
                grid,
                artifacts,
                dynmats,
                patterns,
                "required",
            )
            self.assertEqual(
                {path.name for path in output.iterdir()},
                {
                    "toy.dyn0", "toy.dyn1.xml", "toy.dyn2.xml",
                    *(f"toy.{kind}_q{q}" for q in (1, 2) for kind in ARTIFACT_KINDS),
                    "toy.phsave",
                    "response_manifest.tsv",
                    "README.txt",
                },
            )
            for artifact in (*dynmats, *patterns, *artifacts):
                target = output / artifact.destination
                self.assertEqual(target.stat().st_ino, artifact.source.stat().st_ino)
            self.assertEqual(
                {path.name for path in (output / "toy.phsave").iterdir()},
                {"patterns.1.xml", "patterns.2.xml"},
            )
            manifest = (output / "response_manifest.tsv").read_text(
                encoding="utf-8"
            )
            self.assertIn("2\t5.0000000000000000e-01", manifest)
            self.assertIn("\tdynmat\t", manifest)

    def test_non_hubbard_collection_requires_only_dynmat_and_dvscf(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, dyn0 = self.make_fixture(root, hubbard=False)
            grid = parse_dyn0(dyn0)
            kinds = response_kinds("none")
            artifacts = discover_artifacts(source, "toy", 2, kinds)
            self.assertEqual(
                validate_artifacts(artifacts, 2, kinds), (32, 0, 64)
            )
            dynmats = discover_dynmat_artifacts(dyn0, "toy", 2)
            patterns = discover_pattern_artifacts(source, "toy", 2)

            output = root / "flat_non_u"
            collect(
                output,
                "symlink",
                source,
                dyn0,
                "toy",
                grid,
                artifacts,
                dynmats,
                patterns,
                "none",
            )
            payload_names = {
                "toy.dyn0", "toy.dyn1.xml", "toy.dyn2.xml",
                "toy.dvscf_q1",
                "toy.dvscf_q2",
                "toy.phsave",
            }
            self.assertEqual(
                {path.name for path in output.iterdir()},
                payload_names | {"response_manifest.tsv", "README.txt"},
            )
            for name in payload_names:
                target = output / name
                if name == "toy.phsave":
                    self.assertTrue(target.is_dir())
                else:
                    self.assertTrue(target.is_symlink())
                    self.assertGreater(target.stat().st_size, 0)
            for pattern in patterns:
                target = output / pattern.destination
                self.assertTrue(target.is_symlink())
                self.assertGreater(target.stat().st_size, 0)
            readme = (output / "README.txt").read_text(encoding="utf-8")
            self.assertIn("Hubbard mode: none", readme)
            self.assertNotIn("dnsbare", readme)

    def test_missing_artifact_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, dyn0 = self.make_fixture(root)
            (source / "_ph1" / "toy.q_2" / "toy.dnsscf").unlink()
            grid = parse_dyn0(dyn0)
            with self.assertRaisesRegex(CollectionError, "missing q2 dnsscf"):
                discover_artifacts(source, "toy", len(grid.points))

    def test_missing_dynmat_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _, dyn0 = self.make_fixture(root, hubbard=False)
            (root / "toy.dyn2.xml").unlink()
            with self.assertRaisesRegex(CollectionError, "missing dynmat file"):
                discover_dynmat_artifacts(dyn0, "toy", 2)

    def test_legacy_text_dynmat_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _, dyn0 = self.make_fixture(root)
            (root / "toy.dyn2.xml").unlink()
            (root / "toy.dyn2").write_text("legacy dynmat\n", encoding="utf-8")
            with self.assertRaisesRegex(CollectionError, "legacy text dynmat"):
                discover_dynmat_artifacts(dyn0, "toy", 2)

    def test_non_hubbard_missing_dvscf_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, dyn0 = self.make_fixture(root, hubbard=False)
            (source / "_ph1" / "toy.q_2" / "toy.dvscf1").unlink()
            grid = parse_dyn0(dyn0)
            with self.assertRaisesRegex(CollectionError, "missing q2 dvscf"):
                discover_artifacts(
                    source, "toy", len(grid.points), response_kinds("none")
                )

    def test_missing_pattern_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, _ = self.make_fixture(root, hubbard=False)
            for image in (source / "_ph0", source / "_ph1"):
                (image / "toy.phsave" / "patterns.2.xml").unlink()
            with self.assertRaisesRegex(CollectionError, "missing patterns.2.xml"):
                discover_pattern_artifacts(source, "toy", 2)

    def test_mismatched_replicated_pattern_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, _ = self.make_fixture(root, hubbard=False)
            (source / "_ph1" / "toy.phsave" / "patterns.2.xml").write_text(
                "<Root><q_index>different</q_index></Root>\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                CollectionError, "mismatched replicated pattern q2"
            ):
                discover_pattern_artifacts(source, "toy", 2)


if __name__ == "__main__":
    unittest.main()
