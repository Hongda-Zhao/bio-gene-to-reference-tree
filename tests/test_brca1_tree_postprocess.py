"""Static and optional runtime checks for BRCA1 tree post-processing."""

from __future__ import annotations

import csv
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    REPOSITORY_ROOT
    / "examples"
    / "brca1"
    / "scripts"
    / "postprocess_brca1_tree.R"
)


class Brca1TreePostprocessTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory(prefix="brca1-tree-test-")
        self.root = Path(self.temporary_directory.name)
        self.metadata = self.root / "sequence_metadata.tsv"
        with self.metadata.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=("tip_id", "analysis_role", "inclusion_status"),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(
                (
                    {"tip_id": "A", "analysis_role": "study", "inclusion_status": "selected"},
                    {"tip_id": "B", "analysis_role": "expanded", "inclusion_status": "selected"},
                    {"tip_id": "C", "analysis_role": "expanded", "inclusion_status": "selected"},
                    {"tip_id": "D", "analysis_role": "expanded", "inclusion_status": "selected"},
                    {"tip_id": "OUT1", "analysis_role": "outgroup", "inclusion_status": "selected"},
                    {"tip_id": "OUT2", "analysis_role": "outgroup", "inclusion_status": "selected"},
                    {"tip_id": "REJECTED", "analysis_role": "expanded", "inclusion_status": "rejected"},
                )
            )
        self.primary = self.root / "primary.treefile"
        self.primary.write_text(
            "((OUT1:1,OUT2:1):1,((A:1,B:1):1,C:1):1,D:1);\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    @staticmethod
    def runtime_available() -> bool:
        rscript = shutil.which("Rscript")
        if not rscript:
            return False
        return (
            subprocess.run(
                [rscript, "-e", 'quit(status=if(requireNamespace("ape", quietly=TRUE)) 0 else 1)'],
                check=False,
                capture_output=True,
            ).returncode
            == 0
        )

    def run_script(
        self, *extra: str, prefix_name: str = "result"
    ) -> tuple[subprocess.CompletedProcess[str], Path]:
        prefix = self.root / prefix_name
        completed = subprocess.run(
            [
                shutil.which("Rscript") or "Rscript",
                str(SCRIPT),
                "--tree",
                str(self.primary),
                "--metadata",
                str(self.metadata),
                "--out-prefix",
                str(prefix),
                *extra,
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        return completed, prefix

    def test_script_is_local_ape_only_and_contains_strict_checks(self) -> None:
        content = SCRIPT.read_text(encoding="utf-8")
        for required in (
            "ape::read.tree",
            "ape::write.tree",
            "ape::root",
            "ape::dist.topo",
            'method = "PH85"',
            "Exactly two selected outgroup-role tips",
            "not monophyletic in the unrooted topology",
            "structural root split does not exactly isolate",
            "2L * (tip_count - 3L)",
            "outgroup_isolating_split",
            "root_screen_status",
            "does not contain the approved outgroups as an isolating unrooted split",
            "distance 0: identical unrooted topology",
            "capture_unrooted_split_labels",
            "restore_rooted_split_labels",
            "Rerooting did not preserve the exact support label",
            "single unlabeled duplicate of the root split",
        ):
            self.assertIn(required, content)
        for forbidden in (
            "install.packages",
            "BiocManager::install",
            "download.file",
            "system(",
            "system2(",
            "http://",
            "https://",
        ):
            self.assertNotIn(forbidden, content)

    def test_r_syntax_when_rscript_is_available(self) -> None:
        rscript = shutil.which("Rscript")
        if not rscript:
            self.skipTest("Rscript is not installed")
        completed = subprocess.run(
            [rscript, "-e", "parse(file=commandArgs(TRUE)[1])", str(SCRIPT)],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_runtime_roots_exact_outgroup_and_reports_rf_distances(self) -> None:
        if not self.runtime_available():
            self.skipTest("Rscript plus ape is not installed")
        identical = self.root / "identical.nwk"
        identical.write_text(self.primary.read_text(encoding="utf-8"), encoding="utf-8")
        changed = self.root / "changed.nwk"
        changed.write_text(
            "((OUT1:1,OUT2:1):1,((A:1,C:1):1,B:1):1,D:1);\n",
            encoding="utf-8",
        )
        completed, prefix = self.run_script(
            "--profile-tree",
            f"identical={identical}",
            "--profile-tree",
            f"changed={changed}",
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        for suffix in (
            ".unrooted.nwk",
            ".outgroup-rooted.nwk",
            ".topology_sensitivity.tsv",
        ):
            self.assertTrue(Path(f"{prefix}{suffix}").is_file())
        with Path(f"{prefix}.topology_sensitivity.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            rows = {row["comparison_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(rows["primary"]["rf_distance"], "0")
        self.assertEqual(rows["primary"]["outgroup_isolating_split"], "true")
        self.assertEqual(rows["primary"]["root_screen_status"], "pass")
        self.assertIn("distance 0", rows["primary"]["interpretation"])
        self.assertEqual(rows["identical"]["rf_normalized"], "0.000000")
        self.assertEqual(rows["identical"]["topology_identical"], "true")
        self.assertEqual(rows["changed"]["rf_distance"], "2")
        self.assertEqual(rows["changed"]["rf_maximum"], "6")
        self.assertEqual(rows["changed"]["rf_normalized"], "0.333333")
        self.assertEqual(rows["changed"]["topology_identical"], "false")
        self.assertEqual(rows["changed"]["outgroup_isolating_split"], "true")
        self.assertEqual(rows["changed"]["root_screen_status"], "pass")

        rooted_check = subprocess.run(
            [
                shutil.which("Rscript") or "Rscript",
                "-e",
                (
                    'x<-ape::read.tree(commandArgs(TRUE)[1]); n<-length(x$tip.label); '
                    'r<-setdiff(unique(x$edge[,1]),unique(x$edge[,2])); '
                    'kids<-x$edge[x$edge[,1]==r,2]; desc<-function(z) if(z<=n) '
                    'x$tip.label[z] else unlist(lapply(x$edge[x$edge[,1]==z,2],desc)); '
                    'sets<-lapply(kids,function(z) sort(desc(z))); '
                    'quit(status=if(ape::is.rooted(x) && any(vapply(sets,identical,logical(1),'
                    'c("OUT1","OUT2")))) 0 else 1)'
                ),
                f"{prefix}.outgroup-rooted.nwk",
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(rooted_check.returncode, 0, rooted_check.stderr)

    def test_runtime_reroot_preserves_each_unrooted_split_support_once(self) -> None:
        if not self.runtime_available():
            self.skipTest("Rscript plus ape is not installed")

        # Put the arbitrary Newick root three internal edges away from the
        # requested outgroup edge.  ape::root() alone moves these labels along
        # that path, so this fixture detects label-by-node rather than
        # label-by-unrooted-bipartition handling.
        self.primary.write_text(
            "(A:1,B:1,(C:1,(D:1,(OUT1:1,OUT2:1)11/91:1)33/93:1)22/92:1);\n",
            encoding="utf-8",
        )
        completed, prefix = self.run_script(prefix_name="supported")
        self.assertEqual(completed.returncode, 0, completed.stderr)

        rooted_check = subprocess.run(
            [
                shutil.which("Rscript") or "Rscript",
                "-e",
                r"""
x <- ape::read.tree(commandArgs(TRUE)[1])
n <- length(x$tip.label)
root <- setdiff(unique(x$edge[, 1]), unique(x$edge[, 2]))
stopifnot(length(root) == 1L)
labels <- x$node.label
if (is.null(labels)) labels <- rep(NA_character_, x$Nnode)
labels[labels == ""] <- NA_character_
desc <- function(node) {
  if (node <= n) return(x$tip.label[[node]])
  unlist(lapply(x$edge[x$edge[, 1] == node, 2], desc), use.names = FALSE)
}
find_node <- function(expected) {
  internal <- seq.int(n + 1L, n + x$Nnode)
  hits <- internal[vapply(
    internal,
    function(node) identical(sort(desc(node)), sort(expected)),
    logical(1)
  )]
  stopifnot(length(hits) == 1L)
  hits[[1L]]
}
node_label <- function(node) labels[[node - n]]

outgroup <- find_node(c("OUT1", "OUT2"))
complement <- find_node(c("A", "B", "C", "D"))
ab <- find_node(c("A", "B"))
abc <- find_node(c("A", "B", "C"))
root_children <- x$edge[x$edge[, 1] == root, 2]

stopifnot(
  isTRUE(ape::is.rooted(x)),
  identical(sort(root_children), sort(c(outgroup, complement))),
  identical(node_label(root), "Root"),
  identical(node_label(outgroup), "11/91"),
  is.na(node_label(complement)),
  identical(node_label(ab), "22/92"),
  identical(node_label(abc), "33/93"),
  identical(sort(stats::na.omit(labels)), sort(c("Root", "11/91", "22/92", "33/93")))
)
""",
                f"{prefix}.outgroup-rooted.nwk",
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(rooted_check.returncode, 0, rooted_check.stderr)

    def test_nonmonophyletic_outgroups_fail_without_outputs(self) -> None:
        if not self.runtime_available():
            self.skipTest("Rscript plus ape is not installed")
        self.primary.write_text(
            "((OUT1:1,A:1):1,((OUT2:1,B:1):1,C:1):1,D:1);\n",
            encoding="utf-8",
        )
        completed, prefix = self.run_script()
        self.assertEqual(completed.returncode, 2)
        self.assertIn("not monophyletic", completed.stderr)
        self.assertFalse(Path(f"{prefix}.unrooted.nwk").exists())
        self.assertFalse(Path(f"{prefix}.outgroup-rooted.nwk").exists())

    def test_nonmonophyletic_profile_blocks_rooted_output(self) -> None:
        if not self.runtime_available():
            self.skipTest("Rscript plus ape is not installed")
        unstable = self.root / "unstable-profile.nwk"
        unstable.write_text(
            "((OUT1:1,A:1):1,((OUT2:1,B:1):1,C:1):1,D:1);\n",
            encoding="utf-8",
        )
        completed, prefix = self.run_script(
            "--profile-tree", f"unstable={unstable}", prefix_name="profile-root-fail"
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("does not contain the approved outgroups", completed.stderr)
        self.assertFalse(Path(f"{prefix}.unrooted.nwk").exists())
        self.assertFalse(Path(f"{prefix}.outgroup-rooted.nwk").exists())
        self.assertFalse(Path(f"{prefix}.topology_sensitivity.tsv").exists())

    def test_tip_mismatch_and_negative_branch_fail_closed(self) -> None:
        if not self.runtime_available():
            self.skipTest("Rscript plus ape is not installed")
        self.primary.write_text(
            "((OUT1:1,OUT2:1):1,((A:1,B:1):1,C:1):1,EXTRA:1);\n",
            encoding="utf-8",
        )
        mismatch, prefix = self.run_script(prefix_name="mismatch")
        self.assertEqual(mismatch.returncode, 2)
        self.assertIn("tip set does not exactly equal", mismatch.stderr)
        self.assertFalse(Path(f"{prefix}.unrooted.nwk").exists())

        self.primary.write_text(
            "((OUT1:1,OUT2:1):1,((A:-0.1,B:1):1,C:1):1,D:1);\n",
            encoding="utf-8",
        )
        negative, prefix = self.run_script(prefix_name="negative")
        self.assertEqual(negative.returncode, 2)
        self.assertIn("non-negative branch length", negative.stderr)
        self.assertFalse(Path(f"{prefix}.unrooted.nwk").exists())


if __name__ == "__main__":
    unittest.main()
