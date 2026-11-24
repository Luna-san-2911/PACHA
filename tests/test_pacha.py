import gzip
import importlib.machinery
import importlib.util
import json
import pathlib
import tempfile
import unittest
from unittest.mock import patch

ROOT = pathlib.Path(__file__).resolve().parents[1]
loader = importlib.machinery.SourceFileLoader("pacha", str(ROOT / "PACHA"))
spec = importlib.util.spec_from_loader(loader.name, loader)
p = importlib.util.module_from_spec(spec)
loader.exec_module(p)


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = pathlib.Path(self.tmp.name)
        self.reads = self.root / "reads"
        self.genomes = self.root / "genomes"
        self.reads.mkdir()
        self.genomes.mkdir()
        for name in ("A", "B", "C"):
            for mate in (1, 2):
                (self.reads / f"{name}_R{mate}.fastq").write_text("@r\nACGT\n+\nIIII\n")
        for name in ("a.fa", "b.fa"):
            (self.genomes / name).write_text(">chr1\nACGT\n>chr2\nTGCA\n")
        self.argv = ["-m", "graph", "-r", str(self.reads), "-g", str(self.genomes),
                     "-t", "32", "-a", "2", "-at", "8", "-o", str(self.root / "out")]

    def test_dry_run_no_writes(self):
        self.assertEqual(p.main(self.argv + ["--dry-run"]), 0)
        self.assertFalse((self.root / "out").exists())

    def test_arrays_and_dependencies(self):
        args = p.parser().parse_args(self.argv)
        cfg = p.configuration(args)
        jobs = p.jobs(args, cfg)
        self.assertEqual([j[0] for j in jobs], list(p.STAGES))
        for j in jobs[3:7]:
            cmd = p.submit_command(args, cfg, j, ROOT / "PACHA", self.root / "manifest", "123")
            self.assertIn("--array=0-2%2", cmd)
            self.assertIn("--cpus-per-task=8", cmd)
            self.assertIn("--dependency=afterok:123", cmd)
        self.assertEqual(jobs[1][1:4], (2, 1, 32))

    def test_ambiguous_mates(self):
        (self.reads / "A_1.fastq").write_text("read")
        with self.assertRaises(ValueError):
            p.discover_reads(self.reads)

    def test_missing_mate(self):
        (self.reads / "A_R2.fastq").unlink()
        with self.assertRaises(ValueError):
            p.discover_reads(self.reads)

    def test_legacy_pair(self):
        for mate in (1, 3):
            (self.reads / f"D_{mate}.trimmed.fastq.gz").write_text("read")
        self.assertEqual(len(p.discover_reads(self.reads)), 4)

    def test_chromosome_mismatch(self):
        (self.genomes / "b.fa").write_text(">chr3\nACGT\n")
        with self.assertRaises(ValueError):
            p.discover_genomes(self.genomes, None)

    def test_reference_selection(self):
        args = p.parser().parse_args(self.argv + ["--reference", "b.fa"])
        self.assertEqual(p.configuration(args)["reference_index"], 1)

    def test_vcf_mode_rejected(self):
        argv = list(self.argv)
        argv[1] = "vcf"
        self.assertEqual(p.main(argv + ["--dry-run"]), 1)

    def test_consensus_genotypes(self):
        source = self.root / "in.vcf.gz"
        dest = self.root / "out.vcf"
        with gzip.open(source, "wt") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
            for i, gt in enumerate(("0|1", "1/1", "0/0", "0/1", "./.", "1"), 1):
                f.write(f"chr1\t{i}\t.\tA\tC\t30\tPASS\t.\tGT\t{gt}\n")
        counts = p.select_consensus_variants(source, dest)
        self.assertEqual(counts, dict(records=6, phased_heterozygous=1, homozygous=2,
                                     unphased_heterozygous=1, missing_or_nondiploid=2))
        self.assertEqual(len([l for l in dest.read_text().splitlines() if not l.startswith("#")]), 3)

    def test_receipt_invalidated(self):
        d = self.root / "receipt"
        d.mkdir()
        (d / "out").write_text("ok")
        p.atomic_json(d / ".complete.json", {"files": {"out": 2}})
        self.assertTrue(p.receipt_ok(d))
        (d / "out").write_text("broken")
        self.assertFalse(p.receipt_ok(d))

    def test_launch_mock_scheduler(self):
        args = p.parser().parse_args(self.argv)
        with patch.object(p.shutil, "which", return_value="/bin/tool"), \
             patch.object(p, "active_jobs", return_value=[]), \
             patch.object(p.subprocess, "check_output", side_effect=[f"{n}\n" for n in range(100, 108)]) as submit:
            p.launch(args)
        self.assertEqual(submit.call_count, 8)
        journal = json.loads((self.root / "out" / ".pacha" / "jobs.json").read_text())
        self.assertEqual(journal[-1]["stage"], "finish")
        self.assertEqual(journal[-1]["job_id"], "107")
        self.assertIn("--dependency=afterok:106", journal[-1]["command"])


if __name__ == "__main__":
    unittest.main()
