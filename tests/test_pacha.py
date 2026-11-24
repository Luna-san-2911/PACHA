import gzip
import importlib.machinery
import importlib.util
import json
import os
import shlex
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
        self.argv = ["--executor", "slurm", "-m", "graph", "-r", str(self.reads), "-g", str(self.genomes),
                     "-ref", str(self.genomes / "a.fa"),
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
        self.assertEqual(jobs[0][1:4], (2, 8, 1))
        self.assertTrue(jobs[0][-1])
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
            p.discover_genomes(self.genomes, self.genomes / "a.fa", validate_all=True)

    def test_reference_selection(self):
        args = p.parser().parse_args(self.argv + ["--reference", str(self.genomes / "b.fa")])
        cfg = p.configuration(args)
        self.assertEqual(cfg["reference_index"], 1)
        self.assertEqual(len(cfg["genomes"]), 2)
        self.assertEqual(cfg["genomes"][1]["path"], str(self.genomes / "b.fa"))

    def test_reference_required(self):
        argv = self.argv.copy()
        idx = argv.index("-ref")
        del argv[idx:idx+2]
        with self.assertRaises(SystemExit):
            p.parser().parse_args(argv)

    def test_external_reference(self):
        ref = self.root / "reference.fa"
        ref.write_text(">chr2\nAAAA\n>chr1\nCCCC\n")
        args = p.parser().parse_args(self.argv + ["-ref", str(ref)])
        cfg = p.configuration(args)
        self.assertEqual(len(cfg["genomes"]), 3)
        self.assertEqual(cfg["genomes"][cfg["reference_index"]]["path"], str(ref))
        self.assertEqual(cfg["chromosomes"], ["chr2", "chr1"])
        self.assertEqual(p.jobs(args, cfg)[0][1], 3)

    def test_external_reference_alias_not_duplicated(self):
        ref = self.root / "reference.fa"
        os.link(self.genomes / "a.fa", ref)
        records, index, _ = p.discover_genomes(self.genomes, ref)
        self.assertEqual(len(records), 2)
        self.assertEqual(records[index]["path"], str(ref))

    def test_large_panel_and_array_batches(self):
        for i in range(1001):
            (self.genomes / f"genome{i:04d}.fa").write_text(">chr1\nACGT\n>chr2\nTGCA\n")
        args = p.parser().parse_args(self.argv)
        cfg = p.configuration(args)
        self.assertEqual(len(cfg["genomes"]), 1003)
        self.assertEqual(len({g["id"] for g in cfg["genomes"]}), 1003)
        batches = [(j, offset) for j, offset in p.job_batches(args, cfg) if j[0] == "prepare"]
        self.assertEqual([(j[1], offset) for j, offset in batches], [(1000, 0), (3, 1000)])
        cmd = p.submit_command(args, cfg, batches[1][0], ROOT / "PACHA",
                               self.root / "manifest", "100", batches[1][1])
        self.assertIn("--array=0-2%3", cmd)
        self.assertEqual(shlex.split(cmd[-1])[-1], "1000")
        self.assertIn("--dependency=afterok:100", cmd)

    def test_prepare_uses_selected_reference_only(self):
        args = p.parser().parse_args(self.argv + ["-ref", str(self.genomes / "b.fa")])
        cfg = p.configuration(args)
        def fake_index(cmd, **kwargs):
            self.assertEqual(cmd[:2], ["samtools", "faidx"])
            pathlib.Path(str(cmd[2]) + ".fai").write_text("chr1\t4\t6\t4\t5\nchr2\t4\t17\t4\t5\n")
        with patch.object(p, "run", side_effect=fake_index):
            for i in range(len(cfg["genomes"])):
                work = p.stage_dir(cfg, "prepare", i)
                work.mkdir(parents=True)
                files = p.prepare(cfg, work, i, 1)
                p.atomic_json(work / ".complete.json", {"files": {f: (work / f).stat().st_size for f in files}})
                self.assertEqual((work / "reference.fa").exists(), i == cfg["reference_index"])
        reference = p.reference_inputs(cfg)
        self.assertEqual(reference, p.stage_dir(cfg, "prepare", cfg["reference_index"]))
        self.assertEqual((reference / "reference.fa").read_text(), (self.genomes / "b.fa").read_text())

    def test_finish_two_haplotypes_per_individual(self):
        cfg = p.configuration(p.parser().parse_args(self.argv))
        for i, sample in enumerate(cfg["samples"]):
            work = p.stage_dir(cfg, "consensus", i)
            work.mkdir(parents=True)
            p.atomic_json(work / "summary.json", {
                "sample": sample["name"], "records": 1, "phased_heterozygous": 1,
                "homozygous": 0, "unphased_heterozygous": 0,
                "missing_or_nondiploid": 0, "selected_fraction": 1.0})
            for h in (1, 2):
                for suffix in ("", ".fai", ".gzi"):
                    (work / f"{sample['name']}_hap{h}.fa.gz{suffix}").write_text("test")
            p.atomic_json(work / ".complete.json", {"files": {
                f.name: f.stat().st_size for f in work.iterdir() if f.is_file()}})
        final = self.root / "final"
        final.mkdir()
        p.finish(cfg, final, 0, 1)
        self.assertEqual(sorted(f.name for f in final.glob("*.fa.gz")), [
            f"{s['name']}_hap{h}.fa.gz" for s in cfg["samples"] for h in (1, 2)])
        self.assertTrue((final / "_SUCCESS").is_file())

    def test_vcf_mode_rejected(self):
        argv = list(self.argv)
        argv[argv.index("-m") + 1] = "vcf"
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


    def test_auto_executor_detection(self):
        for env, sbatch, expected in [
            ({}, None, "local"), ({}, "/bin/sbatch", "slurm"),
            ({"SLURM_JOB_ID": "123", "SLURM_CPUS_PER_TASK": "64"}, "/bin/sbatch", "local"),
        ]:
            args = p.parser().parse_args(self.argv + ["--executor", "auto"])
            with patch.dict(os.environ, env, clear=True), \
                 patch.object(p.shutil, "which", return_value=sbatch), \
                 patch.object(p, "available_cpus", return_value=64):
                executor, budget = p.configure_execution(args)
                self.assertEqual(executor, expected)
                self.assertEqual(budget, 64 if expected == "local" else None)

    def test_allocation_cpu_limit(self):
        with patch.dict(os.environ, {"SLURM_JOB_ID": "99", "SLURM_CPUS_PER_TASK": "12"}, clear=True), \
             patch.object(p.os, "sched_getaffinity", return_value=set(range(64))):
            self.assertEqual(p.available_cpus(), 12)
            args = p.parser().parse_args(self.argv + ["--executor", "local"])
            with self.assertRaisesRegex(ValueError, "local CPU budget"):
                p.configure_execution(args)

    def test_local_default_threads(self):
        args = p.parser().parse_args(["--executor", "local", "-m", "graph", "-r", str(self.reads),
            "-g", str(self.genomes), "-ref", str(self.genomes / "a.fa"), "-o", str(self.root / "out")])
        with patch.object(p, "available_cpus", return_value=4):
            self.assertEqual(p.configure_execution(args), ("local", 4))
        self.assertEqual((args.threads, args.array_threads, args.array), (4, 4, 1))

    def test_local_concurrency_cpu_and_memory(self):
        args = p.parser().parse_args(self.argv + ["--executor", "local", "--local-cpus", "32",
            "--local-memory", "64G", "--mem", "16G", "-a", "10"])
        with patch.object(p, "available_cpus", return_value=64):
            _, budget = p.configure_execution(args)
        cfg = p.configuration(args)
        align = next(j for j in p.jobs(args, cfg) if j[0] == "align")
        self.assertEqual(p.local_concurrency(args, align, budget), 3)
        self.assertEqual(p.local_concurrency(args, align, 16), 2)
        args.local_memory = "16G"
        self.assertEqual(p.local_concurrency(args, align, budget), 1)
        args.local_memory = "8G"
        with self.assertRaisesRegex(ValueError, "memory estimate"):
            p.local_concurrency(args, align, budget)

    def test_slurm_resource_overrides(self):
        args = p.parser().parse_args(self.argv + ["--partition", "normal", "--account", "lab",
            "--qos", "standard", "--time", "24:00:00", "--mem", "64G", "--align-mem", "96G",
            "--stage-mem", "align=120G", "--stage-time", "align=2-00:00:00",
            "--stage-partition", "align=bigmem", "--constraint", "amd",
            "--mail-type", "FAIL", "--slurm-option=--exclusive"])
        p.configure_execution(args)
        cfg = p.configuration(args)
        jobs = p.jobs(args, cfg)
        align = next(j for j in jobs if j[0] == "align")
        call = next(j for j in jobs if j[0] == "call")
        self.assertEqual(align[4:6], ("120G", "2-00:00:00"))
        self.assertEqual(call[4:6], ("64G", "24:00:00"))
        cmd = p.submit_command(args, cfg, align, ROOT / "PACHA", self.root / "manifest")
        for item in ("--partition=bigmem", "--account=lab", "--qos=standard",
                     "--constraint=amd", "--exclusive", "--mail-type=FAIL"):
            self.assertIn(item, cmd)

    def test_local_rejects_slurm_submission_flags(self):
        args = p.parser().parse_args(self.argv + ["--executor", "local", "--partition", "normal"])
        with self.assertRaisesRegex(ValueError, "require --executor slurm"):
            p.configure_execution(args)

    def test_unsafe_slurm_override_rejected(self):
        args = p.parser().parse_args(self.argv + ["--slurm-option=--array=0-4"])
        with self.assertRaisesRegex(ValueError, "Unsupported"):
            p.configure_execution(args)

    def test_local_launch_does_not_require_slurm(self):
        args = p.parser().parse_args(self.argv + ["--executor", "local"])
        def which(name):
            return "/bin/tool" if name in p.TOOLS else None
        with patch.object(p.shutil, "which", side_effect=which), \
             patch.object(p, "available_cpus", return_value=64), \
             patch.object(p, "run_local") as run, \
             patch.object(p, "active_jobs") as queue:
            p.launch(args)
        run.assert_called_once()
        queue.assert_not_called()
        self.assertFalse((self.root / "out" / ".pacha" / "jobs.json").exists())

    def fake_local_worker(self, cfg, fail=False):
        out = pathlib.Path(cfg["output"])
        out.mkdir()
        (out / "logs").mkdir()
        manifest = out / "manifest.json"
        p.atomic_json(manifest, cfg)
        script = out / "fake_worker.py"
        script.write_text(
            "import json, pathlib, sys\n"
            "_, _, manifest, stage, threads, index = sys.argv\n"
            "cfg = json.loads(pathlib.Path(manifest).read_text())\n"
            "i = int(index)\n"
            + ("if i == 0: sys.exit(7)\n" if fail else "") +
            "d = pathlib.Path(cfg['output']) / '00_inputs' / cfg['genomes'][i]['id']\n"
            "d.mkdir(parents=True)\n"
            "(d / 'data').write_text('ok')\n"
            "(d / '.complete.json').write_text(json.dumps({'files': {'data': 2}}))\n"
        )
        return script, manifest

    def test_local_real_subprocesses_and_resume(self):
        args = p.parser().parse_args(self.argv)
        cfg = p.configuration(args)
        script, manifest = self.fake_local_worker(cfg)
        job = p.jobs(args, cfg)[0]
        p.local_stage(cfg, job, 2, script, manifest)
        for i in range(2):
            self.assertTrue(p.receipt_ok(p.stage_dir(cfg, "prepare", i)))
        log = pathlib.Path(cfg["output"]) / "logs" / "prepare-0.local.log"
        before = log.read_bytes()
        p.local_stage(cfg, job, 2, script, manifest)
        self.assertEqual(log.read_bytes(), before)

    def test_local_failure_stops_new_tasks(self):
        args = p.parser().parse_args(self.argv)
        cfg = p.configuration(args)
        script, manifest = self.fake_local_worker(cfg, fail=True)
        with self.assertRaisesRegex(ValueError, "downstream stages not started"):
            p.local_stage(cfg, p.jobs(args, cfg)[0], 1, script, manifest)
        self.assertFalse(p.stage_dir(cfg, "prepare", 1).exists())
        self.assertFalse((pathlib.Path(cfg["output"]) / "logs" / "prepare-1.local.log").exists())

    def test_local_stage_order(self):
        args = p.parser().parse_args(self.argv + ["--executor", "local"])
        cfg = p.configuration(args)
        with patch.object(p, "local_stage") as runner, patch.object(p, "completed"):
            p.run_local(args, cfg, 32, ROOT / "PACHA", self.root / "manifest")
        self.assertEqual([call.args[1][0] for call in runner.call_args_list], list(p.STAGES))

    def test_local_worker_ignores_inherited_array_index(self):
        cfg = p.configuration(p.parser().parse_args(self.argv))
        manifest = self.root / "manifest"
        p.atomic_json(manifest, cfg)
        with patch.dict(os.environ, {"SLURM_ARRAY_TASK_ID": "99"}, clear=True), \
             patch.object(p, "available_cpus", return_value=8), \
             patch.object(p, "execute_task") as task:
            p.worker(manifest, "prepare", 1, task_index=1)
        task.assert_called_once_with(cfg, "prepare", 1, 1)


    def sheet_args(self, sheet):
        argv = self.argv.copy()
        idx = argv.index("-r")
        argv[idx:idx+2] = ["--sample-sheet", str(sheet)]
        return argv

    def test_sample_sheet_arbitrary_filenames(self):
        (self.reads / "a forward.fastq").write_text("@r\nACGT\n+\nIIII\n")
        (self.reads / "arbitrary-reverse.fq").write_text("@r\nACGT\n+\nIIII\n")
        sheet = self.root / "samples.tsv"
        sheet.write_text("sample\tr1\tr2\nPlant-X.7\treads/a forward.fastq\treads/arbitrary-reverse.fq\n")
        cfg = p.configuration(p.parser().parse_args(self.sheet_args(sheet)))
        self.assertEqual(len(cfg["samples"]), 1)
        sample = cfg["samples"][0]
        self.assertEqual(sample["name"], "Plant-X.7")
        self.assertEqual(sample["r1"]["path"], str(self.reads / "a forward.fastq"))
        self.assertEqual(p.stage_dir(cfg, "align", 0).name, "Plant-X.7")

    def test_sample_sheet_rejects_duplicate_reads(self):
        sheet = self.root / "samples.tsv"
        sheet.write_text("sample\tr1\tr2\nA\treads/A_R1.fastq\treads/A_R2.fastq\n"
                         "B\treads/A_R1.fastq\treads/B_R2.fastq\n")
        with self.assertRaisesRegex(ValueError, "assigned more than once"):
            p.configuration(p.parser().parse_args(self.sheet_args(sheet)))

    def test_sample_sheet_rejects_unsafe_id_and_bad_header(self):
        sheet = self.root / "samples.tsv"
        sheet.write_text("name\tforward\treverse\nA\treads/A_R1.fastq\treads/A_R2.fastq\n")
        with self.assertRaisesRegex(ValueError, "header"):
            p.read_sample_sheet(sheet)
        sheet.write_text("sample\tr1\tr2\n../A\treads/A_R1.fastq\treads/A_R2.fastq\n")
        with self.assertRaisesRegex(ValueError, "Invalid sample ID"):
            p.read_sample_sheet(sheet)

    def test_concurrent_samples_aliases_and_native_throttle(self):
        args = p.parser().parse_args(self.argv + ["--concurrent-samples", "1",
            "--sample-threads", "4", "--array-batch-size", "2"])
        self.assertEqual((args.array, args.array_threads), (1, 4))
        cfg = p.configuration(args)
        for stage in ("align", "call", "phase", "consensus"):
            batches = [(j, offset) for j, offset in p.job_batches(args, cfg) if j[0] == stage]
            self.assertEqual([(j[1], offset) for j, offset in batches], [(2, 0), (1, 2)])
            for j, offset in batches:
                cmd = p.submit_command(args, cfg, j, ROOT / "PACHA", self.root / "manifest", "123", offset)
                self.assertIn(f"--array=0-{j[1]-1}%1", cmd)
                self.assertIn("--cpus-per-task=4", cmd)
                self.assertIn("--dependency=afterok:123", cmd)

    def test_existing_allocation_auto_never_submits_children(self):
        args = p.parser().parse_args(self.argv + ["--executor", "auto",
            "--concurrent-samples", "3", "--sample-threads", "8"])
        env = {"SLURM_JOB_ID": "123", "SLURM_CPUS_PER_TASK": "32",
               "SLURM_NTASKS": "1", "SLURM_JOB_NUM_NODES": "1"}
        with patch.dict(os.environ, env, clear=True), \
             patch.object(p, "available_cpus", return_value=32), \
             patch.object(p.shutil, "which", return_value="/bin/tool"), \
             patch.object(p, "run_slurm") as slurm, \
             patch.object(p, "run_local") as local:
            p.launch(args)
        slurm.assert_not_called()
        local.assert_called_once()
        self.assertFalse((self.root / "out" / ".pacha" / "jobs.json").exists())

    def test_outer_slurm_array_rejected(self):
        args = p.parser().parse_args(self.argv + ["--executor", "local"])
        with patch.dict(os.environ, {"SLURM_JOB_ID": "123", "SLURM_ARRAY_TASK_ID": "0"}, clear=True):
            with self.assertRaisesRegex(ValueError, "one sbatch job"):
                p.configure_execution(args)

    def test_repository_no_archived_cluster_scripts(self):
        self.assertFalse((ROOT / "legacy").exists())
        for relative in ("PACHA", "examples/run_pacha.sbatch", "README.md"):
            content = (ROOT / relative).read_text()
            for prefix in ("/home/", "/group/", "/blue/", "/lustre/", "/scratch/"):
                self.assertNotIn(prefix, content, relative)


if __name__ == "__main__":
    unittest.main()
