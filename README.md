# PACHA

**Pangenome Aware Chromosome-level Haplotype Assembler**

PACHA infers two reference-guided haplotype consensus sequences per diploid individual from paired-end short reads and chromosome-level haplotype assemblies. Version 0.3.1 supports local execution, automatic SLURM submission and execution inside a single existing SLURM allocation, with explicit reference selection and scalable source panels. Graph mode is implemented; VCF mode is reserved and explicitly rejected.

## Installation (Linux)

```bash
git clone https://github.com/Luna-san-2911/PACHA.git
cd PACHA
conda env create --file environment.yml
conda activate pacha
export PATH="$PWD:$PATH"
PACHA --help
```

Add the absolute clone directory to PATH in your shell startup file if desired. Linux with sufficient storage, RAM and CPUs is required; SLURM is optional. The environment pins PGGB, VG and ODGI to the versions used together in the PGGB Bioconda recipe. The remaining dependencies use compatible version ranges; the environment is not a platform-independent lockfile.

## Quick start: choose how to run

The scientific workflow is identical in every execution mode. The default --executor auto selects:

| Environment | Automatic behavior |
|---|---|
| Existing SLURM allocation (sbatch/salloc) | Execute locally within the current allocation |
| Outside an allocation, sbatch available | Submit dependent SLURM jobs/arrays |
| No sbatch available | Run locally in the terminal |

Specify --executor local or --executor slurm to override detection. The executor is independent of -m graph (the scientific method).

**Sample concurrency:** use --concurrent-samples N (short form -a) to cap the number of individuals running simultaneously in each alignment, variant calling, phasing and consensus stage. Use --sample-threads T (short form -at) for CPUs per individual. Stages do not overlap. These are the same settings in every executor; users do not need to calculate or specify array indices.

### A. Local workstation or server without SLURM

```bash
PACHA --executor local -m graph \
  -r /path/to/reads \
  -g /path/to/genomes \
  -ref /path/to/reference.fa \
  -t 16 --concurrent-samples 2 --sample-threads 8 \
  -o /path/to/output
```

This runs in the foreground until completion or failure. No scheduler is required. Use a persistent terminal session (such as tmux) if disconnecting from a remote server. Local execution uses one node and caps concurrent tasks according to the CPU affinity/available CPU budget. If -t or -at exceeds the available budget, PACHA stops before starting work. For instance, on a 16-CPU machine, --concurrent-samples 10 --sample-threads 8 permits at most two concurrent sample tasks. PACHA starts the next waiting individual as a slot becomes free.

For explicit memory-aware scheduling, add --local-memory 128G together with realistic estimates for every stage using --mem and --stage-mem. For example, if justified by a pilot dataset, use --mem 32G --stage-mem graph=96G --stage-mem index=96G. These values are admission-control estimates, not operating-system memory limits. Without --local-memory, PACHA does not enforce memory-based concurrency. Do not lower estimates simply to bypass resource checks. Large genomes may need far more RAM than a workstation has.

### B. Submit the complete pipeline to SLURM

Run on your cluster's submission/login node, with the Conda environment activated:

```bash
PACHA --executor slurm -m graph \
  -r /path/to/reads \
  -g /path/to/genomes \
  -ref /path/to/reference.fa \
  -t 32 --concurrent-samples 10 --sample-threads 16 \
  -o /path/to/output \
  --partition normal --account mylab --qos normal \
  --time 48:00:00 --mem 128G \
  --stage-mem graph=200G --stage-mem index=256G \
  --stage-mem align=300G \
  --stage-time align=72:00:00
```

Replace the example partition, account and QoS with actual site values; omit those that are unnecessary. PACHA submits dependent jobs and arrays, then returns. Jobs continue independently of the terminal. Each stage starts only after its preceding stage succeeds. CPU and memory requests are per task, not per array. Adjust the time limits to comply with your partition policy.

For 500 individuals with --concurrent-samples 10, PACHA uses a native SLURM array such as --array=0-499%10 for alignment and then another for variant calling. At most 10 sample tasks run at once across these sequential stages, but other array elements may remain pending. The limit applies to this PACHA run only, not to other jobs or total queued-job quotas. If your cluster limits the size of submitted arrays, set --array-batch-size (e.g. 100); batches are dependency-linked so that they cannot overlap. All stage jobs are submitted up front in this mode.

### C. Run PACHA inside your own sbatch script

**Option C submits exactly one SLURM job: your outer sbatch job. PACHA does NOT submit any child sbatch scripts or native SLURM arrays.** Inside that job, PACHA starts ordinary local worker processes as slots become available. The existing allocation supplies all CPU, memory and walltime resources. The sample-concurrency flag caps those processes just as it caps running array elements in option B.

A ready-to-use generic wrapper is included in examples/run_pacha.sbatch. Activate the environment and add the repository to PATH as shown in Installation, then from the repository directory run:

```bash
sbatch --partition=normal --account=mylab \
  --cpus-per-task=32 --mem=256G --time=48:00:00 \
  examples/run_pacha.sbatch \
  -m graph -r /path/to/reads -g /path/to/genomes \
  -ref /path/to/reference.fa -o /path/to/output \
  -t 32 --concurrent-samples 4 --sample-threads 8
```

The wrapper forwards every user argument to PACHA and forces --executor local. No paths, sample names, accounts or partitions are hard-coded in it. Its default resource directives can be overridden on the sbatch command line. If needed, add your cluster-specific Conda initialization to the wrapper.

Alternatively, write your own sbatch script:

```bash
#!/bin/bash
#SBATCH --job-name=PACHA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --output=pacha-%j.log

set -euo pipefail
# Activate the pacha environment here if not already activated.
PACHA --executor local -m graph \
  -r /path/to/reads -g /path/to/genomes \
  -ref /path/to/reference.fa -o /path/to/output \
  -t 32 --concurrent-samples 4 --sample-threads 8
```

Submit with sbatch run_pacha.sbatch. Add your site's partition, account, QoS and other directives as needed.

With 32 allocated CPUs, four sample tasks using eight CPUs each can run concurrently. PACHA will reduce concurrency if CPUs are insufficient. To also cap concurrency based on memory, supply --local-memory 256G and realistic per-task memory estimates with --mem/--stage-mem (e.g. --stage-mem align=64G permits at most four such alignment tasks within that budget, subject to CPU limits). The estimates are not hard limits; SLURM enforces the actual job allocation. These figures are examples, not universal resource recommendations.

Use exactly one node and one task, and do not put the outer script in a SLURM array. The complete workflow must fit the allocation's walltime and memory; if it times out, resubmit the same command to reuse completed stages. An individual stage task that times out must restart.

For multiple nodes or separately scheduled sample arrays, use option B. Explicit --executor slurm inside your own batch script is an opt-in to submitting independent child jobs, not to using the parent allocation; supply child resource flags separately and ensure nested submission is permitted. It is not needed for option C.

## Command-line options

| Flag | Meaning |
|---|---|
| -m | graph (implemented); vcf (future) |
| -r | Directory containing paired-end FASTQ files (automatic pairing) |
| --sample-sheet | TSV mapping arbitrary FASTQ filenames to sample IDs; alternative to -r |
| -g | Directory containing one chromosome-level haplotype per FASTA |
| -ref, --reference | **Required** reference FASTA path, inside or outside -g |
| -o | Output directory (new or existing matching PACHA run) |
| --executor | auto (default), local, slurm |
| -t | CPUs per graph/index task (SLURM default 16; local default up to 16) |
| --concurrent-samples, -a | Maximum running individuals per sample stage (SLURM default 4; local default 1) |
| --sample-threads, -at | CPUs per individual task (SLURM default 8; local default up to 8) |
| --prepare-parallel | Maximum concurrent source FASTA preparation tasks (default 8) |
| --graph-parallel | Maximum concurrent chromosome PGGB tasks (default 1) |
| --local-cpus | Total local CPU budget, at most available/allocated CPUs |
| --local-memory | Optional total local RAM admission budget |
| --array-batch-size | Maximum tasks per SLURM array (default 1000) |
| --dry-run | Print execution plan without creating output or starting jobs |
| --validate-genomes | Exhaustively validate source FASTA chromosome headers before starting |

The old --array and --array-threads spellings remain supported as compatibility aliases. The recommended names describe sample concurrency without implying that local execution requires a scheduler.

The reference is included in the pangenome exactly once by file identity: if it is outside -g and not already represented by the same file, it is added to the source panel. Symlink and hardlink aliases are recognized; distinct files containing identical sequences are not collapsed. No alphabetical fallback is used.

### SLURM resources

| Option | Purpose |
|---|---|
| --partition, --account, --qos | Site-specific allocation and scheduling |
| --mem, --time | Default memory and walltime for each task/job |
| --stage-mem STAGE=VALUE | Per-stage memory override; repeat as needed |
| --stage-time STAGE=VALUE | Per-stage walltime override |
| --stage-partition STAGE=VALUE | Per-stage partition |
| --stage-account STAGE=VALUE | Per-stage account |
| --stage-qos STAGE=VALUE | Per-stage QoS |
| --constraint, --reservation | Node features or reservation |
| --mail-user, --mail-type | Optional scheduler email notifications |
| --slurm-option=OPTION | Additional allowed sbatch options (e.g. --exclusive, --licenses=foo:1) |

Stage names: prepare, graph, index, align, call, phase, consensus, finish. Memory precedence is --stage-mem, then --STAGE-mem (where provided), then --mem, then built-in defaults. Built-in per-task memory estimates/requests are prepare 16G, graph 200G, index 256G, align 300G, call/phase 128G, consensus 16G, finish 2G. Time precedence is --stage-time, then --graph-time (prepare/graph/index) or --sample-time (align/call/phase/consensus), then --time, then built-in defaults (48h graph group, 72h sample group, 1h finish). Local mode uses memory estimates but rejects SLURM-only submission flags; set those in your outer sbatch script instead. Additional allowed sbatch options are exclusive, licenses, nice, tmp, comment, nodelist, exclude, prefer, wckey, gres and begin. PACHA owns dependency, array, CPU, memory, logging, export and command options to prevent conflicting overrides. Some tools (notably WhatsHap phase and BCFtools consensus) are not fully multithreaded.

Always start with --dry-run to review task counts, effective concurrency, resources and the chosen reference. It reads reference headers and discovers inputs without writing or submitting anything; full genome validation occurs during preparation unless --validate-genomes is supplied.

## Large pangenome panels

There is no fixed number of input haplotype assemblies or target individuals. PGGB receives the actual number of source haplotypes (including the chosen reference); each chromosome is built from all source FASTAs, with one source path per FASTA. The preparation array normalizes and indexes each FASTA independently. FASTA extraction for chromosome graphs streams one source at a time. Large arrays are split into sequential, dependency-linked batches, preserving global task indices and concurrency limits. Adjust --array-batch-size to your cluster's MaxArraySize, and --prepare-parallel, --graph-parallel, --concurrent-samples and memory/time options to the resources available.

For example, 200 source FASTAs in -g plus a separate reference not present in that directory give 201 source haplotypes. With 500 paired-end target individuals, a successful run publishes exactly 1,000 FASTA files: hap1 and hap2 for each individual. The number of output haplotypes is determined by the diploid target samples, not the size of the source panel. PGGB graph construction and VG whole-genome indexing still have panel-dependent memory, CPU, disk and walltime costs. There is no claim that arbitrarily large graphs fit the default resources.

## Input contract

* One FASTA file per source haplotype, including a single haploid reference FASTA specified by -ref. For phased diploid source genomes, provide each haplotype as a separate FASTA. All must have the same set of chromosome identifiers as the reference. Chromosome count and sample count are inferred, not hard-coded. Use consistent, unique chromosome names (e.g. Chr01); no whitespace or # in the identifier. Descriptions following the first header token are ignored. FASTA extensions: .fa, .fasta, .fna, optionally .gz.
* Paired reads directly in -r, one pair per sample. Supported patterns: SAMPLE_R1.fastq.gz / SAMPLE_R2.fastq.gz, SAMPLE_R1_001.fastq.gz / SAMPLE_R2_001.fastq.gz, SAMPLE_1.fq.gz / SAMPLE_2.fq.gz, and legacy SAMPLE_1.trimmed.fastq.gz / SAMPLE_3.trimmed.fastq.gz. Uncompressed files and .fq are supported. Merge sequencing lanes beforehand. Mixed or ambiguous pairs fail validation.
* If FASTQ filenames use nonstandard naming or files are in different directories, use --sample-sheet instead of -r. See the example below.
* Inputs must be visible at identical absolute paths on compute nodes; use a shared filesystem. Do not edit inputs during a run. PACHA checks file sizes and modification times, not full input content hashes.
* Assemblies must use homologous chromosome labels and comparable orientation. PACHA does not discover orthology or rename noncorresponding chromosomes. Diploid target samples only.

### Arbitrary sample names and FASTQ filenames

For standard filenames, PACHA discovers each individual's ID from the paired FASTQ names. For nonstandard filenames, create a tab-separated samples.tsv with exactly three columns:

```text
sample	r1	r2
individual_A	reads/forward_file.fq.gz	reads/reverse_file.fq.gz
individual_B	reads/another_forward.fastq.gz	reads/another_reverse.fastq.gz
```

Relative read paths are resolved relative to the sample sheet. Absolute paths also work. Use unique sample IDs consisting of letters, digits, dots, underscores or hyphens, beginning with a letter or digit. Each row must identify one paired-end individual; merge sequencing lanes beforehand. Duplicate sample IDs and reused read files are rejected.

```bash
PACHA --executor local -m graph \
  --sample-sheet samples.tsv \
  -g /path/to/genomes -ref /path/to/reference.fa \
  -t 16 --concurrent-samples 2 --sample-threads 8 \
  -o /path/to/output
```

All sample paths are determined from the run manifest. For example, individual_A receives results/individual_A_hap1.fa.gz and results/individual_A_hap2.fa.gz. PACHA has no dataset-specific paths or sample IDs.

## Workflow

1. Normalize and index input FASTAs in parallel, validate chromosome names, assign stable PanSN source identifiers, and retain the original filename mapping and explicit reference choice in the manifest. Only the selected reference preparation task publishes the shared reference FASTA and reference-path list.
2. Build one PGGB graph per chromosome with the original defaults (segment 20 kb, 90% mapping identity, asm20 POA, minimum match 47, POA lengths 10000,15000, padding 0.03).
3. Merge exactly one final ODGI graph per chromosome; copy the selected reference paths, deconstruct using VG, verify REF alleles, and rebuild Giraffe indexes from the reference FASTA and graph-derived VCF. This is **graph mode**, even though the indexed graph is reconstructed from VCF.
4. Align reads with VG Giraffe, filter at 0.90 normalized substitution-based score and MAPQ 20, and generate packs using the same GBZ.
5. Call diploid variants with VG call (including reference calls); filter QUAL >=20, GQ >=20, DP >=5 and missing genotypes.
6. Surject into the selected reference coordinates, validate BAM chromosome names and lengths, phase with WhatsHap and report phase statistics.
7. Reconstruct two FASTAs with BCFtools consensus, preserving high-confidence homozygous calls (including unphased homozygous ALT) and phased heterozygotes. Unphased heterozygotes and missing calls are excluded. Unsupported symbolic selected alleles and reported overlapping/skipped variants fail explicitly.
8. Verify all outputs and publish final FASTAs and summary.

## Output layout

```text
output/
  .pacha/manifest.json        # input fingerprint, reference and scientific settings
  .pacha/worker.py            # immutable worker snapshot
  .pacha/jobs.json            # latest submitted SLURM job IDs
  logs/
  00_inputs/<source_id>/     # prepared genome.fa and index; chosen source also has reference.fa
  01_graph/<chromosome>/
  02_index/
  03_alignment/<sample>/
  04_variants/<sample>/
  05_phasing/<sample>/
  06_consensus/<sample>/
  results/
    <sample>_hap1.fa.gz
    <sample>_hap2.fa.gz
    summary.tsv
    _SUCCESS
```

For SLURM submission mode, monitor with `squeue -u "$USER"` and `sacct` where enabled. Local mode prints progress in the terminal and writes per-task logs under logs/. Inspect logs and hidden .incomplete-* directories for failed tasks. _SUCCESS is created only after all sample consensus outputs pass validation. A stage is published via atomic directory rename and a receipt records expected output sizes. These receipts detect missing or size-changed files, not same-size corruption.

### Restart

After all prior tasks are terminal (cancel pending SLURM jobs with scancel if necessary), rerun the same command with the same output path. Local mode exits nonzero on task failure and stops launching new tasks in that stage; already running tasks may finish. Ctrl+C/SIGTERM stops active local worker groups. A hard kill or node failure may leave orphan workers or incomplete temporary directories; verify prior tasks are stopped before restarting. Completed stages are reused; failed unpublished stages restart. You may change scheduling resources but not inputs, scientific settings or the executable version. Do not manually delete individual completed intermediates; use a fresh run directory for changes (including a different reference or executable version) or corrupted outputs. A partial submission leaves already submitted jobs running; the launcher refuses a duplicate submission while matching jobs remain active.

## Scientific limitations

These are **reference-guided chromosome-length consensus sequences**, not validated chromosome-wide biological haplotypes. Short reads generally phase local blocks only; phase orientation may switch between blocks and chromosomes. No parental assignment or phase-block stitching is currently implemented. Uncalled, filtered and unphased heterozygous positions retain reference bases; they are not evidence of homozygous reference genotypes. Graph-to-VCF reconstruction can lose complex variation, and variants not represented in the graph may be missed. The summary reports genotype selection counts, not genome-wide accuracy or completeness. Use the phased VCF and WhatsHap statistics to assess phase evidence.

## Tests

```bash
python -m unittest discover -s tests -v
```

Unit tests cover required reference selection, external references, reference aliases, a synthetic 1,003-source panel, array batching, selected reference preparation, two outputs per individual, input discovery, command construction, dry run, receipts, consensus genotype selection and mocked SLURM submission. Additional tests cover executor detection, allocated CPU limits, memory-aware concurrency, resource overrides, actual local subprocess orchestration using fake workers, restart, failure propagation and stage ordering. Real scientific end-to-end validation requires the installed bioinformatics tools; SLURM integration also requires a real cluster. The original dataset-specific research scripts have been removed from the current tree. They remain recoverable from [the historical revision](https://github.com/Luna-san-2911/PACHA/tree/be57159017ab6460ff99e86b0d6640014d70ee0c) for provenance. PACHA is the only active pipeline implementation. This release also tests arbitrary FASTQ filenames, sample-sheet validation, native SLURM throttling, no child submissions in an existing allocation and absence of legacy scripts.

