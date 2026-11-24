# PACHA

**Pangenome Aware Chromosome-level Haplotype Assembler**

PACHA infers two reference-guided haplotype consensus sequences per diploid individual from paired-end short reads and chromosome-level haplotype assemblies. Version 0.3.0 supports local execution, automatic SLURM submission and execution inside a single existing SLURM allocation, with explicit reference selection and scalable source panels. Graph mode is implemented; VCF mode is reserved and explicitly rejected.

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

### A. Local workstation or server without SLURM

```bash
PACHA --executor local -m graph \
  -r /path/to/reads \
  -g /path/to/genomes \
  -ref /path/to/reference.fa \
  -t 16 -a 2 -at 8 \
  -o /path/to/output
```

This runs in the foreground until completion or failure. No scheduler is required. Use a persistent terminal session (such as tmux) if disconnecting from a remote server. Local execution uses one node and caps concurrent tasks according to the CPU affinity/available CPU budget. If -t or -at exceeds the available budget, PACHA stops before starting work. For instance, on a 16-CPU machine, -a 10 -at 8 permits at most two concurrent sample tasks.

For explicit memory-aware scheduling, add --local-memory 128G together with realistic estimates for every stage using --mem and --stage-mem. For example, if justified by a pilot dataset, use --mem 32G --stage-mem graph=96G --stage-mem index=96G. These values are admission-control estimates, not operating-system memory limits. Without --local-memory, PACHA does not enforce memory-based concurrency. Do not lower estimates simply to bypass resource checks. Large genomes may need far more RAM than a workstation has.

### B. Submit the complete pipeline to SLURM

Run on your cluster's submission/login node, with the Conda environment activated:

```bash
PACHA --executor slurm -m graph \
  -r /path/to/reads \
  -g /path/to/genomes \
  -ref /path/to/reference.fa \
  -t 32 -a 10 -at 16 \
  -o /path/to/output \
  --partition normal --account mylab --qos normal \
  --time 48:00:00 --mem 128G \
  --stage-mem graph=200G --stage-mem index=256G \
  --stage-mem align=300G \
  --stage-time align=72:00:00
```

Replace the example partition, account and QoS with actual site values; omit those that are unnecessary. PACHA submits dependent jobs and arrays, then returns. Jobs continue independently of the terminal. Each stage starts only after its preceding stage succeeds. CPU and memory requests are per task, not per array. Adjust the time limits to comply with your partition policy.

### C. Run PACHA inside your own sbatch script

Use this option if you want the entire pipeline contained within one SLURM allocation, using your own scheduler directives. Save as run_pacha.sbatch and adjust the paths and resource requests to suit your cluster and data.

```bash
#!/bin/bash
#SBATCH --job-name=PACHA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=4-00:00:00
#SBATCH --partition=normal
#SBATCH --account=mylab
#SBATCH --output=pacha-%j.log

set -euo pipefail

source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate pacha

/path/to/PACHA/PACHA --executor local -m graph \
  -r /path/to/reads \
  -g /path/to/genomes \
  -ref /path/to/reference.fa \
  -t 32 -a 4 -at 8 \
  -o /path/to/output \
  --local-cpus 32 --local-memory 256G \
  --mem 128G --stage-mem graph=200G \
  --stage-mem index=256G --stage-mem align=128G
```

```bash
sbatch run_pacha.sbatch
```

This does not submit child jobs. CPU budget is limited by SLURM_CPUS_PER_TASK and CPU affinity; the allocation must cover the full workflow runtime and peak memory. The example memory estimates must be validated for your dataset. SLURM enforces the allocation's actual memory and time limits; PACHA's local-memory option only limits simultaneous tasks according to estimates. Request one node and one task. Do not use srun with multiple tasks or launch PACHA as a SLURM array itself. If a site's maximum walltime is too short for the full pipeline, prefer option B; alternatively resubmit the same script to reuse completed tasks. A stage task that exceeds the walltime must be restarted.

If you instead explicitly request --executor slurm inside a batch script, PACHA submits independent child jobs and exits after submission. The parent allocation does not supply resources to those jobs and its #SBATCH settings are not inherited by PACHA; pass the desired child resources on the command line. Avoid this approach unless the cluster permits nested submissions.

## Command-line options

| Flag | Meaning |
|---|---|
| -m | graph (implemented); vcf (future) |
| -r | Directory containing paired-end FASTQ files |
| -g | Directory containing one chromosome-level haplotype per FASTA |
| -ref, --reference | **Required** reference FASTA path, inside or outside -g |
| -o | Output directory (new or existing matching PACHA run) |
| --executor | auto (default), local, slurm |
| -t | CPUs per graph/index task (SLURM default 16; local default up to 16) |
| -a | Maximum concurrent sample tasks (SLURM default 4; local default 1) |
| -at | CPUs per sample task (SLURM default 8; local default up to 8) |
| --prepare-parallel | Maximum concurrent source FASTA preparation tasks (default 8) |
| --graph-parallel | Maximum concurrent chromosome PGGB tasks (default 1) |
| --local-cpus | Total local CPU budget, at most available/allocated CPUs |
| --local-memory | Optional total local RAM admission budget |
| --array-batch-size | Maximum tasks per SLURM array (default 1000) |
| --dry-run | Print execution plan without creating output or starting jobs |
| --validate-genomes | Exhaustively validate source FASTA chromosome headers before starting |

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

There is no fixed number of input haplotype assemblies or target individuals. PGGB receives the actual number of source haplotypes (including the chosen reference); each chromosome is built from all source FASTAs, with one source path per FASTA. The preparation array normalizes and indexes each FASTA independently. FASTA extraction for chromosome graphs streams one source at a time. Large arrays are split into sequential, dependency-linked batches, preserving global task indices and concurrency limits. Adjust --array-batch-size to your cluster's MaxArraySize, and --prepare-parallel, --graph-parallel, -a and memory/time options to the resources available.

For example, 200 source FASTAs in -g plus a separate reference not present in that directory give 201 source haplotypes. With 500 paired-end target individuals, a successful run publishes exactly 1,000 FASTA files: hap1 and hap2 for each individual. The number of output haplotypes is determined by the diploid target samples, not the size of the source panel. PGGB graph construction and VG whole-genome indexing still have panel-dependent memory, CPU, disk and walltime costs. There is no claim that arbitrarily large graphs fit the default resources.

## Input contract

* One FASTA file per source haplotype, including a single haploid reference FASTA specified by -ref. For phased diploid source genomes, provide each haplotype as a separate FASTA. All must have the same set of chromosome identifiers as the reference. Chromosome count and sample count are inferred, not hard-coded. Use consistent, unique chromosome names (e.g. Chr01); no whitespace or # in the identifier. Descriptions following the first header token are ignored. FASTA extensions: .fa, .fasta, .fna, optionally .gz.
* Paired reads directly in -r, one pair per sample. Supported patterns: SAMPLE_R1.fastq.gz / SAMPLE_R2.fastq.gz, SAMPLE_R1_001.fastq.gz / SAMPLE_R2_001.fastq.gz, SAMPLE_1.fq.gz / SAMPLE_2.fq.gz, and legacy SAMPLE_1.trimmed.fastq.gz / SAMPLE_3.trimmed.fastq.gz. Uncompressed files and .fq are supported. Merge sequencing lanes beforehand. Mixed or ambiguous pairs fail validation.
* Inputs must be visible at identical absolute paths on compute nodes; use a shared filesystem. Do not edit inputs during a run. PACHA checks file sizes and modification times, not full input content hashes.
* Assemblies must use homologous chromosome labels and comparable orientation. PACHA does not discover orthology or rename noncorresponding chromosomes. Diploid target samples only.

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

Unit tests cover required reference selection, external references, reference aliases, a synthetic 1,003-source panel, array batching, selected reference preparation, two outputs per individual, input discovery, command construction, dry run, receipts, consensus genotype selection and mocked SLURM submission. Additional tests cover executor detection, allocated CPU limits, memory-aware concurrency, resource overrides, actual local subprocess orchestration using fake workers, restart, failure propagation and stage ordering. Real scientific end-to-end validation requires the installed bioinformatics tools; SLURM integration also requires a real cluster. The original six research scripts are retained under legacy/ for provenance; they are not invoked by PACHA.

