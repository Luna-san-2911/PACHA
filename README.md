# PACHA

**Pangenome Aware Chromosome-level Haplotype Assembler**

PACHA infers two reference-guided haplotype consensus sequences per diploid individual from paired-end short reads and chromosome-level haplotype assemblies. Version 0.2.1 requires an explicit reference FASTA and supports large input panels through parallel genome preparation and batched SLURM arrays. Graph mode is implemented; VCF mode is reserved and explicitly rejected.

## Installation (Linux, SLURM)

```bash
git clone https://github.com/Luna-san-2911/PACHA.git
cd PACHA
conda env create --file environment.yml
conda activate pacha
export PATH="$PWD:$PATH"
PACHA --help
```

Add the absolute clone directory to PATH in your shell startup file if desired. SLURM (sbatch and squeue) must already be configured by your HPC administrator. The environment pins PGGB, VG and ODGI to the versions used together in the PGGB Bioconda recipe. The remaining dependencies use compatible version ranges; the environment is not a platform-independent lockfile.

## One command

```bash
PACHA -m graph \
  -r /path/to/reads \
  -g /path/to/genomes \
  -t 32 -a 10 -at 16 \
  -ref /path/to/reference/Patl.hap2.fa \
  -o /path/to/output
```

The command discovers the inputs, validates the reference headers, creates a manifest and worker snapshot, and submits eight dependent stages (one or more SLURM jobs/arrays per stage). It returns after submission; jobs continue after logout. Use --dry-run to inspect the input layout and job plan without creating outputs or submitting jobs. Add --validate-genomes for an exhaustive chromosome-header preflight; otherwise source FASTA compatibility is validated on compute nodes during preparation, before any PGGB job can start.

| Flag | Meaning |
|---|---|
| -m | graph (implemented); vcf (future) |
| -r | Directory containing paired-end FASTQ files |
| -g | Directory containing one chromosome-level haplotype per FASTA |
| -t | CPUs per chromosome graph job and graph indexing job |
| -a | Maximum concurrent sample tasks in each sample array |
| -at | CPUs requested per sample task |
| -o | Output directory (new or existing matching PACHA run) |
| -ref, --reference | **Required** reference FASTA path, inside or outside -g; used as the graph, VCF, phasing and consensus backbone |
| --prepare-parallel | Concurrent source FASTA preparation tasks (default 8) |
| --array-batch-size | Maximum tasks per individual SLURM array submission (default 1000) |
| --graph-parallel | Concurrent chromosome PGGB tasks (default 1) |

The reference is included in the pangenome exactly once by file identity: if it is outside -g and not already represented by the same file, it is added to the source panel. Symlink and hardlink aliases of the same file are recognized. Distinct files containing identical sequences are not collapsed, since they may represent different source haplotypes. No alphabetical reference fallback is used. Optional --partition, --account and --qos select site resources. Defaults use the cluster's own SLURM settings. Per-job memory defaults: preparation 16G, graph 200G, index 256G, align 300G, call/phase 128G, consensus 16G. Override with --prepare-mem, --graph-mem, --index-mem, --align-mem, --call-mem, --phase-mem, --consensus-mem. --graph-time and --sample-time adjust walltime. Sample stages run sequentially as arrays; e.g., -a 10 -at 16 permits at most 160 allocated CPUs during each sample stage. Some tools (notably WhatsHap phase and BCFtools consensus) are not fully multithreaded.

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

Monitor with `squeue -u "$USER"` and `sacct` where enabled. Inspect logs and hidden .incomplete-* directories for failed tasks. _SUCCESS is created only after all sample consensus outputs pass validation. A stage is published via atomic directory rename and a receipt records expected output sizes. These receipts detect missing or size-changed files, not same-size corruption.

### Restart

After all prior jobs are terminal (cancel pending jobs with scancel if necessary), rerun the same command with the same output path. Completed stages are reused; failed unpublished stages restart. You may change scheduling resources but not inputs, scientific settings or the executable version. Do not manually delete individual completed intermediates; use a fresh run directory for changes (including a different reference or upgrade from v0.2.0) or corrupted outputs. A partial submission leaves already submitted jobs running; the launcher refuses a duplicate submission while matching jobs remain active.

## Scientific limitations

These are **reference-guided chromosome-length consensus sequences**, not validated chromosome-wide biological haplotypes. Short reads generally phase local blocks only; phase orientation may switch between blocks and chromosomes. No parental assignment or phase-block stitching is currently implemented. Uncalled, filtered and unphased heterozygous positions retain reference bases; they are not evidence of homozygous reference genotypes. Graph-to-VCF reconstruction can lose complex variation, and variants not represented in the graph may be missed. The summary reports genotype selection counts, not genome-wide accuracy or completeness. Use the phased VCF and WhatsHap statistics to assess phase evidence.

## Tests

```bash
python -m unittest discover -s tests -v
```

Unit tests cover required reference selection, external references, reference aliases, a synthetic 1,003-source panel, array batching, selected reference preparation, two outputs per individual, input discovery, command construction, dry run, receipts, consensus genotype selection and mocked SLURM submission. Real end-to-end validation requires the installed bioinformatics tools and a SLURM cluster. The original six research scripts are retained under legacy/ for provenance; they are not invoked by PACHA.

