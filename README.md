# bolt

`bolt` implements the UMCCR WGS post-processing logic and supporting functionality, and is used extensively throughout the
UMCCR workflow [`sash`](https://github.com/scwatts/sash).

The post-processing logic is organised into individual processing steps, which are collectively grouped according to the
type of data processed. Each processing step is accessible through specific `bolt` commands (see [Available
commands](#available-commands) for the complete listing) and can be run in isolation. Processing steps are integrated
into a production workflow in the `sash` pipeline.

## Table of contents

* [Available commands](#available-commands)
* [Docker images](#docker-images)
* [Usage](#usage)

## Available commands

### Germline small variants

| Command                      | Purpose                                     |
| ---                          | ---                                         |
| `bolt smlv_germline prepare` | Select variants given a genome regions file |
| `bolt smlv_germline report`  | Generate summary statistics and CPSR report |

### Somatic small variants

| Command                      | Purpose                                                       |
| ---                          | ---                                                           |
| `bolt smlv_somatic rescue`   | Supplement DRAGEN calls with SAGE calls and info              |
| `bolt smlv_somatic annotate` | Annotate variants with PCGR, PON, and defined genomic regions |
| `bolt smlv_somatic filter`   | Set and apply filters for variants                            |
| `bolt smlv_somatic report`   | Generate summary statistics and PCGR report                   |

### Somatic structural variants

| Command                      | Purpose                                                      |
| ---                          | ---                                                          |
| `bolt sv_somatic annotate`   | Annotate SVs (and CNVs) with SnpEff                          |
| `bolt sv_somatic prioritise` | Prioritise variants with rank ordering using `prioritize_sv` |

### Other

| Command                      | Purpose                                                 |
| ---                          | ---                                                     |
| `bolt other cancer_report`   | Generate the UMCCR Cancer Report                        |
| `bolt other multiqc_report`  | Generate the MultiQC report                             |
| `bolt other purple_baf_plot` | Create a modified PURPLE plot with β-allele frequencies |

## Docker images

For convenience, several Docker images are provided to run `bolt` aimed at balancing image size with the complexity of
the software environment and dependencies. Consequently, dependencies are split across six images:

| Name    | Docker image URI                    | Commands                       |
| ---     | ---                                 | ---                            |
| pcgr    | ghcr.io/scwatts/bolt:0.2.14-pcgr    | • `bolt smlv_germline report`<br />• `bolt smlv_somatic annotate`<br />• `bolt smlv_somatic report`<br /> |
| gpgr    | ghcr.io/scwatts/bolt:0.2.14-gpgr    | • `bolt other cancer_report`   |
| snpeff  | ghcr.io/scwatts/bolt:0.2.14-snpeff  | • `bolt sv_somatic annotate`   |
| circos  | ghcr.io/scwatts/bolt:0.2.14-circos  | • `bolt other purple_baf_plot` |
| multiqc | ghcr.io/scwatts/bolt:0.2.14-multiqc | • `bolt other multiqc_report`  |
| base    | ghcr.io/scwatts/bolt:0.2.14         | • `bolt smlv_germline prepare`<br />• `bolt smlv_somatic rescue`<br />• `bolt smlv_somatic filter`<br />• `bolt sv_somatic prioritise`<br /> |

## Usage

Given the nature of software dependencies required, it is strongly recommended to run `bolt` commands via the existing
[Docker images](#docker-images):

```bash
docker run -ti -v $(pwd):$(pwd) -w $(pwd) ghcr.io/scwatts/bolt:0.2.14 \
  bolt smlv_somatic filter \
    --tumor_name tumor_sample \
    --vcf_fp tumor_sample.vcf.gz \
    --output_dir ./
```
