# SediQuest

This pipeline is a friendly-user implementation of the Vernot et al. (2021) workflow. It aims to extract human reads from capture data to perform further analysis.

---

## Required Files

- **`config.yaml`** – Configuration file with paths and parameters.
- **`probeset.csv`** – Table describing the nuclear probes each sample was captured with. This table should include:
  - A BED file with target regions
  - A modified reference genome
  - An optional BED file with burden scores for filtering
- **`samples.csv`** – Table with sample information. Must include at least:
  - `LibraryID`
  - `Probeset`
  - Any other metadata you want to track

---

## Quick Start

### Step 0 – Verify Mapping
Ensure your genome is mapped to the modified reference genome corresponding, it should be a genome where the site captured where modified to avoid reference bias.

```bash
python check_reference_mapping.py  
```

### Step 1 – Basic processing 
Process your reads, duplication removal, quality filtering, length filtering, and deam filtering.

```bash
snakemake process_all --snakefile pipeline.local.v1.smk --configfile config/config.yaml --cores 25
```

### Step 2 – Kraken step 1 
Produce fasta file for Kraken 

```bash
snakemake kraken_step_1 --snakefile pipeline.local.v1.smk --configfile config/config.yaml --cores 25
```

### Step 3 – Kraken step 2
Run Kraken and produce summaries.

```bash
snakemake kraken_step_2 --snakefile pipeline.local.v1.smk --configfile config/config.yaml --cores 25
```

### Step 4 – Summaries
Produce a variety of plots and tables as a summary for your data

```bash
snakemake summaries --snakefile pipeline.local.v1.smk --configfile config/config.yaml --cores 25
```

### Step 5 – Look at your data!
You can have a look at the XXX plot to have an idea about the faunal contamination in you sample and decide if you want to use a kraken or a mammalian diversity score filtering. 


