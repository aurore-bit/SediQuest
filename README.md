# SediQuest

This pipeline is a friendly-user implementation of the Vernot et al. (2021) workflow. It aims to extract human reads from capture data to perform further analysis. The set of sites used to capture the DNA is called a probeset here.

---

## Required Files

- **`config.yaml`** – Configuration file with parameters, see config_example.yaml.

- **`probeset.csv`** – Table with information about the nuclear probes each sample was captured with. 
This table should include:
  - A BED file with target regions
  - A modified reference genome where you replaced the target sites by a "third" base (see Vernot et al. 2021).
  - A BED file with target sites.
  - An optional BED file with mammalian diversity scores for filtering.

| probeset   | path_to_ref | path_to_bed    | path_to_control     |
|--------|-----|-------------|-------------|
| 1240k  | whole_genome_modified_1240k.fa  | mammalian_diversity_score_1240k.bed      | 1240k.bed    |

- **`samples.csv`** – Table with sample information. Must include at least:
  - LibraryID
  - Probeset
  
| probeset_to_lib   | indexlibid  | 
|--------|-----|
| 1240k  | Lib.1608   |

---

## Start

### Step 0 – Verify Mapping
Ensure your genome is mapped to the modified reference genome corresponding.

```bash
python check_reference_mapping.py  
```

### Step 1 – Basic processing 
Process your reads, quality filtering, removal of duplicated, and deam filtering.

```bash
snakemake process_all --snakefile  SediQuest.smk --configfile config.yaml --cores 25
```

### Step 2 – Kraken step 1 
Produce fasta file for Kraken.

```bash
snakemake kraken_step_1 --snakefile  SediQuest.smk --configfile config.yaml --cores 25
```

### Step 3 – Kraken step 2
Run Kraken.

```bash
snakemake kraken_step_2 --snakefile  SediQuest.smk --configfile config.yaml --cores 25
```

### Step 4 – Summaries
Produce a variety of plots and tables as a summary for your data

```bash
snakemake summaries --snakefile  SediQuest.smk --configfile config.yaml --cores 25
```

### Step 5 – Look at your data!
You can have a look at the coverage plot to have an idea about the faunal contamination in your sample and decide if you want to use a kraken or a mammalian diversity score filtering. 


