# OSC Genome Heterozygous TE Insertion Pipeline

Pipeline for re-inserting heterozygous transposable elements (TEs) that were removed during haplotig purging but are present in the OSC cell genome.

Part of the **Handler et al., 2025** publication:

**The Drosophila OSC Genome: A Resource for Studies of Transposon and piRNA Biology**

## Overview

This repository contains the workflow for identifying and re-inserting heterozygous transposable element insertions that were lost during the haplotig purging step. During purging, one haplotype is selected as the primary assembly, which can result in the loss of TE insertions present only in the alternative haplotype. This pipeline recovers these heterozygous TEs to create a more complete representation of the transposon landscape in the OSC cell population.

## Repository Structure

```
├── script-files/        # Core TE insertion scripts
├── utility-files/       # Helper scripts and configuration files
└── insert_het-TEs.sh    # Main submission script for TE insertion
```

## Pipeline Components

### Heterozygous TE Detection
Scripts for identifying TE insertions present in purged haplotigs but absent from the primary assembly.

### Insertion Site Validation
Tools for validating heterozygous TE insertions using read coverage and structural variant analysis.

### TE Re-insertion
Methods for incorporating validated heterozygous TEs back into the primary assembly at appropriate genomic locations.

### Quality Control
Scripts to verify successful TE insertions and assess their impact on assembly quality.

## Requirements

- Apptainer
  
## Usage

Run the main TE insertion pipeline using:

```bash
bash insert_het-TEs.sh
```

Adjust parameters in the script files based on your TE detection thresholds and validation criteria.

## Output

The pipeline produces:
- Updated genome assembly with heterozygous TEs re-inserted
- Annotation of heterozygous TE insertion sites
- Statistics on recovered TE insertions

## Related Resources

### Main Publication Repository
https://github.com/BrenneckeLab/Handler_2025-OSC-genome


### UCSC Genome Browser Hub
https://genome-euro.ucsc.edu/s/Brennecke%2DLab/OSC_r1.01_Handler_et.al._2025

## Citation

Please find the proper citation in https://github.com/BrenneckeLab/Handler_2025-OSC-genome

## Contact

For questions or additional information, please contact:
dominik.handler@imba.oeaw.ac.at

## License

MIT License
