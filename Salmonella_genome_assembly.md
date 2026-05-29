# *Salmonella* Genome Assembly

### Notes
* When using Nextflow on Atlas, you need to enable remote proxies to download containers or data from the internet. The `nextflow.config` file must contain the line:

```
singularity {     
     envWhitelist = "http_proxy, https_proxy"}
```

### Downloading Bakta database
* The Bakta database downloaded, but I received the error:
Command error:
  AMRFinderPlus failed! amrfinder-error-code=1
  ERROR: AMRFinderPlus failed! command: 'amrfinder_update --force_update --database db/amrfinderplus-db', error code: 1
* The Bakta database downloaded to my work directory instead of `${projectDir}'.  I moved it. 
* I manually executed ARMFinderPlus in the project directory:
```
apptainer run apptainer_cache/quay.io-biocontainers-bakta-1.12.0--pyhdfd78af_0.img
amrfinder_update --force_update --database bakta_db/amrfinderplus-db
```
* This database only needs to be downloaded once so I will remove the `BAKTA_DOWNLOAD` process from the Nextflow pipeline. The script was:

```
process BAKTA_DOWNLOAD {
    container 'quay.io-biocontainers-bakta-1.12.0--pyhdfd78af_0'
    storeDir params.bakta_db_path // Persistent storage across runs

    output:
    path "db-full", type: 'dir'

    script:
    """
    bakta_db download --output . --type full
    """
}
```

### HiFi Adapter Filtration
* Attempting to run this in the Nextflow pipeline led to a variety of errors, primarily: "Missing output file(s)" that seems to be related to BusyBox `dirname FILENAME`.
* I ran `HiFiAdapterFilt` manually in interactive node.  Using no options seemed to be the only way to execute without errors.
  - Make soft links to the .bam files in the working directory.
  - Run the container in apptainer and execute `hifiadapterfilt.sh`.
```
ln -s /project/gbru_fy21_salm_compgen/GBRU_Transfers/SalmonellaPool/*.bam /project/gbru_fy21_salm_compgen/annette/salmonella_genome_assembly/work/bam_dat

cd /project/gbru_fy21_salm_compgen/annette/salmonella_genome_assembly/work/bam_dat

module load apptainer

apptainer run /project/gbru_fy21_salm_compgen/annette/salmonella_genome_assembly/apptainer_cache/quay.io-biocontainers-hifiadapterfilt-3.0.0--hdfd78af_0.img

hifiadapterfilt.sh 
```
* Each file prints the comments:
  * Converting [filename].bam to [filename].fastq on Mon Feb  9 18:40:01 CST 2026.
  * Converting [filename].bam to [filename].fasta on Mon Feb  9 18:40:01 CST 2026.
  * Identifying reads in [filename].bam with adapter contamination on Mon Feb  9 18:41:46 CST 2026.
  * BLAST Database error: Database memory map file error.
  * Creating blocklist for [filename].bam on Mon Feb  9 18:41:46 CST 2026.
  * Removing adapter contaminated reads from [filename].bam on Mon Feb  9 18:41:46 CST 2026.

## Authors:  Annette M. Hynes and Adam R. Rivers, USDA-ARS