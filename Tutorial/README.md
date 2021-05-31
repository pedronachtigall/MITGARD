Tutorial
========

A quick tutorial to use MITGARD.

The testing set contained in this repository was designed using the mitochondrial genome of *Crotalus adamanteus* available at NCBI (accession number NC_041524.1).
We used [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq) to simulate a NGS data with the following command `iss generate --genomes reference.fasta --n_reads 8000 --seed 5 --model hiseq --output AD5`.

Running Paired-End
==================
- Create the output directory and change to the output directory: `mkdir test_PE && cd test_PE`

- Run MITGARD in paired-end mode:
  ```
  MITGARD.py -s test -1 ../AD5_R1.fastq.gz -2 ../AD5_R2.fastq.gz -R ../reference.fasta -c 6 -M 8G
  ```

- Check the mitogenome assembly in the file ```test_mitogenome.fa```.

Running Single-End
==================
- If the user wants to use/test single-end mode, reads can concatenated together: `cat AD5_R1.fastq.gz AD5_R2.fastq.gz > AD5_merge.fastq.gz`

- Create the output directory and change to the output directory: `mkdir test_SE && cd test_SE`

- Run MITGARD in single-end mode:
  ```
  MITGARD.py -s test -S ../AD5_merge.fastq.gz -R ../reference.fasta -c 6 -M 8G
  ```

- Check the mitogenome assembly in the file ```test_mitogenome.fa```.

Checking for rearrangements in the gene order
=============================================

- Organelle genomes may present some rearrangements in the gene order, which can be checked by using the option ```--rearrangement True```.
- Run MITGARD in single-end mode with the rearrangement parameter activated:
  ```
  MITGARD.py -s test -S ../AD5_merge.fastq.gz -R reference.fasta -c 6 -M 8G --rearrangement True
  ```
- Check the mitogenome assembly in ```test_mitogenome.fa``` and the rearrangement check in the folder ```RearrangementCheck/```.

Checking the contigs file
=========================

- If the user wants to check all contigs generated and used to perform the mitogenome assembly, they are available at the file ```test_contigs.fasta```.
- The contigs file may be used to any downstream analysis.
