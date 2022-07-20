# ![](/Figures/matedb_dna4.png) MATEdb: Metazoan Assemblies from Transcriptomic Ensembles v1

[![License](https://img.shields.io/badge/license-GPLv3-blue.svg)](http://www.gnu.org/licenses/gpl.html)

MATEdb (Metazoan Assemblies from Transcriptomic Ensembles) is a curated database comprising 335 high-quality transcriptome assemblies from different animal phyla analyzed following the same pipeline. The repository is composed, for each species, of a *de novo* transcriptome assembly, its candidate coding regions within transcripts (both at the level of nucleotide and amino acid sequences), the longest isoform of the amino acid candidate coding regions, the gene content completeness score as assessed against the BUSCO database, the candidate coding regions filtered using their contamination profiles and an orthology-based gene annotation. We complement the repository with gene annotations from high-quality genomes totalling 423 high quality genomic and transcriptomic datasets. We invite the community to provide suggestions for new data sets and new annotation features to be included in subsequent versions, that will be permanently stored in public repositories. 

![Tree showing the diversity comprised in MATEdb](/Figures/MATEdb_Tree_v1.png)

Database DOI: [FigShare](https://doi.org/10.6084/m9.figshare.20178800.v1)

Description DOI: [bioRxiv](https://doi.org/10.1101/2022.07.18.500182)

## How to cite the data repository

Fernández, Rosa; Tonzo, Vanina; Simón Guerrero, Carolina; Lozano-Fernandez, Jesus ; Martínez-Redondo, Gemma I.; Balart-García, Pau; Aristide, Leandro; Eleftheriadi, Klara; Vargas-Chávez, Carlos (2022). MATEdb, a data repository of high quality metazoan transcriptome assemblies to accelerate phylogenomic studies. bioRxiv https://doi.org/10.1101/2022.07.18.500182.

## Data

MATEdb version 1 comprises 423 species, 322 arthropods (57 genomes and 265 transcriptomes) and 101 molluscs (31 genomes and 70 transcriptomes). These species belong to at least 394 different genra which are included in at least 325 families which are grouped in 110 orders which are further grouped in 27 different classes. 8,909,233 proteins are included in the final dataset with eggNOG annotation for 6,308,716 of them. In the Data folder you can find the **Table_S1.txt** file which contains the metadata for all species in the database. 
- Columns Phylum, Subphylum, Class, Subclass, Superorder, Order, Family, Genus and Species_or_Subspecies use the describe the taxonomy of each species using the names from NCBI Taxonomy.
- NCBI_Taxonomy_ID refers to the id used by [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy).
- Code refers to a four letter code that was selected for each species followed by a number.	Genome_or_Transcriptome	refers to whether the data available for that species comes from a genome or a transcriptome.
- Sequencing_Technology refers to the sequencing technology used to generate the data.
- Database and Accession_Number refer to the database used to obtain the raw data and the accession number under which it is stored.
- Reference refers to the publication associated with the data used.
- Alternative_repository refers to additional information sources. 
- BUSCO_Database refers to the BUSCO dataset that was used to estimate the genome completeness.
- Trinity_C_plus_F, Trinity_C, Trinity_S, Trinity_D, Trinity_F and Trinity_M refer to the % of Complete plus Fragmented, Complete, Complete Single-copy, Complete Duplicated, Fragmented and Missing genes from the Trinity assembly.
- Filtered_C_plus_F, Filtered_C, Filtered_S, Filtered_D, Filtered_F and Filtered_M refer to the % of Complete plus Fragmented, Complete, Complete Single-copy, Complete Duplicated, Fragmented and Missing genes after filtering the Trinity assembly using BlobTools.
- Protein_Number refers to the final number of proteins after keeping only the longest isoform for each gene.
- C_plus_F, C, S, D, F and M refer to the % of Complete plus Fragmented, Complete, Complete Single-copy, Complete Duplicated, Fragmented and Missing genes after conserving only the longest isoform for each gene.
- Eggnog refers to the number of proteins that have a functional annotation obtained with eggNOG-mapper. 


## Scripts

We provide all the commands and scripts needed to download the raw data and process it to obtain all files in MATEdb. Here you can see a diagram of the pipeline used to generate the data in MATEdb followed by a detailed description of every step. Software versions and parameters are summarized in the **Table_S2.xlsx** table in the Data folder.

![Pipeline used to generate the data in MATEdb](/Figures/MATEdb_Pipeline.png)

1.- Downloading the raw data

Either directly from a browser, going to a data repository:
- NCBI: https://www.ncbi.nlm.nih.gov
- Figshare: https://figshare.com
- etc...

Or using the SRA Toolkit version 2.10.7.
``` 
prefetch SRR1157986
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files $STORE/Craterostigmus_tasmanianus/SRR1157986 -O $STORE/Craterostigmus_tasmanianus/SRR1157986
```
For more information see https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump

2.- Quality check

To remove the adapters, clean and filter the data we use fastp 0.20.1. 
```
fastp --detect_adapter_for_pe -i $STORE/Craterostigmus_tasmanianus/SRR1157986/SRR1157986_1.fastq -o $STORE/Craterostigmus_tasmanianus/SRR1157986/SRR1157986_1_trimmed.fastq -I $STORE/Craterostigmus_tasmanianus/SRR1157986/SRR1157986_2.fastq -O $STORE/Craterostigmus_tasmanianus/SRR1157986/SRR1157986_2_trimmed.fastq
```

3.- Transcriptome assembly

The transcriptomes were assembled using Trinity 2.11.0.
```
Trinity --seqType fq --left $STORE/Craterostigmus_tasmanianus/SRR1157986/SRR1157986_1_trimmed.fastq --right $STORE/Craterostigmus_tasmanianus/SRR1157986/SRR1157986_2_trimmed.fastq --CPU 24 --max_memory 48G --output trinity --monitoring --no_parallel_norm_stats --full_cleanup --no_version_check
```

4.- Change filenames and headers

The name and the headers of the outputs of Trinity are modified using the following commands:
```
mv trinity.Trinity.fasta CTAS1.trinity.fasta
grep ">" CTAS1.trinity.fasta | sed 's/>//g' > original_headings_CTAS1.txt
paste trinity.Trinity.fasta.gene_trans_map original_headings_CTAS1.txt > CTAS1_conversion.txt
sed 's/TRINITY/CTAS1/g; s/ .*//g' CTAS1.trinity.fasta > CTAS1.mod.trinity.fasta
```

5.- Quality assessment of the assembly

The quality of the assembly was assessed with BUSCO 4.1.4 and using the arthropoda_odb10 for Arthropods or the metazoa_odb10 for Molluscs.
```
busco -i CTAS1.mod.trinity.fasta -l arthropoda_odb10 -o Busco -m transcriptome
```

Among the output files we focused on a file called short_summary.specific.*.Busco.txt where we can see the resulting statistics from comparing our assembly against the database used.
It gives the number of Complete BUSCOs (C), Complete and single-copy BUSCOs (S), Complete and duplicated BUSCOs (D), Fragmented BUSCOs (F) and Missing BUSCOs (M) and the fraction out of the total of BUSCO groups searched.
From these values we added that of C + F and if this sum is > 85% we considered that the assembly is of good quality.
Note: the limit value of 85% was decided by our group.

6.- Extract Open Reading Frames (ORFs)

The ORFs from the transcriptome assembly were identified using TransDecoder 5.5.0.
```
TransDecoder.LongOrfs -t CTAS1.mod.trinity.fasta
counts=$(grep -c '>' CTAS1.mod.trinity.fasta.transdecoder_dir/longest_orfs.pep)
TransDecoder.Predict -t CTAS1.mod.trinity.fasta -T $((counts/4))
```

7.- Elimination of foreign contaminant sequences

The taxonomy of the sequences in the TransDecoder output files was determined using BlobTools 2.3.3 and sequences which did not belong to the expected taxonomical group were discarded. First DIAMOND 2.0.8 was used to compare against the nr database which was downloaded in December 2020 from NCBI.
```
diamond blastp --query CTAS1.mod.trinity.fasta.transdecoder.pep --db nr.dmnd --sensitive --max-target-seqs 1 --evalue 1e-10 --threads 24 --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --out CTAS1.diamond.blastp.out
blobtools create --fasta CTAS1.mod.trinity.fasta.transdecoder.cds --hits CTAS1.diamond.blastp.out --hits-cols 1=qseqid,2=staxids,3=bitscore,5=sseqid,10=sstart,11=send,14=evalue --taxrule bestsum --taxdump taxdump BlobDir
python extract_phyla_for_blobtools.py BlobDir/bestsum_phylum.json | sed "s/', '/,/g" | tr -d "[]'" > contaminants.txt
PHYLA=$(cat contaminants.txt)
blobtools filter --param bestsum_phylum--Keys="$PHYLA" --taxrule bestsum --fasta CTAS1.mod.trinity.fasta.transdecoder.cds --summary STDOUT --summary-rank kingdom BlobDir >CTAS1.blobtools.summary
blobtools filter --param bestsum_phylum--Keys="$PHYLA" --taxrule bestsum --fasta CTAS1.mod.trinity.fasta.transdecoder.pep BlobDir
```

8.- BUSCO scores for the filtered sequences

BUSCO scores were calculated, as previously mentioned, to assess the effects of the filtering by taxonomy.

9.- Obtain the longest isoform for each gene

The **fetch_longest_iso.py** script was used to obtain the longest isoform for each gene. 
For the transcriptomes it is executed as follows:
```
python2.7 fetch_longest_iso.py -i CTAS1.mod.trinity.fasta.transdecoder.filtered.pep -o CTAS1.longiso.pep -t -l
```

For the genomes: 
First you need to download the matching annotation file in gff format for each genome. Then the **remove_isoforms_proteome.sh** script should be executed. Keep in mind that each annotation file is unique and the format may not match. In these cases, you have to modify the script accordingly after manually reviewing the gff file, the fasta and it is also useful to look at the geneprot.txt that we get when executing the script. By default, the script considers the structure of most gff files downloaded from NCBI.

```
remove_isoforms_proteome.sh -g species.gff -f species.aa.fasta -s SPEC1
```

10.- BUSCO scores for longest isoforms

BUSCO scores were calculated, as previously mentioned, to assess the effects of keeping only the longest isoform for each gene.

11.- Annotate the longest isoforms using eggNOG-mapper 

To annotate the longest isoforms we used eggNOG-mapper 2.1.6 against the eukaryota DIAMOND database (which can be downloaded using the download_eggnog_data.py script).
```
emapper.py -i CTAS1.longiso.pep -o CTAS1 --itype proteins --matrix BLOSUM62 --dmnd_db Eukaryota.dmnd --cpu 18 --dbmem --go_evidence all --output_dir results
```

## Docker

We also provide a docker container with the software that was used to generate all the files with the required versions.
- [MATEdb](https://hub.docker.com/repository/docker/vargaschavezc/matedb)

## Figures
- **MATEdb_Tree_v1.png** : tree in png format for visualization
- **matedb_dna4.png** : MATEdb logo in png format
- **MATEdb_Pipeline.png** : MATEdb pipeline in png format
