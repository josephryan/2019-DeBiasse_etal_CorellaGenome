# PLANNED ANALYSES FOR ASSEMBLING THE CORELLA INFLATA GENOME AND TRANSCRIPTOME AND INFERRING CHORDATE/TUNICATE PHYLOGENETIC RELATIONSHIPS 
 Principle Investigator: Joseph Ryan, Bradley Davidson  
 Support Personnel: Melissa DeBiasse, William Colgan  
 Draft or Version Number: v.1.0  
 Date: 6 Feb 2019  
 Note: this document will be updated (updates will be tracked through github)
 
## 1 INTRODUCTION: BACKGROUND INFORMATION AND SCIENTIFIC RATIONALE  

### 1.1 _Background Information_  

As the closest living relatives of vertebrates, tunicates are important for understanding vertebrate origins and identifying the genetic changes that drove vertebrate evolution. To date, tunicate phylogenetic relationships remain poorly resolved across taxonomic levels. Several recent studies have published phylogenies of Tunicata, using slightly overlapping datasets (Kocot et al. 2018; Delsuc et al. 2018; Alie et al. 2018).

### 1.2 _Rationale_  

We need to place Corella inflata in the tree of Tunicata to understand the evolution of genomic features. We will also use this opportunity to combine all of the publicly available data from
the studies mentioned in the Background information into a more comprehensive phylogeny of Chordata.

### 1.3 _Objectives_  

We will assemble a draft genome and transcriptome for the tunicate Corella inflata. We will infer phylogenetic relationships among tunicate species using previously published data and data generated in this study for C. inflata. We will investigate Hox gene family evolution in Corella inflata by estimating a gene tree for C. inflata, Ciona robusta, and Homo sapiens. 

## 2 STUDY DESIGN AND ENDPOINTS  

#### 2.1 Assemble the Corella inflata genome from PacBio and Illumina reads  

2.1.1 Error correct PacBio reads

```
SeqChunker -s 30G -o pb-%03d.fq pacBio/pacbioreads.fastq

proovread -t 40 --coverage=240 -l pb-001.fq -s G34.fastq -s G36.fastq --pre pb-001
```

2.1.2 Remove adaptor sequences from Illumina DNA reads using Trimmomatic executed through Galaxy with a sliding window of 4 and a phred score cutoff of 27. All other settings were default. We used jellyfish and Quake to correct substitution errors in the Illumina DNA reads. 

```
jellyfish count -m 18 -t 6 -s 100000000 --both-strands R1.fq R2.fq
jellyfish merge -o mer_counts_merged.jf mer_counts_0 mer_counts_1 mer_counts_2  mer_counts_3 mer_counts_4 mer_counts_5 mer_counts_6 mer_counts_7 mer_counts_8 mer_counts_9 mer_counts_10 mer_counts_11 mer_counts_12 mer_counts_13
```

```
python quake.py -f input_files.txt —no_count -k 17 -p 6 -q 33
```

2.1.3 Assemble genome with trimmed and corrected Illumina reads

```
run_meraculous.sh -c <configFile>
```

configure file:
```
is_diploid 1  
mer_size 65  
num_prefixed_blocks 16  
local_num_procs 8  
genome_size 0.15  
no_read_validation 1  
local_max_memory 2100  
lib_seq /path/to/fastq.R1 /path/to/fastq.R2 LIB1 550 20 100 0 0 1 1 1 0 0
```

2.1.4 Create mate pairs from PacBio data. Matemaker is available here https://github.com/josephryan/matemaker.

```
matemaker --assembly ../../07-PROOVREAD/pb-001/pb-001.trimmed.fa --out=Cinf2kmates_v4 --insertsize=2000 --matesperkb=50
matemaker --assembly ../../07-PROOVREAD/pb-001/pb-001.trimmed.fa --out=Cinf2kmates_v4 --insertsize=2000 --matesperkb=50
matemaker --assembly ../../07-PROOVREAD/pb-001/pb-001.trimmed.fa --out=Cinf5kmates_v4 --insertsize=5000 --matesperkb=50
matemaker --assembly ../../07-PROOVREAD/pb-001/pb-001.trimmed.fa --out=Cinf10kmates_v4 --insertsize=10000 --matesperkb=10
matemaker --assembly ../../07-PROOVREAD/pb-001/pb-001.trimmed.fa --out=Cinf15kmates_v4 --insertsize=15000 --matesperkb=10

echo 'Lib1 bowtie Cinf2kmates_v4.A.fq Cinf2kmates_v4.B.fq 2000 0.25 RF' > libraries.txt
echo 'Lib1 bowtie Cinf5kmates_v4.A.fq Cinf5kmates_v4.B.fq 5000 0.25 RF' >> libraries.txt
echo 'Lib1 bowtie Cinf10kmates_v4.A.fq Cinf10kmates_v4.B.fq 10000 0.25 RF' >> libraries.txt
echo 'Lib1 bowtie Cinf15kmates_v4.A.fq Cinf15kmates_v4.B.fq 15000 0.25 RF' >> libraries.txt
```

2.1.5 Scaffold the meraculous Illumina assembly with PacBio mate pairs

```
perl /usr/local/SSPACE-STANDARD-3.0_linux-x86_64/SSPACE_Standard_v3.0.pl -l libraries.txt -s  meraculous_v3.fasta -T 44 -x 1 -k 3 -b Cinf_prmmss_v4 > sspacev4.out 2> sspacev4.err

rm -rf Cinf_prmmss_v4/intermediate_results/ Cinf_prmmss_v4/pairinfo Cinf_prmmss_v4/reads;
pigz -9 *.fq Cinf_prmmss_v4/*.fasta
```

#### 2.2 Assemble the Corella transcriptome from Illumina RNA-seq reads

2.2.1 Concatenate the R1 reads and the R2 reads

```
cat FGC1291_s_4_1_CTTGTA.fastq FGC1291_s_4_1_GCCAAT.fastq > core_R1.fq

cat FGC1291_s_4_2_CTTGTA.fastq FGC1291_s_4_2_GCCAAT.fastq > core_R2.fq
```

2.2.2 Trim adaptors from RNAseq reads using bl-filter-illumina v0.4.0 from biolite-tools

```
bl-filter-illumina -a -i core_R1.fq -i core_R2.fq -o core_R1_bfl.fq -o core_R2_bfl.fq -u core_unp_bfl.fq > blf.out 2> blf.err &
```

2.2.3 Concatenate trimmed unpaired reads to trimmed R1 reads

```
cat core_R1_bfl.fq core_unp_bfl.fq > core_R1unp_bfl.fq
```

2.2.4 Assemble the transcriptome in Trinity v2.4.0

```
Trinity --seqType fq --max_memory 750G --CPU 12 --left core_R1unp_bfl.fq --right core_R2.fq --full_cleanup --normalize_reads --normalize_max_read_cov 30 > trin.out 2> trin.err
```

2.2.5 Identify isoforms with the most aligned reads and create a new assembly with these isoforms

```
align_and_estimate_abundance.pl --transcripts trinity_out_dir.Trinity.fasta --seqType fq --left core_R1_bfl.fq --rightcore_R2_bfl.fq --output_dir aea --est_method RSEM --aln_method bowtie2 --thread_count 100 --prep_reference > aea.out 2> aea.err &
```

```
rsemgetbestseqs.py ./aea/RSEM.isoforms.results trinity_out_dir.Trinity.fasta > rgbs.out 2> rgbs.err
```

2.2.6 Collapse transcripts with 97% sequence similarity or better in CDHIT v4.7

```
cd-hit-est -i /bwdata1/lincoln/02-RESTART/Corella_inflata_transcriptopme_v1.fa -o Corella_inflata_transcriptome_v1_cdhit_97.fa -c 0.97 -T 80 -M 3100 > cdhit_97.out 2> cdhit_97.err &
```

2.2.7 We used the program [Alien Index](https://github.com/josephryan/alien_index) and a database of representative tunicate and crustacean sequences to filter sequences from parasitic copepods that contaminate the Corella inflata transcriptome

```
blastp -query [infile.pep.fa] -db ai.fa -outfmt 6 -max_target_seqs 1000 -seg yes -evalue 0.001 -out [file.out] > file.std 2> file.err
```

```
./alien_index --blast=[file_ai.out] --alien_pattern=ALIEN [out.alien_index] > ai.out 2> ai.err 
```

```
remove_aliens.pl [out.alien_index] [original_transcriptome.fa] > [filtered_transcriptome.fa] > ra.out 2> ra.err
```

2.2.8 Rename the definition lines

```
replace_deflines.pl --trinity2 --pad=6 --fasta=RSEM.isoforms.results.best.fa --prefix=Core_infl_trans > Corella_inflata_transcriptome_v1.fa
```

#### 2.3 Translate the nucleotide transcriptomes into amino acid sequences with TransDecoder v3.0.1. We set the –m flag to 50 and used the results from BLAST and HMMscan searches to inform the final TransDecoder prediction step. Included in this analysis are the Corella inflata transcriptome we assembled and published data from Delsuc et al 2018 doi: 10.1186/s12915-018-0499-2 and Alie et al. 2018 doi: 10.1093/molbev/msy068  Our downstream analyses that include Kocot et al (2018) transcriptomic data will use the translations included in their Dryad repository. 

```
TransDecoder.LongOrfs -t [infile.fa] -m 50 > td.out 2> td.err
```

```
blastp -query longest_orfs.pep -db swissprot -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 4 > blastp.out 2> blastp.err 
```

```
hmmscan --cpu 1 --domtblout outfile.domtblout Pfam-A.hmm longest_orfs.pep > hs.out 2> hs.err
```

```
TransDecoder.Predict -t [infile.fa] --retain_pfam_hits out.domtblout --retain_blastp_hits out.blastp.out > tdp.out 2> tdp.err
```

#### 2.4 Estimate Corella inflata gene models in Augustus

```
blat -noHead -stepSize=5 -minIdentity=93 Cinf_prmmss_v4.final.scaffolds.fasta Corella_inflata_trans_v1_cdhit97_AI.fa ali.psl
cat ali.psl | sort -k 10,10 > ali_sort.psl
cat ali_sort.psl | /usr/local/augustus-3.2.3/scripts/filterPSL.pl --uniq > ali_sort_filt.psl
sort -k 14,14 ali_sort_filt.psl > ali_sort_filt_sort.psl
/usr/local/augustus-3.3/bin/aln2wig -f ali_sort_filt_sort.psl > cov.wig
cat cov.wig | /usr/local/augustus-3.2.2/scripts/wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep --UCSC=unstranded.track --radius=4.5 --pri=4 --strand="." > hints.ep.gff
/usr/local/augustus-3.2.2/scripts/blat2hints.pl --intronsonly --in=ali_sort_filt_sort.psl --out=hints.introns.gff
cat hints.introns.gff hints.ep.gff > hints.gff
cp /usr/local/augustus-3.2.2/config/extrinsic/extrinsic.M.RM.E.W.cfg .
augustus --species=caenorhabditis --UTR=on --extrinsicCfgFile=extrinsic.M.RM.E.W.cfg --alternatives-from-evidence=true --hintsfile=hints.gff --allow_hinted_splicesites=atac Cinf_prmmss_v4.final.scaffolds.fasta > Cinf_prmmss_v4.final.scaffolds.augustus 2> aug.err &
perl /usr/local/augustus-3.2.3/scripts/getAnnoFasta.pl Cinf_prmmss_v4.final.scaffolds.augustus
perl /usr/local/augustus-3.3/scripts/getAnnoFasta.pl Cinf_prmmss_v4.final.scaffolds.augustus --seqfile=Cinf_prmmss_v4.final.scaffolds.fasta
```

#### 2.5 We used OrthoFinder v2.2.3 to identify orthologs among taxa for two data sets. The first data set, DS_gene_loss, included C. inflata gene models and transcriptome, C. robusta gene models transcriptome, and Homo sapiens gene models to aid in functional annotation. We used this data set to determine patterns of gene loss in C. inflata (see below). The second data set, DS_phylogeny, included transcriptomes from previously published tunicate and outgroup species and was used to infer phylogenetic relationships among tunicates (see below). We used a custom script (ortho_diamond.pl) available in this repository to run the diamond searches. 

```
orthofinder -f [dir_w_protein_fast_files] -op > of.out 2> of.err

```

```
diamond blastp -d BlastDB -q species.fa -o species_out.txt -e 0.001 -p 10 -f 6 > species.stdout 2> species.err &
```

```
orthofinder -b [dir_w_blast_results] -a 16 -M msa -os > ofb.out 2> ofb.err
```

#### 2.6 Generate single copy orthogroups from DS_phylogeny using the script ```filter_ogs_write_scripts.pl``` (available in https://github.com/josephryan/RyanLabPhylogenomicTools). The script allows the user to define the minimum number of taxa and the maximum number of duplicates per taxon allowed per orthogroup (in this study 80% (40 taxa) and 8 duplicates) and after this filtering, automates the following steps:   

2.6.1 sequences within each orthogroup are aligned using Mafft v7.309 

```mafft-linsi --localpair --maxiterate 1000 --thread 20 [infile] >mafft.out 2> mafft.err```

2.6.2 alignments are refined using Gblockswrapper v0.03 (https://goo.gl/fDjan6)

```Gblockswrapper [infile.mafft] > outfile.mafft-gb > gbw.out 2> gbw.err```

2.6.3 Gblockswrapper sometimes leaves blank sequences that cause downstream issues; the ```remove_empty_seqs``` script, available in https://github.com/josephryan/RyanLabPhylogenomicTools, removes empty sequences and spaces from sequence lines. 

```remove_empty_seqs [outfile.mafft-gb] > res.out 2> res.err```

2.6.4 maximum-likelihood orthogroup gene trees are estimated in IQTree v1.5.5 

```iqtree-omp -s [infile.mafft-gb] -nt AUTO -bb 1000 -m LG -pre [output prefix] > iq.out 2> iq.err```

2.6.5 paralogs are pruned in PhyloTreePruner v1.0

```java PhyloTreePruner [infile.tree] 28 [infile.align] 0.5 u > ptp.out 2> ptp.err```

#### 2.7 Concatenate the single-copy loci filtered from step 2.7 to create a matrix and partition file for use in downstream phylogenomic analyses using ```fasta2phylomatrix``` (available in https://github.com/josephryan/RyanLabPhylogenomicTools). Definition lines in each fasta file were edited (```perl -pi.orig -e 's/\|.*$//;' *.fa```) prior to running ```fasta2phylomatrix```  

#### 2.8 Estimate species phylogeny from the concatenated matrix (tunicate_210.fa) using maximum likelihood and Bayesian inference

2.8.1 ML species tree inference 

```
iqtree-omp -s tunicate_210.fa -pre tunicate_210 -spp tunicate_210.nex -nt AUTO -m TEST -bb 1000 > iq.stdout 2> iq.err
``` 

2.8.2 Bayesian species tree inference. We generated 10 random starting trees in IQTREE and then in PhyloBayes we launched 2 chains for every random starting tree for 20 chains total, with 8 processors per chain.

```
iqtree-omp -nt 9 -s ../tunicate_210.fa -r 47 tunicate_210_random_1.tre > random_1.stdout 2> random_1.err &
iqtree-omp -nt 9 -s ../tunicate_210.fa -r 47 tunicate_210_random_2.tre > random_2.stdout 2> random_2.err &
iqtree-omp -nt 9 -s ../tunicate_210.fa -r 47 tunicate_210_random_3.tre > random_3.stdout 2> random_3.err &
etc...

mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_1.tre random_1_chain1 > random_1_chain1.stdout 2> random_1_chain1.err &
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_1.tre random_1_chain2 > random_1_chain2.stdout 2> random_1_chain2.err &
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_2.tre random_2_chain1 > random_2_chain1.stdout 2> random_2_chain1.err &
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_2.tre random_2_chain2 > random_2_chain2.stdout 2> random_2_chain2.err &
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_3.tre random_3_chain1 > random_3_chain1.stdout 2> random_3_chain1.err &
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_3.tre random_3_chain2 > random_3_chain2.stdout 2> random_3_chain2.err &
etc...
```

After 6 weeks (9/11/18-10/24/18) we tested for convergence (i.e., when the discrepancy observed across all bipartitions was < 0.1). Two of the ten pairs had converged and at that point all chains were terminated. We estimated a consensus tree from both chain 1 trees, sampling every 10th tree after a 100 tree burnin.

```
tracecomp -x 100 2 random_1_chain1 random_1_chain2
readpb -x 100 10 random_1_chain
```

#### 2.9 The ML and Bayesian phylogenies had different topologies. To compare these topologies, in IQ-TREE we estimated likelihood scores for phylogenies where the data were constrained to the Bayesian topology and the ML topology. 

```
iqtree-omp -s tunicate_210.fa -m TEST -g ml_topo.tre -pre ml_topo_constraint -spp tunicate_210.nex -nt AUTO > ml_topo_test.out 2> ml_topo_test.err

iqtree-omp -s tunicate_210.fa -m TEST -g bayes_topo.tre -pre bi_topo_constraint -spp tunicate_210.nex -nt AUTO > bi_topo_test.out 2> bi_topo_test.err
```
  
#### 2.10 Use custom scripts to search for orthogroups that were present in some taxa and missing from others in the DS_gain_loss data set. We used a Markov Chain Monte Carlo approach to test if any functional categories were lost more often than expect by chance. The scripts and files required for these analyses are available in this repository. #### ADD TO REPO ####

#### 2.11 We tested our 210-gene data set and the 798-gene original full data set and 50-gene RCVF data set from Kocot et al. 2018 for compositional heterogeneity using chet v0.03 (github.com/josephryan/chet).

```
chet --fasta=tunicate_210.fa --clade1=Ascidia_sp,Phal_mamm,Core_infl,Core_will --clade2=Clav_lepa,Cyst_dell,Dist_occi
chet --fasta=tunicate_210.fa --clade1=Ascidia_sp,Phal_mamm,Core_infl,Core_will --clade2=Cion_inte,Cion_savi
```

```
chet --fasta=FcC_supermatrix_original_full_dataset.fas --clade1=DIST,CDEL --clade2=CWIL,ASCI
chet --fasta=FcC_supermatrix_original_full_dataset.fas --clade2=CWIL,ASCI --clade1=CSAV,CROB,CINT
chet --fasta=50_best_OGs_average_RCFV.fa --clade1=DIST,CDEL --clade2=CWIL,ASCI
chet --fasta=50_best_OGs_average_RCFV.fa --clade2=CWIL,ASCI --clade1=CSAV,CROB,CINT
```

#### 2.12 Identify C. inflata Hox genes using a hidden Markov model and a tree-based approach 

2.12.1 We will search the translated C. inflata transcriptome and the amino acid gene models against the hd60.hmm hidden Markov model from Zwarycz et al. 2015 using the script ```hmm2aln.pl``` which runs hmmsearch, stockholm2fasta, and custom code to remove indels and fill end gaps.

```./hmm2aln.pl --hmm=hd60.hmm --name=<out_prefix> --fasta=<fasta_file_to_search> --threads=40 --nofillcnf=nofill.hox.conf > outfile.fa 2> std_err.txt```
  
2.12.2 We will estimate a gene tree with the C. inflata loci identified in step 2.12.1 and Branchiostoma Hox genes downloaded from HomeoDB.

```iqtree-omp -s [infile] -nt AUTO -bb 1000 -m TEST -pre [output prefix] > iq.out 2> iq.err &```

2.12.3 We will prune non-Hox genes from the gene tree generated in step 2.12.2 using the script ```make_subalignment```. The output of this script is an alignment of only those genes within the clade descended from the most recent common ancestor of all genes with the specified prefix. We will use Branchiostoma as the ingroup prefix. The script is available here https://github.com/josephryan/make_subalignment. 

```
./make_subalignment --tree=<newick_treefile> --aln=<phylip_alignment> --root=<root_taxa> --pre=<prefix>
```

2.12.4 Concatenate the C. inflata Hox genes we identified in step 2.13.3 with the Branchiostoma and C. robusta Hox genes downloaded from ANIDSEED and estimate the Hox gene tree in IQTree. 

```
iqtree-omp -s [infile.fa] -pre hox -nt AUTO -m TEST -bb 1000 > iq.stdout 2> iq.err
``` 

## 3 LITERATURE REFERENCED  

Emms, D. M., & Kelly, S. (2015). OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy. Genome Biology, 16(1), 157. 

Gblockswrapper: http://bit.ly/2svaKcR

Kocot, K. M., Citarella, M. R., Moroz, L. L., & Halanych, K. M. (2013). PhyloTreePruner: a phylogenetic tree-based approach for selection of orthologous sequences for phylogenomics. Evolutionary Bioinformatics Online, 9, 429.

Lartillot, N., Rodrigue, N., Stubbs, D., & Richer, J. (2013). PhyloBayes MPI: phylogenetic reconstruction with infinite mixtures of profiles in a parallel environment. Systematic Biology, 62(4), 611-615.

Nguyen, L. T., Schmidt, H. A., von Haeseler, A., & Minh, B. Q. (2014). IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Molecular Biology and Evolution, 32(1), 268-274.

TransDecoder: https://transdecoder.github.io/

Yamada, K. D., Tomii, K., & Katoh, K. (2016). Application of the MAFFT sequence alignment program to large data—reexamination of the usefulness of chained guide trees. Bioinformatics, 32(21), 3246-3251.

## APPENDIX

Version : Date : Significant Revisions  
1.1 
1.2  
1.3  
1.4 
