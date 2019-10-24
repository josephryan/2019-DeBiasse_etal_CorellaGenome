# Commands used in the analyses for the _Corella inflata_ genome project
### DeBiasse MB, Colgan W, Harris L, Davidson B, Ryan JF

### NEED TO ADD:
### ```remove_numbers_before_fasta2phylomatrix.pl```
### ```rename_augustus.pl```

1. Install remove_numbers_before_fasta2phylomatrix.pl and rename_augustus.pl from the scripts directory in this repo:

```
cd scripts
perl Makefile.PL
make
make install
```

2. Install JFR::Fasta, fasta2phylomatrix, and replace_deflines.pl available at: https://github.com/josephryan/JFR-PerlModules

3. Install ortho_diamond and parapruner.pl available at: https://github.com/josephryan/RyanLabPhylogenomicTools

4. Install chet available at: https://github.com/josephryan/chet

5. Install hmm2aln.pl available at: https://github.com/josephryan/hmm2aln.pl

6. Install matemaker available at: https://github.com/josephryan/matemaker
 
7. Install make_subalignment available at: https://github.com/josephryan/make_subalignment

8. Install remove_short_and_sort available at: https://github.com/josephryan/RyanLabShortReadAssembly

9. Install rsemgetbestseqs.py from Warren Francis’ BitBucket page: https://bitbucket.org/wrf/sequences/src

10.Create a directory named ```01-FILES``` on your local machine. Download the following and place them in ```01-FILES```

_Available from_ https://github.com/josephryan/2019-DeBiasse_etal_CorellaGenome
```meraculous_config_file```
```subclade.txt```
```branchiostoma.fa```
```Cion_inte_trans.pep.fa```
```Cion_inte_gm.pep.fa```

_Available from the European Nucleotide Archive at accession numbers XXXXX_ 
```Core_infl_Illumina_DNA_PE_R1.fastq```
```Core_infl_Illumina_DNA_PE_R2.fastq``` 
```FGC1291_s_4_1_CTTGTA.fastq```
```FGC1291_s_4_1_GCCAAT.fastq```
```FGC1291_s_4_2_CTTGTA.fastq```
```FGC1291_s_4_2_GCCAAT.fastq``` 
```core_infl_pacbioreads.fastq``` 

### Quality control for Illumina DNA reads
```
java -jar trimmomatic-0.36.jar PE -threads 12 -phred33 Core_infl_Illumina_DNA_PE_R1.fastq Core_infl_Illumina_DNA_PE_R2.fastq Core_infl_DNA_R1_trimmed.fq Core_infl_DNA_R1_trimmed_unp.fq Core_infl_DNA_R2_trimmed.fq Core_infl_DNA_R2_trimmed_unp.fq ILLUMINACLIP:/usr/local/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:12:1:true MINLEN:36
```
### Correct Illumina DNA read sequencing errors
```
perl allpathslg-44837/src/ErrorCorrectReads.pl PAIRED_READS_A_IN=Core_infl_DNA_R1_trimmed.fq PAIRED_READS_B_IN=Core_infl_DNA_R1_trimmed.fq PAIRED_SEP=350 THREADS=46 PHRED_ENCODING=33 READS_OUT=Core_infl_DNA_trim_ecr
```
```
perl /allpathslg-44837/src/ErrorCorrectReads.pl UNPAIRED_READS_IN=Core_infl_DNA_R1_trimmed_unp.fq THREADS=46 PHRED_ENCODING=33 READS_OUT=Core_infl_DNA_R1_trimmed_ecr_unp.fq
```
```
perl /allpathslg-44837/src/ErrorCorrectReads.pl UNPAIRED_READS_IN=Core_infl_DNA_R2_trimmed_unp.fq THREADS=46 PHRED_ENCODING=33 READS_OUT=Core_infl_DNA_R2_trimmed_ecr_unp.fq
```
### Genome assembly with Illumina DNA reads
```
run_meraculous.sh -c meraculous_config_file
```
### PacBio read correction
```
SeqChunker -s 30G -o pb-%03d.fq 01-FILES/Core_infl_pacbioreads.fastq
```
```
proovread -t 40 --coverage=240 -l pb-001.fq -s Core_infl_DNA_trim_ecr_R1.fq -s Core_infl_DNA_trim_ecr_R2.fq --pre pb-001
```
### Generate mated seqs
```
matemaker --assembly pb-001.trimmed.fa --out=Cinf2kmates_v4 --insertsize=2000 --matesperkb=50
matemaker --assembly pb-001.trimmed.fa --out=Cinf2kmates_v4 --insertsize=2000 --matesperkb=50
matemaker --assembly pb-001.trimmed.fa --out=Cinf5kmates_v4 --insertsize=5000 --matesperkb=50
matemaker --assembly pb-001.trimmed.fa --out=Cinf10kmates_v4 --insertsize=10000 --matesperkb=10
matemaker --assembly pb-001.trimmed.fa --out=Cinf15kmates_v4 --insertsize=15000 --matesperkb=10
```
### Scaffold Illumina genome with PacBio reads
```
perl SSPACE-STANDARD-3.0_linux-x86_64/SSPACE_Standard_v3.0.pl -l libraries.txt -s Core_infl_meraculous.fasta -T 44 -x 1 -k 3 -b Core_infl
```
### Sort genome, remove scaffolds less than 200bp, rename deflines
```
remove_short_and_sort Core_infl.final.scaffolds.fasta 200 > Core_infl_genome_sorted_no_short.fa
```
```
replace_deflines.pl --pad=5 --fasta=Core_infl_genome_sorted_no_short.fa --prefix=Cinf4 > Core_infl_genome_v2.fa
```
​
#### Quality control for Illumina RNAseq reads
```
cat 01-FILES/FGC1291_s_4_1_CTTGTA.fastq 01-FILES/FGC1291_s_4_1_GCCAAT.fastq > core_R1.fq
```
```
cat 01-FILES/FGC1291_s_4_2_CTTGTA.fastq 01-FILES/FGC1291_s_4_2_GCCAAT.fastq > core_R2.fq
```
```
bl-filter-illumina -a -i core_R1.fq -i core_R2.fq -o core_R1_bfl.fq -o core_R2_bfl.fq -u core_unp_bfl.fq
```
```
cat core_R1_bfl.fq core_unp_bfl.fq > core_R1unp_bfl.fq
```

#### Transcriptome assembly from RNAseq reads
```
trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --max_memory 750G --CPU 12 --left core_R1unp_bfl.fq --right core_R2_bfl.fq --full_cleanup --normalize_reads --normalize_max_read_cov 30
```
```
trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --transcripts trinity_out_dir.Trinity.fasta --seqType fq --left core_R1_bfl.fq --right core_R2_bfl.fq --output_dir aea --est_method RSEM --aln_method bowtie2 --thread_count 100 --prep_reference
```
```
rsemgetbestseqs.py RSEM.isoforms.results trinity_out_dir.Trinity.fasta
```
```
mv trinity_out_dir.Trinity.fasta.RSEM.transcripts.fa Corella_inflata_transcriptopme_v1.fa
```

#### Collapse transcriptome assembly
```
cd-hit-est -i Corella_inflata_transcriptopme_v1.fa -o Corella_inflata_transcriptome_v1_cdhit_97.fa -c 0.97 -T 80 -M 3100
```

### Transcriptome translation
```
TransDecoder-3.0.1/TransDecoder.LongOrfs -t Corella_inflata_transcriptome_v1_cdhit_97.fa -m 50
```
```
blastp -query Corella_inflata_transcriptome_v1_cdhit_97.fa.transdecoder_dir/longest_orfs.pep -db swissprot -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 16
```
```
hmmscan --cpu 30 --domtblout corella.domtblout Pfam-A.hmm Corella_inflata_transcriptome_v1_cdhit_97.fa.transdecoder_dir/longest_orfs.pep
```
```
TransDecoder-3.0.0/TransDecoder.Predict -t Corella_inflata_transcriptome_v1_cdhit_97.fa --retain_pfam_hits outfile.domtblout --retain_blastp_hits outfile.blastp.out
```
```
mv Corella_inflata_transcriptome_v1_cdhit_97.fa.transdecoder.pep Corella_inflata_transcriptome_v1_cdhit_97.pep.fa
```

### Ortholog identification for species phylogeny
```
orthofinder -f 01-FILES -op
```
```
ortho_diamond
```
```
orthofinder -b 01-FILES/Results/WorkingDirectory -a 46 -M msa -os
```
```
filter_ogs_write_scripts.pl --ogdir=Results/WorkingDirectory/Orthologues/Sequences --name=40sp8dups --min_sp=40 --max_sp_occur=8 --num_scripts=46 --threads=1
```
```
perl parapruner.pl 
```

### Generate concatenated matrix
```
remove_numbers_before_fast2phylomatrix.pl 
```
```
fasta2phylomatrix --dir=01-NO_NUMBERS --raxpartition=tunicate_210.rax --nexpartition=tunicate_210.nex > tunicate_210.fa
```
### Maximum likelihood species tree estimation
```
iqtree-omp -s tunicate_210.fa -pre tunicate_210 -spp tunicate_210.nex -nt AUTO -m TEST -bb 1000
```
### Bayesian species tree estimation
```
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_nj_10.tre nj_10_chain1
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_1.tre random_1_chain1
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_2.tre random_2_chain1
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_3.tre random_3_chain1
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_4.tre random_4_chain1
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_5.tre random_5_chain1
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_6.tre random_6_chain1
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_7.tre random_7_chain1
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_8.tre random_8_chain1
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_9.tre random_9_chain1
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_nj_10.tre nj_10_chain2
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_1.tre random_1_chain2
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_2.tre random_2_chain2
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_3.tre random_3_chain2
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_4.tre random_4_chain2
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_5.tre random_5_chain2
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_6.tre random_6_chain2
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_7.tre random_7_chain2
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_8.tre random_8_chain2
mpirun -n 8 pb_mpi -d tunicate_210.phy -cat -gtr -t tunicate_210_random_9.tre random_9_chain2
```
```
readpb -x 100 10 random_6_chain1.err
```
### ML vs Bayesian topology test
```
iqtree-omp -s tunicate_210.fa -m TEST -g ml_topo.tre -pre ml_topo_constraint -spp tunicate_210.nex -nt AUTO
```
```
iqtree-omp -s tunicate_210.fa -m TEST -g bayes_topo.tre -pre bi_topo_constraint -spp tunicate_210.nex -nt AUTO
```

### Tests for compositional heterogeneity
```
chet --fasta=tunicate_210.fa --clade1=Clav_lepa,Cyst_dell,Dist_occi --clade2=Ascidia_sp,Phal_mamm,Core_infl,Core_will
```
```
chet --fasta=tunicate_210.fa --clade1=Ascidia_sp,Phal_mamm,Core_infl,Core_will --clade2=Cion_inte,Cion_savi
```
```
chet --fasta=FcC_supermatrix_original_full_dataset.fas --clade1=DIST,CDEL --clade2=CWIL,ASCI
```
```
chet --fasta=FcC_supermatrix_original_full_dataset.fas --clade2=CWIL,ASCI --clade1=CSAV,CROB,CINT
```
```
chet --fasta=50_best_OGs_average_RCFV.fa --clade1=DIST,CDEL --clade2=CWIL,ASCI
```
```
chet --fasta=50_best_OGs_average_RCFV.fa --clade2=CWIL,ASCI --clade1=CSAV,CROB,CINT
```
### _10. BaCoCa analyses_
```
perl BaCoCa.v1.105.r.pl -i 50_best_OGs_average_RCFV.fa -c subclade.txt -p Kocot_50.partition
```
```
perl BaCoCa.v1.105.r.pl -i FcC_supermatrix_original_full_dataset.fa -c subclade.txt -p Kocot_orig.partition
```
### Predict gene models
```
blat -noHead -stepSize=5 -minIdentity=93 Core_infl_genome_v2.fa Corella_inflata_trans_v1_cdhit97.fa ali.psl
```
```
cat ali.psl | sort -k 10,10 > ali_sort.psl
```
```
cat ali_sort.psl | /augustus-3.2.3/scripts/filterPSL.pl --uniq > ali_sort_filt.psl
```
```
sort -k 14,14 ali_sort_filt.psl > ali_sort_filt_sort.psl
```
```
/augustus-3.3/bin/aln2wig -f ali_sort_filt_sort.psl > cov.wig
```
```
cat cov.wig | /augustus-3.3/scripts/wig2hints.pl --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep --UCSC=unstranded.track --radius=4.5 --pri=4 --strand="." > hints.ep.gff
```
```
/augustus-3.3/scripts/blat2hints.pl --intronsonly --in=ali_sort_filt_sort.psl --out=hints.introns.gff
```
```
cat hints.introns.gff hints.ep.gff > hints.gff
```
```
cp /augustus-3.3/config/extrinsic/extrinsic.M.RM.E.W.cfg .
```
```
augustus --species=ciona --extrinsicCfgFile=extrinsic.M.RM.E.W.cfg --alternatives-from-evidence=true --hintsfile=hints.gff --allow_hinted_splicesites=atac Core_infl_genome_v2.fa > Core_infl_genome_augustus.out
```
```
perl /augustus-3.3/scripts/getAnnoFasta.pl Core_infl_genome_augustus_2.out --seqfile=Core_infl_genome_v2.fa
```
```
rename_augustus.pl --gff=Core_infl_genome_augustus.out --out=Core_infl_gene_models --pre=Core_infl_genome_augustus.out
```
```
mv Core_infl_gene_models.aa Core_infl_gene_models.aa.fa
mv Core_infl_gene_models.codingseq mv Core_infl_gene_models.nuc.fa
```
### Identify Hox genes
```
cat Corella_inflata_transcriptome_v1_cdhit_97.pep.fa Core_infl_gene_models.aa.fa Cion_inte_trans.pep.fa Cion_inte_gm.pep.fa > cinf_crob_trans_gm.fa
```
```
hmm2aln.pl --hmm=hd60.hmm --name=CCTG --fasta=cinf_crob_trans_gm.fa --threads=22 --nofillcnf=nofill.hox.conf
```
```
mv cctg_hmm.out CCTG.fa
cat branchiostoma.fa CCTG.fa > CCTGB.fa
```
```
iqtree-omp -s CCTGB.fa -nt AUTO -bb 1000 -m TEST -pre CCTGB
```
```
pip3 install --user dendropy
make_subalignment --tree=CCTGB.treefile --aln=CCTGB.fa --root=Core_infl.116187 --pre=Branchiostoma > CCTGB_subaligned.fa
```
```
iqtree-omp -s CCTGB_subaligned.fa -nt AUTO -bb 1000 -m TEST -pre CCTGB_subaligned
```
### Hox gene tree AU test
```
iqtree -s CCTGB_subaligned.fa -nt AUTO -m LG+G4 -g hox10_bflor_10-12.constr -pre hox10_bflor_10-12.constr > hox10_bflor_10-12.constr.out 2> hox10_bflor_10-12.constr.err &
iqtree -s CCTGB_subaligned.fa -nt AUTO -m LG+G4 -g hox10.constr -pre hox10.constr > hox10.constr.out 2> hox10.constr.err &
iqtree -s CCTGB_subaligned.fa -nt AUTO -m LG+G4 -g hox12.constr -pre hox12.constr > hox12.constr.out 2> hox12.constr.err &
iqtree -s CCTGB_subaligned.fa -nt AUTO -m LG+G4 -g hox13.constr -pre hox13.constr > hox13.constr.out 2> hox13.constr.err &
iqtree -s CCTGB_subaligned.fa -nt AUTO -m LG+G4 -g hox4.constr -pre hox4.constr > hox4.constr.out 2> hox4.constr.err &
iqtree -s CCTGB_subaligned.fa -nt AUTO -m LG+G4 -g hox5.constr -pre hox5.constr > hox5.constr.out 2> hox5.constr.err &
iqtree -s CCTGB_subaligned.fa -nt AUTO -m LG+G4 -g hox6.constr -pre hox6.constr > hox6.constr.out 2> hox6.constr.err &
iqtree -s CCTGB_subaligned.fa -nt AUTO -m LG+G4 -g tunicate_bflor_posterior.constr -pre tunicate_bflor_posterior.constr > tunicate_bflor_posterior.constr.out 2> tunicate_bflor_posterior.constr.err &
```
```
cat hox.unconstr.treefile hox4.constr.treefile hox5.constr.treefile hox6.constr.treefile hox10.constr.treefile hox12.constr.treefile hox13.constr.treefile hox10_bflor_10-12.constr.treefile tunicate_bflor_posterior.constr.treefile > all_trees
```
```
iqtree -s CCTGB_subaligned.fa -m LG+G4 -z all_trees -pre hox_au_test -nt AUTO -n 0 -zb 1000 -au
```
