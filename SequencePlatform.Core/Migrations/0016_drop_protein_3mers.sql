-- 0016_drop_protein_3mers.sql (2026-11-05): the amino-acid 3-mer counts of 0014 are dropped. CVTree's 4-mer vector needs 3-mers only
-- for its Markov expectation, and the 4-mer counts summed over the last letter give them, short only of each gene's last 3-mer: the
-- hashed 8,192 vector kept within-family Spearman 0.8511 vs 0.8505 and 5-nearest purity 0.6003 vs 0.6013 on 600 GTDB-placed bacteria
-- (results/CvTree_Derived3mers.csv). About 0.30 GB on the local catalogue. Getting them back needs the gene pass's re-download.

delete from protein_kmer where k = 3;

alter table protein_kmer drop constraint protein_kmer_k_check;
alter table protein_kmer add constraint protein_kmer_k_check check (k = 4);
