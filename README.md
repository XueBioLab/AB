# Genus-level phylogenomics reveals the stepwise evolution driving Acinetobacter baumannii ecological success

The Acinetobacter genus encompasses numerous species occupying diverse habitats, with several emerging as human pathogens posing significant threats. However, the mechanisms underlying the emergence remain elusive. Here, through a genus-wide population genomic analysis, we reveal significant phylogenetic clustering of hosts coupled with distinct habitat-switching profiles, delineating two predominant lineages: an animal-associated lineage (C1) and a human-associated lineage (C2). C2 includes major clinically important pathogens, notably the A. calcoaceticus-baumannii (ACB) complex. Genome reduction driven by relaxation of selection in C1 represents a primary driver of this phylogenetic diversification, resulting in significantly smaller genomes compared to C2. This divergent genome evolution remodels metabolic profiles, enabling C2 to utilize a markedly broader spectrum of substrates than C1. An enrichment of the catalase gene katG under strong positive selection is further detected in the ACB complex, conferring enhanced resistance to hydrogen peroxide. This evolutionary advantage facilitates the prevalence of the complex in clinical settings. The epidemic success of A. baumannii is further bolstered by its exceptionally high pangenome openness compared to other taxa of the complex, driving the accumulation of diverse fitness factors. In particular, the overrepresentation of an ant(3")-II gene not only confers aminoglycoside resistance, but also promotes gut colonization and systematic dissemination. Human and non-human A. baumannii strains share a highly similar genetic background, underscoring the species' robust ecological versatility. Collectively, our findings pinpoint fundamental drivers of phylogenetic diversification in Acinetobacter and elucidate the genomic basis propelling the emergence and success of human pathogenic species.
This repository contains a collection of code and scripts used in the paper.


## 1. Software used in this workflow

- [Perl](https://www.perl.org/)
- [Python3](https://www.python.org/)
- [Trimmomatic](https://github.com/timflutre/trimmomatic)
- [SPAdes](https://github.com/ablab/spades)
- [Prokka](https://github.com/tseemann/prokka)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [RAxML](https://evomics.org/learning/phylogenetics/raxml/)
- [Pandas](https://pandas.pydata.org/)
- [Numpy](https://numpy.org/)
- [CheckM](https://ecogenomics.github.io/CheckM/)
- [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk)
- [Abricate](https://github.com/tseemann/abricate)
- [Snippy](https://github.com/tseemann/snippy)
- [Gubbins](https://github.com/nickjcroucher/gubbins)
- [Phytools](https://cran.r-project.org/web/packages/phytools/index.html)
- [Seqkit](https://bioinf.shenwei.me/seqkit/)

>Take the GCA_018883565.1 genome as an example.

## 2. Dataset

## 3. read trimming
### Trimmomatic
```bash
java -jar trimmomatic-0.36.jar PE -threads 5 GCA_018883565.1_raw_1.fq.gz -2 GCA_018883565.1_raw_2.fq.gz ...
```

## 4. Assembly
### SPAdes
```bash
spades.py -1 GCA_018883565.1_clean_1.fq.gz -2 GCA_018883565.1_clean_2.fq.gz --isolate --cov-cutoff auto -o GCA_018883565.1.fasta
```

## 5. Preprocessing
### checkM
```bash
checkm2 predict --threads 1 --input fasta_dir -x fasta --output-directory checkM.out --database_path your_checkM_database/uniref100.KO.1.dmnd
```

## 6. Genome annotation
### Prokka
```bash
prokka ../04.assembly/GCA_018883565.1.fasta --prefix GCA_018883565.1 --outdir prokka.out --compliant
```

## 7. Taxonomy assignment
### GTDB
```bash
gtdbtk classify_wf --genome_dir fasta_dir/ --out_dir fasta_dir.GTDB.out --extension fasta
# fasta_dir, the input directory containing a set of genomic assembly sequences.
# fasta_dir.GTDB.out, output directory
```

## 8. ancestral state reconstruction
### phytools
```R
library(phytools)
library(ape)
tree <- read.tree("gtdbtk.bac120.user_msa_selected.rmdup_root.nwk")
tree <- root(tree, "GCA_002174125.1_ASM217412v1_genomic", resolve.root = TRUE)
x <- read.csv("1251genomes_metadata.CSV", row.names = 1)
x <- as.matrix(x)[,1]
mtree<-make.simmap(tree,x,model="ER")
pd<-summary(mtree,plot=FALSE)

cols <- read.csv("cols.csv",row.names =1)
cols <- as.matrix(cols)[,1]
plot(mtree,cols,type="fan",fsize=0.8,ftype="off")
```

## 9. Panaroo
```bash
panaroo -i GCA_018883565.1.gff -o results --clean-mode strict --alignment pan --aligner mafft --merge_paralogs -f 0.5 --len_dif_percent 0.7 -t 30
```

## 10. Identification of ARGs
### Abricate
```bash
mkdir ARG_dir

for f in `ls fasta_dir`
do 

  abricate -db resfinder --nopath --minid 80 --mincov 80 --quiet fasta_dir/${f} > ARG_dir/${f%%.fasta}.tab

done

abricate --nopath --summary ARG_dir/*tab > ARG.tab
# fasta_dir, the input directory containing a set of genomic assembly sequences.
```

## 11. HyPhy
```bash
for i in *.fna
do 
  
  hyphy ../FitMG94.bf --alignment $i --tree ../02.core_gene_subtree/${i%%.pal2nal.fna}.faa.nwk --lrt Yes --type global --output $i.hyphy.out

done
```

## 12. Roary
```bash
N=19
B=100

mkdir -p subsamples

for i in $(seq 1 $B)
do
  echo "Iteration $i"

  mkdir -p subsamples/run_$i

  ls AB_gff/*.gff | shuf | head -n $N > subsamples/run_$i/list.txt

  roary -p 8 -f subsamples/run_$i $(cat subsamples/run_$i/list.txt)

done
```

## 13. micropan_alpha
```R
library(micropan)

B <- 100
alphas <- numeric(B)

for (i in 1:B) {

  file <- paste0("./subsamples/run_", i, "/gene_presence_absence.csv")
  
  if (file.exists(file)) {
    g <- tryCatch({
      roary <- read.csv(file, header=TRUE, check.names=FALSE)

      start_col <- which(colnames(roary) == "Genome.Fragment")[1] + 1
      if (is.na(start_col)) start_col <- 15

      mat <- roary[, start_col:ncol(roary)]
      mat[mat != ""] <- 1
      mat[mat == ""] <- 0

      mat <- as.matrix(apply(mat, 2, as.numeric))
      pan <- t(mat)

      h <- heaps(pan, n.perm=100)

      if (is.list(h)) h$alpha
      else as.numeric(h["alpha"])
    }, error=function(e) NA)

    if (is.finite(g) && g > 0 && g < 2) {
      alphas[i] <- g
    } else {
      alphas[i] <- NA
    }

    cat("Run", i, "alpha =", g, "\n")
  } else {
    alphas[i] <- NA
  }
}

alphas <- alphas[!is.na(alphas)]

mean_alpha <- mean(alphas)
sd_alpha <- sd(alphas)
ci <- quantile(alphas, c(0.025, 0.975))

cat("\n====================\n")
cat("Mean alpha:", mean_alpha, "\n")
cat("SD:", sd_alpha, "\n")
cat("95% CI:", ci[1], "-", ci[2], "\n")
cat("====================\n")

write.csv(alphas, "alpha_values_roary_resampling.csv", row.names=FALSE)

pdf("alpha_distribution_roary.pdf")
hist(alphas, breaks=20, main="alpha distribution (Roary resampling)", xlab="alpha")
abline(v=mean_alpha, lwd=2)
dev.off()
```

## 14. Use PIRATE to build the accessory gene matrix
### PIRATE
```bash
PIRATE -i gff/ -t 40

构建accessory gene matrix
PIRATE_to_Rtab.pl -i ./PIRATE.*.tsv --low 0 --high 0.95 -o ./binary_presc_absc.tsv
Allele frequency:
     -l|--low               min allele frequency to include in output 
                            [default: 0.05]
     -h|--high              max allele frequency to include in output 
                            [default: 0.95]
```
## 15. Use GraPPLE to calculate the similarity between isolates.
### GraPPLE
```python
python3 pw_similarity.py -i binary_presc_absc.tsv -o GraPPLE -r "isolates" -s "jaccard" -f 0.5 -t 2
```


