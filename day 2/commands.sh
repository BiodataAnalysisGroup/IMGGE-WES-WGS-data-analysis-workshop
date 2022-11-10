#create baits.bed
#uses .bed file and caverage from .bam to infer the targeted regions 
guess_baits.py Tumor.bam -t gene_targets.bed -o baits.bed

#devides larger bins to smaller
cnvkit.py target baits.bed  --split --annotate refFlat.txt -o baits_target.bed

#run antitarget
cnvkit.py antitarget gene_targets.bed -g access-5kb-mappable.hg19_chr5_chr12_chr17.bed -o my_antitargets.bed

#calculate coverage in the target/antitarget regions from BAM read depths
cnvkit.py coverage Tumor.bam baits_target.bed -o Tumor.targetcoverage.cnn
cnvkit.py coverage Tumor.bam my_antitargets.bed -o Tumor.antitargetcoverage.cnn

cnvkit.py coverage Normal.bam baits_target.bed -o Normal.targetcoverage.cnn
cnvkit.py coverage Normal.bam my_antitargets.bed -o Normal.antitargetcoverage.cnn

#copy number reference from Normal.bam
cnvkit.py reference Normal.targetcoverage.cnn --fasta reference.fasta -o my_reference.cnn

#combine the uncorrected target and antitarget coverage tables (.cnn) and correct for biases
cnvkit.py fix Tumor.targetcoverage.cnn Tumor.antitargetcoverage.cnn my_reference.cnn -o Tumor.cnr
