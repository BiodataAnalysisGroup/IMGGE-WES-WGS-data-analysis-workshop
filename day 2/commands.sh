##################################
#     User defined pipeline      #
################################## 

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
cnvkit.py reference *Normal.{,anti}targetcoverage.cnn --fasta reference.fasta -o my_reference.cnn

#combine the uncorrected target and antitarget coverage tables (.cnn) and correct for biases
cnvkit.py fix Tumor.targetcoverage.cnn Tumor.antitargetcoverage.cnn my_reference.cnn -o Tumor.cnr
cnvkit.py fix Normal.targetcoverage.cnn Normal.antitargetcoverage.cnn my_reference.cnn -o Normal.cnr

#infer discrete copy number segments from the given coverage table
cnvkit.py segment Tumor.cnr -o Tumor.cns

#given segmented log2 ratio estimates (.cns), derive each segment’s absolute integer copy number
cnvkit.py call Tumor.cns -o Tumor.call.cns -v VarScan_somatic.vcf

#plot bin-level log2 coverages and segmentation calls together
cnvkit.py scatter Tumor.cnr -s Tumor.cns  -v VarScan_somatic.vcf -o Tumor-scatter.pdf

#we can try genemetrics both with and without the segment files, 
#take the intersection of those as a list of “trusted” genes, and visualize each of them
cnvkit.py genemetrics -y Tumor.cnr -s Tumor.cns | tail -n+2 | cut -f1 | sort > segment-genes.txt
cnvkit.py genemetrics -y Tumor.cnr | tail -n+2 | cut -f1 | sort > ratio-genes.txt
comm -12 ratio-genes.txt segment-genes.txt > trusted-genes.txt
mkdir gene_scatter_plots

for gene in `cat trusted-genes.txt`
do
    cnvkit.py scatter -s Tumor.cn{s,r} -g $gene -o gene_scatter_plots/Tumor-$gene-scatter.pdf
done

##################################
#      Automated pipeline        #
################################## 

#create baits.bed
#uses .bed file and caverage from .bam to infer the targeted regions 
guess_baits.py Tumor.bam -t gene_targets.bed -o baits.bed

#with batch command
cnvkit.py batch Tumor.bam --normal Normal.bam  --targets baits.bed --annotate refFlat.txt  --fasta reference.fasta --access access-5kb-mappable.hg19_chr5_chr12_chr17.bed --output-reference my_reference.cnn --output-dir results/ --diagram --scatter

#we can try genemetrics both with and without the segment files, 
#take the intersection of those as a list of “trusted” genes, and visualize each of them
cnvkit.py genemetrics -y Tumor.cnr -s Tumor.cns | tail -n+2 | cut -f1 | sort > segment-genes.txt
cnvkit.py genemetrics -y Tumor.cnr | tail -n+2 | cut -f1 | sort > ratio-genes.txt
comm -12 ratio-genes.txt segment-genes.txt > trusted-genes.txt
mkdir gene_scatter_plots

for gene in `cat trusted-genes.txt`
do
    cnvkit.py scatter -s Tumor.cn{s,r} -g $gene -o gene_scatter_plots/Tumor-$gene-scatter.pdf
done