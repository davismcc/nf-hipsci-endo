chrom1=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr1.fa
chrom2=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr2.fa
chrom3=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr3.fa
chrom4=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr4.fa
chrom5=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr5.fa
chrom6=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr6.fa
chrom7=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr7.fa
chrom8=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr8.fa
chrom9=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr9.fa
chrom10=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr10.fa
chrom11=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr11.fa
chrom12=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr12.fa
chrom13=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr13.fa
chrom14=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr14.fa
chrom15=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr15.fa
chrom16=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr16.fa
chrom17=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr17.fa
chrom18=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr18.fa
chrom19=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr19.fa
chrom20=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr20.fa
chrom21=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr21.fa
chrom22=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chr22.fa
chromY=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chrY.fa
chromX=/nfs/research/stegle/datasets/references/human/hg19/chromosomes/chrX.fa
ercc=/nfs/research/stegle/datasets/references/ERCC92/ERCC92.fa


cd /homes/stegle/research/datasets/references/human/tophat2-GRHC37_ERCC
bowtie2-build ./GRCh37.p13.genome.fa,./ERCC92.fa GRCh37.p13.genome_ERCC

#bowtie2-build $chrom1,$chrom2,$chrom3,$chrom4,$chrom5,$chrom6,$chrom7,$chrom8,$chrom9,$chrom10,$chrom11,$chrom12,$chrom13,$chrom14,$chrom15,$chrom16,$chrom17,$chrom18,$chrom19,$chrom20,$chrom21,$crom22,$chromX,$chromY,$ercc /homes/stegle/research/datasets/references/human/bowtie2-GRHC37_ERCC/GRCH37_ERCC

#tophat2 -G /nfs/research/stegle/datasets/references/human/GRCh37/Homo_sapiens.GRCh37.75.gtf --transcriptome-index /homes/stegle/research/datasets/references/human/tophat2-GRHC37_ERCC/GRCH37_ERCC /homes/stegle/research/datasets/references/human/bowtie2-GRHC37_ERCC/GRCH37_ERCC

