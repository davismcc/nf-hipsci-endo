genome=/homes/stegle/research/datasets/references/human/GRCh37/GRCh37.p13.genome.fa
ercc=/homes/stegle/research/datasets/references/ERCC92/ERCC92.fa
annotation=/homes/stegle/research/datasets/references/human/GRCh37/gencode.v19.annotation_ERCC.gtf

STAR --runMode genomeGenerate --genomeDir /nfs/research/stegle/datasets/references/human/STAR_GRCh37.75_ERCC --genomeFastaFiles $genome $ercc --runThreadN 5 --sjdbGTFfile /nfs/research/stegle/datasets/references/human/GRCh37/Homo_sapiens.GRCh37.75.gtf --sjdbOverhang 100
