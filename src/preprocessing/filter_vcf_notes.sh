variants='variants.txt'
vcf_in='21554_8#9.filtered.vcf'
vcf_out='variants_of_interest.vcf.gz'
HIPSCI_VCF='/hps/nobackup/stegle/projects/hipsci/hipsci_genotypes/unreleased/hipsci.wec.gtarray.unreleased.REL-2016-09.imputed_phased.INFO_0.4_filtered.20160912.genotypes.allchr.temp_for_endodiff_20170217.vcf.gz'

bcftools view -O v "$vcf_in" | grep -v ^# | awk '{sub(/chr/,""); print $1"\t"$2"\t"$4"\t"$5}' > variants.txt

bcftools view -O z -l 9 -R "$variants" "$HIPSCI_VCF"  "$vcf_out"




