#! /usr/bin/env nextflow

//vim: syntax=groovy -*- mode: groovy;-*-

// Copyright (C) 2022 IRB Barcelona

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

log.info ""
log.info "--------------------------------------------------------------------------"
log.info "  Rloop_mutations: nextflow pipeline to extract SNVs found in/near Rloops "
log.info "--------------------------------------------------------------------------"
log.info "Copyright (C) IRB Barcelona"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------------------------"
log.info ""

params.help = null
if (params.help) {
    log.info ''
    log.info '--------------------------------------------------'
    log.info '  USAGE              '
    log.info '--------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run main.nf --input_file Hartwig_all_samples.csv --closer_value 1000 --close_value 10000 --output_folder /g/*/*cancer1*/SV* --fasta_file /g/*/*cancer1*/SV*/hg19.fasta --hg19 /g/*/*cancer1*/SV*/hg19.genome'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --input_file                   FILE           Input .csv file containing 2 columns: sample name and SNV vcf path.'
    log.info '    --output_folder                FOLDER         Output folder.'
    log.info '    --Rloops_bed                   FILE           Bed containing the coordinates of Rloops (strand agnostic).'
    log.info '    --Surrounding_bed              FILE           Bed containing the coordinates of regions surrounding Rloops (strand agnostic and up/downstream agnostic.'
    log.info '    --Unclustered_bed              FILE           Bed containing the cooridnates of all other genomic regions not covered by the Rloops and Surrounding beds.'
    log.info ''
    log.info 'Flags:'
    log.info '    --help                                        Display this message'
    log.info ''
    exit 0
}

params.input_file = "/g/strcombio/fsupek_cancer1/Rloop_clusters_project/Hartwig_Lung+Esophagus_samples.csv"
params.Rloops_bed = "/g/strcombio/fsupek_cancer1/Rloop_clusters_project/Rloops.hg19.bed"
params.Surrounding_bed = "/g/strcombio/fsupek_cancer1/Rloop_clusters_project/Surrounding.hg19.bed"
params.Unclustered_bed = "/g/strcombio/fsupek_cancer1/Rloop_clusters_project/Unclustered.hg19.bed"
params.output_folder = "/g/strcombio/fsupek_cancer1/Rloop_clusters_project/SNV_clusters_VCFs/"

pairs_list = Channel.fromPath(params.input_file, checkIfExists: true).splitCsv(header: true, sep: '\t', strip: true)
                   .map{ row -> [ row.sample, file(row.snv) ] }.view()
                   
process get_vcfs {

      publishDir params.output_folder, mode: 'copy', pattern: '*.vcf.gz'
      tag {sample}

       input:
       set val(sample), file(snv) from pairs_list




}
                   
process make_sv_beds {

       publishDir params.output_folder+"/SV_BEDs/", mode: 'copy', pattern: '*.bed'

       tag {sample}

       input:
       set val(sample), file(sv), file(snv) from pairs_list
       file hg19

       output:
       set val(sample), file(sv), file(snv), file("*unclustered.bed"), file("*close.bed"), file("*closer.bed"), file("*cluster.bed") into beds
       

       shell:
       '''
       close_bp=!{params.close_value}
       bp_per_kb=1000
       close=$((close_bp / bp_per_kb))
       echo $close
    
       bcftools view -f 'PASS' !{sv} -Oz > !{sample}.sv.filt.vcf.gz
       
       Rscript !{baseDir}/vcf_to_bed.R --VCF !{sample}.sv.filt.vcf.gz --close !{params.close_value} --closer !{params.closer_value}
       bedtools complement -i !{sample}_0_${close}kb_cluster.bed -g !{hg19} > !{sample}_unclustered.bed
       '''
  }


process make_vcfs {
    
    publishDir params.output_folder+"/SNV_clusters_VCFs/", mode: 'move', pattern: '*.snv.vcf*'    
    publishDir params.output_folder+"/MNV_clusters_VCFs/", mode: 'move', pattern: '*.mnv.vcf*'
    publishDir params.output_folder+"/INDEL_clusters_VCFs/", mode: 'move', pattern: '*.indel.vcf*'
    
    tag {sample}

    input:
    set val(sample), file(sv), file(snv), file("*unclustered.bed"), file("*close.bed"), file("*closer.bed"), file("*cluster.bed") from beds
    file fasta_ref
    
    output:
    set val(sample), file(sv), file(snv), file("*snv*"), file("*mnv*"), file("*indel*") into vcfs
    
    shell:
    '''
    close_bp=!{params.close_value}
    closer_bp=!{params.closer_value}
    bp_per_kb=1000
    close=$((close_bp / bp_per_kb))
    closer=$((closer_bp / bp_per_kb))
    echo $close $closer
   
    bcftools view -f 'PASS' !{snv} -Oz > !{sample}.filt.vcf.gz
    tabix -p vcf !{sample}.filt.vcf.gz
    
    bcftools view -i '%QUAL>500' -f PASS --types snps --regions-file unclustered.bed !{sample}.filt.vcf.gz | bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_unclustered_highc.snv.vcf.gz
    bcftools view -i '%QUAL>500' -f PASS --types snps --regions-file closer.bed !{sample}.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_clustered_0_${closer}kb_highc.snv.vcf.gz
    bcftools view -i '%QUAL>500' -f PASS --types snps --regions-file close.bed !{sample}.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_clustered_${closer}kb_${close}kb_highc.snv.vcf.gz
    bcftools view -i '%QUAL<=500' -f PASS --types snps --regions-file unclustered.bed !{sample}.filt.vcf.gz | bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_unclustered_lowc.snv.vcf.gz
    bcftools view -i '%QUAL<=500' -f PASS --types snps --regions-file closer.bed !{sample}.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_clustered_0_${closer}kb_lowc.snv.vcf.gz
    bcftools view -i '%QUAL<=500' -f PASS --types snps --regions-file close.bed !{sample}.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_clustered_${closer}kb_${close}kb_lowc.snv.vcf.gz
    gunzip *.snv.vcf.gz
    
    bcftools view -i '%QUAL>500' -f PASS --types mnps --regions-file unclustered.bed !{sample}.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_unclustered_highc.mnv.vcf.gz
    bcftools view -i '%QUAL>500' -f PASS --types mnps --regions-file closer.bed !{sample}.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_clustered_0_${closer}kb_highc.mnv.vcf.gz
    bcftools view -i '%QUAL>500' -f PASS --types mnps --regions-file close.bed !{sample}.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_clustered_${closer}kb_${close}kb_highc.mnv.vcf.gz
    bcftools view -i '%QUAL<=500' -f PASS --types mnps --regions-file unclustered.bed !{sample}.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_unclustered_lowc.mnv.vcf.gz
    bcftools view -i '%QUAL<=500' -f PASS --types mnps --regions-file closer.bed !{sample}.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_clustered_0_${closer}kb_lowc.mnv.vcf.gz
    bcftools view -i '%QUAL<=500' -f PASS --types mnps --regions-file close.bed !{sample}.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_clustered_${closer}kb_${close}kb_lowc.mnv.vcf.gz
    gunzip *.mnv.vcf.gz
    
    bcftools view -i '%QUAL>500' -f PASS --types indels --regions-file unclustered.bed !{sample}.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_unclustered_highc.indel.vcf.gz
    bcftools view -i '%QUAL>500' -f PASS --types indels --regions-file closer.bed !{sample}.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_clustered_0_${closer}kb_highc.indel.vcf.gz
    bcftools view -i '%QUAL>500' -f PASS --types indels --regions-file close.bed !{sample}.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_clustered_${closer}kb_${close}kb_highc.indel.vcf.gz
    bcftools view -i '%QUAL<=500' -f PASS --types indels --regions-file unclustered.bed !{sample}.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_unclustered_lowc.indel.vcf.gz
    bcftools view -i '%QUAL<=500' -f PASS --types indels --regions-file closer.bed !{sample}.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_clustered_0_${closer}kb_lowc.indel.vcf.gz
    bcftools view -i '%QUAL<=500' -f PASS --types indels --regions-file close.bed !{sample}.filt.vcf.gz |  bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_clustered_${closer}kb_${close}kb_lowc.indel.vcf.gz
    gunzip *.indel.vcf.gz
    '''
    
} 
