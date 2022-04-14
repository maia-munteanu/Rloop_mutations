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
    log.info '    --fasta_ref                    FILE           Fasta reference file.'
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
params.fasta_ref = "/g/strcombio/fsupek_cancer1/Rloop_clusters_project/hg19.fa"

fasta_ref=file(params.fasta_ref)
Unclustered_bed=file(params.Unclustered_bed)
Surrounding_bed=file(params.Surrounding_bed)
Rloops_bed=file(params.Rloops_bed)

pairs_list = Channel.fromPath(params.input_file, checkIfExists: true).splitCsv(header: true, sep: '\t', strip: true)
                   .map{ row -> [ row.sample, file(row.snv) ] }.view()
 

process get_vcfs {

       publishDir params.output_folder, mode: 'move', pattern: '*.snv.vcf*'
       tag {sample}

       input:
       set val(sample), file(snv) from pairs_list
       file fasta_ref
       file Rloops_bed
       file Surrounding_bed
       file Unclustered_bed
    
    
       output:
       set val(sample), file(snv), file("*snv*") into vcfs

       shell:
       '''
       bcftools view -f 'PASS' !{snv} -Oz > !{sample}.filt.vcf.gz
       tabix -p vcf !{sample}.filt.vcf.gz
       bcftools view --types snps --regions-file !{Unclustered_bed} !{sample}.filt.vcf.gz | bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_unclustered.snv.vcf.gz
       bcftools view --types snps --regions-file !{Surrounding_bed} !{sample}.filt.vcf.gz | bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_surrounding.snv.vcf.gz
       bcftools view --types snps --regions-file !{Rloops_bed} !{sample}.filt.vcf.gz | bcftools norm -d all -f !{fasta_ref} | bcftools sort -Oz > !{sample}_rloops.snv.vcf.gz
       gunzip *.snv.vcf.gz
       '''
}
             
