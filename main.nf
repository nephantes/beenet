$HOSTNAME = ""
params.outdir = 'results'  


if (!params.custom_fastq_file){params.custom_fastq_file = ""} 
if (!params.Ref){params.Ref = ""} 
if (!params.reads){params.reads = ""} 
if (!params.barcode_num){params.barcode_num = ""} 
if (!params.params){params.params = ""} 
if (!params.star_args){params.star_args = ""} 
if (!params.mate){params.mate = ""} 
// Stage empty file to be used as an optional input where required
ch_empty_file_1 = file("$baseDir/.emptyfiles/NO_FILE_1", hidden:true)
ch_empty_file_2 = file("$baseDir/.emptyfiles/NO_FILE_2", hidden:true)

g_3_custom_fasta2_g_2 = params.custom_fastq_file && file(params.custom_fastq_file, type: 'any').exists() ? file(params.custom_fastq_file, type: 'any') : ch_empty_file_1
Channel.value(params.Ref).set{g_4_ref_flat0_g_1}
if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_6_reads2_g_5}
 } else {  
	g_6_reads2_g_5 = Channel.empty()
 }

Channel.value(params.barcode_num).set{g_7_barcode3_g_5}
Channel.value(params.params).set{g_8_params4_g_5}
Channel.value(params.star_args).set{g_10_params2_g_0}
Channel.value(params.mate).set{g_11_mate5_g_5}



process Add_custom_seq_to_genome_gtf {

input:
 file custom_fasta from g_3_custom_fasta2_g_2

output:
 file "${genomeName}_custom.fa"  into g_2_genome00_g_0
 file "${gtfName}_custom_sorted.gtf"  into g_2_gtfFile11_g_0

when:
params.add_sequences_to_reference == "yes"

script:
genomeName = genome.baseName
gtfName = gtf.baseName
is_custom_genome_exists = custom_fasta.name.startsWith('NO_FILE') ? "False" : "True" 

"""
#!/usr/bin/env python 
import requests
import os
import pandas as pd
import re
import urllib
from Bio import SeqIO

def add_to_fasta(seq, sqid, out_name):
	new_line = '>' + sqid + '\\n' + seq + '\\n'
	with open(out_name + '.fa', 'a') as f:
		f.write(new_line)

def createCustomGtfFromFasta(fastaFile, outCustomGtfFile):

    fasta_sequences = SeqIO.parse(open(fastaFile),'fasta')
    with open(outCustomGtfFile, "w") as out_file:
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            last = len(sequence)
            line1 = "{gene}\\tKNOWN\\tgene\\t{first}\\t{last}\\t.\\t+\\t.\\tgene_id \\"{gene}\\"; gene_version \\"1\\"; gene_type \\"protein_coding\\"; gene_source \\"KNOWN\\"; gene_name \\"{gene}\\"; gene_biotype \\"protein_coding\\"; gene_status \\"KNOWN\\"; level 1;".format(gene=name, first="1", last=last)
            line2 = "{gene}\\tKNOWN\\ttranscript\\t{first}\\t{last}\\t.\\t+\\t.\\tgene_id \\"{gene}\\"; gene_version \\"1\\"; transcript_id \\"{gene}_trans\\"; transcript_version \\"1\\"; gene_type \\"protein_coding\\"; gene_source \\"KNOWN\\"; transcript_source \\"KNOWN\\"; gene_status \\"KNOWN\\"; gene_name \\"{gene}\\"; gene_biotype \\"protein_coding\\"; transcript_type \\"protein_coding\\"; transcript_status \\"KNOWN\\"; transcript_name \\"{gene}_1\\"; level 1; tag \\"basic\\"; transcript_biotype \\"protein_coding\\"; transcript_support_level \\"1\\";".format(gene=name, first="1", last=last)
            line3 = "{gene}\\tKNOWN\\texon\\t{first}\\t{last}\\t.\\t+\\t.\\tgene_id \\"{gene}\\"; gene_version \\"1\\"; transcript_id \\"{gene}_trans\\"; transcript_version \\"1\\"; exon_number 1; gene_type \\"protein_coding\\"; gene_source \\"KNOWN\\"; transcript_source \\"KNOWN\\"; gene_status \\"KNOWN\\"; gene_name \\"{gene}\\"; gene_biotype \\"protein_coding\\"; transcript_type \\"protein_coding\\"; transcript_status \\"KNOWN\\"; transcript_biotype \\"protein_coding\\"; transcript_name \\"{gene}_1\\"; exon_number 1; exon_id \\"{gene}.1\\"; level 1; tag \\"basic\\"; transcript_support_level \\"1\\";".format(gene=name, first="1", last=last)
            out_file.write("{}\\n{}\\n{}\\n".format(line1, line2, line3))

	
os.system('cp ${genomeName}.fa ${genomeName}_custom.fa')  
os.system('cp ${gtfName}.gtf ${gtfName}_custom.gtf')  

if ${is_custom_genome_exists}:
	os.system('cat ${custom_fasta} >> ${genomeName}_custom.fa')
	createCustomGtfFromFasta("${custom_fasta}", "${custom_fasta}.gtf")
	os.system('cat ${custom_fasta}.gtf >> ${gtfName}_custom.gtf')
	
os.system('samtools faidx ${genomeName}_custom.fa')
os.system('igvtools sort ${gtfName}_custom.gtf ${gtfName}_custom_sorted.gtf')
os.system('igvtools index ${gtfName}_custom_sorted.gtf')

"""
}

g_10_params2_g_0= g_10_params2_g_0.ifEmpty("") 


process make_ref {

input:
 file genome from g_2_genome00_g_0
 file gtf from g_2_gtfFile11_g_0
 val star_args from g_10_params2_g_0

output:
 val "ref"  into g_0_ref_flat00_g_5

when:
params.add_sequences_to_reference == "yes"

script:
"""
#!/bin/sh 
beenet make-ref --ref=${genome} --gtf=${gtf} --star-args='${star_args}' --out=ref
"""
}


if (!(params.run_Download_Genomic_Sources == "yes")){
g_4_ref_flat0_g_1.set{g_1_ref_flat01_g_5}
} else {

process Check_Genome_GTF {

input:
 val ref from g_4_ref_flat0_g_1

output:
 val "${ref}"  into g_1_ref_flat01_g_5

when:
params.run_Download_Genomic_Sources == "yes"

script:
genome_dir  = params.genome_dir
"""
[ ! -d "${genome_dir}/${ref}" ] && beenet download-ref ${ref} && mv ${ref} ${genome_dir} || echo ${ref}
"""




}
}


g_0_ref_flat00_g_5= g_0_ref_flat00_g_5.ifEmpty("") 
g_8_params4_g_5= g_8_params4_g_5.ifEmpty("") 


process Analyze {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /out\/.*$/) "Report/$filename"}
input:
 val custom_ref from g_0_ref_flat00_g_5
 val ref from g_1_ref_flat01_g_5
 set val(name),file(reads) from g_6_reads2_g_5
 val num_barcodes from g_7_barcode3_g_5
 val other_params from g_8_params4_g_5
 val mate from g_11_mate5_g_5

output:
 file "out/*"  into g_5_outputDir00

script:
ref = custom_ref ? custom_ref : ref
"""
mkdir -p out
#gzip -f ${reads}
beenet analyze --sample-name=${name} --ref=${params.genome_dir}/${ref}/${params.star_version} --num-barcodes=${num_barcodes} ${other_params} --out=out *.gz
"""
}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
