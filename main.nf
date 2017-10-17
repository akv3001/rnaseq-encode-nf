/*
 * Copyright (c) 2013-2017, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'RNASEQ-NF'.
 *
 *   RNASEQ-NF is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   RNASEQ-NF is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with RNASEQ-NF.  If not, see <http://www.gnu.org/licenses/>.
 */
 
 
/* 
 * Proof of concept of a RNAseq pipeline for ENCODE data implemented with Nextflow
 * 
 * Authors:
 * - Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * - Emilio Palumbo <emiliopalumbo@gmail.com> 
 * - Evan Floden <evanfloden@gmail.com> 
 * - Francesco Strozzi <francesco.strozzi@gmail.com>
*/ 

 
/*
 * Default pipeline parameters. They can be overriden on the command line eg. 
 * given `params.foo` specify on the run command line `--foo some_value`.  
 */
 
params.transcriptome = "ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
params.metadata = "metadata.tsv"
params.output = "results"

log.info """\
         R N A S E Q - N F   P I P E L I N E   E N C O D E   
         =================================================
         transcriptome: ${params.transcriptome}
         metadata     : ${params.metadata}
         output       : ${params.output}
         """
         .stripIndent()


process index {
    
    tag "$transcriptome"
   
    cpus 4

    memory '30 GB'
 
    input:
    file(transcriptome) from Channel.fromPath(params.transcriptome)
     
    output:
    file 'index' into index_ch

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}
 
 
process parseEncode {

    tag "$params.metadata"

    cpus 2

    memory '4 GB'

    container 'job-definition://python3'

    input:
    file(metadata) from Channel.fromPath(params.metadata)

	output:
	stdout encode_files_ch_1
    stdout encode_files_ch_2
    
    """
    #!/usr/bin/env python

    from __future__ import print_function
    import sys
    from collections import defaultdict

    pairs = defaultdict(list)
    for line in open("$metadata",encoding='utf-8'):
        data = line.rstrip().split(\"\t\")
        file_type = data[1]
        sample_type = data[6].replace('\\'','').replace(' ','_')
        file_url = data[41]
        dbxref = data[40].replace(':','-')
        seq_type = data[33]
        strand_specific = data[23]
        if file_type == "fastq" and seq_type == "paired-ended":
            pairs[dbxref].append(file_url)
            if len(pairs[dbxref]) == 2:
                print(",".join([dbxref,sample_type,strand_specific]+pairs[dbxref]+list(map(lambda x: x.split("/")[-1],pairs[dbxref])) ))
                break
    """

}

process quant {
    
    tag "$dbxref"
    
    cpus 8

    memory '8 GB' 
 
    input:
    file index from index_ch
    set dbxref,sample_type,strand_specific,fastq_url_1,fastq_url_2,fastq_1,fastq_2 from encode_files_ch_1.splitCsv()
 
    output:
    file("${sample_type}-${dbxref}") into quant_ch

    script:
    def libType = strand_specific == "True" ? "S" : "U"
    """
    wget $fastq_url_1
    wget $fastq_url_2
    salmon quant --threads $task.cpus --libType=${libType} -i index -1 ${fastq_1} -2 ${fastq_2} -o ${sample_type}-${dbxref}
    """
}
  
process fastqc {
    
    tag "FASTQC on $dbxref"

    cpus 2

    memory '8 GB'
 
    input:
    set dbxref,sample_type,strand_specific,fastq_url_1,fastq_url_2,fastq_1,fastq_2 from encode_files_ch_2.splitCsv()

    output:
    file("fastqc_${dbxref}_logs") into fastqc_ch


    script:
    """
    wget $fastq_url_1
    wget $fastq_url_2
    mkdir fastqc_${dbxref}_logs
    fastqc -o fastqc_${dbxref}_logs -f fastq -q {$fastq_1} ${fastq_2}
    """  
} 
  
  
process multiqc {
    
    cpus 2

    memory '8 GB'
    
    publishDir params.output, mode:'copy'
       
    input:
    file('*') from quant_ch.mix(fastqc_ch).collect()
    
    output:
    file('multiqc_report.html')  
     
    script:
    """
    multiqc . 
    """
}
 
workflow.onComplete { 
	println ( workflow.success ? "\nDone! Open the following report in your browser --> $params.output/multiqc_report.html\n" : "Oops .. something went wrong" )
}
