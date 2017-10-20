/*
 *
 *   RNASEQ-ENCODE-NF is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with RNASEQ-ENCODE-NF.  If not, see <http://www.gnu.org/licenses/>.
 */
 
 
/* 
 * Example of RNAseq pipeline for ENCODE data implemented with Nextflow
 * This example was built based on the Nextflow RNA-seq pipeline https://github.com/nextflow-io/rnaseq-nf
 *  
 * Authors:
 * - Francesco Strozzi <francesco.strozzi@gmail.com>
 * - Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * - Evan Floden <evanfloden@gmail.com>
*/ 

 
/*
 * Default pipeline parameters. They can be overriden on the command line eg. 
 * given `params.foo` specify on the run command line `--foo some_value`.  
 */
 
params.transcriptome = "ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
params.metadata = "$baseDir/data/metadata.small.tsv"
params.output = "results"

log.info """\
         R N A S E Q - N F   E N C O D E   P I P E L I N E
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
    file transcriptome from file(params.transcriptome)   
 
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

    input:
    file(metadata) from Channel.fromPath(params.metadata)

	output:
    stdout into (encode_files_ch_1, encode_files_ch_2)

    """
    #!/usr/bin/env python

    from __future__ import print_function
    import sys
    from collections import defaultdict

    sra = []
    with open("$metadata") as f:
        next(f)
        for line in f:
            data = line.rstrip().split(\"\t\")
            file_type = data[1]
            sample_type = data[6].replace('\\'','').replace(' ','_')
            if data[40]:
                dbxref = data[40].split(':')[1]
                seq_type = data[33]
                strand_specific = data[23]
                if file_type == "fastq" and seq_type == "paired-ended":
                    if not dbxref in sra:
                        sra_id = dbxref.split("SRR")[1]
                        url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{0}/".format(dbxref[0:6])
                        if len(sra_id) == 6: url += "{0}".format(dbxref)
                        if len(sra_id) == 7: url += "00{0}/{1}".format(sra_id[-1],dbxref)
                        if len(sra_id) == 8: url += "0{0}/{1}".format(sra_id[-2:-1],dbxref)
                        if len(sra_id) == 9: url += "{0}/{1}".format(sra_id[-3:-1],dbxref)
                        print(",".join([dbxref,sample_type,strand_specific,url]))
                        sra.append(dbxref)
    """
}

process quant {
    
    tag "$dbxref"
    
    cpus 8

    memory '16 GB' 
 
    input:
    file index from index_ch
    set dbxref,sample_type,strand_specific,url from encode_files_ch_1.splitCsv()
 
    output:
    file("${sample_type}-${dbxref}") into quant_ch

    script:
    def libType = strand_specific == "True" ? "SF" : "U"
    """
    wget -q ${url}/${dbxref}_1.fastq.gz
    wget -q ${url}/${dbxref}_2.fastq.gz
    salmon quant --threads $task.cpus --libType=${libType} -i index -1 ${dbxref}_1.fastq.gz -2 ${dbxref}_2.fastq.gz -o ${sample_type}-${dbxref}
    """
}
  
process fastqc {
    
    tag "FASTQC on $dbxref"

    cpus 2

    memory '8 GB'
 
    input:
    set dbxref,sample_type,strand_specific,url from encode_files_ch_2.splitCsv()

    output:
    file("fastqc_${dbxref}_logs") into fastqc_ch


    script:
    """
    wget -q ${url}/${dbxref}_1.fastq.gz
    wget -q ${url}/${dbxref}_2.fastq.gz
    mkdir fastqc_${dbxref}_logs
    fastqc -t $task.cpus -o fastqc_${dbxref}_logs -f fastq -q ${dbxref}_1.fastq.gz ${dbxref}_2.fastq.gz
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
