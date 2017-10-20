# RNAseq-NF ENCODE pipeline 

A basic pipeline for quantification of genomic features from short read data coming from ENCODE project
implemented with Nextflow. 

This example can also run locally, however instructions are also given to test it specifically on AWS Batch service

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg)](http://nextflow.io)

## Requirements 

* Unix-like operating system (Linux, macOS, etc)
* Java 8 

## Quickstart 

1. If you don't have it already install Docker in your computer. Read more [here](https://docs.docker.com/).

2. Install Nextflow (version 0.26.x or higher):
 
        export NXF_VER=0.26.0-SNAPSHOT
      
        curl -s https://get.nextflow.io | bash

3. Launch the pipeline execution: 

        ./nextflow run fstrozzi/rnaseq-encode-nf -with-docker
        
4. When the execution completes open in your browser the report generated at the following path:

        results/multiqc_report.html 
	
You can see an example report at the following [link](http://multiqc.info/examples/rna-seq/multiqc_report.html).	
	
Note: the very first time you execute it, it will take a few minutes to download the pipeline 
from this GitHub repository and the the associated Docker images needed to execute the pipeline.  


## AWS Batch support

To run this pipeline on AWS using the Batch service, you need to:

1. Create an AWS Batch computing environment and queue, you can skip the job definition since Nextflow will do that for you. Follow the instructions on the [AWS website](http://docs.aws.amazon.com/batch/latest/userguide/Batch_GetStarted.html).

2. Follow the indications on the Nextflow [docs](https://github.com/nextflow-io/nextflow/blob/master/docs/awscloud.rst#allows-batch) to prepare an AMI with enough disk space to run the workflow.

3. Install the AWS keys on the machine where you will execute Nextflow:

        pip install awscli
        aws configure

4. Create a local nextflow.config file to specify the AWS Batch executor and parameters, plus the path of the AWS CLI on the AMI:

        executor {
            name = 'awsbatch'
            awscli = '/scratch/miniconda/bin/aws'
        }

        process {
            queue = 'my-aws-batch-queue'
        }

5. Run the pipeline

        ./nextflow fstrozzi/rnaseq-encode-nf -w s3://bucket/prefix

## Cluster support

RNASeq-NF execution relies on [Nextflow](http://www.nextflow.io) framework which provides an 
abstraction between the pipeline functional logic and the underlying processing system.

This allows the execution of the pipeline in a single computer or in a HPC cluster without modifying it.

Currently the following resource manager platforms are supported:

  + Univa Grid Engine (UGE)
  + Platform LSF
  + SLURM
  + PBS/Torque

By default the pipeline is parallelized by spawning multiple threads in the machine where the script is launched.

To submit the execution to a UGE cluster create a file named `nextflow.config` in the directory
where the pipeline is going to be executed with the following content:

    process {
      executor='uge'
      queue='<queue name>'
    }

To lean more about the avaible settings and the configuration file read the 
Nextflow [documentation](http://www.nextflow.io/docs/latest/config.html).


## Components 

RNASeq-NF uses the following software components and tools: 

* [Salmon](https://combine-lab.github.io/salmon/) 0.8.2
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 0.11.5
* [Multiqc](https://multiqc.info) 1.2

