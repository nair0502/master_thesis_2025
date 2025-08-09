#!/bin/bash
PATH=$PATH:....../sratoolkit.3.1.1-centos_linux64/bin/fasterq-dump
find vaishnavi/vaishnavi_masterthesis/sra_id -name "*.sra" -exec fasterq-dump {} -O vaishnavi/vaishnavi_masterthesis/fastq \;
