#!/bin/bash
PATH=$PATH:/dss/dsshome1/08/ge85jek2/sratoolkit.3.1.1-centos_linux64/bin/fasterq-dump
find /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/vaishnavi/vaishnavi_masterthesis/sra_id -name "*.sra" -exec fasterq-dump {} -O /dss/dssfs03/pn57ba/pn57ba-dss-0001/computational-plant-biology/vaishnavi/vaishnavi_masterthesis/fastq \;
