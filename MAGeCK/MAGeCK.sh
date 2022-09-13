mageck count -l 'Brie library(mageck)+non-target.csv' -n demo --sample-label lib,input1,input2,outpu1,output2  --fastq Brie.top20.fasta.gz Input1.top20.fasta.gz  Input2.top20.fasta.gz Output1.top20.fasta.gz Output2.top20.fasta.gz 

mageck test -k demo.count.txt -c input1,input2 -t outpu1,output2 -n demo
