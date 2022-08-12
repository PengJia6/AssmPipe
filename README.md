# AssmPipe
Pipeline for genome assembly


## Quick start
### 1. change the snakemake file according to you environment.

   


### 2. run with snakemake           

       nohup snakemake -s pipeline.smk -j 10 -k --ri >sublog 2>&1 &
       nohup snakemake -s pipeline.smk -j 10 -k --ri --cluster "qsub -l nodes=1:ppn=20 -l walltime=999:00:00" >sublog 2>&1 & 
      

## Support tools 
  - hifiasm
  - hicanu
  - flye
  - shasta
  - canu
  - ...
## Contribution 
   If you want to apply other tools to this pipeline, we encourage you to pull a request or email us.
   

## Contact  
 * Peng Jia at Xi'an Jiaotong University (pengjia@stu.xjtu.edu.cn)
