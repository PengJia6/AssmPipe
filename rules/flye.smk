import os

flye = "/data/home/pengjia//miniconda3/envs/flye/bin/flye"
dir_work = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/data/HiFi/"
dir_assm = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/assm_res/"


rule all:
    input:
        expand(dir_assm + "02_flye_{para}/{assm}.{hap}.flye_{para}.ok",para=["default"],
            hap=["hap1", "hap2", "collapsed"],assm="ChineseQuartet")

wildcard_constraints:
    assm="ChineseQuartet"


def get_fqs(wildcards):
    fqs = []
    if wildcards.hap in ["hap1", "hap2"]:
        fqs_tag = expand(dir_work + "{sample}/{sample}.T2T.minimap2.HiFi.{{hap}}.fq.gz",sample=["LCL5",
                                                                                                "LCL6"],hap=wildcards.hap)
        fqs_untag = expand(dir_work + "{sample}/{sample}.T2T.minimap2.HiFi.untag.{{hap}}.fq.gz",sample=["LCL5",
                                                                                                        "LCL6"],hap=wildcards.hap)
        fqs.extend(fqs_untag)
        fqs.extend(fqs_tag)
    elif wildcards.hap == "collapsed":
        for sample in ["LCL5", "LCL6"]:
            files = os.listdir("/data/DATA/ChineseQuartet/RAWDATA/HiFi/fastqs/{}".format(sample))
            files =[i for i in files if i.endswith("fq.gz")]
            for i in files:
                fqs.append("/data/DATA/ChineseQuartet/RAWDATA/HiFi/fastqs/{s}/{i}".format(s=sample,i=i))
    return fqs


rule flye:
    input:
        fqs=get_fqs
    output:
        dir_assm + "02_flye_{para}/{assm}.{hap}.flye_{para}.ok"
    threads: 60
    log:
        dir_assm + "logs/02_flye_{para}/{assm}.{hap}.flye_{para}.log"
    benchmark:
        dir_assm + "logs/02_flye_{para}/{assm}.{hap}.flye_{para}.tsv"
    run:
        output_pre = str(output)[:-3]
        if str(wildcards.para) == "default":
            shell("{flye} --pacbio-hifi {input.fqs} --genome-size 3.1g --out-dir {output_pre} --threads {threads}  2>{log} 1>{log}")
        shell("touch {output}")