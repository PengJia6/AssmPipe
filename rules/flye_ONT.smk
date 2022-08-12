import os


flye = "/data/home/pengjia//miniconda3/envs/flye/bin/flye"
dir_assm = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/assm_res/"
dir_assm_data_hap = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/hap_assm/01_assm_data/"
dir_assm_data_collapsed = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/collapse_assm/01_assm_data/"

rule all:
    input:
        expand(dir_assm + "02_flye_ONT_{para}/{assm}.{hap}.flye_ONT_{para}.ok",para=["default"],
            hap=["hap1", "hap2", "collapsed"],assm="ChineseQuartet")

wildcard_constraints:
    assm="ChineseQuartet"
def get_fqs(wildcards):
    fqs = []
    if wildcards.hap in ["hap1", "hap2"]:
        fqs_ont = expand(dir_assm_data_hap + "{assm}.{sample}.{hap}.ONT.long.fq.gz",
            sample=["LCL56", ],hap=wildcards.hap,assm=wildcards.assm)
        fqs_ulont = expand(dir_assm_data_hap + "{assm}.LCL5.{hap}.ONTUL.long.fq.gz",
            sample=["LCL5"],hap=wildcards.hap,assm=wildcards.assm)
        fqs.extend(fqs_ont)
        fqs.extend(fqs_ulont)
    elif wildcards.hap == "collapsed":
        fqs_ont = expand(dir_assm_data_collapsed + "{assm}.{sample}.ONT.long2.fq.gz",
            sample=["LCL56", ],hap=wildcards.hap,assm=wildcards.assm)
        fqs_ulont = expand(dir_assm_data_collapsed + "{assm}.LCL5.ONTUL.long2.fq.gz",
            sample=["LCL5"],hap=wildcards.hap,assm=wildcards.assm)
        fqs.extend(fqs_ont)
        fqs.extend(fqs_ulont)

    return fqs


rule flye:
    input:
        fqs=get_fqs
    output:
        dir_assm + "02_flye_ONT_{para}/{assm}.{hap}.flye_ONT_{para}.ok"
    threads: 60
    log:
        dir_assm + "logs/02_flye_ONT_{para}/{assm}.{hap}.flye_ONT_{para}.log"
    benchmark:
        dir_assm + "logs/02_flye_ONT_{para}/{assm}.{hap}.flye_ONT_{para}.tsv"
    run:
        output_pre = str(output)[:-3]
        if str(wildcards.para) == "default":
            shell("{flye} --nano-raw {input.fqs} --genome-size 3.1g --out-dir {output_pre} --threads {threads}  2>{log} 1>{log}")
        shell("touch {output}")