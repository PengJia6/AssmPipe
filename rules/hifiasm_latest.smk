# hifiasm = "/data/home/pengjia/mybin/hifiasm"
hifiasm = "/data/home/pengjia//mysoftware/assm/hifiasm_latest/hifiasm/hifiasm"
dir_work = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/data/HiFi/"
dir_assm = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/assm_res/"


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
            for i in files:
                if not i.endswith("fq.gz"): continue
                fqs.append("/data/DATA/ChineseQuartet/RAWDATA/HiFi/fastqs/{s}/{i}".format(s=sample,i=i))
    return fqs


wildcard_constraints:
    assm="ChineseQuartet",
    hap="|".join(["hap1", "hap2", "collapsed"])

rule all:
    input:
        expand(dir_assm + "02_hifiasm_latest_{para}/{assm}.{hap}.hifiasm_latest_{para}/{assm}.{hap}.hifiasm_latest_{para}.ok",para=[
            "default"],
            hap=["hap1", "hap2", "collapsed"],assm="ChineseQuartet"),
        expand(dir_assm + "02_hifiasm_latest_trio/{assm}.hifiasm_latest_trio/{assm}.hifiasm_latest_trio.ok",assm="ChineseQuartet")

rule hifiasm0_15_5:
    input:
        get_fqs
    output:
        dir_assm + "02_hifiasm_latest_{para}/{assm}.{hap}.hifiasm_latest_{para}/{assm}.{hap}.hifiasm_latest_{para}.ok"
    threads: 120
    log:
        dir_assm + "logs/02_hifiasm_latest_{para}/{assm}.{hap}.hifiasm_latest_{para}.log"
    benchmark:
        dir_assm + "logs/02_hifiasm_latest_{para}/{assm}.{hap}.hifiasm_latest_{para}.tsv"
    run:
        output_pre = str(output)[:-3]
        if str(wildcards.para) == "default":
            shell("{hifiasm} -o {output_pre}  -t {threads} {input} 2>{log} 1>{log}")
        elif str(wildcards.para) == "l0":
            shell("{hifiasm} -o {output_pre}  -l0 -t {threads} {input} 2>{log} 1>{log}")
        elif str(wildcards.para) == "l1":
            shell("{hifiasm} -o {output_pre}  -l1 -t {threads} {input} 2>{log} 1>{log}")
        elif str(wildcards.para) == "l2":
            shell("{hifiasm} -o {output_pre}  -l2 -t {threads} {input} 2>{log} 1>{log}")
        elif str(wildcards.para) == "l3":
            shell("{hifiasm} -o {output_pre}  -l3 -t {threads} {input} 2>{log} 1>{log}")
        shell("touch {output}")


samples = ["LCL5", "LCL6", "LCL7", "LCL8"]
samples_info = {}
for sample in samples:
    info = open("/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/rawdata/HiFi/fofn/{}.HiFi.fofn".format(sample),"r").readlines()
    # thisinfo = {i.rstrip(".bam\n").rstrip(".fq.gz").rstrip(".fastq.gz").split("/")[-1]: i.rstrip("\n") for i in info}
    samples_info[sample] = [i.rstrip("\n") for i in info]


rule hifiasm_trio_binning:
    input:
        LCL5_fq_list=samples_info["LCL5"],
        LCL6_fq_list=samples_info["LCL6"],
        father_LCL7_yak="/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/common/yak/LCL7.yak",
        mother_LCL8_yak="/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/common/yak/LCL8.yak",
    output:
        dir_assm + "02_hifiasm_latest_trio/{assm}.hifiasm_latest_trio/{assm}.hifiasm_latest_trio.ok"
    priority: 9

    log:
        dir_assm + "logs/{assm}..hifiasm_latest_trio/{assm}.hifiasm_latest_trio.log"
    benchmark:
        dir_assm + "logs/{assm}.hifiasm_latest_trio/{assm}.hifiasm_latest_trio.tsv"
    threads: 120
    run:
        output_pre = str(output)[:-3]
        shell("{hifiasm} -o {output_pre} -t {threads}  -1 {input.father_LCL7_yak} -2 {input.mother_LCL8_yak} "
              "{input.LCL5_fq_list} {input.LCL6_fq_list} 2>{log} 1>{log}")
        touch(output.tag)
