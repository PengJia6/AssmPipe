falcon = "/data/home/pengjia/mybin/falcon"
dir_work = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/data/HiFi/"
dir_assm = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/assm_res/"
import os
seqtk = "/data/home/pengjia/miniconda3/envs/assm/bin/seqtk"
falcon_env_python = "/data/home/pengjia/miniconda3/envs/falcon/bin"
falcon_src = "/data/home/pengjia/mysoftware/assm/falcon/falcon-pg1.6.3/py/scripts/pg_run.py"
ctg_falcon = "/data/home/pengjia/mysoftware/assm/Falcon/fc_run_human.cfg"

rule all:
    input:
        expand(dir_assm + "02_falcon_{para}/{assm}.{hap}.falcon_{para}/{assm}.{hap}.falcon_{para}.ok",para=[
            "default"],
            hap=["collapsed", ],assm="ChineseQuartet")


# expand(   dir_assm + "}",
#     hap=["hap1", "hap2"]),
# expand(dir_assm + "ChineseQuartet/falcon_r4/{hap}/ChineseQuartet.HiFi.falcon_r4.{hap}", hap=["hap1", "hap2"])


def get_fqs(wildcards):
    fqs = []
    if wildcards.hap in ["hap1", "hap2"]:
        fqs_tag = expand(dir_work + "{sample}/{sample}.T2T.minimap2.HiFi.{{hap}}.fa",sample=["LCL5",
                                                                                             "LCL6"],hap=wildcards.hap)
        fqs_untag = expand(dir_work + "{sample}/{sample}.T2T.minimap2.HiFi.untag.{{hap}}.fa",sample=["LCL5",
                                                                                                     "LCL6"],hap=wildcards.hap)
        fqs.extend(fqs_untag)
        fqs.extend(fqs_tag)
    elif wildcards.hap == "collapsed":
        for sample in ["LCL5", "LCL6"]:
            files = os.listdir("/data/DATA/ChineseQuartet/RAWDATA/HiFi/fastqs/{}".format(sample))
            files = [i[:-5] + "fa" for i in files if i.endswith("fq.gz")]

            for i in files:
                fqs.append("/data/DATA/ChineseQuartet/RAWDATA/HiFi/fastqs/{s}/{i}".format(s=sample,i=i))
    return fqs


rule falcon:
    input:
        fqs=get_fqs
    output:
        dir_assm + "02_falcon_{para}/{assm}.{hap}.falcon_{para}/{assm}.{hap}.falcon_{para}.ok"
    threads: 48
    log:
        dir_assm + "logs/02_falcon_{para}/{assm}.{hap}.falcon_{para}.log"
    benchmark:
        dir_assm + "logs/02_falcon_{para}/{assm}.{hap}.falcon_{para}.tsv"
    run:
        input_file = "/".join(str(output).split("/")[:-1]) + "/input.fofn"
        input_dir = "/".join(str(output).split("/")[:-1])
        output_file = str(output)[:-3]
        if str(wildcards.para) == "default":
            dir_input = "/".join(input.fqs[0].split("/")[:-2])
            dir_output = "/".join(str(output).split("/")[:-1])
            shell("mkdir -p {dir_output}")
            file = open(str(input_file),"w")
            for i in input.fqs:
                file.write(i + "\n")
            file.close()

        shell("export PATH={falcon_env_python}:$PATH && "
              "export PATH=/data/home/pengjia/mybin:$PATH && "
              "cd {input_dir} && "
              "/usr/bin/cp {ctg_falcon} ./ && "
              "fc_run {input_dir}/fc_run_human.cfg 2>{log} 1>{log}")
        shell("touch {output}")


rule fq2fa:
    input:
        "{prefix}.fq.gz"
    output:
        "{prefix}.fa"
    run:
        shell("{seqtk} seq -A {input} > {output}")

# input_file = str(output)[:-2] + "fofn"
# docker_input_file = "/output/" + input_file.split("/")[-1]
# docker_output_file = "/output/" + str(output).split("/")[-1][:-2]
# if str(wildcards.para) == "default":
#     dir_input = "/".join(input.fqs[0].split("/")[:-2])
#     dir_output = "/".join(str(output).split("/")[:-1])
#     shell("mkdir -p {dir_output}")
#     file = open(str(input_file),"w")
#     for i in input.fqs:
#         file.write("/input/" + "/".join(i.split("/")[-2:]) + "\n")
#     file.close()
#     shell("docker run  -v {dir_input}:/input -v {dir_output}:/output "
#           "cschin/falcon:1.6.3 asm {docker_input_file} 24 24 24 24 24 24 24 24 24 "
#           "--with-consensus --shimmer-r 3 --best_n_ovlp 8 "
#           "--output {docker_output_file} <<< 'yes' ")
# shell("touch {output}")
#
