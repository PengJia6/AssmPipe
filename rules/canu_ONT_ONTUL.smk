# ======================================================================================================================
# Project: Project_Human_iharbor
# Script : 02_canu_HiFi_ONT_ONTUL.smk TODO check
# Author : Peng Jia
# Date   : 2021.05.31
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
# ======================================================================================================================
canu = "/data/home/pengjia/miniconda3/envs/assm/bin/canu"
# dir_raw_data = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/rawdata/HiFi/fofn/"
dir_work_dir = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/hap_assm/"
dir_assm_data = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/hap_assm/01_assm_data/"

# samples = ["LCL5", "LCL6", "LCL7", "LCL8"]
# samples_info = {}
# for sample in samples:
#     info = open(dir_raw_data + "{}.HiFi.fofn".format(sample), "r").readlines()
#     thisinfo = {i.rstrip(".bam\n").split("/")[-1]: i.rstrip("\n") for i in info}
#     samples_info[sample] = thisinfo
# print(samples_info)
#
# rule all:
#     input:
#          dir_work_dir + "CHT/CHT_girl/CHT_girl.canu_hifi.ok",
#          expand(dir_work_dir + "CHT/CHT_girl_hap{hap}/CHT_girl_hap{hap}.canu_hifi.ok", hap=["Maternal", "Paternal"])
# # expand(dir_work_dir + "CHT/girl_LCL5_LCL6/{hap}_asm/{hap}.fq.gz", hap=["hap1", "hap2"]),
# # expand(dir_work_dir + "CHT/girl_LCL5_LCL6/{hap}_asm/asm_res/{hap}", hap=["hap1", "hap2"])
#
rule all:
    input:
        expand(dir_work_dir + "02_canu_ONT_ONTUL/{cohort}.{sample}.{hap}.ONT_ONTUL.canu_default/{cohort}.{sample}.{hap}.ONT_ONTUL.canu_default.ok",cohort="ChineseQuartet",
            sample="LCL56",hap=["hap1", "hap2"])

rule cano:
    input:
        # ont=dir_assm_data + "{cohort}.{sample}.ONT.filter.fq.gz",
        ont=dir_assm_data + "{cohort}.{sample}.{hap}.ONT.long.fq.gz",
        ont_ul=dir_assm_data + "{cohort}.LCL5.{hap}.ONTUL.long.fq.gz"
    output:
        tag=dir_work_dir + "02_canu_ONT_ONTUL/{cohort}.{sample}.{hap}.ONT_ONTUL.canu_default/{cohort}.{sample}.{hap}.ONT_ONTUL.canu_default.ok"
    log:
        dir_work_dir + "02_canu_ONT_ONTUL/logs/{cohort}.{sample}.{hap}.ONT_ONTUL.canu_default.log"
    benchmark:
        dir_work_dir + "02_canu_ONT_ONTUL/logs/{cohort}.{hap}.{sample}.ONT_ONTUL.canu_default.tsv"
    run:
        assm_prefix = str(output).split("/")[-1].rstrip(".ok")
        assm_dir = "/".join(str(output).split("/")[:-1])
        # shell(" ulimit -Su 100000 && {canu} useGrid=false -haplotype "
        #       "-p {assm_prefix} -d {assm_dir} genomeSize=3.1g  -haplotypePaternal {input.LCL7_fq_list} "
        #       "-haplotypeMaternal {input.LCL8_fq_list} "
        #       "-pacbio-hifi {input.LCL5_fq_list} {input.LCL6_fq_list} 2>{log} 1>{log}")
        shell(" ulimit -Su 100000 && {canu} useGrid=false "
              "-p {assm_prefix} -d {assm_dir} genomeSize=3.1g  "
              "-nanopore {input.ont} {input.ont_ul} 2>{log} 1>{log}")
        shell(" touch {output.tag}")

#
# minimap2 = "/data/home/pengjia/mybin/minimap2"
# rule purge_dups_minimap2:
#     input:
#          fa_tag=dir_work_dir + "02_canu_ONT_ONTUL/{cohort}.{sample}.ONT_ONTUL.ok",
#          ont=dir_assm_data + "{cohort}.{sample}.ONT.filter.fq.gz",
#          ont_ul=dir_assm_data + "{cohort}.{sample}.ONTUL.filter.fq.gz"
#     output:
#           pdf=dir_work_dir + "02_canu_ONT_ONTUL/purge_dups/{cohort}.{sample}.ONT_ONTUL.pdf.gz",
#     threads: 48
#     run:
#         fa = str(input.fa_tag).rstrip(".ok") + ".contigs.fasta"
#         shell("{minimap2} -t {threads} -x map-ont {fa} {input.ont} {input.ont_ul}|gzip -c - > {output.paf}")


# rule cano_trio:
#     input:
#          LCL5_fq_list=expand(dir_raw_data + "fastqs/LCL5/LCL5_{subsample}.fq.gz",
#                              subsample=samples_info["LCL5"]),
#          LCL6_fq_list=expand(dir_raw_data + "fastqs/LCL6/LCL6_{subsample}.fq.gz",
#                              subsample=samples_info["LCL6"]),
#          LCL7_fq_list=expand(dir_raw_data + "fastqs/LCL7/LCL7_{subsample}.fq.gz",
#                              subsample=samples_info["LCL7"]),
#          LCL8_fq_list=expand(dir_raw_data + "fastqs/LCL8/LCL8_{subsample}.fq.gz",
#                              subsample=samples_info["LCL8"]),
#     output:
#           tag=dir_work_dir + "CHT/CHT_girl/CHT_girl.canu_hifi.ok"
#     log:
#        dir_work_dir + "logs/CHT_girl.canu_hifi.logs"
#     benchmark:
#              dir_work_dir + "logs/CHT_girl.canu_hifi.tsv"
#     threads: 288
#     run:
#         assm_prefix = str(output).split("/")[-1].rstrip(".ok")
#         assm_dir = "/".join(str(output).split("/")[:-1])
#         shell(" ulimit -Su 100000 && {02_canu_HiFi_ONT_ONTUL} useGrid=false -haplotype "
#               "-p {assm_prefix} -d {assm_dir} genomeSize=3.1g  -haplotypePaternal {input.LCL7_fq_list} "
#               "-haplotypeMaternal {input.LCL8_fq_list} "
#               "-pacbio-hifi {input.LCL5_fq_list} {input.LCL6_fq_list} 2>{log} 1>{log}")
#         shell(" touch {output.tag}")
#
# rule canu_hap:
#     input:
#          fq=dir_work_dir + "CHT/CHT_girl/CHT_girl.canu_hifi.ok"
#     output:
#           tag=dir_work_dir + "CHT/CHT_girl_hap{hap}/CHT_girl_hap{hap}.canu_hifi.ok"
#     log:
#        dir_work_dir + "logs/CHT_girl_hap{hap}.canu_hifi.logs"
#     benchmark:
#              dir_work_dir + "logs/CHT_girl_hap{hap}.canu_hifi.tsv"
#     threads: 60
#     run:
#         assm_prefix = str(output).split("/")[-1].rstrip(".ok")
#         assm_dir = "/".join(str(output).split("/")[:-1])
#         fa_this_hap = "/".join(str(input).split("/")[:-1]) + "/haplotype/haplotype-" + wildcards.hap + ".fasta.gz"
#         fa_unknown = "/".join(str(input).split("/")[:-1]) + "/haplotype/haplotype-unknown.fasta.gz"
#         shell(" ulimit -Su 100000 && {02_canu_HiFi_ONT_ONTUL} useGrid=false "
#               "-p {assm_prefix} -d {assm_dir} -untrimmed genomeSize=3.1g  "
#               "-pacbio-hifi {fa_unknown} {fa_this_hap} 2>{log} 1>{log}")
#         shell(" touch {output.tag}")
# rule purge_dups:
#     input:
#          ref="",
#          fq_list=""
#     output:
#           tag=""
#     shell:
#          """
#
#          """

#
# def get_one_sample_fq(wildcards):
#     return samples_info[wildcards.sample][wildcards.subsample]
#
#
# rule bam2fastqgz2:
#     input:
#          get_one_sample_fq
#     output:
#           dir_raw_data + "raw_data/fastqs/{sample}/{sample}_{subsample}.fq.gz"
#     threads: 2
#     run:
#         shell("{samtools} fastq -0 {output} -@ {threads} {input}")
