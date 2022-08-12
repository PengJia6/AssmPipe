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
minimap2 = "/data/home/pengjia/mybin/minimap2"

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
        expand(dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.hicanu_default/{cohort}.{sample}.{hap}.HiFi.hicanu_default.ok",cohort="ChineseQuartet",
            sample="LCL56",hap=["hap1", "hap2"]),
        expand(dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.Hicanu_default_purge_dups/{cohort}.{sample}.{hap}.purge_dups.fa",cohort="ChineseQuartet",
            sample="LCL56",hap=["hap1", "hap2"]),


rule hicano:
    input:
        hifi=dir_assm_data + "{cohort}.{sample}.HiFi.{hap}.filter.fq.gz",
    # ont_h1=dir_assm_data + "{cohort}.{sample}.ONT.hap1.filter.fq.gz",
    # ont_h2=dir_assm_data + "{cohort}.{sample}.ONT.hap2.filter.fq.gz",
    # ont_ul=dir_assm_data + "{cohort}.{sample}.ONTUL.filter.fq.gz"
    output:
        tag=dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.hicanu_default/{cohort}.{sample}.{hap}.HiFi.hicanu_default.ok"
    log:
        dir_work_dir + "02_canu_HiFi/logs/{cohort}.{sample}.{hap}.hicanu_default.log"
    benchmark:
        dir_work_dir + "02_canu_HiFi/logs/{cohort}.{sample}.{hap}.hicanu_default.tsv"
    threads: 48
    run:
        assm_prefix = str(output).split("/")[-1].rstrip(".ok")
        assm_dir = "/".join(str(output).split("/")[:-1])
        # shell(" ulimit -Su 100000 && {canu} useGrid=false -haplotype "
        #       "-p {assm_prefix} -d {assm_dir} genomeSize=3.1g  -haplotypePaternal {input.LCL7_fq_list} "
        #       "-haplotypeMaternal {input.LCL8_fq_list} "
        #       "-pacbio-hifi {input.LCL5_fq_list} {input.LCL6_fq_list} 2>{log} 1>{log}")
        shell(" ulimit -Su 100000 && {canu} useGrid=false maxThreads={threads} "
              "-p {assm_prefix} -d {assm_dir} genomeSize=3.1g  "
              "-pacbio-hifi {input.hifi} 2>{log} 1>{log}")
        shell(" touch {output.tag}")

# minimap2 = "/data/home/pengjia/mybin/minimap2"
purge_dup_prefix = "/data/home/pengjia/mysoftware/assm/purge_dups/purge_dups/bin/"

rule purge_dups_minimap2:
    input:
        fa_tag=dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.hicanu_default/{cohort}.{sample}.{hap}.HiFi.hicanu_default.ok",
        hifi=dir_assm_data + "{cohort}.{sample}.HiFi.{hap}.filter.fq.gz",

    output:
        paf=dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.Hicanu_default_purge_dups/{cohort}.{sample}.{hap}.hifi.paf.gz",
    threads: 48
    run:
        fa = str(input.fa_tag)[:-3] + ".contigs.fasta"
        shell("{minimap2} -t {threads} -x map-hifi {fa} {input.hifi}|gzip -c - > {output.paf}")

rule cutoff:
    input:
        paf=dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.Hicanu_default_purge_dups/{cohort}.{sample}.{hap}.hifi.paf.gz",
    output:
        cutoff=dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.Hicanu_default_purge_dups/{cohort}.{sample}.{hap}.hifi.cutoff",
        basecov=dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.Hicanu_default_purge_dups/{cohort}.{sample}.{hap}.hifi.basecov",
    log: dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.Hicanu_default_purge_dups/{cohort}.{sample}.{hap}.cutoff.log",
    run:
        output_prefix = "/".join(str(output.cutoff).split("/")[:-1])
        shell("cd {output_prefix} && "
              "{purge_dup_prefix}pbcstat {input.paf} && "
              "{purge_dup_prefix}calcuts PB.stat >{output.cutoff} 2> {log}")
        shell("cd {output_prefix} && "
              "ln -s PB.base.cov {output.basecov}")

rule split_fa:
    input:
        fa_tag=dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.hicanu_default/{cohort}.{sample}.{hap}.HiFi.hicanu_default.ok",
    output:
        split_fa=dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.Hicanu_default_purge_dups/{cohort}.{sample}.{hap}.split.fa",
    run:
        fa = str(input.fa_tag)[:-3] + ".contigs.fasta"
        shell("{purge_dup_prefix}split_fa {fa} > {output.split_fa}")

rule selfalign:
    input:
        split_fa=dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.Hicanu_default_purge_dups/{cohort}.{sample}.{hap}.split.fa",
    output:
        paf=dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.Hicanu_default_purge_dups/{cohort}.{sample}.{hap}.split.self.paf.gz",
    threads: 48
    run:
        shell("{minimap2} -t {threads} -xasm5 -DP {input.split_fa} {input.split_fa} | gzip -c > {output.paf}")

rule purge_dups:
    input:
        cutoff=dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.Hicanu_default_purge_dups/{cohort}.{sample}.{hap}.hifi.cutoff",
        basecov=dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.Hicanu_default_purge_dups/{cohort}.{sample}.{hap}.hifi.basecov",
        paf=dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.Hicanu_default_purge_dups/{cohort}.{sample}.{hap}.split.self.paf.gz",
    output:
        dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.Hicanu_default_purge_dups/{cohort}.{sample}.{hap}.split.dups.bed",
    log: dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.Hicanu_default_purge_dups/{cohort}.{sample}.{hap}.purge_dups.log",
    run:
        shell("{purge_dup_prefix}purge_dups -2 -T {input.cutoff} -c  {input.basecov}  {input.paf} > {output} 2> {log}")

rule get_seq:
    input:
        fa_tag=dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.hicanu_default/{cohort}.{sample}.{hap}.HiFi.hicanu_default.ok",
        bed=dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.Hicanu_default_purge_dups/{cohort}.{sample}.{hap}.split.dups.bed",

    output:
        fa=dir_work_dir + "02_canu_HiFi/{cohort}.{sample}.{hap}.HiFi.Hicanu_default_purge_dups/{cohort}.{sample}.{hap}.purge_dups.fa",
    run:
        output_prefix = "/".join(str(output).split("/")[:-1])
        fa = str(input.fa_tag)[:-3] + ".contigs.fasta"
        shell("cd {output_prefix} && "
              "{purge_dup_prefix}get_seqs {input.bed} {fa} && "
              "ln -s purged.fa {output}")
