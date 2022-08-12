falcon = "/data/home/pengjia/mybin/falcon"
dir_work = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/data/HiFi/"
dir_assm = "/data/home/pengjia/Project_Human_DATA/ChineseQuartet/Assembly/assm_res/"
import os

seqtk = "/data/home/pengjia/miniconda3/envs/assm/bin/seqtk"
falcon_env_python = "/data/home/pengjia/miniconda3/envs/falcon/bin"
falcon_src = "/data/home/pengjia/mysoftware/assm/falcon/falcon-pg1.6.3/py/scripts/pg_run.py"
ctg_falcon = "/data/home/pengjia/mysoftware/assm/Falcon/fc_run_human.cfg"
minimap2 = "/data/home/pengjia/mybin/minimap2"

rule all:
    input:
        expand(dir_assm + "02_falcon_{para}/{assm}.{hap}.falcon_{para}/{assm}.{hap}.falcon_{para}.ok",
            para=["default"],hap=["hap1", "hap2", ],assm="ChineseQuartet"),
        expand(dir_assm + "02_falcon_{para}/{assm}.{hap}.falcon_{para}/{assm}.{hap}.falcon_{para}/purge_dups/{assm}.{hap}.purge_dups.fa",
            para=["default"],hap=["hap1", "hap2", ],assm="ChineseQuartet")


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


purge_dup_prefix = "/data/home/pengjia/mysoftware/assm/purge_dups/purge_dups/bin/"


def get_fa(wildcards):
    return dir_assm + "02_falcon_{para}/{assm}.{hap}.falcon_{para}/2-asm-falcon/p_ctg.fa".format(para=wildcards.para,
        assm=wildcards.assm,hap=wildcards.hap)

    pass


rule purge_dups_minimap2:
    input:
        fa_tag=get_fa,
        hifi=get_fqs,
    output:
        paf=dir_assm + "02_falcon_{para}/{assm}.{hap}.falcon_{para}/{assm}.{hap}.falcon_{para}/purge_dups/{assm}.{hap}.paf.gz",
    threads: 48
    run:
        shell("{minimap2} -t {threads} -x map-hifi {input.fa_tag} {input.hifi}|gzip -c - > {output.paf}")

rule cutoff:
    input:
        paf=rules.purge_dups_minimap2.output.paf,
    output:
        cutoff=dir_assm + "02_falcon_{para}/{assm}.{hap}.falcon_{para}/{assm}.{hap}.falcon_{para}/purge_dups/{assm}.{hap}.hifi.cutoff",
        basecov=dir_assm + "02_falcon_{para}/{assm}.{hap}.falcon_{para}/{assm}.{hap}.falcon_{para}/purge_dups/{assm}.{hap}.hifi.basecov",
    log: dir_assm + "02_falcon_{para}/{assm}.{hap}.falcon_{para}/{assm}.{hap}.falcon_{para}/purge_dups/{assm}.{hap}.cutoff.log",
    run:
        output_prefix = "/".join(str(output.cutoff).split("/")[:-1])
        shell("cd {output_prefix} && "
              "{purge_dup_prefix}pbcstat {input.paf} && "
              "{purge_dup_prefix}calcuts PB.stat >{output.cutoff} 2> {log}")
        shell("cd {output_prefix} && "
              "ln -s PB.base.cov {output.basecov}")

rule split_fa:
    input:
        fa_tag=get_fa,

    output:
        split_fa=dir_assm + "02_falcon_{para}/{assm}.{hap}.falcon_{para}/{assm}.{hap}.falcon_{para}/purge_dups/{assm}.{hap}.split.fa",
    run:
        # fa = str(input.fa_tag)[:-3] + ".contigs.fasta"
        shell("{purge_dup_prefix}split_fa {input.fa_tag} > {output.split_fa}")

rule selfalign:
    input:
        split_fa=rules.split_fa.output.split_fa
    output:
        paf=dir_assm + "02_falcon_{para}/{assm}.{hap}.falcon_{para}/{assm}.{hap}.falcon_{para}/purge_dups/{assm}.{hap}.split.self.paf.gz",
    threads: 48
    run:
        shell("{minimap2} -t {threads} -xasm5 -DP {input.split_fa} {input.split_fa} | gzip -c > {output.paf}")

rule purge_dups:
    input:
        cutoff=rules.cutoff.output.cutoff,
        basecov=rules.cutoff.output.basecov,
        paf=rules.selfalign.output.paf,
    output:
        dir_assm + "02_falcon_{para}/{assm}.{hap}.falcon_{para}/{assm}.{hap}.falcon_{para}/purge_dups/{assm}.{hap}.split.dups.bed",
    log: dir_assm + "02_falcon_{para}/{assm}.{hap}.falcon_{para}/{assm}.{hap}.falcon_{para}/purge_dups/{assm}.{hap}.purge_dups.log",
    run:
        shell("{purge_dup_prefix}purge_dups -2 -T {input.cutoff} -c  {input.basecov}  {input.paf} > {output} 2> {log}")

rule get_seq:
    input:
        fa_tag=get_fa,
        bed=rules.purge_dups.output,

    output:
        fa=dir_assm + "02_falcon_{para}/{assm}.{hap}.falcon_{para}/{assm}.{hap}.falcon_{para}/purge_dups/{assm}.{hap}.purge_dups.fa",
    run:
        output_prefix = "/".join(str(output).split("/")[:-1])
        # fa = str(input.fa_tag)[:-3] + ".contigs.fasta"
        shell("cd {output_prefix} && "
              "{purge_dup_prefix}get_seqs {input.bed} {input.fa_tag} && "
              "ln -s purged.fa {output}")
