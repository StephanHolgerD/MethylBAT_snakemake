import pandas as pd
import os
from glob import glob

annotation_file = '/mnt/raid/users/pilgram/methylome/analysis/Homo_sapiens.GRCh38.108.chr.chr.gff3'
annotation_file_db = '/mnt/raid/users/pilgram/methylome/analysis/Homo_sapiens.GRCh38.108.chr.chr.db'
minAlignments = [10,20]

bed = 'data/all_Gene_and_CpG_ID_annotations_covered_targets_Twist_Methylome_V1_TE-96341190_hg38.bed'
df=pd.read_csv("samples.tsv",sep="\t",dtype=object)

SAMPLES=list(df["ID"])
SAMPLES_PATH=list(df["Bedgraph"])


case = list(df[df['type']=='case']['ID'])
control = list(df[df['type']=='control']['ID'])

casestr=','.join(case)
controlstr=','.join(control)

print(casestr)
wkd= os.getcwd()
outdir = '/'.join(wkd.split('/')[:-1])

def check_symlink(file1, file2):
    try:
        os.symlink(file1, file2)
    except FileExistsError:
        print("Link for " + file1 + " is already present in 01_raw")

for a,b in zip(SAMPLES,SAMPLES_PATH):
    a=str(a)
    os.makedirs("../01_bedgraphs/"+a, exist_ok=True)
    check_symlink(b,f'../01_bedgraphs/{a}/{a}.bedgraph')



rule all:
    input:
        expand('../02_FilteredBedgraphs/{sample}/{sample}_min{filter}.bedgraph',sample=SAMPLES,filter=minAlignments),
        expand('{outdir}/03_metilene_call_{filter}/dmr_metilene_qval.0.05.bed',filter=minAlignments,outdir=outdir),
        expand('{outdir}/03_bat_summarize_{filter}/BAT_summarize_summary_case_control.bedgraph.gz',filter=minAlignments,outdir=outdir)
        #f"{outdir}/03_bat/out.txt"
        #'../03_ConcatFilteredBedgraphs/ConcatFilteredBedgraphs.bed',
        #expand('../04_FilterSplitUnionBed/{sample}/{sample}.bed',sample=SAMPLES)


rule FilterBedgraphs:
    input:
        bedgraph = "../01_bedgraphs/{sample}/{sample}.bedgraph"
    threads: 1
    conda: 'envs/pybedtools.yaml'
    output:
        bedgraph_filtered = "../02_FilteredBedgraphs/{sample}/{sample}_min{filter}.bedgraph"
    shell:
        "awk '$5+$6>={wildcards.filter}' {input.bedgraph} | cut -f1-4 | bedtools sort  -header > {output.bedgraph_filtered}"





rule BAT_rule_summrarize:
    input:
          cases=expand('../02_FilteredBedgraphs/{sample}/{sample}_min{{filter}}.bedgraph',sample=case),
          controls=expand('../02_FilteredBedgraphs/{sample}/{sample}_min{{filter}}.bedgraph',sample=control)

    threads: 1
    container: "docker://christianbioinf/bat:latest"
    output:
        bat_sum = "{outdir}/03_bat_summarize_{filter}/BAT_summarize_metilene_case_control.txt",
        bat_sum_bg = "{outdir}/03_bat_summarize_{filter}/BAT_summarize_summary_case_control.bedgraph"

    shell:
        "cases=$(echo {input.cases}|sed 's/ /,/g');\
         controls=$(echo {input.controls}|sed 's/ /,/g');\
         BAT_summarize --in1 $cases --in2 $controls --groups case,control --h1 {casestr} --h2 {controlstr} --out {outdir}/03_bat_summarize_{wildcards.filter}/BAT_summarize --cs hg38_chromSize.txt"



rule BAT_rule_DMRcalling:
    input:
          bat_sum = "{outdir}/03_bat_summarize_{filter}/BAT_summarize_metilene_case_control.txt"
    threads: 1
    container: "docker://christianbioinf/bat:latest"
    output:
        metilene_call = "{outdir}/03_metilene_call_{filter}/dmr_metilene_qval.0.05.bed"
    shell:
        "BAT_DMRcalling -q {input.bat_sum} -o {outdir}/03_metilene_call_{wildcards.filter}/dmr_metilene -a case -b control"

rule bgzip_tabix:
    input:
        bat_sum_bg = "{outdir}/03_bat_summarize_{filter}/BAT_summarize_summary_case_control.bedgraph"
    threads: 1
    conda: 'envs/bedtools.yaml'
    output:
        bat_sum_bg_bgzip = "{outdir}/03_bat_summarize_{filter}/BAT_summarize_summary_case_control.bedgraph.gz"
    shell:
        "bgzip {input.bat_sum_bg}; tabix -p bed {output.bat_sum_bg_bgzip}"