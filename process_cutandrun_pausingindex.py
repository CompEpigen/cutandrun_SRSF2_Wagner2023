import os
import numpy as np
import pandas as pd
import deeptools.countReadsPerBin as crpb
import matplotlib
import matplotlib.pyplot as plt
import pysam
import scipy.stats
from statsmodels.stats.multitest import multipletests
from pathlib import Path



working_dir = "/omics/groups/OE0219/internal/Etienne/Projects/SRSF2/cutandrun/analysis/2022-11-09"
Path(working_dir).mkdir(parents=True, exist_ok=True)
windowsize_TSS = 500
chromosomes = ["chr"+str(x) for x in range(1,23)] + ["chrX","chrY"]

# Select downregulated genes (or upregulated, or background)
chromosomes = [str(x) for x in range(1,23)] + ["X","Y"]
#df = pd.read_csv("/omics/groups/OE0219/internal/Etienne/Projects/SRSF2/cutandrun/analysis/final_48kd_ctr_annotation.txt",sep="\t")
df = pd.read_csv("/omics/groups/OE0219/internal/Etienne/Projects/SRSF2/cutandrun/analysis/res_TcRC_60min_kd_ctr.txt",sep=",")
df.columns = ["gene_name","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"]
df = df.loc[~df["pvalue"].isnull(),]
df = df.loc[df["baseMean"]>5,:]
df["padj"] = multipletests(df["pvalue"],alpha=0.05,method="fdr_bh")[1]

df_downregulated = df.loc[ (df["padj"]<0.05) & (df["log2FoldChange"]<0),:] 
genes_selected_downregulated = list(df_downregulated["gene_name"])

df_upregulated = df.loc[ (df["padj"]<0.05) & (df["log2FoldChange"]>0),:]
genes_selected_upregulated = list(df_upregulated["gene_name"])

df_background = df.loc[ (df["baseMean"]>5),:] # (df["padj"]>0.20) & 
genes_selected_background = list(df_background["gene_name"])

# Find canonical transcripts for the selected genes
#transcripts_ids = set()
canonicalTranscripts2gene_downregulated = {}
canonicalTranscripts2gene_upregulated = {}
canonicalTranscripts2gene_background = {}
df_GRCh38 = pd.read_csv("/omics/groups/OE0219/internal/Etienne/data/reference/RNA/Transcripts_GRCh38.tsv",sep="\t")
for i in range(df_GRCh38.shape[0]):
    if df_GRCh38.loc[i,"Ensembl Canonical"] == 1.0:
        if df_GRCh38.loc[i,"Gene name"] in genes_selected_downregulated and not df_GRCh38.loc[i,"Gene name"] in canonicalTranscripts2gene_downregulated.values():
            canonicalTranscripts2gene_downregulated[df_GRCh38.loc[i,"Transcript stable ID"]] = df_GRCh38.loc[i,"Gene name"]
        if df_GRCh38.loc[i,"Gene name"] in genes_selected_upregulated and not df_GRCh38.loc[i,"Gene name"] in canonicalTranscripts2gene_upregulated.values():
            canonicalTranscripts2gene_upregulated[df_GRCh38.loc[i,"Transcript stable ID"]] = df_GRCh38.loc[i,"Gene name"]
        if df_GRCh38.loc[i,"Gene name"] in genes_selected_background and not df_GRCh38.loc[i,"Gene name"] in canonicalTranscripts2gene_background.values():
            canonicalTranscripts2gene_background[df_GRCh38.loc[i,"Transcript stable ID"]] = df_GRCh38.loc[i,"Gene name"]

df_GRCh37 = pd.read_csv("/omics/groups/OE0219/internal/Etienne/data/reference/RNA/Transcripts_GRCh37.tsv",sep="\t")
df_GRCh37.sort_values(by=["Chromosome/scaffold name","Transcript start (bp)"],inplace=True)
df_GRCh37.reset_index(inplace=True)

with open(os.path.join(working_dir,"selected_genes_downregulated.bed"),"w") as outfile:
    for i in range(df_GRCh37.shape[0]):
        if df_GRCh37.loc[i,"Transcript stable ID"] in canonicalTranscripts2gene_downregulated and df_GRCh37.loc[i,"Chromosome/scaffold name"] in chromosomes:
            strand = "+" if df_GRCh37.loc[i,"Strand"] >0 else "-"
            tmp = outfile.write("\t".join(["chr"+str(df_GRCh37.loc[i,"Chromosome/scaffold name"]),str(df_GRCh37.loc[i,"Transcript start (bp)"]),str(df_GRCh37.loc[i,"Transcript end (bp)"]),canonicalTranscripts2gene_downregulated[df_GRCh37.loc[i,"Transcript stable ID"]],"0",strand])+"\n")

with open(os.path.join(working_dir,"selected_genes_upregulated.bed"),"w") as outfile:
    for i in range(df_GRCh37.shape[0]):
        if df_GRCh37.loc[i,"Transcript stable ID"] in canonicalTranscripts2gene_upregulated and df_GRCh37.loc[i,"Chromosome/scaffold name"] in chromosomes:
            strand = "+" if df_GRCh37.loc[i,"Strand"] >0 else "-"
            tmp = outfile.write("\t".join(["chr"+str(df_GRCh37.loc[i,"Chromosome/scaffold name"]),str(df_GRCh37.loc[i,"Transcript start (bp)"]),str(df_GRCh37.loc[i,"Transcript end (bp)"]),canonicalTranscripts2gene_upregulated[df_GRCh37.loc[i,"Transcript stable ID"]],"0",strand])+"\n")

with open(os.path.join(working_dir,"selected_genes_background.bed"),"w") as outfile:
    for i in range(df_GRCh37.shape[0]):
        if df_GRCh37.loc[i,"Transcript stable ID"] in canonicalTranscripts2gene_background and df_GRCh37.loc[i,"Chromosome/scaffold name"] in chromosomes:
            strand = "+" if df_GRCh37.loc[i,"Strand"] >0 else "-"
            tmp = outfile.write("\t".join(["chr"+str(df_GRCh37.loc[i,"Chromosome/scaffold name"]),str(df_GRCh37.loc[i,"Transcript start (bp)"]),str(df_GRCh37.loc[i,"Transcript end (bp)"]),canonicalTranscripts2gene_background[df_GRCh37.loc[i,"Transcript stable ID"]],"0",strand])+"\n")




# Get dataframe of pausing indices transcripts x sample
antibody = "phospho-Rpb1-ser5"
replicates = ["1","2","3"]
pipeline_dir = "/omics/groups/OE0219/internal/Etienne/Projects/SRSF2/cutandrun/pipeline2_RPKM/"

#antibody = "RPB1"
#replicates = ["3","4","5"]
#pipeline_dir = "/omics/groups/OE0219/internal/Etienne/Projects/SRSF2/cutandrun/pipeline_RPKM/"

sample2bam={}
sample2pysam={}
sample2counter={}
for condition in ["KD48h","control"]:
    for replicate in replicates:
        sample = antibody+"_"+condition+"_R"+replicate
        sample2bam[sample] = pipeline_dir+"02_alignment/bowtie2/target/"+condition+"-"+antibody+"_R"+replicate+".target.markdup.bam"
        sample2pysam[sample] = pysam.AlignmentFile(sample2bam[sample])
        sample2counter[sample] = crpb.CountReadsPerBin([sample2bam[sample]], binLength=100, stepSize=100)

chromosomes = ["chr"+str(x) for x in range(1,23)] + ["chrX","chrY"]

def compute_pausing_indices(genes_file,outfile):
    # We define the pausing index as the number of reads aligned to the TSS (+/- 250bp) divided by the number of reads aligned to the gene body. 
    d={"gene":[]}
    for sample in sample2bam:
        d[sample] = []
    with open(genes_file,"r") as infile:
        for line in infile:
            linesplit= line.rstrip("\n").split("\t")
            chr = linesplit[0]
            if not chr in chromosomes: continue
            start = int(linesplit[1])
            end = int(linesplit[2])
            gene = linesplit[3]
            strand = linesplit[5]
            if strand=="+":
                TSS = start
            else:
                TSS = end
            d["gene"].append(gene)
            for condition in ["KD48h","control"]:
                for replicate in replicates:
                    sample =  antibody+"_"+condition+"_R"+replicate
                    coverage_TSS = 1+sample2counter[sample].get_coverage_of_region(sample2pysam[sample],chr, [(TSS-windowsize_TSS//2, TSS+windowsize_TSS//2)])[0]
                    coverage_genebody = 1+sample2counter[sample].get_coverage_of_region(sample2pysam[sample],chr, [(start -windowsize_TSS//2 , end +windowsize_TSS//2)])[0]
                    ratio = coverage_TSS / coverage_genebody
                    d[sample].append(ratio)
    df_pausing_indices = pd.DataFrame(d)
    df_pausing_indices.to_csv(outfile,sep="\t",index=False)

compute_pausing_indices(os.path.join(working_dir,"selected_genes_downregulated.bed"),os.path.join(working_dir,antibody+"_pausing-indices_downregulated.tsv"))
compute_pausing_indices(os.path.join(working_dir,"selected_genes_upregulated.bed"),os.path.join(working_dir,antibody+"_pausing-indices_upregulated.tsv"))
compute_pausing_indices(os.path.join(working_dir,"selected_genes_background.bed"),os.path.join(working_dir,antibody+"_pausing-indices_background.tsv"))



#######################################################
# Compute pvalues for differences in pausing index between the two conditions

for type in ["downregulated","upregulated","background"]:
    print(type)
    df_pausing_indices = pd.read_csv(os.path.join(working_dir,antibody+"_pausing-indices_"+type+".tsv"),sep="\t",index_col = "gene")
    conditions = ["KD48h","control"]
    

    # For each gene, perform a t-test between the pausing index of 2 conditions (requires several replicates per condition).
    # Then, combine the pvalues of all the genes into a single pvalue

    # Alternatively, for each gene, compute the log ratio of the pausing index between the two conditions, and then perform a t-test for all genes
    pvalues = []
    log_ratios_conditions=[]
    diff_ratios_conditions=[]

    for gene in df_pausing_indices.index:
        # Collect the pausing indices for this gene, for each of the 2 conditions
        pausing_indices={}
        for condition in conditions:
            pausing_indices[condition]=[]
            for replicate in replicates: 
                pausing_indices[condition].append(df_pausing_indices.loc[gene,antibody+"_"+condition+"_R"+str(replicate)])
        # t-test between the two conditions
        pvalue = scipy.stats.ttest_ind(pausing_indices[conditions[0]],pausing_indices[conditions[1]])[1] # ,alternative="greater"
        if pvalue == pvalue:
            pvalues.append(pvalue)

        # log ratio between the 2 conditions
        log_ratios_conditions.append(np.log(np.mean(pausing_indices[conditions[0]]) / np.mean(pausing_indices[conditions[1]])))
        diff_ratios_conditions.append(np.mean(pausing_indices[conditions[0]]) - np.mean(pausing_indices[conditions[1]]))

    ratios_condition = {}
    for condition in conditions:
        ratios_condition[condition] = []
        for replicate in replicates:
            ratios_replicate = []
            for gene in df_pausing_indices.index:
                ratios_replicate.append(df_pausing_indices.loc[gene,antibody+"_"+condition+"_R"+str(replicate)])
            ratios_condition[condition].append(np.mean(ratios_replicate))

    pvalue_samples = scipy.stats.ttest_ind(ratios_condition[conditions[0]],ratios_condition[conditions[1]])[1]
    print(ratios_condition)
    print("p value:" +str(pvalue_samples))

    # use Fisher's method to combine the pvalues.
    print(scipy.stats.combine_pvalues(pvalues))
    pvalue_combined = scipy.stats.combine_pvalues(pvalues)[1]
    print("pvalue combined: "+str(pvalue_combined))


    # t-test between the two conditions
    pvalue2 = scipy.stats.ttest_1samp(diff_ratios_conditions,0)[1] # alternative="greater"
    print("pvalue diff: "+str(pvalue2))
    pvalue3 = scipy.stats.ttest_1samp(log_ratios_conditions,0)[1] # alternative="greater"
    print("pvalue log ratio: "+str(pvalue3))

    matplotlib.rcParams.update({'font.size': 20})
    #plt.hist(log_ratios_conditions,bins=np.arange(-1.6,1.6,0.1),density=False)
    plt.hist(diff_ratios_conditions,bins=np.arange(-0.2,0.2,0.01),density=False)
    plt.xlabel("Difference of pausing index between KD and control")
    plt.ylabel("Number of transcripts")
    plt.title("pvalue: "+ "{:.2E}".format(pvalue2))
    plt.savefig(os.path.join(working_dir,antibody+"_ttest_"+type+".svg"),bbox_inches="tight")
    plt.cla() 
    plt.clf() 
    plt.close('all')
    series_diff_ratio = pd.Series(diff_ratios_conditions,index=df_pausing_indices.index)
    series_diff_ratio.to_csv(os.path.join(working_dir,antibody+"_diff_pausing-indices_"+type+".tsv"),sep="\t")

    matplotlib.rcParams.update({'font.size': 20})
    #plt.hist(log_ratios_conditions,bins=np.arange(-1.6,1.6,0.1),density=False)
    plt.hist(log_ratios_conditions,bins=np.arange(-0.2,0.2,0.01),density=False)
    plt.xlabel("Log ratio of pausing index between KD and control")
    plt.ylabel("Number of transcripts")
    plt.title("pvalue: "+ "{:.2E}".format(pvalue3))
    plt.savefig(os.path.join(working_dir,antibody+"_ttestLogRatio_"+type+".svg"),bbox_inches="tight")
    plt.cla() 
    plt.clf() 
    plt.close('all')


    plt.hist(pvalues)
    plt.xlabel("pvalue")
    plt.ylabel("Number of transcripts")
    plt.title("Combined pvalue:" + "{:.2E}".format(pvalue_combined))
    plt.savefig(os.path.join(working_dir,antibody+"_pvalues_"+type+".svg"),bbox_inches="tight")
    plt.cla() 
    plt.clf() 
    plt.close('all')