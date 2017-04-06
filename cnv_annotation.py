import sys
import pandas as pd
import math
import numpy as np
import pybedtools as pybed
import os


def annotated_gene( cnv_file ):
    '''
    find the genes covered by the cnvs.
    '''
    gene_file = "resources/refgene.bed"
    cnv = pybed.BedTool(cnv_file)
    gene = pybed.BedTool(gene_file)
    all = cnv.intersect(gene, wa=True, wb=True)
    return all


def get_gene_score( data ):
    '''
    generate the dictionary that contains the genes and their
    coresponding  scores
    '''
    gene_score = pd.read_csv(data, sep = "\t", usecols=[1,17,19])
    z_score = gene_score.set_index('gene')['mis_z'].to_dict()
    pLI = gene_score.set_index('gene')['pLI'].to_dict()
    return z_score, pLI

def get_max_score( list ):
    '''
    get the lagest score of the genes covered by the CNVs
    '''
    if len(list) == 0:
        return [0]
    else:
        array = np.array(list)
        max = np.amax(array)
        return max

def get_num_largePLI( list ):
    '''
    Count the number of genes covered by the cnvs with
    high pLI score (indicate the gene is important)
    '''
    i = 0
    for pLI in list:
        if pLI >= 0.9:
            i += 1
    return i

def get_num_largeZmis( list ):
    '''
    Count the number of genes covered by the cnvs with
    high pLI score (indicate the gene is important)
    '''
    j = 0
    for z_miss in list:
        if z_miss > 0:
            j += 1
    return j

def get_disease_gene( file ):
    '''
    generate the gene list and the disease gene dictionary
    '''
    disease_gene_association = {}
    df = pd.read_table( file, sep = "\t" )
    gene_all_list = df['AssociatedGenes'].tolist()
    disease_set = df['DiseaseName'].tolist()
    for i in range(len(gene_all_list)):
        disease_gene_association[gene_all_list[i]] = disease_set[i]
    return set(gene_all_list), disease_gene_association

def find_disease_genes( all_genes, test_gene_set ):
    '''
    get the disease genes covered by the CNVs and
    the disease associated.
    '''
    disease_gene_list = []
    for gene in test_gene_set:
        if gene in all_genes:
            disease_gene_list.append(gene)
    return disease_gene_list, len(disease_gene_list)

def get_CNV_gene( file ):
    '''
    merge the genes covered by the same cnvs
    '''
    site_dic = {}
    for line in open(file, 'r'):
        if line[0] == "#":
            continue
        data = line.strip().split("\t")
        keys = ":".join(data[:-4])
        if keys not in site_dic.keys():
            site_dic[keys] = [data[-1]]

        else:
            if data[-1] not in site_dic[keys]:
                site_dic[keys].append(data[-1])

    return site_dic


gene_zscore, gene_pLI = get_gene_score('resources/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt')

disease_genes, disease_gene_dict = get_disease_gene( 'resources/gene_condition_source_id'  )

# take the input file
input = sys.argv[1]

# generate the temp file with gene annotated
cnv_bed_gene = annotated_gene( input )
c = cnv_bed_gene.saveas('temp.bed', trackline='#cnv_anno_temp')
cnv_gene_dic = get_CNV_gene('temp.bed')

f = open('cnv_w_genes.txt', 'w')
for key, value in cnv_gene_dic.items():
    cnv_info = key.strip().split(":")
    num_gene = len(value)
    cnv_info.append(str(num_gene))
    f.write("\t".join(cnv_info) + "\t" + ",".join(value) + "\n")
f.close()

out_file = sys.argv[2]

header = ['Chr', 'Start', 'End', 'VarType', 'SampleID', 'Length', 'num_gene', 'genes',
              'mis_z_max', 'pli_max', 'num_largePLI', 'num_largeZmis',
            'gene_breaked_start', 'gene_breaked_end', 'start_break_gene_mis_z', 'start_break_gene_pli',
            'end_break_gene_mis_z', 'end_break_gene_pli', 'disease_genes', 'num_disease_genes']
f = open(out_file, 'w')
f.write("\t".join(header) + "\n")

#annotated the genes
for line in open('cnv_w_genes.txt', 'r'):
    z_scores = []
    pLIs = []
    data = line.strip().split("\t")
    genes = data[7].strip().split(",")
    cnv_disease_gene, num_disease_gene = find_disease_genes( disease_genes, genes )
    start_gene = genes[0]
    end_gene = genes[-1]
    if start_gene in gene_zscore.keys():
        start_z_score = gene_zscore[start_gene]
        start_pli = gene_pLI[start_gene]
    else:
        # replace the gene without the score find
        start_z_score = -99999.0
        start_pli = 0.0
    if end_gene in gene_zscore.keys():
        end_z_score = gene_zscore[end_gene]
        end_pli = gene_pLI[end_gene]
    else:
        end_z_score = -99999.0
        end_pli = 0.0
    for i in range(len(genes)):
        if genes[i] in gene_zscore.keys():
            z_scores.append(gene_zscore[genes[i]])
            pLIs.append(gene_pLI[genes[i]])
        else:
            z_scores.append(-99999.0)
            pLIs.append(0.0)
    z_score_max = get_max_score(z_scores)
    pLI_max = get_max_score(pLIs)
    num_largeP = get_num_largePLI(pLIs)
    num_largeZ = get_num_largeZmis(z_scores)
    if cnv_disease_gene != []:
        f.write( line.strip() + "\t"+ "\t".join([str(z_score_max), str(pLI_max), str(num_largeP), str(num_largeZ),
            start_gene, end_gene, str(start_z_score), str(start_pli),
            str(end_z_score), str(end_pli), ",".join(cnv_disease_gene),  str(num_disease_gene)])+"\n")
    elif cnv_disease_gene == []:
        f.write( line.strip()+"\t"+ "\t".join([str(z_score_max), str(pLI_max), str(num_largeP), str(num_largeZ),
            start_gene, end_gene, str(start_z_score), str(start_pli),
            str(end_z_score), str(end_pli), "NaN",  str(num_disease_gene)])+"\n")
    else:
        print("Warning: There is something wrong with the disease gene anotation")
f.close()

# remove temp files
os.remove('temp.bed')
os.remove('cnv_w_genes.txt')
