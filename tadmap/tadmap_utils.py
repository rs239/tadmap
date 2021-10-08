#!/usr/bin/env python

import pandas as pd
import numpy as np
import scipy, sklearn, os, sys, string, fileinput, glob, re, math, itertools, functools, copy, multiprocessing, traceback, tarfile, gzip, csv
import scipy.stats, sklearn.decomposition, sklearn.preprocessing, sklearn.covariance
from scipy.stats import describe
from scipy import sparse
import os.path
import scipy.sparse
from scipy.sparse import csr_matrix, csc_matrix
from sklearn.preprocessing import normalize
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns

from base_config import *


def read_Ensembl_v102_refdata(sp, do_prot_coding):
    assert sp in ['hs','mm']
    assert do_prot_coding #only implemented for protein_coding genes. non-PC genes was deprecated early on

    if sp=='mm':
        df = pd.read_csv("~/work/vartad/data/raw2/genes/mm10_ens102_refdata.tsv", delimiter='\t')
        df = df[~df['MGI symbol'].isnull()].reset_index(drop=True)
        df = df[df['Chromosome/scaffold name'].apply(lambda s: len(s)<=2 and s!="MT")].reset_index(drop=True)
        df["gene1"] = np.where(~df["Gene name"].isnull(), df["Gene name"], df["Gene stable ID"])
        df1 = pd.DataFrame({ "gene": df['gene1'] + "|Ensembl:" + df['Gene stable ID'] + "|MGI:" + df['MGI symbol'],
                             "chr": df['Chromosome/scaffold name'].apply(lambda s: "chr"+s if len(s)<=2 and s!="MT" else s.lower()),
                             "strand": df["Strand"].apply(lambda a: "+" if int(a)>0 else "-"),
                             "txstart": df['Gene start (bp)'], "txend": df['Gene end (bp)']})
        return df1

    if sp=='hs':
        df = pd.read_csv("~/work/vartad/data/raw2/genes/hg38_ens102_refdata.tsv", delimiter='\t')
        df = df[~df['HGNC symbol'].isnull()].reset_index(drop=True)
        df = df[df['Chromosome/scaffold name'].apply(lambda s: len(s)<=2 and s!="MT")].reset_index(drop=True)
        df["gene1"] = np.where(~df["Gene name"].isnull(), df["Gene name"], df["Gene stable ID"])
        df1 = pd.DataFrame({ "gene": df['Gene name'] + "|Ensembl:" + df['Gene stable ID'] + "|HGNC:" + df['HGNC symbol'],
                             "chr": df['Chromosome/scaffold name'].apply(lambda s: "chr"+s if len(s)<=2 and s!="MT" else s.lower()),
                             "strand": df["Strand"].apply(lambda a: "+" if int(a)>0 else "-"),
                             "txstart": df['Gene start (bp)'], "txend": df['Gene end (bp)']})

        return df1
        


def do_scanpy_processing_and_filtering(adata, sample_frac=1.0, do_log1p=True):
    import scanpy as sc

    if sample_frac > 0.999:
        sample_n = int(adata.shape[0]* float(sample_frac))
        sample_idxs = np.random.choice(adata.shape[0], sample_n, replace=False)
        adata2 = adata[sample_idxs].copy()
    else:
        adata2 = adata.copy()
    dbg_print("Flag 123.10 ", adata.shape, adata2.shape, sample_frac)
    
    cinfo, _ = sc.pp.calculate_qc_metrics(adata2, inplace=False)
    c2depth_thresh = sorted(cinfo["total_counts"])[max(10,int(0.001*adata2.shape[0]))]
    c2gcnt_thresh = sorted(cinfo["n_genes_by_counts"])[max(10,int(0.001*adata2.shape[0]))]
    dbg_print("Flag 123.20 ", c2depth_thresh, c2gcnt_thresh)
    
    sc.pp.filter_cells(adata2, min_counts = c2depth_thresh+1)
    sc.pp.filter_cells(adata2, min_genes = c2gcnt_thresh+1)
    sc.pp.normalize_total(adata2, target_sum=1e6)
    if do_log1p:
        sc.pp.log1p(adata2)
    dbg_print("Flag 123.40 ", adata2.shape)
    return adata2


def convert_to_counts(adata):
    import scanpy as sc
    max1 = np.max(np.max(adata.X))
    if max1 < 15: # likely log values
        adata.X = np.expm1(adata.X)
        sc.pp.normalize_total(adata, target_sum=1e6)
        
    else:
        sc.pp.normalize_total(adata, target_sum=1e6)


        
def parse_geneids(s):
    l = s.split("|")
    sym = l[0]
    try:
        hgncid = [a.split(':')[2] for a in l if a.startswith("HGNC:")][0]
    except:
        hgncid = sym
    try:
        ensemblid = [a.split(':')[1] for a in l if a.startswith("Ensembl:")][0]
    except:
        ensemblid = sym
    return (sym, ensemblid, hgncid)


def readTADs(tad_file):
    df = pd.read_csv(tad_file)
    dbg_print("Flag 652.01 ", df.shape)
    df = df[~df["genelist"].isnull()]
    
    df["gene_count"] = df["genelist"].apply(lambda v: len(v.split(';')))
    dbg_print("Flag 652.02 ", df.shape, df['gene_count'].describe())

    
    geneset = set()
    tad2genelist = {}
    for i in range(df.shape[0]):
        tad = df["tad"].iat[i]
        genelist = [parse_geneids(s) for s in df['genelist'].iat[i].split(';')]
        
        #dbg_print("Flag 652.10 ", tad, len(genelist), genelist[0])        
        tad2genelist[tad] = [a[0] for a in genelist]
        geneset.update(genelist)

    return geneset, tad2genelist



def match_adata_to_geneset(adata, geneset, adata_genes=None):
    if adata_genes is not None:
        genes = adata_genes
    else:
        genes = list(adata.var_names)
        
    symbol_matchcnt = len( set(genes) & set(a[0] for a in geneset))
    ensembl_matchcnt = len( set(genes) & set(a[1] for a in geneset))
    print ("Flag 342.24 ", symbol_matchcnt, ensembl_matchcnt, len(geneset), len(genes))
    
    if ensembl_matchcnt > symbol_matchcnt:
        ensembl2gene = {a[1]:a[0] for a in geneset}
        genes1 = copy.deepcopy(genes)
        genes = [ensembl2gene.get(g,g) for g in genes1]

    geneset_symbols = set(a[0] for a in geneset)
    valid_genes = np.array([(g in geneset_symbols) for g in adata.var_names])

    print ("Flag 342.29 ", adata.shape)
    return adata[:, valid_genes].copy()



def read_trajectorama_data(rnaseq_data_path, geneset):
    for pat in ["/*tab.npz", "*tab.npz", "*.npz"]:
        l = glob.glob('{0}{1}'.format(rnaseq_data_path, pat))
        print ("Flag 342.10 ", l)
        if l:
            npz_file = l[0]
            break
    else:
        raise "Could not find tab file in "+rnaseq_data_path

    
    for pat in ["/*genes.txt", "*genes.txt"]:
        l = glob.glob('{0}{1}'.format(rnaseq_data_path, pat))
        print ("Flag 342.15 ", l)
        if l:
            genes_file = l[0]
            break
    else:
        raise "Could not find genes.txt file in "+rnaseq_data_path
    

        
    X = scipy.sparse.load_npz(npz_file)
    genes = np.array(open(genes_file,'r').read().rstrip().split('\n'))

    from anndata import AnnData
    adata = AnnData(X)
    print ("Flag 342.20 ", adata.shape, len(genes))
    adata.var_names = genes
    adata.var_names_make_unique() # because of the matching code below, all but one of the dupes will be dropped

    adata = match_adata_to_geneset(adata, geneset, genes)
    print ("Flag 342.30 ", adata.shape)
        
    return adata


def read_scvelo_data(rnaseq_data_path, geneset):
    from anndata import AnnData
    import scanpy as sc
    
    adata = sc.read(rnaseq_data_path)
    print ("Flag 342.20 ", adata.shape, len(geneset))
    adata.var_names = [c.upper() for c in adata.var_names]
    adata.var_names_make_unique() # because of the matching code below, all but one of the dupes will be dropped

    adata = match_adata_to_geneset(adata, geneset)
    
    print ("Flag 342.30 ", adata.shape)
    if 'cluster' not in adata.obs.columns:
        if 'clusters_coarse' in adata.obs.columns:
            adata.obs['cluster'] = adata.obs['clusters_coarse']
            
        if 'clusters' in adata.obs.columns:
            adata.obs['cluster'] = adata.obs['clusters']
            
        
    return adata
    
