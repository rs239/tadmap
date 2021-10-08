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

from .tadmap_base_config import *
import .tadmap_utils as U


def score_tads_by_alphas_poisson(tUMIs, tCapacity, lambda1, lambda2):
    #https://www.statlect.com/fundamentals-of-statistics/Poisson-distribution-maximum-likelihood
    #https://people.stat.sc.edu/Hitchcock/slides535day5spr2014.pdf
    
    f_ll = lambda L:   -tCapacity*L + np.log(L)*tUMIs
    
    dbg_print("Flag 616.10 ", tUMIs.shape, tCapacity.shape, lambda1, lambda2)
    
    l1, l2 = f_ll(lambda1), f_ll(lambda2)
    dbg_print("Flag 616.30 ", l1.shape, l2.shape, l1[:3,:3], l2[:3,:3])
    
    df_mle_prob = 1.0/ (np.exp(l2-l1) + 1)  # a/(a+b) = 1/ (np.exp(log(b)-log(a)) + 1)
    dbg_print("Flag 616.40 ", df_mle_prob.shape, df_mle_prob[:3,:3])
    
    return df_mle_prob




def map_sc_to_tad_poisson(adata, tad2genelist, skip_singleton_tads=False):
    dbg_print("Flag 874.10 ", adata.shape, len(tad2genelist))
    l1 = sorted([ (t, len(v)) for t,v in tad2genelist.items()], key=lambda s:int(s[0].split('|')[0]))
    
    if skip_singleton_tads:
        l0 = [(a,b) for a,b in l1 if b>1]
    else:
        l0 = [(a,b) for a,b in l1 if b>0]
        
    dbg_print("Flag 874.15 ", len(l0), l0[:3])
    
    tadCnt = pd.Series( [a[1] for a in l0], index=[a[0] for a in l0])
    nT = len(tadCnt)
    tadNames = list(tadCnt.index)
    dbg_print("Flag 874.20 ", nT, tadNames[:3], tadCnt.head(3))
    
    gene2tadidx = defaultdict(list)
    for i, (t,_) in enumerate(l0):
        for g in tad2genelist[t]:
            if g in gene2tadidx:
                dbg_print("Flag 874.22 INFO: gene {1} seems to span multiple tads. Already seen {0}, now seeing {2}".format([tadNames[a] for a in gene2tadidx[g]], g, t))
            gene2tadidx[g].append(i)

    dbg_print("Flag 874.25 ", len(gene2tadidx))
    
    adata_geneidx2tadidx = [[] for i in range(adata.shape[1])] #-1*np.ones(adata.shape[1])
    for i, g in enumerate(adata.var_names):
        adata_geneidx2tadidx[i] = gene2tadidx[g]

    adata_tadidx2geneidx = [[] for i in range(nT)]
    for i, tidxlist in enumerate(adata_geneidx2tadidx):
        for tidx in tidxlist:
            adata_tadidx2geneidx[tidx].append(i)
    
        
    dbg_print("Flag 874.30 ", len(adata_geneidx2tadidx), len(adata_tadidx2geneidx), adata_geneidx2tadidx[:5], adata_tadidx2geneidx[:5])
        
    if scipy.sparse.isspmatrix_csr(adata.X):
        X1 = adata.X.todense()
    else:
        X1 = adata.X
        
    df_tCapacity = np.tile( tadCnt.values, (adata.shape[0],1))

    df_tUMIs = np.full((adata.shape[0], nT),0.0)
    for tidx in range(nT):
        dbg_print("Flag 874.75 ", tidx, df_tCapacity.shape, df_tUMIs.shape, df_tUMIs.sum(), X1.shape)
        if len(adata_tadidx2geneidx[tidx]) == 0: continue
        gidxlist = adata_tadidx2geneidx[tidx]
        dbg_print("Flag 874.76 ", df_tUMIs[:,tidx].shape, X1[:, gidxlist].sum(axis=1).shape)
        df_tUMIs[:,tidx] += np.ravel(X1[:, gidxlist].sum(axis=1)) 
            
    dbg_print("Flag 874.80 ", df_tUMIs.shape, df_tCapacity.shape, df_tUMIs[:3,:3], df_tCapacity[:3,:3])
        
    return df_tUMIs, df_tCapacity, tadCnt, tadNames, gene2tadidx, adata_geneidx2tadidx





#https://people.stat.sc.edu/Hitchcock/slides535day5spr2014.pdf
def score_tads_by_EM_poisson(tUMIs, tCapacity, gammaDistHyperParams=None):
    
    #mean_umi_rate = np.expm1(np.mean(np.log1p(tUMIs/tCapacity)) # geometric mean instead of arithmetic mean
    mean_umi_rate = np.mean(tUMIs/tCapacity)
                             
    dbg_print("Flag 635.20 ", tUMIs.shape, tCapacity.shape, mean_umi_rate)
    
    lambda1 = 2*mean_umi_rate; lambda2 = 0.2*mean_umi_rate

    gammaA1, gammaA2, gammaB = 0.5*mean_umi_rate, 0, 0.5*mean_umi_rate #3, 0, 3

    dbg_print("Flag 635.22 ", gammaA1, gammaA2, gammaB)
    
    if gammaDistHyperParams is not None:
        gammaA1, gammaA2, gammaB = gammaDistHyperParams
    
    newlambda1 = newlambda2 = -1
    itercnt = 0; EPS=0.00001
    while itercnt < 100 and (abs(newlambda1-lambda1) > EPS  or abs(newlambda2-lambda2) > EPS):
        if itercnt>0:
            lambda1 = newlambda1
            lambda2 = newlambda2
            
        df_tad_probs = score_tads_by_alphas_poisson(tUMIs, tCapacity, lambda1, lambda2)

        #tU1, tU2 = (tUMIs + gammaA1), (tUMIs + gammaA2)
        tU1, tU2 = (tUMIs + gammaA1*tCapacity), (tUMIs + gammaA2*tCapacity)
        tC = tCapacity + gammaB
        
        newlambda1 = np.sum(np.ravel(df_tad_probs*tU1))/np.sum(np.ravel(df_tad_probs*tC))
        newlambda2 = np.sum(np.ravel(df_tad_probs*tU2))/np.sum(np.ravel(df_tad_probs*tC))


        itercnt += 1
        dbg_print("Flag 635.40 ",df_tad_probs.shape, lambda1, lambda2, newlambda1, newlambda2, EPS, itercnt)
        
    return lambda1, lambda2, df_tad_probs




def compute_tad_occupancy_by_EM_poisson( adata, tad2genelist, extra_args):

    if int(extra_args.get("adata_is_logcpm_normalized",0))  < 0.5:
        U.convert_to_counts(adata)
        
        #did the log transform to limit the variability of the data and the int conversion since we're assuming poisson
        adata.X = np.log1p(adata.X).astype(int).astype(float)

        
    dbg_print("Flag 935.20 ", adata.X.shape, adata.X[:7,:7].todense())
    
    gammaDistHyperParams = None
    if "alpha2_gammaprior_hyperparams" in extra_args:
        l = [float(a) for a in extra_args.get("alpha2_gammaprior_hyperparams").split(",")]
        gammaDistHyperParams = (l[0], l[1], l[2])
        dbg_print("Flag 935.625 ", gammaDistHyperParams)
        
    # convert each cell to TAD probs. save these probability vectors
    nsample = min(adata.shape[0], max(5000, int(0.1*adata.shape[0])))
    
    # adata_small = adata[np.random.choice(adata.shape[0], nsample, False), :]
    
    # print("Flag 935.63 ", adata_small.shape)
    # adata_small_tadUMIs, adata_small_tadCapacity = map_sc_to_tad_poisson(adata_small, tad2genelist)[:2]
    
    # print("Flag 935.64 ", adata_small_tadUMIs.shape, adata_small_tadCapacity.shape)
    
    adata_tadUMIs, adata_tadCapacity, tadCnt, tadNames, _, _ = map_sc_to_tad_poisson(adata, tad2genelist)    
    print("Flag 935.65 ", adata_tadUMIs.shape, adata_tadCapacity.shape)
    
    
    # lambda1, lambda2, _ = score_tads_by_EM_poisson(adata_small_tadUMIs, adata_small_tadCapacity, gammaDistHyperParams)
    
    lambda1, lambda2, _ = score_tads_by_EM_poisson(adata_tadUMIs, adata_tadCapacity, gammaDistHyperParams)
    
    print("Flag 935.67 ",  lambda1, lambda2)
    tad_occupancy_mtx = score_tads_by_alphas_poisson(adata_tadUMIs, adata_tadCapacity, lambda1, lambda2)

    return lambda1, lambda2, tadCnt, tadNames,  tad_occupancy_mtx    





def compute_tad_representation(sp, outdir, outsfx, tad_file, rnaseq_data_path, rnaseq_data_type, extra_args):
    dbg_print("Flag 934.01 ", sp, outdir, outsfx, tad_file, rnaseq_data_path, rnaseq_data_type, extra_args)
    
    # read tads. also get a list of valid genes
    geneset, tad2genelist = U.readTADs(tad_file)
    dbg_print("Flag 934.05 ", len(geneset), len(tad2genelist))

    
    # read rnaseq data to scanpy. supply list of valid genes to subset
    if rnaseq_data_type == 'trajectorama':
        adata = U.read_trajectorama_data(rnaseq_data_path, geneset)
    elif rnaseq_data_path == 'scvelo':
        adata = U.read_scvelo_data(rnaseq_data_path, geneset)
    else:
        assert 1==0

    dbg_print("Flag 934.20 ", adata.shape, adata.obs_names[:5])
    
    alpha1, alpha2, tadCnt, tadNames,  tad_occupancy_mtx = compute_tad_occupancy_by_EM_poisson( adata, tad2genelist, extra_args)

        
    print("Flag 934.67 ",  tad_occupancy_mtx.shape, alpha1, alpha2)
    tad_occupancy_df = pd.DataFrame(tad_occupancy_mtx, columns = tadNames)
    print("Flag 934.70 ", tad_occupancy_df.shape, alpha1, alpha2)
    
    
    # print the variation in TAD values
    mu_tad = tad_occupancy_df.mean(axis=0)
    z_tad = np.sqrt(tad_occupancy_df.shape[0])*tad_occupancy_df.mean(axis=0)/tad_occupancy_df.std(axis=0)

    n1 = tad_occupancy_df.shape[1]
    l = [['alpha1']+[alpha1]*n1, ['alpha2'] +[alpha2]*n1, ['tadmemb_count'] + list(tadCnt.values),  ['mu']+list(mu_tad.values), ['z']+list(z_tad.values)]
    tad_aux = pd.DataFrame(l, columns = ['name']+list(tad_occupancy_df.columns))

    print("Flag 934.80 ", tad_aux.shape)

    outfile1 = "{0}/tad_occupancy_{1}.h5".format(outdir, outsfx)
    outfile2 = "{0}/tad_auxinfo_{1}.csv".format(outdir, outsfx)
    tad_occupancy_df.to_hdf(outfile1, key="df", index=False)
    tad_aux.to_csv(outfile2, index=False)
    
    return tad_occupancy_df, tad_aux


    
        
############################################################################

if __name__ == "__main__":
    sys.path.append(os.path.join(sys.path[0],PROJ_DIR+'/src'))


    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", help="2-letter code of species (hs or mm)", type=str, default='hs')
    parser.add_argument("--outdir", help="output directory (can set to '.')", type=str, default=PROJ_DIR+"/data/processed/")
    parser.add_argument("--outsfx", help="suffix to use when producing output files")
    parser.add_argument("--tad_file", help="the path to the TAD file", type=str)
    parser.add_argument("--rnaseq_data_path", help="the path to the rnaseq_data_path", type=str)
    parser.add_argument("--rnaseq_data_type", help="one of trajectorama|cytotrace", type=str, choices=("trajectorama","cytotrace"), default="trajectorama")
    parser.add_argument("--extra", help="put this as the LAST option and arbitrary space-separated key=val pairs after that", type=str, nargs='*')


    args = parser.parse_args()
    extra_args = dict([a.split("=") for a in args.extra]) if args.extra else {}

    tad_occupancy_df, tad_aux = run(args.species, args.outdir, args.outsfx, args.tad_file, args.rnaseq_data_path, args.rnaseq_data_type, extra_args)

