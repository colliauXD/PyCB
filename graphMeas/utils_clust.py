import numpy as np
from NeuroTools.signals import spikes as spk

def get_N_non_silent(fname,s=None,N=0):
   if s==None: s=spk.load_spikelist(fname)
   idl=s.id_list
   ids=idl[np.transpose(s.mean_rates(1))>0]
   print len(ids)
   if N==0:
      return s.id_slice(ids)
   elif len(ids)<N: return s.id_slice(ids)
   else: 
      np.random.shuffle(ids)
      return s.id_slice(ids[0:N])  

def get_all_pairs(s,NS):
   p=NS*(NS-1)/2
   p_idx=np.zeros([p,2])
   k=0
   for i in range(NS):
     for j in range(i):
        p_idx[k]=[s.id_list[i],s.id_list[j]]
        k+=1     
   NMAX=np.max(s.id_list)
   rev_idx=np.zeros(NMAX+1)
   for i in range(NS):
      rev_idx[s.id_list[i]]=i
   return spk.CustomPairs(s,s,p_idx),rev_idx

def get_cc_mat(s,NS,tb=25):
   p,rev_idx=get_all_pairs(s,NS)
   cc=s.pairwise_pearson_corrcoeff(NS*(NS-1)/2,p,time_bin=tb,all_coef=True)
   cc_mat=np.zeros([NS,NS])
   for i in range(len(cc)):
        cc_mat[rev_idx[p.pairs[i][0]],rev_idx[p.pairs[i][1]]]=cc[i]
        cc_mat[rev_idx[p.pairs[i][1]],rev_idx[p.pairs[i][0]]]=cc[i]
   return cc_mat 

def get_dk_mat(s,NS,tb=25):
   p,rev_idx=get_all_pairs(s,NS)

   cc=s.pairwise_pearson_corrcoeff(NS*(NS-1)/2,p,time_bin=tb,all_coef=True)

   cc_mat=np.zeros([NS,NS])
   for i in range(len(cc)):
        cc_mat[rev_idx[p.pairs[i][0]],rev_idx[p.pairs[i][1]]]=cc[i]
        cc_mat[rev_idx[p.pairs[i][1]],rev_idx[p.pairs[i][0]]]=cc[i]
   return cc_mat 
