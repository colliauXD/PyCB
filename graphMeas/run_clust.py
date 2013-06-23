import numpy as np
import pylab as pl
import cPickle as cp
import utils_clust as uc
import graph_utils as gu
import measures as mea
from sklearn.cluster import AffinityPropagation
from NeuroTools.signals import spikes as spk
from NeuroTools.signals import PyNNNumpyBinaryFile

def aff_prop_scan(cc,s,pref_range=None,N=150,print_figs=False):
   """
   Scan various (homogeneous) preference values
   return preferences, #clusters  
   """
   print s.id_list
   if not pref_range: pref_range=[np.min(cc),0.1]
   ps=np.linspace(pref_range[0],pref_range[1],N)
   ls=np.zeros(N)
   af=AffinityPropagation(max_iter=100,copy=True)
   print "scanning preferences for affinity propagation"
   for i in range(N):
      f=af.fit(cc,ps[i])
      ls[i]=len(f.cluster_centers_indices_)
      print s.id_list
      if print_figs:
         plot_clust_ord(s,f.labels_,f.cluster_centers_indices_,fname="clust_%i.png"%i)
         pl.clf()
      print ps[i],"done:",ls[i],"centers"
   return ps,ls

def plot_scan(p,l,fig_name="scan_p.eps"):
    """
    plot the number of clusters 
    """
    pl.plot(p,l,"ko")
    pl.title(r"$N_{cluster}$")
    pl.xlabel("preference")
    pl.savefig(fig_name)
    return 0

def plot_clust_ord(s,labs,cents,fname="clust_ord.eps"):
   plot_idx=0
   size_Cs=np.zeros(len(cents))
   cols=np.random.random([len(cents),3])
   k=0
   for i in range(len(cents)):
      cl_idx=s.id_list[np.where(labs==i)[0]]
      size_Cs[i]=len(cl_idx)
      for j in cl_idx:
         st=s.spiketrains[j].spike_times
         pl.plot(st,plot_idx*np.ones(len(st)),",",color=cols[i])
         plot_idx+=1
   pl.xlabel("Time (ms)")
   pl.ylabel("Neurons ID")
   a = pl.axes([.65, .6, .2, .2], axisbg="0.5")
   pl.bar(range(len(size_Cs)),size_Cs,color=cols)
   pl.title('Cluster sizes')
   pl.setp(a, xticks=[], yticks=[0,np.max(size_Cs)])
   pl.savefig(fname)
   return 0

def plot_clust_ord2(s,labs,cents,fname="clust_ord.eps"):
   #Sort idx before plotting clusters
   plot_idx=0
   size_Cs=np.zeros(len(cents))
   cols=np.random.random([len(cents),3])
   k=0
   idx=s.sort_by('cell.mean_rate()')
   cs=np.ones(800)
   
   for i in range(len(cents)):
      cl_idx=s.id_list[np.where(labs==i)[0]]
      size_Cs[i]=len(cl_idx)
   for k in range(len(idx)):
      print labs[k]
      st=s.spiketrains[idx[k]].spike_times
      pl.plot(st,plot_idx*np.ones(len(st)),",",color=cols[labs[k]])
      plot_idx+=1
   pl.xlabel("Time (ms)")
   pl.ylabel("Neurons ID")
   a = pl.axes([.65, .6, .2, .2], axisbg="0.5")
   pl.bar(range(len(size_Cs)),size_Cs,color=cols)
   pl.title('Cluster sizes')
   pl.setp(a, xticks=[], yticks=[0,np.max(size_Cs)])
   pl.savefig(fname)
   return 0
   
def get_clust(params,fname,i,j):
   NS=300
   f=PyNNNumpyBinaryFile(fname)
   sp=f.read_spikes({'t_start':0,'t_stop': 7200})
   s=uc.get_N_non_silent(fname,sp.time_slice(500,1500),NS)   
   cc=uc.get_cc_mat(s,len(s.id_list),tb=2)
   cc+=np.random.normal(0,0.001*np.abs(np.min(cc-1)),len(cc))
   pref=1.3*np.min(cc-1)
   af=AffinityPropagation(max_iter=100,copy=True).fit(cc-1,pref)
   cents=af.cluster_centers_indices_
   labs=af.labels_
   plot_clust_ord2(s,labs,cents,fname="coldiag/clust_%s_%s2.png"%(i,j))
   

def diag_clusts():
   params=cp.load(open("coldiag/col_v_ab/col_v_ab.par"))
   for i in range(10):
      for j in range(10):
         fname="coldiag/col_v_ab/10x10/E_Bloc_%s_%s.gdf"%(i,j)
         get_clust(params,fname,i,j) 
         pl.clf()

def plot_distrib_gath():
   params=cp.load(open("coldiag/col_v_ab/col_v_ab.par"))
   for i in range(10):
      fname="coldiag/col_v_ab/10x10/E_Bloc_%s_%s.gdf"%(0,i)
      f=PyNNNumpyBinaryFile(fname)
      sp=f.read_spikes({'t_start':0,'t_stop': 7200})
      s=sp.time_slice(1000,2000) 
      r=np.array(s.mean_rates()) 
      r.sort()
      pl.plot(r,color=[1-i/15.,0.,0.])
   pl.savefig("0i_sortrates.png")   

def plot_distrib():
   params=cp.load(open("coldiag/col_v_ab/col_v_ab.par"))
   for i in range(10):
      for j in range(10):
         fname="coldiag/col_v_ab/10x10/E_Bloc_%s_%s.gdf"%(i,j)
         f=PyNNNumpyBinaryFile(fname)
         sp=f.read_spikes({'t_start':0,'t_stop': 7200})
         s=sp.time_slice(1000,2000) 
         r=np.array(s.mean_rates())
         #a,b=np.histogram(r,1000)
         r.sort() 
         pl.plot(r)
         #pl.plot(b[:-1],a)
         pl.savefig("hist_rates%s_%s.png"%(i,j))
         pl.clf()

"""
i=9
j=0
fname="Results/col_v_ab/11x11/E_Bloc_%s_%s.gdf"%(i,j)
params=cp.load(open("Results/col_v_ab/col_v_ab.par"))
f=PyNNNumpyBinaryFile(fname)
sp=f.read_spikes({'t_start': 7000,'t_stop': 8000})

s_idx=sp.sort_by('cell.mean_rate()')
for i in range(len(s_idx)): pl.plot(sp.spiketrains[s_idx[i]].spike_times,i*np.ones(len(sp.spiketrains[s_idx[i]])),",")
"""
#ns.raster_plot()
#pl.savefig("ras_ab.png")
diag_clusts()
#plot_distrib_gath()

"""
NS=500
f=PyNNNumpyBinaryFile("Results/col_v_ab/11x11/E_Bloc_0_0.gdf")
sp=f.read_spikes({'t_start':7000,'t_stop': 7200})
s=uc.get_N_non_silent("Results/col_v_ab/11x11/E_Bloc_0_0.gdf",sp,NS)
cc=uc.get_cc_mat(s,NS,tb=5)

path="Results/col_v_ab"
s=mea.get_spikes(params["diag"]["N"], 2000, 8000, path,"E")
time_r=mea.get_t_rates(s, 5)
mr,stdr=mea.get_rates_m_std(s,5)
cv=mea.get_cvs(s)
cc=mea.pw_pearson_corcoef(s, 100, 15)

np.save(path+"_time_r",time_r)
np.save(path+"_mr",mr)
np.save(path+"_stdr",stdr)
np.save(path+"_cv",cv)
np.save(path+"_cc",cc)

gu.plot_mr(mr, params,path)
gu.plot_stdr(stdr, params,path)
gu.plot_cc(cc, params,path)
gu.plot_cv(cv, params,path)
"""
#get_clust(params,fname) 
#pl.clf()


#p,l=aff_prop_scan(cc-1,s,[3.*np.min(cc-1),0.01],N=50,print_figs=True)
#plot_scan(p,l)
#print cc

"""
fname="col_ab.spk"
spk=mea.get_spikes(10,0,2000,fname,"E")
for i in range(99):
   s=uc.get_N_non_silent(fname,spk[i],NS)

   NS=len(s.id_list)
   print NS
   if (NS>100):
      cc=uc.get_cc_mat(s,NS,tb=1)
      cc+=np.random.normal(0,0.01*np.abs(np.min(cc-1)),len(cc))

      #p,l=aff_prop_scan(cc-1,[1.3*np.min(cc-1),0.01],N=150)
      #plot_scan(p,l,"scan_trans.eps")
      #pl.clf()
      pref=1.3*np.min(cc-1)
      af=AffinityPropagation(max_iter=100,copy=True).fit(cc-1,pref)
      cents=af.cluster_centers_indices_
      labs=af.labels_
      plot_clust_ord(s,labs,fname="clust_ord%s.png"%i)
      pl.clf()
"""

