# -*- coding: utf-8 -*-
import numpy as np
from NeuroTools.signals import *
from pyNN.recording import files
from scipy import optimize

def get_spikes(N=10, t_start=500, t_stop = 2000, path='Result',pop=""):
  layer = []
  path += "/%dx%d" %(N, N)
  print "Loading data for all columns..."
  for i in xrange(N):
    for j in xrange(N):
      file = PyNNNumpyBinaryFile("%s/%s_Bloc_%d_%d.gdf" %(path,pop, i, j))
      layer.append(file.read_spikes(params={'t_start':t_start,'t_stop':t_stop}))   
      file.fileobj.close()
  print "Shifting all the id..."
  for idx in xrange(len(layer)):
    layer[idx].id_offset(N*idx)
  return layer

def get_t_rates(spikes, time_bin=10):
  rates = numpy.array([numpy.mean(s.firing_rate(time_bin,display=False),axis=0) for s in spikes])
  N, L   = rates.shape
  rates  = rates.reshape((numpy.sqrt(N),numpy.sqrt(N), L))
  return rates

def get_rates_m_std(spikes, time_bin=10):
  m_rates = numpy.array([numpy.mean(numpy.mean(s.firing_rate(time_bin,display=False),axis=1)) for s in spikes])
  std_rates = numpy.array([numpy.std(numpy.mean(s.firing_rate(time_bin,display=False),axis=1)) for s in spikes])
  N,    = m_rates.shape
  m_rates  = m_rates.reshape((numpy.sqrt(N),numpy.sqrt(N)))
  std_rates  = std_rates.reshape((numpy.sqrt(N),numpy.sqrt(N)))
  return m_rates, std_rates

def get_cvs(spks):
  spkf    = [s.id_slice(s.select_ids("len(cell) > 0")) for s in spks]
  cvs=numpy.array([s.cv_isi(True).mean() for s in spkf])
  N=len(spkf)
  cvs  = cvs.reshape((numpy.sqrt(N), numpy.sqrt(N)))
  return cvs

def pw_pearson_corcoef(spks, nb_pairs=100, time_bin=5):
   N2      = int(len(spks))
   N       = int(numpy.sqrt(len(spks)))
   print "Filtering only the non silent cells...."
   spkf    = [s.id_slice(s.select_ids("len(cell) > 0")) for s in spks]
   m_cc=np.zeros(N2)
   
   for i in xrange(N2):
     if len(spkf[i])<np.sqrt(nb_pairs):
        m_cc[i]=0
     else:
       m_cc[i]=spkf[i].pairwise_pearson_corrcoeff(nb_pairs, RandomPairs(spkf[i], spkf[i], True, True),5)[0]
       """
        cells = spkf[i].id_list
        Nc     = len(cells)
        id1    = cells[numpy.random.uniform(0, nb_pairs)]
        id2    = cells[numpy.random.uniform(0, nb_pairs)]
        doublons=(numpy.array(id1)==numpy.array(id2))
        print cells
        while (len(id2[doublons])>0):
          id2[doublons]=cells[numpy.random.uniform(0, len(id2[doublons]))]
          doublons=(numpy.array(id1)==numpy.array(id2))
        timehist1 = spkf[id1].time_histogram(time_bin)
        timehist2 = spkf[id2].time_histogram(time_bin)
        m_cc[i] = numpy.corrcoef(timehist1,timehist2)[0,1]
       """
   return m_cc

def get_MeanVm(N=30,dt=0.1, t_start=500, t_stop = 2000, path='Results'):
  layer = []
  print "Loading data for all columns..."
  for i in xrange(N):
    for j in xrange(N):
      layer.append(AnalogSignal(pylab.load("%s/Bloc_%d_%d_exc.lfp" %(path, i, j))[-1+int(t_start/dt):int(t_stop/dt)], dt=dt,t_start=t_start, t_stop=t_stop))
  return layer


