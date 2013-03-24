# -*- coding: utf-8 -*-
from NeuroTools.stgen import StGen
from pyNN.nest import *
from pyNN.random import NumpyRNG
from pyNN.recording import files
import os,time, gc
import cPickle
from pyNN.nest.simulator import state
state.optimize = True
from myconnectors2 import *
import utils_pinwheel as up
import shortrangeGen as srg
import numpy as np
nest.Install("mymodule")

global timestep, myseed
timestep = 0.1
rng      = NumpyRNG(seed=249856, parallel_safe=False)
seed     = 512468


def invert_rf_fields(connections,layer):
    res = {}
    for neuron in connections:
      idx = neuron[0]
      if up.is_in_range(idx,layer.N,layer.N):
         for target in neuron[1:]:
       	    if up.is_in_range(target[1],layer.N,layer.N):
               if abs(target[0]) > 0.001:
	          src = target[1]
 	          if res.has_key(src):
	             res[src].append((target[0], idx))
	          else:
	             res[src] = [(target[0], idx)]
    new_list = []
    for key, value in res.items():
      new_list.append((key, value))
    return new_list


class Column(object):  

  def __init__(self, cell_params, col_params):
    global seed    
    self.exc = Population(col_params['EI_ratio']*col_params['size'],native_cell_type("aeif_beuler"), cell_params)
    self.inh = Population(col_params['size'], native_cell_type("aeif_beuler"), cell_params)
    self.inh.set('b', 0.0)
    uniformDistr = RandomDistribution('uniform', [cell_params['V_reset'],cell_params['V_th']], rng=NumpyRNG(seed))
    seed += 1
    self.all.initialize('v',uniformDistr)
    self.w_E = col_params['wE']
    self.w_I = col_params['wI']
  
  @property
  def all(self):
      return self.exc + self.inh

  def _record(self, variable='spikes', to_file=True):
    if variable is 'spikes':
      self.all.record(to_file=to_file)
    if variable is 'spikes_E':
      self.exc.record(to_file=to_file)
    if variable is 'spikes_I':
      self.inh.record(to_file=to_file)
    elif variable is 'lfp':
      if num_processes() > 2:
        raise Exception("Mean voltmeters can't be used in a distributed context .... Sorry for that")
      res = nest.Create('multimeter', params={'record_from': ["V_m"],
                                              'withgid'    : False, 
                                              'withtime'   : False, 
                                              #'to_file'    : True, 
                                              #'to_memory'  : True, 
                                              'interval'   : 1.,
                                              'to_accumulator': True
                                              })
      nest.DivergentConnect(res, self.exc.all_cells.tolist(), model='static_synapse')
      self.mean_recorder = [res]
    elif variable is 'vm':
      self.exc.record_v(to_file=to_file)
    elif variable is 'g':
      self.exc.record_gsyn(to_file=to_file)

  def _print(self, variable, file, gather=True):
    if variable is 'spikes':
      self.all.printSpikes(files.NumpyBinaryFile(file, 'w'), gather=gather)
    if variable is 'spikes_E':
      self.exc.printSpikes(files.NumpyBinaryFile(file, 'w'), gather=gather)
    if variable is 'spikes_I':
      self.inh.printSpikes(files.NumpyBinaryFile(file, 'w'), gather=gather)
    elif variable is 'lfp':
      #vm_name = nest.GetStatus(self.mean_recorder)[0]['filename']
      #os.system("mv %s %s" %(vm_name, file))
      lfp = nest.GetStatus(self.mean_recorder)[0]["events"]["V_m"]
      numpy.save(file,lfp/len(self.exc))
    if variable is 'vm':
      self.exc.print_v(files.NumpyBinaryFile(file, 'w'), gather=gather)      
    if variable is 'g':
      self.exc.print_gsyn(files.NumpyBinaryFile(file, 'w'), gather=gather)


class LGN(object):
  def __init__(self, params,N_layer):
    self.N   = params["N"]
    self.LGN = Population((self.N,self.N), SpikeSourceArray, {'spike_times':[]})
    self.wE=params["wE"]
    self.wI=params["wI"]
    self.w_thresh=params["w_thresh"]
    srg.export_RF(self.N,N_layer)  

  def build_static_bar(self, f0=100,stim_start=0, stim_duration=4000, phi=0, save=None):
      if rank()==0: print "Creating or loading the LGN input: static bar at %s rad...."%phi
      sx=.7*self.N
      sy=0.1*self.N
      n0=(self.N-1)/2.
      spk  = StGen(seed=5423689)
      st   = []
      for i in range(self.N):
         for j in range(self.N):
           bar_fr=f0*np.exp(-((i-n0)*np.cos(phi)+(j-n0)*np.sin(phi))**2/sx**2-(-(i-n0)*np.sin(phi)+(j-n0)*np.cos(phi))**2/sy**2)
           spk_times = spk.poisson_generator(bar_fr, 0, stim_duration, array=True)
	   st.append(spk_times)
      for cell,spikes in zip(self.LGN, st):
	cell.spike_times = spikes+stim_start
      if save:
	cPickle.dump(st, open(save,'w'))

  def build_moving_bar(self, stim_duration=4000, omega=5, save=None):
      if rank()==0: print "Creating or loading the LGN input ...."
      f0   = 100
      time = numpy.arange(0, stim_duration, 1)
      spk  = StGen(seed=5423689)
      st   = []
      for i in range(self.N):
         for j in range(self.N):
	   signal    = f0/(1+numpy.exp(-100*(numpy.cos(omega*.001*time+numpy.pi*i/(self.N*1.))-.95)))
           spk_times = spk.inh_poisson_generator(signal, time, stim_duration, array=True)
	   st.append(spk_times)
      for cell,spikes in zip(self.LGN, st):
	cell.spike_times = spikes
      if save:
	cPickle.dump(st, open(save,'w'))

  def _connect_cell2pop(self, source_cell, targets, p_conn, layer,factor_inh):
     # !! ON/OFF regions of the RF can be modelled in various ways !!
     #Here we suppose LGN cells from ON (OFF) excite (inhibit) both excitatory and inhibitory cells
     tei=np.where(p_conn>0)[0]
     target_exc_idx=[targets[tei[i]] for i in range(len(tei))]
     tii=np.where(p_conn<0)[0]
     target_inh_idx=[targets[tii[i]] for i in range(len(tii))]

     d      = timestep
     n_e, n_i=0, 0
     lgn_cell   = self.LGN[source_cell]

     
     if len(target_exc_idx)>0: 
        #p_conn[tei]=1.
        p_exc = Projection(lgn_cell, layer.get_exc_from_blocks(target_exc_idx), MyConnector(np.repeat(p_conn[tei], layer.Ne), weights=self.wE, delays=d), target='excitatory', rng=rng)
        print p_exc.size(), "exc. connections from LGN"
     if len(target_inh_idx)>0: 
        #p_conn[tii]=1
        p_inh = Projection(lgn_cell, layer.get_exc_from_blocks(target_inh_idx), MyConnector(np.repeat(-factor_inh*p_conn[tii], layer.Ne), weights=self.wI, delays=d), target='inhibitory', rng=rng)
        print p_inh.size(), "inh. connections from LGN"
     return 0


  def project(self, layer, factor_inh=1, cl_val=[-1.,1.]):
    if rank()==0: print "Building the connections from the LGN ...."
    LGNOn  = cPickle.load(open('RFs.dat','r'))
    LGNOn  = invert_rf_fields(LGNOn,layer)
    cPickle.dump(LGNOn,open('invertedRF.dat','w'))  
    n_conn = 0
    for neuron in LGNOn:
      idx = neuron[0]
      in_range_targets = []
      in_range_weights=[]
      if up.is_in_range(idx,layer.N,layer.N):
        for target in neuron[1]:
          if (up.is_in_range(target[1],layer.N,layer.N) and np.abs(target[0])>self.w_thresh):
             in_range_targets.append(target[1])
             in_range_weights.append(target[0])
      self._connect_cell2pop(idx, in_range_targets, np.array(in_range_weights).clip(cl_val[0],cl_val[1]), layer,factor_inh)
    print "--> connections were established from LGN to network..." 

class Layer(object):

  def __init__(self, N, cell_params, col_params, le=1., li=1., self_conn=1., cross_conn=1., velocity=None):
    print "Creation of the layer ..."
    self.velocity = velocity
    if velocity is None:
      self.max_delay = timestep
    else:
      self.max_delay = sqrt(2)*N/self.velocity
    if rank()==0: print "Max delay is set to", self.max_delay
    setup(timestep, max_delay=self.max_delay)
    self.N          = N
    self.blocs      = numpy.array([Column(cell_params, col_params) for i in xrange(N**2)])
    self.blocs      = self.blocs.reshape(N,N)
    self.recordings = []
    self.self_conn  = self_conn
    self.cross_conn = cross_conn
    self.le         = le
    self.li         = li
    self.Ne         = len(self.blocs[0,0].exc)
    self.Ni         = len(self.blocs[0,0].inh)


  def __getitem__(self, i,j):
    return self.blocs[i,j]

  def _indices(self):
    res = []
    for i in xrange(self.N):
        for j in xrange(self.N):
            res += [(i, j)]
    return res

  def get_all_from_blocks(self, indices):
    res = self.blocs[indices[0]].all
    for i in xrange(1, len(indices)):
        res += self.blocs[indices[i]].all
    return res

  def get_exc_from_blocks(self, indices):
    res = self.blocs[indices[0]].exc
    for i in xrange(1, len(indices)):
        res += self.blocs[indices[i]].exc
    return res

  def get_inh_from_blocks(self, indices):
    res = self.blocs[indices[0]].inh
    for i in xrange(1, len(indices)):
        res += self.blocs[indices[i]].inh
    return res

  def mpi_print(self, string):
    if simulator.state.mpi_rank == 0: print string  

  def stimulate(self, stimulus):
    if stimulus['type'] == 'Poisson_E':
      print stimulus["rate_E"], stimulus["N_E"], stimulus["duration"],stimulus['ee']
      self.ext_stim = Population(stimulus["N_E"], native_cell_type('poisson_generator'), {'rate': stimulus["rate_E"], 'stop': float(stimulus["duration"])})
      N_target=self.Ne*self.N**2
      target_cells=self.get_exc_from_blocks(self._indices())
      print target_cells.size, N_target
      conn    = MyConnector(numpy.repeat(stimulus['ee'],N_target), weights = stimulus['wE'], delays=self.max_delay)
      prj     = Projection(self.ext_stim, target_cells, conn, target='excitatory', rng=rng)
      self.mpi_print("--> %d excitatory connections were established from external sources ..." %prj.size())

  def record(self, variable='spikes', vline=False, RefSpk=[], to_file=True):
      self.mpi_print("Setting up the recorders ...")
      if variable is 'spikes':
        if len(RefSpk)==0:
          for bloc in self.blocs.flatten():
            bloc._record(variable, to_file)
        else:
            self.blocs[tuple(RefSpk)]._record(variable, to_file)
      if variable is 'spikes_E':
        if len(RefSpk)==0:
          for bloc in self.blocs.flatten():
            bloc._record(variable, to_file)
        else:
            self.blocs[tuple(RefSpk)]._record(variable, to_file)
      if variable is 'spikes_I':
        if len(RefSpk)==0:
          for bloc in self.blocs.flatten():
            bloc._record(variable, to_file)
        else:
            self.blocs[tuple(RefSpk)]._record(variable, to_file)
      elif variable is 'lfp':
        middle = int(self.N)/2
        RBlocs = self.blocs.flatten()
        if vline==True:
            RBlocs=self.blocs[middle]
        for bloc in RBlocs:
            #for bloc in self.blocs.flatten():
            bloc._record(variable)
      elif variable is 'vm':
        middle = int(self.N**2)/2
        self.blocs.flatten()[middle]._record(variable, to_file)
      elif variable is 'g':
        middle = int(self.N**2)/2
        self.blocs.flatten()[middle]._record(variable, to_file)
      self.recordings.append(variable)

  def simulate(self, simtime, print_rates=True):
    if rank()==0: 
         print "Launching simulation ..."
         nest.SetStatus([0],{'print_time':True})
    os.system('rm -rf g.gdf') 
    run(simtime)
    if print_rates:
       firing_rate = 0
       for bloc in self.blocs.flatten():
	 firing_rate += bloc.exc.meanSpikeCount()
       firing_rate = (firing_rate/(self.N**2))*1000./simtime
       print "Mean Firing rate: %f Hz" % firing_rate

  def save(self, path, vline=False, RefSpk=[], gather=True):
    self.mpi_print("Saving simulation results ...")
    if not os.path.exists(path):
      os.system('mkdir %s' % path)
      print path
    path = "%s/%dx%d" %(path, self.N, self.N)
    if not os.path.exists(path):
      os.system('mkdir %s' % path)
    if 'spikes' in self.recordings:
      if len(RefSpk)==0: 
        for i in xrange(self.N):
          for j in xrange(self.N):
            file = "%s/Bloc_%d_%d.gdf" %(path, i, j)
            self.blocs[i,j]._print('spikes', file, gather=gather)
      else:
        file = "%s/RefSpk.gdf" %path
        self.blocs[tuple(RefSpk)]._print('spikes', file, gather=gather)

    if 'spikes_E' in self.recordings:
      if len(RefSpk)==0: 
        for i in xrange(self.N):
          for j in xrange(self.N):
            file = "%s/E_Bloc_%d_%d.gdf" %(path, i, j)
            self.blocs[i,j]._print('spikes_E', file, gather=gather)
      else:
        file = "%s/E_RefSpk.gdf" %path
        self.blocs[tuple(RefSpk)]._print('spikes_E', file, gather=gather)

    if 'spikes_I' in self.recordings:
      if len(RefSpk)==0: 
        for i in xrange(self.N):
          for j in xrange(self.N):
            file = "%s/I_Bloc_%d_%d.gdf" %(path, i, j)
            self.blocs[i,j]._print('spikes_I', file, gather=gather)
      else:
        file = "%s/I_RefSpk.gdf" %path
        self.blocs[tuple(RefSpk)]._print('spikes_I', file, gather=gather)

    if 'lfp' in self.recordings:
      middle = int(self.N)/2
      if vline==True:
         for j in xrange(self.N):
           file    = "%s/Bloc_%d_%d_exc.lfp" %(path, middle, j)
           self.blocs[middle, j]._print('lfp', file)
      else:
        for i in xrange(self.N):
          for j in xrange(self.N):
            file    = "%s/Bloc_%d_%d_exc.lfp" %(path, i, j)
            self.blocs[i, j]._print('lfp', file)
    if 'vm' in self.recordings:
      middle = int(self.N**2)/2
      file   = "%s/vm.dat" %(path)
      self.blocs.flatten()[middle]._print('vm', file, gather=gather)
    if 'g' in self.recordings:
      middle = int(self.N**2)/2
      file   = "%s/g.dat" %(path)
      self.blocs.flatten()[middle]._print('g', file, gather=gather)

  def connect_nnb(self, idx, w, p_self=0.2, p_nnb=0.1, synapse='excitatory'):
     #Get sources
     if synapse=="excitatory": source_cells=self.get_exc_from_blocks([idx])
     elif synapse=="inhibitory": source_cells=self.get_inh_from_blocks([idx])

     #Get targets
     nnb_idx=up.get_nn_torus(idx,self.blocs.shape[0],self.blocs.shape[1])
     target_nnb=self.get_all_from_blocks(nnb_idx)
     target_self=self.get_all_from_blocks([idx])
     #Connect!

     pr_nnb = Projection(source_cells, target_nnb, MyConnector(numpy.repeat(p_nnb,target_nnb.size), weights = w, delays=self.max_delay), target=synapse, rng=rng)
     pr_self = Projection(source_cells, target_self, MyConnector(numpy.repeat(p_self,target_self.size), weights = w, delays=self.max_delay), target=synapse, rng=rng)
     self.mpi_print("%s rec. conn. and %s lat. conn. etablished from %s %s population"%(pr_self.size(),pr_nnb.size(),idx,synapse) )

def generic_launch_network(path,params, conn_path='./Connections'):
  if not os.path.exists('Results'):
      os.system('mkdir Results')
  save_path="Results/"+path 
  if not os.path.exists(path):
      os.system('mkdir %s' % save_path)
  print "Bulding hypercolumn:"
  print " %sx%s columns"%(params["diag"]["N"],params["diag"]["N"])
  network = Layer(params["diag"]["N"], params["cell_params"],params["col"],le=0,li=0,velocity=None)
  print "Background stim ON..." 
  #network.stimulate(params["stim"])
  if params["sim"]["lgn"]:
    print "LGN loading ..."
    lgn = LGN(params["lgn"],params["diag"]["N"])
    lgn.build_static_bar(params["lgn"]["rate"],params["lgn"]["t_start"], params["lgn"]["duration"], phi=params["lgn"]["angle"],save="inp.dat")
    lgn.project(network,params["lgn"]["factor_inh"], params["lgn"]["clip_val"])
  #for idx in network._indices(): 
     #network.connect_nnb( idx, params["col"]["wE"], params["col"]["p_self"], params["col"]["p_nnb"], "excitatory" )
     #network.connect_nnb( idx, params["col"]["wI"], params["col"]["p_self"], params["col"]["p_nnb"], "inhibitory" )
  
  if "spikes_E" in params["rec"]["rec_var"]: network.record('spikes_E',to_file=False)
  if "spikes_I" in params["rec"]["rec_var"]:  network.record('spikes_I')
  if "lfp" in params["rec"]["rec_var"]:  network.record('lfp', to_file=False)
  network.simulate(params["sim"]["duration"])
  network.save(path=save_path, gather=True)
  
