# -*- coding: utf-8 -*-
from myconnectors2 import *
from NeuroTools.stgen import StGen
import cPickle, time, gc
import numpy 
from pyNN.nest import *
from pyNN.recording import files
import measures as mea
import graph_utils as gu

nest.Install("mymodule")

global  seed
rng      = NumpyRNG(seed=249856, parallel_safe=False)
seed     = 512468

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
  def __init__(self, N):
    self.N   = N
    self.LGN = Population((N, N), SpikeSourceArray, {'spike_times':[]})

  def mpi_print(self, string):
    if simulator.state.mpi_rank == 0: print string    

  def build_input(self, stim_duration=4000, omega=20, angle=0):
      self.mpi_print("Creating or loading the LGN input ....")
      f0    = 0
      f1    = 60
      angle = (angle/360.) * 2 * numpy.pi 
      time  = numpy.arange(0, stim_duration, 1)
      spk   = StGen(seed=5423689)
      st    = []
      for i in range(self.N):
         for j in range(self.N):            
            signal    = f0 + f1*(numpy.cos(omega*.001*time+numpy.pi*(i*numpy.cos(angle) + j*numpy.sin(angle))/(self.N*1.))+1)
            spk_times = spk.inh_poisson_generator(signal, time, stim_duration, array=True)
            st.append(spk_times)
      for cell, spikes in zip(self.LGN, st):
        cell.spike_times = spikes      

  def _connect_cell2pop(self, cell, targets, layer, factor=10):
    weights   = numpy.array([], float)
    Ne        = len(layer.blocs[0, 0].exc)
    Ni        = len(layer.blocs[0, 0].inh)  
    indices   = []    
    for target in targets:
      indices  += [target[1]]
      weights   = numpy.concatenate((weights, [target[0]]))
    d      = RandomDistribution('uniform', [timestep, 30*timestep], rng=rng)
    weight = RandomDistribution('normal', (factor*layer.w_E, layer.w_E/3), boundaries=(0, layer.w_E*100), constrain='redraw', rng=rng)
    cell   = self.N*cell[0]+cell[1]
    c1     = MyConnector(numpy.repeat(weights, layer.Ne+layer.Ni), weights=weight, delays=d)
    prj1   = Projection(self.LGN[[cell]], layer.get_all_from_blocks(indices), c1, target='excitatory', rng=rng)
    return prj1.size()

  def project_on(self, layer):
    self.mpi_print("Building the connections from the LGN ....")
    LGNOn  = cPickle.load(open('Connections/learn/longlearn/V1/LGNOnAfferent.dat','r'))
    n_conn = 0
    for neuron in LGNOn:
      idx = neuron[0]
      in_range_targets = []
      if layer.is_in_range(idx):
        for target in neuron[1:]:
          if layer.is_in_range(target[1]):
            if target[0] > 0:
                in_range_targets.append(target)
        n_conn += self._connect_cell2pop(idx, in_range_targets, layer)
    self.mpi_print("--> %d excitatory connections were established from LGN to network..." %n_conn)

  def record(self, to_file=True):
      self.LGN.record(to_file=to_file)

  def save(self, file, gather=True):
      self.LGN.printSpikes(files.NumpyBinaryFile(file + "/lgn.gdf", 'w'), gather=gather)        


class Layer(object):

  def __init__(self, N, cell_params,col_params,sim, le=1., li=1., velocity=None):
    self.mpi_print("Creation of the layer ...")
    self.velocity = velocity
    self.scaling  = 0.06 # width, in mm, of one column
    self.timestep=sim['dt']
    if velocity is None:
      max_delay = self.timestep
    else:      
      max_delay    = numpy.sqrt(2)*N*self.scaling/self.velocity
    self.mpi_print("Timestep is %g ms, max delay has been set to %g ms" %(self.timestep, max_delay))

    self.node_id    = setup(self.timestep, max_delay=max_delay)
    self.N          = N
    self.blocs      = numpy.array([Column(cell_params,col_params) for i in xrange(N**2)])
    self.blocs      = self.blocs.reshape(N,N)
    self.recordings = []
    self.le         = le
    self.li         = li
    self.w_E        = self.blocs[0,0].w_E
    self.w_I        = self.blocs[0,0].w_I
    self.Ne         = len(self.blocs[0,0].exc)
    self.Ni         = len(self.blocs[0,0].inh)
    self.mpi_print("A layer of size %d x %d x %d = %d cells has been built" %(self.N, self.N, self.Ne+self.Ni, (self.Ne+self.Ni) * self.N**2))
    

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

  def make_stim_pop(self, stimulus):
     if (stimulus['type'] == 'Poisson_E'):
        self.ext_stim = Population(stimulus['N_E'], native_cell_type('poisson_generator'), {'rate': stimulus['rate_E'], 'stop' : float(stimulus['duration'])})

  def stim_n_sources(self, stimulus, idx, n_source):
    if stimulus["target"]=="E":
       N_target=self.Ne
       target_cells=self.get_exc_from_blocks(idx)
    if stimulus["target"]=="I":
       N_target=self.Ni
       target_cells=self.get_inh_from_blocks(idx)
    if stimulus["target"]=="EI":
       N_target=self.Ne+self.Ni
       target_cells=self.get_all_from_blocks(idx)
    conn    = MyConnector(numpy.repeat(stimulus['ee'],N_target), weights = stimulus['wE'], delays=self.timestep)
    prj     = Projection(self.ext_stim.sample(n_source), target_cells, conn, target='excitatory', rng=rng)
    
  def stimulate_on_target(self, stimulus, idx):
    if (stimulus['type'] == 'Poisson_E' or stimulus['type'] == 'Poisson_EI'):
      self.ext_stim = Population(stimulus['N_E'], native_cell_type('poisson_generator'), {'rate': stimulus['rate_E'], 'stop' : float(stimulus['duration'])})
    #weight  = RandomDistribution('normal', (self.w_E, self.w_E/3), boundaries=(0, self.w_E*100), constrain='redraw', rng=rng)

    if stimulus["target"]=="E":
       N_target=self.Ne
       target_cells=self.get_exc_from_blocks(idx)
    if stimulus["target"]=="I":
       N_target=self.Ni
       target_cells=self.get_inh_from_blocks(idx)
    if stimulus["target"]=="EI":
       N_target=self.Ne+self.Ni
       target_cells=self.get_all_from_blocks(idx)
    conn    = MyConnector(numpy.repeat(stimulus['ee'],N_target), weights = stimulus['wE'], delays=self.timestep)
    prj     = Projection(self.ext_stim, target_cells, conn, target='excitatory', rng=rng)	     	
    self.mpi_print("-->o %d excitatory connections were established from external sources ..." %prj.size())

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
    self.mpi_print("Launching simulation ...")
    os.system('rm -rf *.gdf')
    if simulator.state.mpi_rank == 0:
      nest.SetStatus([0], {'print_time':True})
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

  def _connect_blocs(self, source, targets, synapse='exc'):
    weights   = numpy.array([], float)
    distances = numpy.array([], float)
    indices   = []
    for target in targets:      
      weights   = numpy.concatenate((weights, [target[0]]))
      distances = numpy.concatenate((distances, [distance_no_torus(source, target[1])]))
      indices  += [target[1]]
    
    if not self.velocity is None:
      d = self.timestep + distances*self.scaling/self.velocity
      d = numpy.repeat(d, self.Ne+self.Ni)
    else:
      d = RandomDistribution('uniform', [self.timestep, 30*self.timestep], rng=rng)
      
    if synapse is 'exc':
      source = self.blocs[source].exc
      factor = self.le * weights
      weight = RandomDistribution('normal', (self.w_E, self.w_E/3), boundaries=(0, self.w_E*100), constrain='redraw', rng=rng)
      c1     = MyConnector(numpy.repeat(factor, self.Ne+self.Ni), weights=weight, delays=d) 
      prj1   = Projection(source, self.get_all_from_blocks(indices), c1, target='excitatory', rng=rng)
    elif synapse is 'inh':
      source = self.blocs[source].inh
      factor = self.li * weights
      weight = RandomDistribution('normal', (self.w_I, self.w_I/3), boundaries=(0, self.w_I*100), constrain='redraw', rng=rng)
      c1     = MyConnector(numpy.repeat(factor, self.Ne+self.Ni), weights=weight, delays=d)
      prj1   = Projection(source, self.get_all_from_blocks(indices), c1, target='inhibitory', rng=rng)
    result = prj1.size()
    return result


  def connect_single_col(self,idx, col_params):
     prjEE=Projection(self.blocs[idx].exc,self.blocs[idx].exc,MyConnector(numpy.repeat(col_params['ee'],self.Ne),weights=col_params['wE'],delays=self.timestep),target="excitatory",rng=rng)
     prjEI=Projection(self.blocs[idx].exc,self.blocs[idx].inh,MyConnector(numpy.repeat(col_params['ei'],self.Ni),weights=col_params['wE'],delays=self.timestep),target="excitatory",rng=rng)
     prjIE=Projection(self.blocs[idx].inh,self.blocs[idx].exc,MyConnector(numpy.repeat(col_params['ie'],self.Ne),weights=col_params['wI'],delays=self.timestep),target="inhibitory",rng=rng)
     prjII=Projection(self.blocs[idx].inh,self.blocs[idx].inh,MyConnector(numpy.repeat(col_params['ii'],self.Ni),weights=col_params['wI'],delays=self.timestep),target="inhibitory",rng=rng)
     return prjEE.size(), prjEI.size(), prjIE.size(), prjII.size()
     #return prjEE.size()

  def is_in_range(self, position):
    x,y  = position 
    x_ok = (0 <= x < self.blocs.shape[0])
    y_ok = (0 <= y < self.blocs.shape[1])
    return x_ok and y_ok

  def _import_map(self, conn_map, synapse='exc'):  
    projections = cPickle.load(open(conn_map,'r'))
    nb_conn     = 0
    for projection in projections:   
      source = projection[0]
      in_range_targets = []
      if self.is_in_range(source):  
        for target in projection[1:]:
            if self.is_in_range(target[1]):
                if target[0] > 0:
                    in_range_targets.append(target)
        nb_conn += self._connect_blocs(source, in_range_targets, synapse)
    return nb_conn

  def topographica_connect(self, exc_map, inh_map):
    self.mpi_print("Building local excitatory AND inhibitory connections ...")
    nb_conn = self._import_map(exc_map, 'exc')
    self.mpi_print("--> %d local excitatory connections were established within columns ..." %nb_conn)
    nb_conn = self._import_map(exc_map, 'inh')
    self.mpi_print("--> %d local inhibitory connections were established within columns ..." %nb_conn)
    self.mpi_print("Building long-range excitatory connections ...")
    #self.le = 2.5
    #nb_conn = self._import_map(inh_map, 'exc')
    #self.mpi_print("--> %d long range excitatory connections were established within columns ..." %nb_conn)


def generic_launch_network(path,params,verbose=True):
  if not os.path.exists('Results'):
      os.system('mkdir Results')

  if not os.path.exists(path):
      os.system('mkdir %s' % path)
  #Build Layer 
  network = Layer(params["diag"]["N"], params['cell_params'],params['col'], params['sim'],le=0, li=0, velocity=None)
  
  #Stimulate with Poisson and setup internal connections
  #n_source=numpy.arange(params["diag"]["n_source"][0],params["diag"]["n_source"][1],step=(params["diag"]["n_source"][1]-params["diag"]["n_source"][0])/params["diag"]["N"],dtype=numpy.int)
  sr=numpy.linspace( params["diag"]["stim_rates"][0],params["diag"]["stim_rates"][1],params["diag"]["N"])
  network.make_stim_pop(params["stim"])
  wis=numpy.linspace(params["diag"]["wis"][0],params["diag"]["wis"][1],params["diag"]["N"])
  #As=numpy.linspace( params["diag"]["as"][0],params["diag"]["as"][1],params["diag"]["N"])
  #Bs=numpy.linspace( params["diag"]["bs"][0],params["diag"]["bs"][1],params["diag"]["N"])

  for idx in network._indices():
    print "Setting connections in column", idx
    print "with", sr[idx[0]], "Hz sources"
    print "and internal inh conductance", wis[idx[1]]
    #network.stim_n_sources(params['stim'],[idx],n_source[idx[0]])
    params["stim"]["rate_E"]=sr[idx[0]]
    #params["cell_params"]["a"]=As[idx[1]]
    #params["cell_params"]["b"]=Bs[idx[1]]
    #set(network.get_exc_from_blocks([idx]),param={'a': params["cell_params"]["a"]})
    #,'b': params["cell_params"]["b"]})
    network.stimulate_on_target(params["stim"],[idx])
    params["col"]["wI"]=wis[idx[1]]
    r=network.connect_single_col(idx, params['col'])
    print "+->o"
    print "|__|",r, "(ee,ei,ie,ii) internal connections"
  network.record('spikes_E', to_file=False)
  network.record('spikes_I', to_file=False)
  network.record('lfp', to_file=False)
  #network.record('vm', to_file=False)
  #network.record('g', to_file=False)
  network.simulate(params['sim']['duration'], print_rates=True)
  network.save(path=path, gather=True)
  end()

    
file_name="col_stimE_wI"
path="Results/"+file_name
params=cPickle.load(open(file_name+".par"))
generic_launch_network(path,params)

s=mea.get_spikes(params["diag"]["N"], 500, 10000, path,"E")
time_r=mea.get_t_rates(s, 5)
mr,stdr=mea.get_rates_m_std(s,5)
cv=mea.get_cvs(s)
cc=mea.pw_pearson_corcoef(s, 100, 5)

numpy.save(path+"_time_r",time_r)
numpy.save(path+"_mr",mr)
numpy.save(path+"_stdr",stdr)
numpy.save(path+"_cv",cv)
numpy.save(path+"_cc",cc)

gu.plot_mr(mr, params,path)
gu.plot_stdr(stdr, params,path)
gu.plot_cc(cc, params,path)
gu.plot_cv(cv, params,path)
