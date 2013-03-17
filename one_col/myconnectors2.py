# -*- coding: utf-8 -*-
from pyNN import common, connectors, random
import nest, numpy
import pyNN.connectors
from pyNN.random import RandomDistribution

class MyProbabilisticConnector(pyNN.connectors.Connector):

    def __init__(self, projection, weights=0.0, delays=None, allow_self_connections=True, space=pyNN.space.Space(), safe=True):

        pyNN.connectors.Connector.__init__(self, weights, delays, space, safe)
        if isinstance(projection.rng, random.NativeRNG):
            raise Exception("Use of NativeRNG not implemented.")
        else:
            self.rng = projection.rng
        if self.delays is None:
            self.delays = projection._simulator.state.min_delay
        self.local             = projection.post._mask_local
        self.N                 = projection.post.size
        self.weights_generator = pyNN.connectors.WeightGenerator(weights, self.local, projection, safe)
        self.delays_generator  = pyNN.connectors.DelayGenerator(delays, self.local, kernel=projection._simulator.state, safe=safe)
        self.probas_generator  = pyNN.connectors.ProbaGenerator(RandomDistribution('uniform', (0,1), rng=self.rng), self.local)
        self._distance_matrix  = None
        self.projection        = projection
        self.candidates        = projection.post.local_cells
        self.size              = self.local.sum()
        self.allow_self_connections = allow_self_connections

    def _probabilistic_connect(self, src, p):
        rarr       = self.probas_generator.get(self.N)
        create     = rarr < p[self.local]
    	mytargets  = self.candidates[create].tolist()	
        weights    = self.weights_generator.get(self.N)[create]
        delays     = self.delays_generator.get(self.N)[create]
        if not self.allow_self_connections and src in mytargets:
            assert len(mytargets) == len(weights) == len(delays)
            i       = mytargets.index(src)
            weights = numpy.delete(weights, i)
            delays  = numpy.delete(delays, i)
            mytargets.remove(src)
        if len(mytargets) > 0:
            self.projection._divergent_connect(src, mytargets, weights, delays)


class MyConnector(pyNN.connectors.Connector):
    
    def __init__(self, p_connect, allow_self_connections=True, weights=0.0, delays=None, space=pyNN.space.Space(), safe=True, verbose=False):
        pyNN.connectors.Connector.__init__(self, weights, delays, space, safe, verbose)
        assert isinstance(allow_self_connections, bool)
        self.allow_self_connections = allow_self_connections
        self.p_connect              = numpy.minimum(p_connect, 1)    
	
    def connect(self, projection):
        connector = MyProbabilisticConnector(projection, self.weights, self.delays, self.allow_self_connections, self.space, safe=self.safe)  
        self.progressbar(len(projection.pre))
        for count, src in enumerate(projection.pre.all()):      
            connector._probabilistic_connect(src, self.p_connect)
            self.progression(count, projection._simulator.state.mpi_rank)

