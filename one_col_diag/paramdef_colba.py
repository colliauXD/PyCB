#Parameter set for the Brunel99 network

from pyNN.nest import *
import cPickle

cp_aeif={'C_m': 280,
         'Delta_T': 2.0,
         'E_L': -70.0,
         'E_ex': 0.0,
         'E_in': -80.0,
         'I_e': 0.0,
         'V_peak': -30.0,
         'V_reset': -70.0,
         'V_th': -50.0,
         'a': 0.0,
         'b': 0.0,
         'g_L': 30.0,
         't_ref': 5.0,
         'tau_syn_ex': 2.0,
         'tau_syn_in': 13.0,
         'tau_w': 100.0}

params={
  'cell_params': cp_aeif,
                  
  'col':{'subpop'  : {"exc","inh"}, 
         'scale'   : 1,
         'size'    : 200,
         'EI_ratio': 4,
         #'model'   : EIF_cond_alpha_isfa_ista,
         'model': {'name': 'aeif_beuler', 'type': 'Native'},
         'wE'      : 6e-4,
         'wI'      : 5e-3,
         'ee'      : 0.01,
         'ii'      : 0.01,
         'ei'      : 0.1,
         'ie'      : 0.1
        },

  'stim':{'type'     : "Poisson_E",
          'target'   : "EI", 
          'N_E'      : 200,
          'wE'       : 6e-4,
          'ee'       : 0.01,
          'ei'       : 0.01, 
          'rate_E'   : 2700,
          'duration' : 2000
          },

  'rec':{'spikes_inp': True,
         'spikes_E'  : 800,
         'Vm_E'      : 1,
         'Mv_E'      : 0,
         'spikes_I'  : 200,
         'Vm_I'      : 1,
         'Mv_I'      : False
         },

  'sim':{'duration': 2000,
         'dt': 0.1,
         'integ': "default"
        },

  'diag':{'N'         : 3,
          'stim_rates':[200,10000],
          'wis'       :[1e-3,1e-2]
        }
}

cPickle.dump(params,open("Big_col.par",'w'))
params["cell_params"]["a"]=2.
cPickle.dump(params,open("Big_col_a.par",'w'))
params["cell_params"]["a"]=0.
params["cell_params"]["b"]=3.
cPickle.dump(params,open("Big_col_b.par",'w'))
params["cell_params"]["a"]=2.
params["cell_params"]["b"]=3.
cPickle.dump(params,open("Big_col_ab.par",'w'))
