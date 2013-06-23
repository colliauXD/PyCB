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
         'tau_syn_in': 3.0,
         'tau_w': 100.0}

params={
  'cell_params': cp_aeif,
                  
  'col':{'subpop'  : {"exc","inh"}, 
         'scale'   : 1,
         'size'    : 200,
         'EI_ratio': 4,
         #'model'   : EIF_cond_alpha_isfa_ista,
         'model': {'name': 'aeif_beuler', 'type': 'Native'},
         'wE'      : 1e-3,
         'wI'      : 6e-3,
         'ee'      : 0.05,
         'ii'      : 0.05,
         'ei'      : 0.05,
         'ie'      : 0.05
        },

  'stim':{'type'     : "Poisson_E",
          'target'   : "E", 
          'N_E'      : 200,
          'wE'       : 1e-3,
          'ee'       : 0.01,
          'ei'       : 0.1, 
          'rate_E'   : 850,
          'duration' : 800
          },

  'rec':{'spikes_inp': True,
         'spikes_E'  : 800,
         'Vm_E'      : 1,
         'lfp_E'      : 1,
         'spikes_I'  : 200,
         'Vm_I'      : 1,
         'Mv_I'      : False
         },

  'sim':{'duration': 1000,
         'dt': 0.1,
         'integ': "default"
        },

  'diag':{'N'         : 5,
          'n_source'  : [20,200],  
          'stim_rates': [20,2000],
          'wes'       : [1e-4,5e-3],
          'wis'       : [1e-4,1.5e-2],
          'as'        : [0,200],
          'bs'        : [0,200]
        }
}

cPickle.dump(params,open("col_v_ab.par",'w'))

