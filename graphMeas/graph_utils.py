import pylab as pl

def plot_mr(mr,params,path):
   pl.imshow(mr,origin="lowest",interpolation="nearest")
   pl.title("Pop. exc.: firing rate")
   pl.colorbar()
   pl.yticks([0,params["diag"]["N"]],(params["diag"]["stim_rates"][0],params["diag"]["stim_rates"][1]))
   pl.xticks([0,params["diag"]["N"]],(params["diag"]["wis"][0],params["diag"]["wis"][1]))
   pl.xlabel(r"inh. max cond. $g_{max}^I$")
   pl.ylabel(r"Input rate : $r_{ext}$")
   pl.savefig(path+"_mr.png")
   pl.clf()

def plot_stdr(stdr,params,path):
   pl.imshow(stdr,origin="lowest",interpolation="nearest")
   pl.title("Pop. exc.: std firing rate")
   pl.colorbar()
   pl.yticks([0,params["diag"]["N"]],(params["diag"]["stim_rates"][0],params["diag"]["stim_rates"][1]))
   pl.xticks([0,params["diag"]["N"]],(params["diag"]["wis"][0],params["diag"]["wis"][1]))
   pl.xlabel(r"inh. max cond. $g_{max}^I$")
   pl.ylabel(r"Input rate : $r_{ext}$")
   pl.savefig(path+"_stdr.png")
   pl.clf()

def plot_cv(cv,params,path):
   pl.imshow(cv,origin="lowest",interpolation="nearest")
   pl.title("Pop. exc.: CV")
   pl.colorbar()
   pl.yticks([0,params["diag"]["N"]],(params["diag"]["stim_rates"][0],params["diag"]["stim_rates"][1]))
   pl.xticks([0,params["diag"]["N"]],(params["diag"]["wis"][0],params["diag"]["wis"][1]))
   pl.xlabel(r"inh. max cond. $g_{max}^I$")
   pl.ylabel(r"Input rate : $r_{ext}$")
   pl.savefig(path+"_cv.png")
   pl.clf()

def plot_cc(cc,params,path):
   pl.imshow(cc.reshape(params["diag"]["N"],params["diag"]["N"]),origin="lowest",interpolation="nearest")
   pl.title("Pop. exc.:  pw. cc.")
   pl.colorbar()
   pl.yticks([0,params["diag"]["N"]],(params["diag"]["stim_rates"][0],params["diag"]["stim_rates"][1]))
   pl.xticks([0,params["diag"]["N"]],(params["diag"]["wis"][0],params["diag"]["wis"][1]))
   pl.xlabel(r"inh. max cond. $g_{max}^I$")
   pl.ylabel(r"Input rate : $r_{ext}$")
   pl.savefig(path+"_cc.png")
   pl.clf()
