import measures as mea
import pylab as pl

filename="col_stimE_wI"
s=mea.get_spikes(5,500,2000,filename,"E")
for i in range(len(s)):
   s[i].raster_plot()
   pl.savefig("rs%s.png"%i)
   pl.clf()
   pl.plot(s[i].firing_rate(5,average=True))
   pl.savefig("rate%s.png"%i)
   pl.clf()
