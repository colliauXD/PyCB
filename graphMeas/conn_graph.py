import networkx as nx
import pylab as pl

class OneColGraph(object):
    def __init__(self, idx):
       self.G = nx.DiGraph()
       self.col_idx = idx
       self.stim_ids = []
       self.e_ids = []
       self.i_ids = []

    def my_in_degree_iter(self, sources, targets):
       nodes_nbrs=((n,self.G.pred[n]) for n in self.G.nbunch_iter(targets))
       for n,nbrs in nodes_nbrs:
          yield (n,len(set(nbrs.keys()).intersection(sources)))

    def my_in_degree(self, sources, targets):
       return dict(self.my_in_degree_iter(sources, targets))

    def my_out_degree_iter(self, sources, targets):
       nodes_nbrs=((n,self.G.succ[n]) for n in self.G.nbunch_iter(sources))
       for n,nbrs in nodes_nbrs:
          yield (n,len(set(nbrs.keys()).intersection(targets)))

    def my_out_degree(self, sources, targets):
       return dict(self.my_out_degree_iter(sources, targets))


    def make_stim_nodes(self, stim_ids, edg):
    	self.stim_ids=stim_ids
        self.G.add_nodes_from( stim_ids, key={"cell": "stim"})
        self.G.add_edges_from(edg)

    def make_layer_nodes(self, e_ids, i_ids):
        self.e_ids = e_ids
        self.i_ids = i_ids
       	self.G.add_nodes_from( e_ids, key={"cell": "E"})
        self.G.add_nodes_from( i_ids, key={"cell": "I"})

    def conn_col(self,conn_EE,conn_EI,conn_IE,conn_II):
        self.G.add_edges_from(conn_EE)
        self.G.add_edges_from(conn_EI)
        self.G.add_edges_from(conn_IE)
        self.G.add_edges_from(conn_II)

    def get_deg(self, in_out="in"):
        if in_out=="in":
           self.se=self.my_in_degree(self.stim_ids,self.e_ids)
           self.ee=self.my_in_degree(self.e_ids,self.e_ids)
           self.ei=self.my_in_degree(self.e_ids,self.i_ids)
           self.ie=self.my_in_degree(self.i_ids,self.e_ids)
           self.ii=self.my_in_degree(self.i_ids,self.i_ids)
        else:
           self.se=self.my_out_degree(self.stim_ids,self.e_ids)
           self.ee=self.my_out_degree(self.e_ids,self.e_ids)
           self.ei=self.my_out_degree(self.e_ids,self.i_ids)
           self.ie=self.my_out_degree(self.i_ids,self.e_ids)
           self.ii=self.my_out_degree(self.i_ids,self.i_ids)

        
    def init_plot(self):
        f=pl.figure(facecolor="k")
    	ax=f.add_axes([0,0,1,1],axisbg="k")
        stim_g=self.G.subgraph(self.stim_ids)
        self.ps=nx.spring_layout(stim_g,iterations=0)
        for pos in self.ps.values():
           pos+=[0.,1]
        e_g=self.G.subgraph(self.e_ids)
        self.pe=nx.spring_layout(e_g,iterations=15)
        for pos in self.pe.values():
           pos+=1
        i_g=self.G.subgraph(self.i_ids)
        self.pi=nx.spring_layout(i_g,iterations=15)
        for pos in self.pi.values():
           pos-=[-1,1]
        return ax

    def plot_stim(self,my_ax):
       stim_g=self.G.subgraph(self.stim_ids)
       nx.draw_networkx_nodes(self.G, nodelist=self.stim_ids, ax=my_ax, pos=self.ps, node_size=20,
       	                              node_color="lightsalmon", alpha=.02)

    def plot_col(self,my_ax):
       e_g=self.G.subgraph(self.e_ids)
       i_g=self.G.subgraph(self.i_ids)
       nx.draw_networkx_nodes(e_g, ax=my_ax, pos=self.pe, node_size=20,
                                   node_color="tomato", alpha=.02)
       nx.draw_networkx_nodes(i_g, ax=my_ax, pos=self.pi, node_size=20,
                                   node_color="steelblue", alpha=.02)

    def plot_stim2colE(self,my_ax):
       stim2E_p=dict(self.ps, **self.pe) 
       edg=self.G.edges(self.stim_ids+self.e_ids)
       edg=[ei for ei in edg if (not(ei[0] in self.i_ids) and (ei[0] not in self.e_ids))]
       nx.draw_networkx_edges(self.G, stim2E_p, edgelist=edg, ax=my_ax, 
                              edge_color="lightsalmon", alpha=.01, arrows=False)
                     
    def plot_col2col(self,my_ax):
       col_p=dict(self.pe, **self.pi)
       edg=self.G.edges(self.e_ids+self.i_ids)
       edg=[ei for ei in edg if ((ei[0] in self.e_ids) and (ei[1] not in self.stim_ids))]
       nx.draw_networkx_edges(self.G, col_p,edgelist=edg, ax=my_ax, 
                              edge_color="tomato", alpha=.01, arrows=False)
       edg=self.G.edges(self.e_ids+self.i_ids)
       edg=[ei for ei in edg if ((ei[0] in self.i_ids) and (ei[1] not in self.stim_ids))]
       nx.draw_networkx_edges(self.G, col_p,edgelist=edg, ax=my_ax, 
                                      edge_color="steelblue", alpha=.01,arrows=False)
       
       
    def plot_target(self,my_ax,idx):
       stim2E_p=dict(self.ps, **self.pe)
       target_pred=self.G.predecessors(idx)
       s_pred=[]
       e_pred=[]
       i_pred=[]
       for t in target_pred:
           if t in self.stim_ids: s_pred.append(t)
           if t in self.e_ids: e_pred.append(t)
           if t in self.i_ids: i_pred.append(t)
           
       nx.draw_networkx_nodes(self.G, nodelist=[idx],ax=my_ax, pos=self.pe, node_size=20, node_color="red", alpha=1)
       nx.draw_networkx_nodes(self.G, nodelist=e_pred,ax=my_ax, pos=self.pe, node_size=20, node_color="tomato", alpha=0.8)
       nx.draw_networkx_nodes(self.G, nodelist=s_pred,ax=my_ax, pos=self.ps, node_size=20, node_color="lightsalmon", alpha=0.8)
       nx.draw_networkx_nodes(self.G, nodelist=i_pred,ax=my_ax, pos=self.pi, node_size=20, node_color="steelblue", alpha=0.8)

       nx.draw_networkx_edges(self.G, stim2E_p, edgelist=zip(s_pred, [idx]*len(s_pred)) ,ax=my_ax, edge_color="lightsalmon", alpha=0.5, arrows=False)
       nx.draw_networkx_edges(self.G, self.pe, edgelist=zip(e_pred, [idx]*len(e_pred)) ,ax=my_ax, edge_color="tomato", alpha=0.5, arrows=False)
       nx.draw_networkx_edges(self.G, dict(self.pe, **self.pi), edgelist=zip(i_pred, [idx]*len(e_pred)) ,ax=my_ax, edge_color="steelblue", alpha=0.3, arrows=False)

    def make_text(self,my_ax):
       my_ax.text(0.3,1.4,"%s Stim"%len(self.stim_ids),{"color":"salmon","fontsize":30}) 
       my_ax.text(1.2,1.4,"%s E"%len(self.e_ids),{"color":"red","fontsize":30}) 
       my_ax.text(1.2,-0.5,"%s I"%len(self.i_ids),{"color":"blue","fontsize":30}) 

    def plot_deg(self, my_ax):
       deg=self.se.values() 
       a = pl.axes([.25, .88, .2, .1], axisbg='none') 
       h=a.hist(deg,max(deg)-min(deg),color="salmon")
       pl.plot([sum(deg)/(1.*len(deg)),sum(deg)/(1.*len(deg))],[0,max(h[0])],"red",linewidth=3)
       pl.setp(a,xticks=[min(deg),round(sum(deg)/(1.*len(deg)),2),max(deg)],
               xticklabels=[min(deg),round(sum(deg)/(1.*len(deg)),2),int(max(deg))],
               yticks=[max(h[0])],
               ylabel=r"$Stim \rightarrow E$")
       pl.setp(a.xaxis.get_ticklabels(),color="salmon")
       pl.setp(a.yaxis.get_ticklabels(),color="salmon")
       pl.setp(a.yaxis.label,color="salmon")

       deg=self.ee.values() 
       a = pl.axes([0.6, .88, .2, .1], axisbg='none') 
       h=a.hist(deg,max(deg)-min(deg),color="tomato")
       pl.plot([sum(deg)/(1.*len(deg)),sum(deg)/(1.*len(deg))],[0,max(h[0])],"red",linewidth=3)
       pl.setp(a,xticks=[min(deg),round(sum(deg)/(1.*len(deg)),2),max(deg)],
               xticklabels=[min(deg),round(sum(deg)/(1.*len(deg)),2),int(max(deg))],
               yticks=[max(h[0])],
               ylabel=r"$E \rightarrow E$")
       pl.setp(a.xaxis.get_ticklabels(),color="tomato")
       pl.setp(a.yaxis.get_ticklabels(),color="tomato")
       pl.setp(a.yaxis.label,color="tomato")

       deg=self.ii.values() 
       a = pl.axes([0.6, 0.1, .2, .1], axisbg='none') 
       h=a.hist(deg,max(deg)-min(deg),color="steelblue")
       pl.plot([sum(deg)/(1.*len(deg)),sum(deg)/(1.*len(deg))],[0,max(h[0])],"blue",linewidth=3)
       pl.setp(a,xticks=[min(deg),round(sum(deg)/(1.*len(deg)),2),max(deg)],
               xticklabels=[min(deg),round(sum(deg)/(1.*len(deg)),2),int(max(deg))],
               yticks=[max(h[0])],
               ylabel=r"$I \rightarrow I$")
       pl.setp(a.xaxis.get_ticklabels(),color="steelblue")
       pl.setp(a.yaxis.get_ticklabels(),color="steelblue")
       pl.setp(a.yaxis.label,color="steelblue")

       deg=self.ie.values() 
       a = pl.axes([0.35, 0.4, .2, .1], axisbg='none') 
       h=a.hist(deg,max(deg)-min(deg),color="steelblue")
       pl.plot([sum(deg)/(1.*len(deg)),sum(deg)/(1.*len(deg))],[0,max(h[0])],"blue",linewidth=3)
       pl.setp(a,xticks=[min(deg),round(sum(deg)/(1.*len(deg)),2),max(deg)],
               xticklabels=[min(deg),round(sum(deg)/(1.*len(deg)),2),int(max(deg))],
               yticks=[max(h[0])],
               ylabel=r"$I \rightarrow E$")
       pl.setp(a.xaxis.get_ticklabels(),color="steelblue")
       pl.setp(a.yaxis.get_ticklabels(),color="steelblue")
       pl.setp(a.yaxis.label,color="steelblue")

       deg=self.ei.values() 
       a = pl.axes([0.75, 0.55, .2, .1], axisbg='none') 
       h=a.hist(deg,max(deg)-min(deg),color="tomato")
       pl.plot([sum(deg)/(1.*len(deg)),sum(deg)/(1.*len(deg))],[0,max(h[0])],"red",linewidth=3)
       pl.setp(a,xticks=[min(deg),round(sum(deg)/(1.*len(deg)),2),max(deg)],
               xticklabels=[min(deg),round(sum(deg)/(1.*len(deg)),2),int(max(deg))],
               yticks=[max(h[0])],
               ylabel=r"$E \rightarrow I$")
       pl.setp(a.xaxis.get_ticklabels(),color="tomato")
       pl.setp(a.yaxis.get_ticklabels(),color="tomato")
       pl.setp(a.yaxis.label,color="tomato")

       #ax.bar(h[1],h[0])

    def plot_it(self,fname="lala.png"):
       ax=self.init_plot()
       self.plot_stim2colE(ax)
       self.plot_col2col(ax)
       self.plot_stim(ax)
       self.plot_col(ax)
       self.plot_target(ax,30)
       self.make_text(ax)
       self.plot_deg(ax)
       pl.savefig(fname)
       
