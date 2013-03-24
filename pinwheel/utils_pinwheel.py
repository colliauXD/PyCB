import numpy as np

#indices coping and distances
def is_in_range(position,Nx,Ny):
    x,y  = position 
    x_ok = (0 <= x < Nx)
    y_ok = (0 <= y < Ny)
    return x_ok and y_ok

def rescale(position,Nx,Ny):
    x,y=position
    if x >= Nx: x-= Nx
    if y >= Ny: y-= Ny
    return (x,y) 

def distance_eucl(pos_1, pos_2):
    dx = pos_1[0]-pos_2[0]
    dy = pos_1[1]-pos_2[1]
    return np.sqrt(dx*dx + dy*dy)

def distance_torus(pos_1, pos_2,N):
    #!Should be tested!
    dx = min(np.abs(pos_1[0]-pos_2[0]),np.abs(N-pos_2[0]-pos_1[0]))
    dy = min(np.abs(pos_1[1]-pos_2[1]),np.abs(N-pos_2[1]-pos_1[1]))
    return np.sqrt(dx*dx + dy*dy)


#Nearest neighbors
def get_nn_torus(idx,Nx,Ny):
    nnb_idx=[]
    for i in range(-1,1):
        for j in range(-1,1):
           x,y=(idx[0]+i,idx[1]+j) 
           if x >= Nx: x-= Nx
           if y >= Ny: y-= Ny
           nnb_idx.append((x,y))
    return nnb_idx

def get_nn_square(idx,Nx,Ny):
    nnb_idx=[]
    for i in range(-1,1):
        for j in range(-1,1):
           x,y=(idx[0]+i,idx[1]+j) 
           if (0 <= x < Nx) and  (0 <= y < Ny): nnb_idx.append((x,y))
    return nnb_idx 

#Compute angles and pinwheel map 
def angleP(ig,jg,N,H):
   ih=ig/N
   jh=jg/N
   il=ig%N-1.*N/2.+1./2.
   jl=jg%N-1.*N/2.+1./2.
   if (ih%2==0 and jh%2==0): a=atan2(il,jl)/2.
   if (ih%2==1 and jh%2==0): a=atan2(-il,jl)/2.
   if (ih%2==0 and jh%2==1): a=atan2(il,-jl)/2.
   if (ih%2==1 and jh%2==1): a=atan2(-il,-jl)/2.
   return a

def buildMap(N=10,H=3):
   Nt=N*H
   phi=zeros([Nt,Nt])
   for ig in xrange(Nt):
      for jg in xrange(Nt):
         t=angleP(ig,jg,N,H)
         #print t 
         phi[ig,jg]=t   
   return phi
   

def exc_proj(N,lgr):
  E=[]
  for i in xrange(N):
     for j in xrange(N):
        neuron=[]
        neuron.append((i,j))
        for k in xrange(3):
           for l in xrange(3):
              x=i-1+k
              y=j-1+l
	      w=1./(8+lgr)
              if (i==x and j==y): w=w*lgr 
              if is_in_range(x,y,N): neuron.append((w,(x,y)))
        E.append(neuron)
  f=open('./Connections/LateralExcitatory.dat','w')
  cPickle.dump(E,f)
  f.close()

def inh_proj(N,lgr):
  
  I=[]
  for i in xrange(N):
     for j in xrange(N):
        neuron=[]
        neuron.append((i,j))
        for k in xrange(3):
           for l in xrange(3):
              x=i-1+k
              y=j-1+l
              w=1./(8+lgr)
              if (i==x and j==y): w=w*lgr 
              if is_in_range(x,y,N): neuron.append((w,(x,y)))
        I.append(neuron)
  f=open('./Connections/LateralInhibitory.dat','w')
  cPickle.dump(I,f)
  f.close()

def RF(size=11,sigma=.4,alpha=0.):
  RF=zeros([size,size])
  x=linspace(-1,1,size)
  y=linspace(-1,1,size)
  for i,xi in enumerate(x):
    for j,yi in enumerate(y):
      RF[i,j]=(1 - 2.*(xi*cos(alpha)+yi*sin(alpha))**2/sigma**2) *exp(- ( (xi*cos(alpha)+yi*sin(alpha))**2 + (-xi*sin(alpha)+yi*cos(alpha))**2 )/sigma**2)
  return RF
  
def export_RF():
   N=11
   om=buildMap(N,1)
   LGN=[]
   for i in range(N):
     for j in range(N):
       conn=[]
       inp=RF(alpha=om[i,j])
       print i, j
       conn.append((i,j))
       for k in range(len(inp)):
          for l in range(len(inp[:,0])):
             if abs(inp[k,l])>0.001:
	        conn.append((inp[k,l],(k,l)))
       LGN.append(conn)       
   f=open('RFs.dat','w')
   cPickle.dump(LGN,f)
   f.close()
   
def test():
   N=11
   N=11
   om=buildMap(N,1)
   LGN=[]
   for i in range(N):
     for j in range(N):
       
       conn=[]
       inp=RF(alpha=om[i,j]*pi/2.)
       print i, j
