from pylab import *
from math import *
import cPickle 

def is_in_range(i,j,N):
   a=True
   #if i<0: a=False
   #if j<0: a=False
   #if i>N-1: a=False
   #if j>N-1: a=False
   return a 
   

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
  
def export_RF(NLGN,N):
   om=buildMap(N,1)
   LGN=[]
   for i in range(N):
     for j in range(N):
       conn=[]
       inp=RF(NLGN,alpha=om[i,j])
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
  
