# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 17:59:38 2023

@author: yanbw
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 15:56:49 2023

@author: yanbw
"""

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from sys import argv
import P
####MODELO######
vw = 0.074671545
phi = 0.30975736
mi_w = 1
mi_g = 0.0172
MRF = 1.0


def PermEff(S):
        lamb = 5.0 #mobilidade total*
        krg = 1.0 #permeabilidade efetiva (gás)
        krw = 0.75 #permeabilidade efetiva (água)
        Swc = 0.99 #Saturação da água*
        Sgr = 0.01 #Saturação do gás* 
        Swe = (S-Swc)/(1-S-Sgr)  #equação 7 
        print(Swe)
        k_w = krw*Swe**lamb #equação 8
        k_g = krg*(1-Swe)**(3+(2/5.0)) #equação 9
        #if(S!=0):
            #print(S,k_w,k_g)
        return k_w,k_g

def f(Sw):
   krw1,krg1=PermEff(Sw)
   lw=  (krw1 * 0.2 / mi_w) 
   lg=       ((krg1 * 0.00114667) / (mi_g * MRF))      
   lt=lw+lg
   print(lw)
   return lw/lt
####\MODELO######

def init(L,dl):
    nl=int(L/dl)
    r= np.zeros((nl,nl))
    return r



####SOLVER######
def Solve2d(L,dl,t,dt,Sw0):  #dx=dy=dl, X=Y=L
  nt,nl=int(t/dt),int(L/dl)

  rw = (vw*dt)/(dl*phi)
  W,U=P.poission(dl,L)

  Sw=np.zeros((nt+1,nl,nl))
  Sw[0]=init(L,dl)
  for k in range(1,nt):
    ##Borda
    #NEUMAN NA MAIORIA DO DOMINIO -> DSW=0
    Sw[k,0,:]=Sw[k-1,0,:]
    Sw[k,:,0]=Sw[k-1,:,0]
    
    ##SOURCE E SINK
    Sw[k,nl//2-2:nl//2 + 2,0]=0.99
    #Sw[k,:,0]=0.99
    
    
    
    for i in range (1,nl-1):
      for j in range (1,nl-1):
        v=0.5
        w=0.5
        Fy=(1/dl)*((Sw[k-1,i,j])-(Sw[k-1,i-1,j]))
        Fx=(1/dl)*((Sw[k-1,i,j])-(Sw[k-1,i,j-1]))
        
        
       
        Sw[k,i,j]= Sw[k-1,i,j]- dt*(v*Fy + w*Fx)
        
    plt.imshow(Sw[k], cmap='viridis', origin='lower', vmin=0, vmax=1)
    plt.colorbar(label='Sw')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Matrix Sw')
    plt.show()
  # border


  return Sw

s=Solve2d(0.5,0.01,10,0.001,0)
print(s)