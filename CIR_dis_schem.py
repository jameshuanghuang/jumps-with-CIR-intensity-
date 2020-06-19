#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 24 20:58:22 2020

@author: jameshuang
"""
import numpy as np
import random 
import matplotlib.pyplot as plt
import itertools
import math


#Descretize Scheme
alpha =0.9 ;lam0=0.9 ;delt=1 ;sigma=1; MT=2; beta=1.2;
mu=1/beta;
k=np.sqrt(delt**2+2*sigma**2);
elp =delt-mu
#exact sol
#EN_T = alpha*delt*MT/elp+(lam0-alpha*delt/elp)*(1-np.exp(-elp*MT))/elp;


grid=10;
h=MT/grid;
mesh= np.arange(0,MT+h,h);
J= len(mesh);
#initialize Lam
Lam=[]
for i in range(J):
    if i==0:
        Lam.append(lam0);
    else:
        temp= Lam[i-1] + delt*(alpha-Lam[i-1])*h + sigma*np.sqrt(max(Lam[i-1],0)*h)*np.random.normal(0,1);
        Lam.append(temp);
    

j=0; #j 0~J-1
T=[0];N=[0];

while mesh[j]<MT:
    Comp_lam = 0 if j==0 else sum(Lam[:j])*h;
    temp = np.random.exponential(1);
    prev=j;flag=0;
    ar_temp=[];
    
    #find the set
    for i in range(j,J):
        if sum(Lam[:i])*h >= Comp_lam + temp:
            ar_temp.append(1);
        else:
            ar_temp.append(0);
       
    #find the inf of the set       
    key=prev;
    for i in range(len(ar_temp)-1):
        if ar_temp[i]==0 and ar_temp[i+1]==1 and flag==0:
            flag=1; 
            key=prev+i+1;
            temp2 = np.random.exponential(beta);
            #self excite jump
            Lam[key]+=temp2;
            #update lam
            for e in range(prev,key-1):
                Lam[e+1]=Lam[e] + delt*(alpha-Lam[e])*h + sigma*np.sqrt(max(Lam[i-1],0)*h)*np.random.normal(0,1);
            
    if flag==0:
        j=J-1;
    else:
        T.append(mesh[key]);
        N.append(1);
        j=key;
    

N= list(itertools.accumulate(N));
plt.plot(T,N)


        
        
        
        
    