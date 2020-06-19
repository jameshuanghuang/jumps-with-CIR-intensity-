#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 21:18:46 2020

@author: jameshuang
"""

#import quant_func
import  numpy as np
import random 
import matplotlib.pyplot as plt
import itertools
#import quant_func as qf

alpha =0.9 ;lam0=0.9 ;delt=1 ;sigma=1; MT=2; beta=1.2;
mu=1/beta;
k=np.sqrt(delt**2+2*sigma**2);
elp =delt-mu;
#exact sol
EN_T = alpha*delt*MT/elp+(lam0-alpha*delt/elp)*(1-np.exp(-elp*MT))/elp;

    
Lam=[lam0];
Ti=0;
N=[0];
interS=[];
Y=[np.random.exponential(beta)];



#simulate S_star
trial=0;
while 1:
    U1=np.random.uniform(0.1);
    Wg=2*k*(U1**(-(sigma**2)/(alpha*delt*k*(k-delt)))-1)/(k+delt);
    
    U2=np.random.uniform(0.1);
    p1=((k+delt)/(2*k))**(  ( (k+delt)/(2*k) )*( (2*alpha*delt)/sigma**2 )  );
    p2=((Wg+1)/(2*k/(k+delt)+Wg))**(alpha*delt*(k+delt)/(sigma**2*k));
    p3=Wg/(Wg+1);
    trial+=1;
    if U2<=p1*p2*p3:
        S=(1/k)*np.log(Wg+1);
        #print('n:',trial);
        break;

i=0;
while 1: #Ti<=MT

    #simulate S(i)
    U3=np.random.uniform(0.1);
    lambi=Lam[i];
    d = ( 1+ np.log(U3)/(2*lambi)*(k+delt) )/( 1-np.log(U3)/(2*lambi)*(k-delt) );
    if d>0:
        V = -(1/k)*np.log(d);
        Si=min(S,V);
    else:
        Si=S; 
    interS.append(Si);
    
    Ti=sum(interS);
    if Ti>MT:
        break;
    

    s = Si;
    B = sigma**2*(np.exp(k*s)-1);
    C = (k-delt)+(k+delt)*np.exp(k*s);
    D = 2*alpha*delt/sigma**2;
    E = (k+delt)+(k-delt)*np.exp(k*s);
    F = 2*(np.exp(k*s)-1);

    J=np.random.poisson(lambi*(E/B-F/C));
    G1 = np.random.gamma(J+D+1,B/C);  #shape and rate (C/B)^-1
    G2 = np.random.gamma(J+D+2,B/C);
    w1 = D*B/(D*B+lambi*(E-F*(B/C)));
    
    #simulate lambi i+1
    temp = np.random.choice([G1,G2],p=[w1,1-w1]);
    #self excite jump
    if i>0:
        Y.append(np.random.exponential(beta));
    
    Lam.append(temp+Y[i]);
    N.append(N[i]+1);
    i+=1; #num of iter

T=list(itertools.accumulate(interS))[:-1];
T.insert(0,0);
plt.plot(T,N);
print(max(N));
#plt.plot(np.arange(0,len(Lam)),Lam);








