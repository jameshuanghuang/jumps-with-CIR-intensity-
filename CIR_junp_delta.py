#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 23:20:09 2020

@author: jameshuang
"""

import numpy as np
import matplotlib.pyplot as plt
import itertools
import math
from scipy import integrate
import sympy as sym
from scipy.interpolate import lagrange
from scipy.integrate import quad

alpha =0.9 ;lam0=0.9 ;delt=1 ;sigma=1; MT=1; beta=1.2;
k=np.sqrt(delt**2+2*sigma**2);
Lam=[lam0];

def Fnc(s):
    lam= Lam[-1];
    delt=1 ;sigma=1;
    p1 = (  2*k*np.exp((k+delt)/2*s)/( (k+delt)*(np.exp(k*s)-1)+2*k )  )**(2*alpha*delt/(sigma**2));
    p2 = np.exp( lam*(-2)*(np.exp(k*s)-1)/((k+delt)*(np.exp(k*s)-1)+2*k));

    return 1-p1*p2;

def diff_F(s_v,alpha_v,lam_v,delt_v,sigma_v):
    s = sym.Symbol('s');
    alpha = sym.Symbol('alpha');
    delt = sym.Symbol('delt');
    sigma = sym.Symbol('sigma');
    k=sym.sqrt(delt**2+2*sigma**2);
    lam = sym.Symbol('lam');
    p1 = (  2*k*sym.exp((k+delt)/2*s)/( (k+delt)*(sym.exp(k*s)-1)+2*k )  )**(2*alpha*delt/(sigma**2));
    p2 = sym.exp( lam*(-2)*(sym.exp(k*s)-1)/((k+delt)*(sym.exp(k*s)-1)+2*k));
    Fc = 1-p1*p2;
    f  = sym.diff(Fc, s);
    ans = f.subs({s:s_v,alpha:alpha_v,lam:lam_v,delt:delt_v,sigma:sigma_v});
    #print(sym.simplify(f));
    #diff_F(0.5,0.9,0.9,1,1);
    return float(ans) 


def L(X,j):
    
    y = np.array([0 for i in range(len(X))]);
    y[j]=1;
    #return lagrange(x, y);
    return lagrange(X, y);
    
def H(X,j):
    
    hw = lambda x: np.exp(-x**2)*L(X,j)(x); #Hermite 
    return quad(hw,min(X),max(X))[0] ;

def F_inv(a,b,lam):
    U = np.random.uniform(0,1);
    
    el= 2**(-2);
    p = lambda s: s*np.exp(-(Fnc(s)-U)**2/(4*el))/(2*np.sqrt(np.pi*el));
    fsds = lambda s: diff_F(s,0.9,lam,1,1);
    f = lambda s: p(s)*fsds(s); 
    ans = quad(f,a,b);
    print('s_dirct',ans[0]);
    return ans[0];
    
'''
def GQ_H(a,b,lam):
    n=2**(2);  #exact up to def 2n+1;
    h=(b-a)/n;
    x = np.arange(a,b+h,h);
    el=2**(-2);
    U=np.random.uniform(0,1);
    X = [(Fnc(e)-U)/(2*el) for e in x];
    
    # sum of w_i * X_i *  f(x_i)
    #Xi = np.array( [np.sqrt(2)*el*xi for xi in X] );
    #f_xi = np.array([diff_F(e,0.9,lam,1,1) for e in x]);
    w = np.array( [H(X,j) for j in range(len(X))] );
    
    #ans = sum(w*x)*2*el; 
    
    #ans1 = sum(w*x)*2*el/np.sqrt(np.pi); 
    
    ans2 = sum(w*x)/np.sqrt(np.pi)*np.sqrt(el);
    
    
    #print(ans);
    #print(ans1);
    print('s_GQ:',ans2);
    return ans2;
'''

b=5;a=0;
trunc=b+1;
S=[];

while sum(S)<MT:
    lam=Lam[-1];
    ans1= F_inv(a,trunc,lam); #GQ_H(a,trunc,lam)
    if  ans1>0:
        S.append(ans1);
        
        s = ans1;
        B = sigma**2*(np.exp(k*s)-1);
        C = (k-delt)+(k+delt)*np.exp(k*s);
        D = 2*alpha*delt/sigma**2;
        E = (k+delt)+(k-delt)*np.exp(k*s);
        F = 2*(np.exp(k*s)-1);
    
        J=np.random.poisson(lam*(E/B-F/C));
        G1 = np.random.gamma(J+D+1,B/C);  #shape and rate (C/B)^-1
        G2 = np.random.gamma(J+D+2,B/C);
        w1 = D*B/(D*B+lam*(E-F*(B/C)));
        
        #simulate lambi i+1
        temp2 = np.random.choice([G1,G2],p=[w1,1-w1]);
        #self excite jump    
        Lam.append(temp2+np.random.exponential(beta));
        
    else:
        print('s<0',ans1);
        break;





T= list(itertools.accumulate(S));   
Nw=[1 for e in range(len(T)) ];
Nw= list(itertools.accumulate(Nw));

plt.plot(T,Nw)
 


       
             
    

#F_inv(a,b,lam0);
#GQ_H(a,b,lam0);



















