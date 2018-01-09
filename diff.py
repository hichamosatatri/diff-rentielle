# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 14:09:04 2018

@author: wislamos
"""

import numpy as np
import matplotlib.pyplot as plt

class Diff():
    def __init__(self,n):
        self.n=n
        self.h=1/(n+1)
        self.a=self.uexact(0)
        self.b=self.uexact(1)
        self.A=1/(self.h**2)*(np.diag([-1]*(self.n-1),-1)+np.diag([2]*(self.n),0)+np.diag([-1]*(self.n-1),1))
        self.B=np.zeros((self.n,1))
        self.X=np.zeros((self.n,1))
        self._genererX()
        self._genererB()        
        self.solution=self.calculerU()
    def _genererB(self):
        self.B[0]=self.f(0)+self.a/(self.h**2)
        self.B[self.n-1]=self.f(1)+self.b/(self.h**2)
        for i in range(1,self.n-1):
            self.B[i]=self.f(self.X[i])
    def calculerU(self):
        return np.linalg.inv(self.A)@self.B
        
    def f(self,x):
        return 2*x
    
    def uexact(self,x):
        return -1/3*(x**3)
    def _genererX(self):
        for i in range(self.n):
            self.X[i][0]=self.h*(i+1)
        
    def trace(self):
        v=self.uexact(self.X)
        plt.plot(self.X,self.solution,color='red')
        plt.plot(self.X,v,color='green')
        plt.show()
        return v

u=Diff(21)

#print("Matrice A :\n",u.A)
#print("B : \n",u.b)
#print("Les Xi : \n",u.X)
#print("Solution : \n",u.solution)
v=u.trace()
erreur=np.linalg.norm(u.solution-v)
print(erreur)
