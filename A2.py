# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 18:27:08 2022

@author: Tanmeet singh
"""

import numpy as np
import matplotlib.pyplot as plt

# 2(a)

def RK4(func, X0, tmin, tmax, N):
    h = (tmax-tmin)/N
    t = np.linspace(tmin, tmax, N+1)
    X = np.zeros([N+1, len(X0)])
    X[0] = X0
    for i in range(N):
        k1 = func(t[i], X[i])
        k2 = func(t[i] + h/2, X[i] + (h * k1)/2)
        k3 = func(t[i] + h/2, X[i] + h/2 * k2)
        k4 = func(t[i] + h, X[i] + h * k3)
        X[i+1] = X[i] + h / 6. * (k1 + 2*k2 + 2*k3 + k4)
    return X, t


def simps(a, b, n, y):
    h = (b-a)/n
    integral = (h/3)*(2*np.sum(y[2:-2:2]) + 4*np.sum(y[1:-1:2]) + y[0] + y[-1])

    return integral

# 2(b) : e=8, u'=1

def func(x, Y):
    y, y1 = Y           # y1 is the first order derivative
    f1 = y1
    f2 = -(8)*y  # e=8
    return np.array([f1, f2])


ic = [0, 1]  # initial conditions: u=0, u'=1
X, t = RK4(func, ic, -1/2, 1/2, 50)

wf = X[:, 0]  # wave function
nv = (wf)/np.sqrt(simps(-1/2, 1/2, 50, wf**2))  # normalized function

si = np.cos(np.pi*t)
nsi = si/(np.sqrt(simps(-1/2, 1/2, 50, si**2)))

plt.plot(t, nsi, label='analytical')
plt.plot(t, nv, label='computed')
plt.grid()
plt.legend()
plt.xlabel('x')
plt.ylabel('Ψ(x)')
plt.title("for e=8 and u'=1")

plt.show()

# changing the value of u'
ic1 = [0, 30]  # initial conditions: u=0, u'=30
X1, t1 = RK4(func, ic1, -1/2, 1/2, 50)

wf1 = X1[:, 0]  # wave function
nv1 = (wf1)/np.sqrt(simps(-1/2, 1/2, 50, wf1**2))  # normalized function

si = np.cos(np.pi*t1)
nsi = si/(np.sqrt(simps(-1/2, 1/2, 50, si**2)))

plt.plot(t1, nsi, label='analytical')
plt.plot(t1, nv1, label='computed')
plt.grid()
plt.legend()
plt.xlabel('x')
plt.ylabel('Ψ(x)')
plt.title("for e=8 and u'=30")

plt.show()

# 2(c)----------------------------------------------------------------------


def func1(x, Y):
    y, y1 = Y  # y1 is the first order derivative
    f1 = y1
    f2 = -(11)*y  # e=11
    return np.array([f1, f2])


ic2 = [0, 1]  # initial conditions: u=0, u'=1
X2, t2 = RK4(func1, ic2, -1/2, 1/2, 50)

wf2 = X2[:, 0]  # wave function
nv2 = (wf2)/np.sqrt(simps(-1/2, 1/2, 50, wf2**2))  # normalized function

si = np.cos(np.pi*t2)
nsi = si/(np.sqrt(simps(-1/2, 1/2, 50, si**2)))

plt.plot(t2, nsi, label='analytical')
plt.plot(t2, nv2, label='computed')
plt.grid()
plt.legend()
plt.xlabel('x')
plt.ylabel('Ψ(x)')
plt.title("for e=11 and u'=1")
plt.show()

# changing the value of u'
ic3 = [0, 30]  # initial conditions: u=0, u'=30
X3, t3 = RK4(func1, ic3, -1/2, 1/2, 50)

wf3 = X3[:, 0]  # wave function
nv3 = (wf3)/np.sqrt(simps(-1/2, 1/2, 50, wf3**2))  # normalized function

si = np.cos(np.pi*t3)
nsi = si/(np.sqrt(simps(-1/2, 1/2, 50, si**2)))

plt.plot(t3, nsi, label='analytical')
plt.plot(t3, nv3, label='computed')
plt.grid()
plt.legend()
plt.xlabel('x')
plt.ylabel('Ψ(x)')
plt.title("for e=11 and u'=30")
plt.show()

# 2(d)-----------------------------------------------------------------------

e = 0.9*(np.pi**2)
while e <= 1.1*(np.pi**2):

    def func2(x, Y):
        y, y1 = Y  # y1 is the first order derivative
        f1 = y1
        f2 = -(e)*y
        return np.array([f1, f2])

    ics = [0, 1]
    x, t = RK4(func2, ics, -1/2, 1/2, 50)
    wfn = x[:, 0]
    nwfn = wfn/(np.sqrt(simps(-1/2, 1/2, 50, wfn**2)))

    if abs(wfn[-1]) <= 0.5*10**(-5):
        print('ground state eigen value: ', e)

        si = np.cos(np.pi*t)
        nsi = si/(np.sqrt(simps(-1/2, 1/2, 50, si**2)))
        plt.plot(t, nsi, label='analytical', color='red')
        plt.plot(t, nwfn, label='computed', ls='dotted', color='black')
        plt.grid()
        plt.legend()
        plt.xlabel('x')
        plt.ylabel('Ψ(x)')
        plt.title("the ground state normalised eigen function")
        
        plt.show()
        break

    e += 10**(-4)

# 2(e)---------------------------------------------------------------------

e = 3.9*(np.pi**2)
while e <= 4.1*(np.pi**2):

    def func2(x, Y):
        y, y1 = Y  # y1 is the first order derivative
        f1 = y1
        f2 = -(e)*y
        return np.array([f1, f2])

    ics = [0, -1]
    x, t = RK4(func2, ics, -1/2, 1/2, 50)
    wfn = x[:, 0]
    nwfn = wfn/(np.sqrt(simps(-1/2, 1/2, 50, wfn**2)))

    if abs(wfn[-1]) <= 0.5*10**(-5):
        print('first state eigen value: ', e)

        si = np.sin(2*np.pi*t)
        nsi = si/(np.sqrt(simps(-1/2, 1/2, 50, si**2)))
        plt.plot(t, nwfn)
        plt.plot(t, nsi)
        plt.plot(t, nsi, label='analytical')
        plt.plot(t, nwfn, label='computed', ls='dotted')
        plt.grid()
        plt.legend()
        plt.xlabel('x')
        plt.ylabel('Ψ(x)')
        plt.title("the first state normalised eigen function")
        
        plt.show()
        break

    e += 10**(-4)
