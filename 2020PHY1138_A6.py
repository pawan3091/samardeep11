
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
import pandas as pd
from scipy.special import hermite
import math

om = 5.5e14
h_bar = 1.05e-34
m = 9.31e-31

x_o = np.sqrt(h_bar/(m*om))

def Numerov(IC,a,b,N,alp,inti,arr = None):
    alp = np.vectorize(alp)
    x = np.linspace(a,b,int(N)+1)
    h = x[1] - x[0]
    u = np.zeros(len(x))
    IC[1] = inti(h,arr)
    u[:2] = IC
    c = (1 + ((h**2/12)*alp(x)))
    
    for i in range(2,len(x)):
        
        u[i] = (1/c[i])*(((12-(10*c[i-1]))*u[i-1]) - (c[i-2]*u[i-2]))
   
    return x,u

def inti_4(h,arr= None):
    return -h    

def Analytic(x,n):
    p = hermite(n)
    return (p(x))/(np.exp(x**2/2))

Analytic = np.vectorize(Analytic)

prob_an = lambda x,n : Analytic(x,n)**2
area = lambda x,n : integrate.simps(prob_an(x,n),x)
new_anl = lambda x,n : Analytic(x,n)/np.sqrt(area(x,n))

new_pro = lambda x,n : new_anl(x,n)**2




def solve_eigen(x_min,x_max,N,e):
    def slope(x):
        return 2*e - x**2 
    slope = np.vectorize(slope)

    def chk_turn(x_min,x_max,N) :
        x = np.linspace(x_min,x_max,int(N)+1)  
        f = slope(x)
        k = -1
        for i in range(len(x)-1):
            if f[i]*f[i+1] < 0:
               k = i
            else:
                k = -1
        if k > N-10:
            x_m = x_max+50
        else:
            x_m = x_max
        return x_m
    
    x,si_o = Numerov([0,1],chk_turn(x_min,x_max,N),x_min,N,slope,inti_4 )
    
    a = integrate.simps(np.power(si_o,2),x)
    nom_sk = si_o/np.sqrt(abs(a))   

    return x,nom_sk



del_e = np.logspace(-2,-8,base = 10,num = 4)

e = lambda n,d : n + (1/2) + d
e = np.vectorize(e)

def check_energy(n,x_max,N,d):
    
    x,u= solve_eigen(0,x_max,N,e(n,d))
    
    X = np.concatenate((-x[:-1] , x[::-1]))
    if n%2 ==0:
        U = np.concatenate((u[:-1] , u[::-1]))
       
    else :
        U = np.concatenate((-u[:-1] , u[::-1]))
        
    a = integrate.simps(np.power(U,2),X)
    U_f = U/np.sqrt(abs(a))    
    
    return X,U_f


def plotting(x,y,title,key,key2,key3,func = None,):
    m = np.shape(y)[0]
    if key == 'p':
       fig,ax = plt.subplots()
       ax.set_title(title)
    for i in range(m):
        dic = {0:f'for energies tweaked by {del_e[i]}',1:'numerical'}
        if key == 'w':
           fig,ax = plt.subplots()
           ax.set_title(f'{title} for {i+1} state')
        ax.scatter(x[i],y[i],s = 10,label=f'wave function {dic[key3]}')
        if key2 == 'y':
           ax.plot(x[i],func(x[i],i),label = 'Analytic')
        ax.set_xlabel('x')
        ax.set_ylabel('u')
        ax.legend()
    plt.show()  


check_energy = np.vectorize(check_energy,otypes=[np.ndarray,np.ndarray])

x,u= check_energy(0,4,100,del_e)        

x1,u1 = check_energy([0,1,2,3],5,200,0)

plotting(x, u,'ground state for different energies' ,'p',0,0)  
plotting(x1,u1,f'wavefunction with energies tweaked by {del_e[-2]}' ,'w','y',1,new_anl)    

U_den = np.power(u1,2)

plotting(x1,U_den,'Probablity density','p','y',1,new_pro)


data = {'state':[0,1,2,3],f'E tweaked by {del_e[-2]}':e([0,1,2,3],del_e[-2])*h_bar*om/(1.6e-19),'Analytic energy':e([0,1,2,3],0)*h_bar*om/(1.6e-19)}
df = pd.DataFrame(data)
print(df)


def prob_for(e,n):
    
    X,U = check_energy(n,4,100,0)
    
    def slope(x):
        return 2*e - x**2 
    slope = np.vectorize(slope)
    f = slope(X)
    index = []
    for i in range(len(f)-1):
        if f[i]*f[i+1] < 0:
            index.append(i)
        elif f[i] == 0:
            index.append(i)
    
    x1 = X[:index[0]+1] ; x2 = X[index[1]:]
    U1 = U[:index[0]+1] ; U2 = U[index[1]:]
    
    P = integrate.simps(np.power(U1,2),x1) + integrate.simps(np.power(U2,2),x2)
    return P
    
print('Probablity of finding electron in classically forbidden region: ',prob_for(e(0,0),0))

