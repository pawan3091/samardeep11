import numpy as np
from scipy.integrate import quad


def MyTrap(f,a,b,n):  #input parameter are func,lower,upper limit,no.of panels
    # calculating step size
    h = (b - a) / n
    
    trpzint = f(a) + f(b)
    
    for i in range(1,n):
        k = a + i*h
        trpzint = trpzint + 2 * f(k)
    
    # multiply h/2 with the obtained integration to get trapezoidal integration 
    trpzint =trpzint * h/2
    
    return trpzint


"""function for Simpson with 
upper limit b lower limit a and n intervals"""
def MySimp(f,a,b,n):
    # calculating step size
    h = (b - a) / n
    
    simpint = f(a) + f(b)
    
    for i in range(1,n):
        k = a + i*h
        
        if i%2 == 0:
            simpint = simpint + 2 * f(k)
        else:
            simpint = simpint + 4 * f(k)          
    
    # multiply h/2 with the obtained integration to get Simpson integration 
    simpint =simpint * h/3
    
    return simpint




def MyTrap_tol(f,a,b,N=2,N_max = 1000,tol=0.5e-3):
    trial = 0    #trial is a dummy to check whether tolerance is achieved to give out certain output mssg
    while N <= N_max :
        nlist = np.array([])     #we will append the evaluated values for N and 2*N in nlist array
        for i in range(N,2*N +1,N):    #i have only 2 possible values i.e N and 2*N
        
            h = abs((b - a)/i)   #h would be different for N and 2*N
            trpzint = f(a) + f(b)
            for j in range(1,i):
                k = a + j*h
                trpzint = trpzint + 2 * f(k)
    
    # multiply h/2 with the obtained integration to get trapezoidal integration 
            trpzint =trpzint * h/2
                
                          
            nlist = np.append(nlist,trpzint)         #we will append the final values for N and 2*N in this array 
                                                   
        
        if abs((nlist[1]-nlist[0])/nlist[1]) <= tol:     #break the loop if tolerence is achieved
           
            break
        elif abs((nlist[1]-nlist[0])/nlist[1]) > tol and 2*N <=N_max:    #if tolerence isn't achieved but 2*N <= Nf then continue with N = 2*N
            
            N = 2*N
        else :
            trial =1              #if tolerence is not achieved within N_max,change trial value to 1,and show error! 
            
            break
        
    eval_int = nlist[1]
    err = eval_int - quad(f, a, b)[0]   #quad func gives the value of integration within given limits  
    if trial == 0 :
        msg = "tolerance is obtained at n = "+str(2*N)
        return [eval_int,err,msg]
    else :
        msg = "ERROR !!! tolerance is not obtained,try increasing maximum limit of N" 
        print(msg)
        
        
        
def MySimp_tol(f,a,b,N=2,N_max = 1000,tol=0.5e-3):
    trial = 0      #trial is a dummy to check whether tolerance is achieved to give out certain output mssg
    while N <= N_max :
        nlist = np.array([])    #we will append the evaluated values for N and 2*N in nlist array
        h1 = abs((b-a)/N)      #h1 and h2 are step-sizes for N and 2*N
        h2 = abs((b-a)/(2*N))
        
        S1 = 1/3*(f(a) + f(b))  #S defines summation for all odd terms and T is summation for all even terms
        T1,T2,I1,I2 = 0,0,0,0  
        for i in range(2,N,2):
            S1 += (2/3)*(f(a+i*h1))
        for i in range(1,N,2):
            T1 += (2/3)*(f(a+i*h1)) 
        for i in range(1,2*N,2):
            T2 += (2/3)*(f(a+(i)*h2)) 
            
        I1 = h1*(S1 + 2*T1)
        I2 = h2*(S1 + T1 +2*T2 ) #here,we're calculating S2 from S1 and T1,this is to reduce computation
            
        nlist = np.append(nlist,I1)
        nlist = np.append(nlist,I2)  
        #    nlist = np.append(nlist,simpson_integration)       
        #we will append the final values for N and 2*N in this array 
        if abs((nlist[1]-nlist[0])/nlist[1]) <= tol:       #break the loop if tolerence is achieved
           
            break
        elif abs((nlist[1]-nlist[0])/nlist[1]) > tol and 2*N <=N_max:    #if tolerence isn't achieved but 2*N <= Nf then continue with N = 2*N
            N = 2*N
        else :
            trial =1             #if tolerence is not achieved within N_max,change trial value to 1,and show error! 
            
            break
        
    eval_int = nlist[1]
    err = eval_int - quad(f, a, b)[0]    #quad func gives the value of integration within given limits  
    if trial == 0 :
        msg = "tolerance is obtained at n = "+str(2*N)
        return [eval_int,err,msg]
    else :
        msg = "ERROR !!! tolerance is not obtained,try increasing maximum limit of N" 
        print(msg)
