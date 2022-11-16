
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import odeint as solve

def f(x):
    return -1*(1+x**2)

def cond_b(a,b,N):
    h = (b-a)/N
    y_1 = 1 + (h**2)/2 + (1/8)*h**4
    return y_1


def numerov(a,b,y_0,y_1,N,f):
    
    h = (b-a)/N
    x = np.arange(a,b+h,h)
    y_arr = np.array([y_0,y_1])
    
    for i in range(2,N+1):
        c1 = 2*(1-(5/12)*(h**2)*f(a+(i-1)*h))*y_arr[i-1]
        c2 = (1+((h**2)/12)*f(a+(i-2)*h))*y_arr[i-2]
        c3 = 1+((h**2) /12)*f(a+i*h)
        
        y_arr = np.append(y_arr,(c1-c2)/c3)
        
    return [x,y_arr]

def func(u,x):
    return (u[1],(1+x**2)*u[0])
        
if __name__ == "__main__":
    a = 0
    b = 1
    #for N = 2
    y_0 = 1
    
    N = 2
    y_1 = cond_b(a, b, N)
    
    solution1 = numerov(a, b, y_0, y_1, N, f)
    inbuilt1 = solve(func,(1,0),solution1[0])[:,0]
    print(inbuilt1)
    
    
    N = 4
    y_1 = cond_b(a, b, N)
    
    solution2 = numerov(a, b, y_0, y_1, N, f)
    inbuilt2 = solve(func,(1,0),solution2[0])[:,0]
    print(inbuilt2)
    
    
    #part c
    
    data = {"x":solution1[0],'u_num':solution1[1],'u_inbuilt':inbuilt1,'E':abs(inbuilt1-solution1[1])}
    print("DATA TABLE FOR N = 2 ")
    print(pd.DataFrame(data))
    print('')
    
    data2 = {"x":solution2[0],'u_num':solution2[1],'u_inbuilt':inbuilt2,'E':abs(inbuilt2-solution2[1])}
    print("DATA TABLE FOR N = 4")
    print(pd.DataFrame(data2))
    
    
    #part d
    
    k = np.arange(1,8)
    N_list = 2**k
    for i in N_list:
        y_1 = cond_b(a, b, i)
        solution = numerov(a, b, y_0, y_1, i, f)
        inbuilt = solve(func,(1,0),solution[0])[:,0]
        plt.plot(solution[0],solution[1],label='computed ,N='+str(i))
        plt.plot(solution[0],inbuilt,label='inbuilt ,N='+str(i),ls = '--')
        
    plt.xlabel('X')
    plt.ylabel('U(X)')
    plt.grid()
    plt.legend()
    plt.title("PLOTS OF U VS X FOR DIFFERENT N")
    plt.show()