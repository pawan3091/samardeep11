import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import hbar,pi,epsilon_0,e,m_e

#frequency of revolution of an electron revolving around the nucleus in the nth Bohr orbit of hydrogen atom
def fcln(n):
    return m_e*e**4/(32*(pi**3)*(epsilon_0**2)*(hbar**3)*n**3)

#frequency of radiation emitted by an atom when the electron makes transition from n to n-1
def fqn(n):
    return (m_e*e**4*(2*n-1))/(64*(pi**3)*(epsilon_0**2)*(hbar**3)*(n**2)*(n-1)**2)
 
def fun(tol ):
    p = np.array([0.5])
    n = 10**p

    cln_arr = np.array([fcln(n[-1])])
    qn_arr = np.array([fqn(n[-1])])

    rel_err = np.array([abs(qn_arr[-1] - cln_arr[-1])*100/qn_arr[-1]])
    
    
    while rel_err[-1] >= tol :

        p = np.append(p,p[-1] + 0.5)
        
        n = 10**p[-1]

        cln_arr = np.append(cln_arr,fcln(n))
        qn_arr = np.append(qn_arr,fqn(n))

        rel_err = np.append(rel_err,abs(qn_arr[-1] - cln_arr[-1])*100/qn_arr[-1])
    
    
    n = 10**p
    lists = [n,cln_arr,qn_arr,rel_err]
    mat = np.zeros(len(n)*4).reshape(len(n),4)

    for i in range(len(n)):
        for j in range(4):
            mat[i][j] = lists[j][i]

    print(mat)

    #plot of percentage relative difference as a function of ln(n)
    plt.plot(n,rel_err)
    plt.title("plot between relative error and ln(n)")
    plt.xlabel("rel error")
    plt.ylabel("ln(n)")
    plt.xscale('log')
    plt.grid()
    plt.show()

tol= 10**(-5)
fun(tol)