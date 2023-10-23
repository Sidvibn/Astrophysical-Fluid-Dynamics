import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

#Radius of Planet = 1 and GM = 1 where M is the mass of the planet.

x_arr_athena = np.load('/home/sid/athena/isothermal_hotjupiter/x_arr_512.npy')
print(x_arr_athena)

# x_arr_athena = np.linspace(1,1.0084,2052)

# print(x_arr_athena)

# def f_numerical(rho,pres,i):
#     x_arr = x_arr_athena #np.linspace(0,L,arr_size)
#     phi_arr = -230.14/x_arr # 10*x_arr 
#     # phi_arr[0] = CubicSpline(x_arr[1:],phi_arr[1:],extrapolate=True)(1.)*2
#     phi_x1 = phi_arr[i+1]
#     phi_x = phi_arr[i]
#     return (rho-pres+(pres+rho)/2*(phi_x1-phi_x))

# def f_numerical_prime(rho,pres,i):
#     x_arr = x_arr_athena #np.linspace(0,L,arr_size)
#     phi_arr = -230.14/x_arr # 10*x_arr
#     # phi_arr[0] = CubicSpline(x_arr[1:],phi_arr[1:],extrapolate=True)(1.)*2
#     phi_x1 = phi_arr[i+1]
#     phi_x = phi_arr[i]
#     return (1+(1/2)*(phi_x1-phi_x))

def f_brentq_iso(x,rho,x_x,x_x1,phi_x,phi_x1): 
    return (0.011*x-0.011*rho)/(x_x1-x_x)+((rho+x)*(phi_x1-phi_x))/(2*(x_x1-x_x))

arr_size = len(x_arr_athena)
L = x_arr_athena[-1]-x_arr_athena[0]
rho_arr = np.zeros(arr_size)
rho0 = 1
rho_arr[0] = rho0

for i in range(arr_size-1):
    rho = rho_arr[i]
    x_arr = x_arr_athena
    phi_arr = -3.15/x_arr #230.14
    x_x = x_arr[i]
    x_x1 = x_arr[i+1]
    phi_x = phi_arr[i]
    phi_x1 = phi_arr[i+1]
    # print(f_brentq_iso(rho-0.5,rho,phi_x,phi_x1))
    # print(f_brentq_iso(rho+0.5,rho,phi_x,phi_x1))
    # print(rho_arr[i])
    rho_arr[i+1] = brentq(f_brentq_iso,0,1,args=(rho,x_x,x_x1,phi_x,phi_x1))
    # print(rho_arr[i+1])
    # print(i)
    # print(rho_arr[i])

x_arr = x_arr_athena
np.savetxt("iso.csv",rho_arr,delimiter=",")

plt.plot(x_arr,rho_arr)
plt.xlabel("r")
plt.ylabel(r"$\rho$")
plt.yscale("log")
plt.show()

# print(rho_arr)

def f_analytic_iso(x):
    phi = -3.15/x
    return np.exp(-2.3*(phi-phi[0])/0.02476)

x_arr = x_arr_athena
rho_arr = f_analytic_iso(x_arr)
np.savetxt("ana_iso.csv",rho_arr,delimiter=",")

plt.plot(x_arr,rho_arr)
plt.xlabel("r")
plt.ylabel(r"$\rho$")
plt.yscale("log")
plt.show()

def f_brentq_isen(x,rho,x_x,x_x1,phi_x,phi_x1): 
    return (0.011*x**(1.666666666667)-0.011*rho**(1.666666666667))/(x_x1-x_x)+((rho+x)*(phi_x1-phi_x))/(2*(x_x1-x_x))

arr_size = len(x_arr_athena)
L = x_arr_athena[-1]-x_arr_athena[0]
rho_arr = np.zeros(arr_size)
rho0 = 1
rho_arr[0] = rho0

for i in range(arr_size-1):
    rho = rho_arr[i]
    x_arr = x_arr_athena
    phi_arr = -3.15/x_arr #230.14
    x_x = x_arr[i]
    x_x1 = x_arr[i+1]
    phi_x = phi_arr[i]
    phi_x1 = phi_arr[i+1]
    gamma = 1.666666666667
    # print(np.abs(gamma/(gamma-1)*0.02476/(2.3*phi_arr[0]/x_arr[0])))
    # print(f_brentq_isen(rho-0.5,rho,phi_x,phi_x1))
    # print(f_brentq_isen(rho+0.5,rho,phi_x,phi_x1))
    # print(rho_arr[i])
    # if rho > 0.02638903529502747:
    rho_arr[i+1] = brentq(f_brentq_isen,0,1,args=(rho,x_x,x_x1,phi_x,phi_x1)) #0.0264
    # else: #rho > 0.007581049862317402:
        # rho_arr[i+1] = brentq(f_brentq_isen,rho-7e-3,rho+7e-3,args=(rho,phi_x,phi_x1))
    # print(rho_arr[i])

x_arr = x_arr_athena
np.savetxt("isen.csv",rho_arr,delimiter=",")

plt.plot(x_arr,rho_arr)
plt.xlabel("r")
plt.ylabel(r"$\rho$")
plt.yscale("log")
plt.show()


def f_analytic_isen(x):
    phi = -3.15/x
    gamma = 1.666666666667
    return (1-(gamma-1)*((phi-phi[0])*2.3)/(gamma*0.02476))**(1/(gamma-1))

x_arr = x_arr_athena
rho_arr = f_analytic_isen(x_arr)
np.savetxt("ana_isen.csv",rho_arr,delimiter=",")

plt.plot(x_arr,rho_arr)
plt.xlabel("r")
plt.ylabel(r"$\rho$")
plt.yscale("log")
plt.show()





# print(rho_arr)

# def f_analytical(rho,pres,i,L,arr_size):
#     x_arr = x_arr_athena #np.linspace(0,L,arr_size)
#     phi_arr = -1/x_arr # 10*x_arr
#     # phi_arr[0] = CubicSpline(x_arr[1:],phi_arr[1:],extrapolate=True)(1.)*2
#     phi_x1 = phi_arr[i+1]
#     phi_x = phi_arr[i]
#     return rho-pres*np.exp(-(phi_x1-phi_x))

# def f_analytical_prime(rho,pres,i,L,arr_size):
#     x_arr = x_arr_athena #np.linspace(0,L,arr_size)
#     phi_arr = -1/x_arr # 10*x_arr
#     # phi_arr[0] = CubicSpline(x_arr[1:],phi_arr[1:],extrapolate=True)(1.)*2
#     phi_x1 = phi_arr[i+1]
#     phi_x = phi_arr[i]
#     return 1

# def Newton_Raphson_numerical(arr_size,L):
#     rho_arr = np.zeros(arr_size)
#     pres_arr = np.zeros(arr_size)
#     rho0 = 1
#     pres0 = 0.011
#     rho_arr[0] = rho0
#     pres_arr[0] = pres0
#     for i in range(arr_size-1):
#         print(f_numerical(rho_arr[i],pres_arr[i],i))        
#         rho_arr[i+1] = rho_arr[i] - f_numerical(rho_arr[i],pres_arr[i],i)/f_numerical_prime(rho_arr[i],pres_arr[i],i)
#         pres_arr[i+1] = (0.011)*rho_arr[i+1]
#     return rho_arr,pres_arr

# def Newton_Raphson_analytical(arr_size,L):
#     rho_arr = np.zeros(arr_size)
#     pres_arr = np.zeros(arr_size)
#     rho0 = 100
#     pres0 = 1
#     rho_arr[0] = rho0
#     pres_arr[0] = pres0
#     for i in range(arr_size-1):        
#         rho_arr[i+1] = rho_arr[i] - f_analytical(rho_arr[i],pres_arr[i],i,L,arr_size)/f_analytical_prime(rho_arr[i],pres_arr[i],i,L,arr_size)
#         pres_arr[i+1] = (0.1**2)*rho_arr[i+1]
#     return rho_arr,pres_arr

# arr_size = len(x_arr_athena)
# # print(Newton_Raphson_numerical(arr_size,2))
# x_arr = x_arr_athena
# plt.plot(x_arr,Newton_Raphson_numerical(arr_size,x_arr[-1]-x_arr[0])[0],label='density')
# plt.plot(x_arr,Newton_Raphson_numerical(arr_size,x_arr[-1]-x_arr[0])[1],label='pressure')
# # plt.yscale('log')
# plt.xlabel("r")
# plt.ylabel(r"$\rho$, P")
# plt.title("Numerical (spherical gravity)")
# plt.legend()
# # plt.yscale("log")
# plt.show()

# arr_size = len(x_arr_athena)
# # print(Newton_Raphson_analytical(arr_size,2))
# x_arr = x_arr_athena
# plt.plot(x_arr,Newton_Raphson_analytical(arr_size,2)[0],label='density')
# plt.plot(x_arr,Newton_Raphson_analytical(arr_size,2)[1],label='pressure')
# # plt.yscale('log')
# plt.xlabel("r")
# plt.ylabel(r"$\rho$, P")
# plt.title("Analytical (spherical gravity)")
# plt.legend()
# # plt.yscale("log")
# plt.show()

# arr_size = len(x_arr_athena)
# # print(Newton_Raphson_analytical(arr_size,2))
# x_arr = x_arr_athena
# plt.plot(x_arr,np.abs(Newton_Raphson_numerical(arr_size,2)[0]-Newton_Raphson_analytical(arr_size,2)[0]),label='density')
# plt.plot(x_arr,np.abs(Newton_Raphson_numerical(arr_size,2)[1]-Newton_Raphson_analytical(arr_size,2)[1]),label='pressure')
# # plt.yscale('log')
# plt.xlabel("r")
# plt.ylabel(r"|$\rho$|, |P|")
# plt.title("Residual plot (spherical gravity)")
# plt.legend()
# # plt.yscale("log")
# plt.xlim([x_arr[2],x_arr[-2]])
# plt.show()

# print((Newton_Raphson_numerical(arr_size,0.00354)[0]).tolist())
# np.savetxt("Newton_Raphson_numerical.csv",Newton_Raphson_numerical(arr_size,0.00354)[0],delimiter=",")

# x_arr_brent = x_arr_athena #np.linspace(0,L,arr_size)
# phi_arr_brent = -1/x_arr_brent # 10*x_arr
# arr_size = len(x_arr_athena)
# L = 2
# rho_arr_brent = np.zeros(arr_size)
# pres_arr_brent = np.zeros(arr_size)
# rho0_brent = 1
# pres0_brent = 1
# rho_arr_brent[0] = rho0_brent
# pres_arr_brent[0] = pres0_brent
# analytical_rho = np.exp(phi_arr_brent[0])*np.exp(-phi_arr_brent)
# for i in range(arr_size-1):        
#     rho_arr_brent[i+1] = brentq(f_numerical,a=analytical_rho-2*10**(-2),b=analytical_rho-10**(-2),args=(pres_arr_brent[i],i,L,arr_size))
#     pres_arr_brent[i+1] = rho_arr_brent[i+1]

# print(rho_arr_brent)

# x_arr = x_arr_athena
# phi_arr = -1/x_arr
# rho_sol = np.exp(-phi_arr)
# plt.plot(x_arr,rho_sol)
# plt.show()

# root = brentq(f_numerical,args=(pres,i=0,L,arr_size))

# def f_numerical(rho,pres,i,L,arr_size):
#     x_arr = x_arr_athena #np.linspace(0,L,arr_size)
#     phi_arr = 10*x_arr 
#     # phi_arr[0] = CubicSpline(x_arr[1:],phi_arr[1:],extrapolate=True)(1.)*2
#     phi_x1 = phi_arr[i+1]
#     phi_x = phi_arr[i]
#     return (rho-pres+(pres+rho)/2*(phi_x1-phi_x))/(arr_size/L)

# def f_numerical_prime(rho,pres,i,L,arr_size):
#     x_arr = x_arr_athena #np.linspace(0,L,arr_size)
#     phi_arr = 10*x_arr
#     # phi_arr[0] = CubicSpline(x_arr[1:],phi_arr[1:],extrapolate=True)(1.)*2
#     phi_x1 = phi_arr[i+1]
#     phi_x = phi_arr[i]
#     return (1+(1/2)*(phi_x1-phi_x))/(arr_size/L)

# def f_analytical(rho,pres,i,L,arr_size):
#     x_arr = x_arr_athena #np.linspace(0,L,arr_size)
#     phi_arr = 10*x_arr
#     # phi_arr[0] = CubicSpline(x_arr[1:],phi_arr[1:],extrapolate=True)(1.)*2
#     phi_x1 = phi_arr[i+1]
#     phi_x = phi_arr[i]
#     return rho-pres*np.exp(-(phi_x1-phi_x))

# def f_analytical_prime(rho,pres,i,L,arr_size):
#     x_arr = x_arr_athena #np.linspace(0,L,arr_size)
#     phi_arr = 10*x_arr
#     # phi_arr[0] = CubicSpline(x_arr[1:],phi_arr[1:],extrapolate=True)(1.)*2
#     phi_x1 = phi_arr[i+1]
#     phi_x = phi_arr[i]
#     return 1

# def Newton_Raphson_numerical(arr_size,L):
#     rho_arr = np.zeros(arr_size)
#     pres_arr = np.zeros(arr_size)
#     rho0 = 1
#     pres0 = 1
#     rho_arr[0] = rho0
#     pres_arr[0] = pres0
#     for i in range(arr_size-1):        
#         rho_arr[i+1] = rho_arr[i] - f_numerical(rho_arr[i],pres_arr[i],i,L,arr_size)/f_numerical_prime(rho_arr[i],pres_arr[i],i,L,arr_size)
#         pres_arr[i+1] = rho_arr[i+1]
#     return rho_arr,pres_arr

# def Newton_Raphson_analytical(arr_size,L):
#     rho_arr = np.zeros(arr_size)
#     pres_arr = np.zeros(arr_size)
#     rho0 = 1
#     pres0 = 1
#     rho_arr[0] = rho0
#     pres_arr[0] = pres0
#     for i in range(arr_size-1):        
#         rho_arr[i+1] = rho_arr[i] - f_analytical(rho_arr[i],pres_arr[i],i,L,arr_size)/f_analytical_prime(rho_arr[i],pres_arr[i],i,L,arr_size)
#         pres_arr[i+1] = rho_arr[i+1]
#     return rho_arr,pres_arr

# arr_size = len(x_arr_athena)
# # print(Newton_Raphson_numerical(arr_size,2))
# x_arr = x_arr_athena
# plt.plot(x_arr,Newton_Raphson_numerical(arr_size,2)[0],label='density')
# plt.plot(x_arr,Newton_Raphson_numerical(arr_size,2)[1],label='pressure')
# # plt.yscale('log')
# plt.xlabel("r")
# plt.ylabel(r"$\rho$, P")
# plt.title("Numerical (constant gravity)")
# plt.legend()
# # plt.yscale("log")
# plt.show()

# arr_size = len(x_arr_athena)
# # print(Newton_Raphson_analytical(arr_size,2))
# x_arr = x_arr_athena
# plt.plot(x_arr,Newton_Raphson_analytical(arr_size,2)[0],label='density')
# plt.plot(x_arr,Newton_Raphson_analytical(arr_size,2)[1],label='pressure')
# # plt.yscale('log')
# plt.xlabel("r")
# plt.ylabel(r"$\rho$, P")
# plt.title("Analytical (constant gravity)")
# plt.legend()
# # plt.yscale("log")
# plt.show()

# arr_size = len(x_arr_athena)
# # print(Newton_Raphson_analytical(arr_size,2))
# x_arr = x_arr_athena
# plt.plot(x_arr,np.abs(Newton_Raphson_numerical(arr_size,2)[0]-Newton_Raphson_analytical(arr_size,2)[0]),label='density')
# plt.plot(x_arr,np.abs(Newton_Raphson_numerical(arr_size,2)[1]-Newton_Raphson_analytical(arr_size,2)[1]),label='pressure')
# # plt.yscale('log')
# plt.xlabel("r")
# plt.ylabel(r"|$\rho$|, |P|")
# plt.title("Residual plot (constant gravity)")
# plt.legend()
# # plt.yscale("log")
# plt.show()

# print((Newton_Raphson_numerical(arr_size,2)[0]).tolist())
