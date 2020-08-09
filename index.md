---
layout: default
---

Evanescent electromagnetic waves possess spin-momentum locking, where the direction of propagation (momentum) is locked to the inherent polarization of the wave (transverse spin). We study the optical forces arising from this universal phenomenon and show that the fundamental origin of recently reported non-trivial optical chiral forces is spin-momentum locking. For evanescent waves, we show that the direction of energy flow, the direction of decay, and the direction of spin follow a right hand rule for three different cases of total internal reflection, surface plasmon polaritons, and HE₁₁ mode of an optical fiber. Furthermore, we explain how the recently reported phenomena of lateral optical force on chiral and achiral particles are caused by the transverse spin ofthe evanescent field and the spin-momentum locking phenomenon. Finally, we propose an experiment to identify the unique lateral forces arising from the transverse spin in the optical fiber and point to fundamental differences of the spin density from the well-known orbital angular momentum of light. Our work presents a unified view on spin-momentum locking and how it affects optical forces on chiral and achiral particles.

<p align="center">
    <img src="EField.PNG" alt="EField" width="500" height = "350" align = "middle"/>
</p>
<p align = "center">
    <em>Figure 1: Direction of transverse spin in evanescent fields for (a) total internal reflection</em>
</p>     

<p align = "center">
    <img src="SPP.PNG" alt="SPP" width="500" height = "350" align = "middle"/>
</p>
<p align = "center">
    <em>Figure 2: Direction of transverse spin in evanescent fields for (b) surface plasmon polariton</em>
</p>   

<p align = "center">
    <img src="HEGeneral.PNG" alt="HEGeneral" width="500" height = "350" align = "middle"/>
</p>
<p align = "center">
    <em>Figure 3: Direction of transverse spin in evanescent fields for  (c) optical fiber</em>
</p>

<p align = "center">
    <img src="force1.PNG" alt="force1" width="500" height = "350" align = "middle"/>
</p>
<p align = "center">
    <em>Figure 4: Time averaged optical force in different directions plotted versus radius of the fiber.</em>
</p>   

<p align = "center">
    <img src="force1.PNG" alt="force1" width="500" height = "350" align = "middle"/>
</p>
<p align = "center">
    <em>Figure 5: Time averaged optical force in different directions plotted versus radius of the fiber.</em>
</p>   

<p align = "center">
    <img src="force3.PNG" alt="force2" width="500" height = "350" align = "middle"/>
</p>
<p align = "center">
    <em>Figure 6: Time averaged optical force in different directions plotted versus radius of the fiber.</em>
</p>   

### References

*¹Farid Kalhor,Thomas Thundat,and Zubin Jacob, "Universal spin-momentum locked optical forces" Appl. Phys. Lett. 108, 061102 (2016); https://doi.org/10.1063/1.4941539*
https://aip.scitation.org/doi/10.1063/1.4941539 

### Python Code:

For Figure 1:

```python
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import math
import cmath
from scipy.linalg import expm
```

```python
epsilon_0 = 8.85e-12;
mu_0 = 1.26e-6;
N = 2     # N is the number of layers inclusing the two semi infinite layers so N>1
N_d = 1   # N_d is the resolution for distance of dipole from the surface
N_z = 20  # N_z is the resolution for calculationg the fields
N_x = 20
N_lambda = 1
N_kx = 1
```

```python
# alpha_** are the elements of the polarizability matrix

alpha_ee = 2
alpha_mm = 1
alpha_em = 0.1
alpha_me = alpha_em
```

```python
n = np.ones((N,1),dtype = complex)
# n is the refractive indecies for the N layers

permit = n**2
# permit is a Nx1 matrix of the permittivity of the layers
```


```python

if N-1 == 1:
    z = np.array([np.linspace(0, 0e-9, 2)])
    z = np.array([[z[0][1]]])
else:
    z = np.array([np.linspace(0, 0e-9, N-1)])
    
# Z is the position of the interfaces of the layers
```

```python
z_obs = np.array([np.linspace(-75e-9, 75e-9, N_z)]);  #this matrix shows the z component of the points in which the fields are being calculated
x = np.array([np.linspace(0, 150e-9, N_x)]); 

z_layer = np.ones((z_obs.shape), dtype = int);
```

```python
for i_layer in range(0,z_obs.shape[1]):
    if z_obs[0][i_layer] > 0:
        z_layer[0][i_layer] = z_layer[0][i_layer] + 1
```


```python
c_0 = 3e8;

if N_lambda == 1 :
    lambda_0 = np.array([np.linspace(300e-9, 300e-9, 2)])
    lambda_0 = np.array([[lambda_0[0][1]]])   
else:
    lambda_0 = np.array([np.linspace(300e-9, 300e-9, N_lambda)]) 


w = 2*np.pi*c_0/lambda_0
k_0 = w/c_0;

if N_kx == 1 :
    k_x = np.array([np.linspace(0, (k_0[0]*15*n[0])[0], 2)])
    k_x =  np.array([[k_x[0][1]]])
else:
    k_x = np.array([np.linspace(0, (k_0[0]*15*n[0])[0], N_kx)])
#c_0 is the speed of light in vacuum. w is the angular frequency. 
#lambda_0 is the wavelength in vacuum and k_0 is the wave number in vacuum.
```

```python
# [r, t]= meshgrid(k_x, lambda_0);

r_p = np.zeros((lambda_0.shape[1],k_x.shape[1]), dtype = complex)
t_p = np.zeros((lambda_0.shape[1],k_x.shape[1]), dtype = complex)
R_p = np.zeros((lambda_0.shape[1],k_x.shape[1]), dtype = complex)
T_p = np.zeros((lambda_0.shape[1],k_x.shape[1]), dtype = complex)

#t and r are Transmission and reflection for H and T and R are transmission and reflection for power.
r_s = np.zeros((lambda_0.shape[1],k_x.shape[1]), dtype = complex)
t_s = np.zeros((lambda_0.shape[1],k_x.shape[1]), dtype = complex)
R_s = np.zeros((lambda_0.shape[1],k_x.shape[1]), dtype = complex)
T_s = np.zeros((lambda_0.shape[1],k_x.shape[1]), dtype = complex)
```

```python
H_y = np.zeros((N_z, N_x), dtype = complex);
E_x = np.zeros((N_z, N_x), dtype = complex);
E_z = np.zeros((N_z, N_x), dtype = complex);
S_x = np.zeros((1, N_z), dtype = complex);
S_y = np.zeros((1, N_z), dtype = complex);
S_z = np.zeros((1, N_z), dtype = complex);
Le_x = np.zeros((1, N_z), dtype = complex);
Le_y = np.zeros((1, N_z), dtype = complex);
Le_z = np.zeros((1, N_z), dtype = complex);
Lm_x = np.zeros((1, N_z), dtype = complex);
Lm_y = np.zeros((1, N_z), dtype = complex);
Lm_z = np.zeros((1, N_z), dtype = complex);
U = np.zeros((1, N_z), dtype = complex);
```

```python
F_gf_matrix = np.zeros((N_lambda, N_kx, 3));
F_rp_matrix = np.zeros((N_lambda, N_kx, 3));
F_esd_matrix = np.zeros((N_lambda, N_kx, 3));
F_matrix = np.zeros((N_lambda, N_kx, 3));

ro = r_s;
# ro is W-LDOS defined in "Quantum nanophotonics using HMM" apendix C.

x_plot, z_plot = np.meshgrid(x/lambda_0, z_obs/lambda_0);
#k_x_normalized and lambda_normalized are used for plotting
```

```python
for i_d in range(0,N_d):
    for i_lambda in range(0,N_lambda):
        n[0] = 3.5;
        n[N-1] = 1;
        
        if N_kx == 1:
            k_x = np.array([np.linspace(k_0[i_lambda], k_0[i_lambda]/math.sqrt(2)*n[0], 2)])
            k_x = np.array([k_x[0][1]])
        else:
            k_x = np.array([np.linspace(k_0[i_lambda], k_0[i_lambda]/math.sqrt(2)*n[0], N_kx)])
            
        permit= n**2;
        for i_kx in range(0,N_kx):
            
            k_z = ((np.ones((N, 1),dtype = np.complex_)*(k_0[i_lambda-1]*n))**2 - k_x[i_kx-1]**2)**0.5
            #k_z is a Nx1 matrix of the wavenumbers along the z direction for rach layer
            
            M_p = np.identity(2,dtype = complex);
            for i_layer in range(0,N-1):
                 M_p = np.dot(M_p,np.dot(np.linalg.inv((np.array([[cmath.exp(1j*k_z[i_layer]*z[i_layer]),0],[0,cmath.exp(-1j*k_z[i_layer]*z[i_layer])]])))\
                ,np.dot(np.linalg.inv((np.array([[1, 1],[(k_z[i_layer]/permit[i_layer])[0], (-k_z[i_layer]/permit[i_layer])[0]]])))\
                ,np.dot((np.array([[1, 1],[(k_z[i_layer+1]/permit[i_layer+1])[0], (-k_z[i_layer+1]/permit[i_layer+1])[0]]]))\
                ,(np.array([[cmath.exp(1j*(k_z[i_layer+1]*z[i_layer])[0]), 0],[0, cmath.exp(-1j*(k_z[i_layer+1]*z[i_layer])[0])]]))))))
              
                
            #M is a 2x2 matrix which relates the amplitudes of the waves in the last layer to those of the first one
            M_s = np.identity(2,dtype = complex)
            for i_layer in range(0,N-1):
                 M_s = np.dot(M_s,np.dot(np.linalg.inv((np.array([[cmath.exp(1j*k_z[i_layer]*z[i_layer]),0],[0,cmath.exp(-1j*k_z[i_layer]*z[i_layer])]])))\
                ,np.dot(np.linalg.inv((np.array([[1, 1],[(k_z[i_layer])[0], (-k_z[i_layer])[0]]])))\
                ,np.dot((np.array([[1, 1],[(k_z[i_layer+1])[0], (-k_z[i_layer+1])[0]]]))\
                ,(np.array([[cmath.exp(1j*(k_z[i_layer+1]*z[i_layer])[0]), 0],[0, cmath.exp(-1j*(k_z[i_layer+1]*z[i_layer])[0])]]))))))
            
            
            t_p[i_lambda, i_kx] = (1/M_p[0, 0])*cmath.exp(1j*k_z[i_layer+1]*z[i_layer])
            r_p[i_lambda, i_kx] = M_p[1, 0]/M_p[0, 0]
            T_p[i_lambda, i_kx] = abs(t_p[i_lambda, i_kx])**2/abs(n[N-1])
            R_p[i_lambda, i_kx] = abs(r_p[i_lambda, i_kx])**2 
            
            t_s[i_lambda, i_kx] = (1/M_p[0, 0])*cmath.exp(1j*k_z[i_layer+1]*z[i_layer])
            r_s[i_lambda, i_kx] = M_p[1, 0]/M_p[0, 0]
            T_s[i_lambda, i_kx] = abs(t_s[i_lambda, i_kx])**2/abs(n[N-1])
            R_s[i_lambda, i_kx] = abs(r_s[i_lambda, i_kx])**2
            
            amp = np.zeros((2, N), dtype = complex)
            vamp = np.zeros((2, N), dtype = complex)
            amp[:,0] = np.array([[1 , r_p[0][0]]])
            
            for i_layer in range(1,N):
                amp[:,i_layer] = np.dot(np.dot(np.linalg.inv((np.array([[cmath.exp(1j*k_z[i_layer]*z[i_layer-1]),0],[0,cmath.exp(-1j*k_z[i_layer]*z[i_layer-1])]])))\
                ,np.dot(np.linalg.inv((np.array([[1, 1],[(k_z[i_layer]/permit[i_layer])[0], (-k_z[i_layer]/permit[i_layer])[0]]])))\
                ,np.dot((np.array([[1, 1],[(k_z[i_layer-1]/permit[i_layer-1])[0], (-k_z[i_layer-1]/permit[i_layer-1])[0]]]))\
                ,(np.array([[cmath.exp(1j*(k_z[i_layer-1]*z[i_layer-1])[0]), 0],[0, cmath.exp(-1j*(k_z[i_layer-1]*z[i_layer-1])[0])]]))))), amp[:,i_layer - 1])
                
            amp[1][-1] = 0    
            
            for i_point in range(0,N_z):
                temp2 = np.array([amp[:,(z_layer[0][i_point] - 1)]]).reshape(2,1)
                temp3 = np.array([np.exp(1j*k_x[i_kx]*x[0])])           
                temp = np.array([[cmath.exp(1j * k_z[z_layer[0][i_point] - 1] * z_obs[0][i_point]), cmath.exp(-1j * k_z[z_layer[0][i_point] - 1] * z_obs[0][i_point])]])
                H_y[i_point,:] = np.dot(np.dot(temp,temp2),temp3)
                
                
                temp4 = np.array([[cmath.exp(1j * k_z[z_layer[0][i_point] - 1] * z_obs[0][i_point]), -cmath.exp(-1j * k_z[z_layer[0][i_point] - 1] * z_obs[0][i_point])]])
                temp5 = np.array([amp[:,(z_layer[0][i_point] - 1)]]).reshape(2,1)
                temp6 = np.array([k_z[(z_layer[0][i_point] - 1)]/permit[(z_layer[0][i_point] - 1)]/epsilon_0**0.5*(mu_0**0.5)/k_0[0][i_lambda]*[np.exp(1j*k_x[i_kx]*x[0])]])
                E_x[i_point,:] = np.dot(np.dot(temp4,temp5),temp6)
                
                temp7 = np.array([[-cmath.exp(1j * k_z[z_layer[0][i_point] - 1] * z_obs[0][i_point]), -cmath.exp(-1j * k_z[z_layer[0][i_point] - 1] * z_obs[0][i_point])]])
                temp8 = np.array([amp[:,(z_layer[0][i_point] - 1)]]).reshape(2,1)
                temp9 = np.array([k_x[i_kx]/permit[(z_layer[0][i_point] - 1)]/epsilon_0**0.5*(mu_0**0.5)/k_0[0][i_lambda]*[np.exp(1j*k_x[i_kx]*x[0])]])
                E_z[i_point,:] = np.dot(np.dot(temp7,temp8),temp9)

```

```python
position = [200,200,380*4/3,380]
plt.figure(figsize = (200,200))
fig, ax = plt.subplots()
q = ax.quiver(x_plot,z_plot,E_x.real,E_z.real,color = "Blue")
plt.axhline(color='r', linestyle='-')
plt.title('Figure 1')
plt.xlabel(r"$ \vec E  $")
plt.ylabel(r"$ Z/\lambda_0  $")
plt.show()
```
For Figure 2:

```python
import numpy as np
import cmath
import math
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt
import scipy.interpolate
from scipy.linalg import expm
```

```python
def Refractive_Silver(lambdaa):
    h = 6.626*1e-34
    c = 3e8
    E = h/1.6e-19*c/lambdaa
    Energy=[0.64,
    0.77,
    0.89,
    1.02,
    1.14,
    1.26,
    1.39,
    1.51,
    1.64,
    1.76,
    1.88,
    2.01,
    2.13,
    2.26,
    2.38,
    2.50,
    2.63,
    2.75,
    2.88,
    3,
    3.12,
    3.25,
    3.37,
    3.5,
    3.62,
    3.74,
    3.87,
    3.99,
    4.12,
    4.24,
    4.36,
    4.49,
    4.61,
    4.74,
    4.86,
    4.98,
    5.11,
    5.23,
    5.36,
    5.48,
    5.6,
    5.73,
    5.85,
    5.98,
    6.1,
    6.22,
    6.35,
    6.47,
    6.6]
    
    Silver_eps=[0.24+1j*14.08,
    0.15+1j*11.85,
    0.13+1j*10.10,
    0.09+1j*8.828,
    0.04+1j*7.795,
    0.04+1j*6.992,
    0.04+1j*6.312,
    0.04+1j*5.727,
    0.03+1j*5.242,
    0.04+1j*4.838,
    0.05+1j*4.483,
    0.06+1j*4.152,
    0.05+1j*3.858,
    0.06+1j*3.586,
    0.05+1j*3.324,
    0.05+1j*3.093,
    0.05+1j*2.869,
    0.04+1j*2.657,
    0.04+1j*2.462,
    0.05+1j*2.275,
    0.05+1j*2.07,
    0.05+1j*1.864,
    0.07+1j*1.657,
    0.10+1j*1.419,
    0.14+1j*1.142,
    0.17+1j*0.829,
    0.81+1j*0.392,
    1.13+1j*0.616,
    1.34+1j*0.964,
    1.39+1j*1.161,
    1.41+1j*1.264,
    1.41+1j*1.331,
    1.38+1j*1.372,
    1.35+1j*1.387,
    1.33+1j*1.393,
    1.31+1j*1.389,
    1.30+1j*1.378,
    1.28+1j*1.367,
    1.28+1j*1.357,
    1.26+1j*1.344,
    1.25+1j*1.342,
    1.22+1j*1.336,
    1.20+1j*1.325,
    1.18+1j*1.312,
    1.15+1j*1.296,
    1.14+1j*1.277,
    1.12+1j*1.255,
    1.10+1j*1.232,
    1.07+1j*1.212]
    
    Energy = np.array([Energy])
    Silver_eps = np.array([Silver_eps])
  
    eps = np.interp(E,Energy[0],Silver_eps[0]**2)
    return eps
```

```python
epsilon_0 = 8.85e-12
mu_0 = 1.26e-6

N = 2     #N is the number of layers inclusing the two semi infinite layers so N>1
N_d = 1   #N_d is the resolution for distance of dipole from the surface
N_z = 20  #N_z is the resolution for calculationg the fields
N_x = 20
N_lambda = 1  
N_kx = 6000  
```

```python
n = np.ones((N, 1),dtype = complex) #ok
#n is the refractive indecies for the N layers

permit = n**2

if N-1 == 1:
    z = np.array([np.linspace(0, 0e-9, 2)])
    z = np.array([[z[0][1]]])
else:
    z = np.array([np.linspace(0, 0e-9, N-1)])
    
z_obs = np.array([np.linspace(-300e-9, 300e-9, N_z)])
# this matrix shows the z component of the points in which the fields are being calculated

x = np.array([np.linspace(0,600e-9,N_x)])
```

```python
z_layer = np.ones((1,20), dtype = int)

for i_layer in range(0,z_obs.shape[1]):
    if z_obs[0][i_layer] > 0:
        z_layer[0][i_layer] = z_layer[0][i_layer] + 1
        
 #d is the distance of the dipole from the surface.
if N_d == 1 :
    d = np.array([np.linspace(0,300*1e-9, 2)])
    d = np.array([[d[0][1]]])   
else:
    d = np.array([np.linspace(0,300*1e-9, N_d)]) 
```

```python
c_0 = 3e8;

if N_lambda == 1 :
    lambda_0 = np.array([np.linspace(200e-9, 500e-9, 2)])
    lambda_0 = np.array([[lambda_0[0][1]]])   
else:
    lambda_0 = np.array([np.linspace(200e-9, 500e-9, N_lambda)]) 

w = 2*np.pi*c_0/lambda_0
k_0 = w/c_0;

if N_kx == 1 :
    k_x = np.array([np.linspace(0, (k_0[0]*15*n[0])[0], 2)])
    k_x =  np.array([[k_x[0][1]]])
else:
    k_x = np.array([np.linspace(0, (k_0[0]*15*n[0])[0], N_kx)])
    
# k_x is the wavenumber in x direction and theta is the angle of incidence    
```

```python
r_p = np.zeros((lambda_0.shape[1],k_x.shape[1]), dtype = complex )
t_p = np.zeros((lambda_0.shape[1],k_x.shape[1]), dtype = complex )
R_p = np.zeros((lambda_0.shape[1],k_x.shape[1]), dtype = complex )
T_p = np.zeros((lambda_0.shape[1],k_x.shape[1]), dtype = complex )
# t and r are Transmission and reflection for H and T and R are transmission and reflection for power.

r_s = np.zeros((lambda_0.shape[1],k_x.shape[1]), dtype = complex )
t_s = np.zeros((lambda_0.shape[1],k_x.shape[1]), dtype = complex )
R_s = np.zeros((lambda_0.shape[1],k_x.shape[1]), dtype = complex )
T_s = np.zeros((lambda_0.shape[1],k_x.shape[1]), dtype = complex )

H_y = np.zeros((N_z,N_x),dtype = complex)
E_z = np.zeros((N_z,N_x),dtype = complex)
E_x = np.zeros((N_z,N_x),dtype = complex)
S_x = np.zeros((1,N_z),dtype = complex)
S_y = np.zeros((1,N_z),dtype = complex)
S_z = np.zeros((1,N_z),dtype = complex)
Le_x = np.zeros((1,N_z),dtype = complex)
Le_y = np.zeros((1,N_z),dtype = complex)
Le_z = np.zeros((1,N_z),dtype = complex)
Lm_x = np.zeros((1,N_z),dtype = complex)
Lm_y = np.zeros((1,N_z),dtype = complex)
Lm_z = np.zeros((1,N_z),dtype = complex)
U = np.zeros((1,N_z),dtype = complex)

F_gf_matrix = np.zeros((N_lambda, N_kx, 3), dtype = complex)
F_rp_matrix = np.zeros((N_lambda, N_kx, 3), dtype = complex)
F_esd_matrix = np.zeros((N_lambda, N_kx, 3), dtype = complex)
F_matrix = np.zeros((N_lambda, N_kx, 3), dtype = complex)
```


```python
ro = np.zeros((lambda_0.shape[1],k_x.shape[1]), dtype = complex )
# ro is W-LDOS defined in "Quantum nanophotonics using HMM" apendix C.
x_plot,z_plot = np.meshgrid(x/lambda_0, z_obs/lambda_0)

kx_plot,lambda_plot = np.meshgrid(k_x, lambda_0)
# k_x_normalized and lambda_normalized are used for plotting
```


```python
for i_d in range(0,N_lambda):
    for i_lambda in range(0,N_lambda):
        n[0] = 1;
        n[1] = Refractive_Silver(lambda_0[0][i_lambda])**0.5
        eps_silver = Refractive_Silver(lambda_0[0][i_lambda])
        k_x = np.linspace(k_0[i_lambda]*0,k_0[i_lambda]*6*n[0],N_kx)
        kxplot = k_x
        kx_normal_plot = k_x/k_0[i_lambda]
        permit = n**2
        
        for i_kx in range(0,N_kx):
            k_z = ((np.ones((N,1))*k_0[i_lambda]*n)**2 - k_x[i_kx]**2)**0.5
            M_p = np.identity(2,dtype = complex)
            for i_layer in range(0,N-1):
                 M_p = np.dot(M_p,np.dot(np.linalg.inv((np.array([[cmath.exp(1j*k_z[i_layer]*z[i_layer]),0],[0,cmath.exp(-1j*k_z[i_layer]*z[i_layer])]])))\
                ,np.dot(np.linalg.inv((np.array([[1, 1],[(k_z[i_layer]/permit[i_layer])[0], (-k_z[i_layer]/permit[i_layer])[0]]])))\
                ,np.dot((np.array([[1, 1],[(k_z[i_layer+1]/permit[i_layer+1])[0], (-k_z[i_layer+1]/permit[i_layer+1])[0]]]))\
                ,(np.array([[cmath.exp(1j*(k_z[i_layer+1]*z[i_layer])[0]), 0],[0, cmath.exp(-1j*(k_z[i_layer+1]*z[i_layer])[0])]]))))))
                
            M_s = np.identity(2,dtype = complex)
            for i_layer in range(0,N-1):
                 M_s = np.dot(M_s,np.dot(np.linalg.inv((np.array([[cmath.exp(1j*k_z[i_layer]*z[i_layer]),0],[0,cmath.exp(-1j*k_z[i_layer]*z[i_layer])]])))\
                ,np.dot(np.linalg.inv((np.array([[1, 1],[(k_z[i_layer])[0], (-k_z[i_layer])[0]]])))\
                ,np.dot((np.array([[1, 1],[(k_z[i_layer+1])[0], (-k_z[i_layer+1])[0]]]))\
                ,(np.array([[cmath.exp(1j*(k_z[i_layer+1]*z[i_layer])[0]), 0],[0, cmath.exp(-1j*(k_z[i_layer+1]*z[i_layer])[0])]]))))))
                
            t_p[i_lambda, i_kx]=   1/M_p[0, 0]*cmath.exp(1j*k_z[i_layer+1]*z[i_layer]);
            r_p[i_lambda, i_kx]=   M_p[1, 0]/M_p[0, 0];
            T_p[i_lambda, i_kx]=   abs(t_p[i_lambda, i_kx])**2/abs(n[N-1]);
            R_p[i_lambda, i_kx]=   abs(r_p[i_lambda, i_kx])**2;    
            
            t_s[i_lambda, i_kx]=   1/M_s[0, 0]*cmath.exp(1j*k_z[i_layer+1]*z[i_layer]);
            r_s[i_lambda, i_kx]=   M_s[1, 0]/M_p[0, 0];
            T_s[i_lambda, i_kx]=   abs(t_s[i_lambda, i_kx])**2/abs(n[N-1]);
            R_s[i_lambda, i_kx]=   abs(r_s[i_lambda, i_kx])**2; 
            
            
        k_x = k_x.reshape(k_x.shape[1],k_x.shape[0])    
        max_val = np.amax(T_p[0])
        ind = np.where(T_p == max_val)
        beta = k_x[0][ind[1]]
        
        k_z = ((np.ones((N, 1),dtype = complex)*k_0[i_lambda]*n)**2 - beta**2)**0.5;
        
        M_p = np.identity(2,dtype = complex)
        
        for i_layer in range(0,N-1):
            M_p = np.dot(M_p,np.dot(np.linalg.inv((np.array([[cmath.exp(1j*k_z[i_layer]*z[i_layer]),0],[0,cmath.exp(-1j*k_z[i_layer]*z[i_layer])]])))\
            ,np.dot(np.linalg.inv((np.array([[1, 1],[(k_z[i_layer]/permit[i_layer])[0], (-k_z[i_layer]/permit[i_layer])[0]]])))\
            ,np.dot((np.array([[1, 1],[(k_z[i_layer+1]/permit[i_layer+1])[0], (-k_z[i_layer+1]/permit[i_layer+1])[0]]]))\
            ,(np.array([[cmath.exp(1j*(k_z[i_layer+1]*z[i_layer])[0]), 0],[0, cmath.exp(-1j*(k_z[i_layer+1]*z[i_layer])[0])]]))))))
                        
        
        amp = np.zeros((2,N),dtype = complex)
        amp[:,0]= [1,r_p[i_lambda,ind[1]]]
        
        for i_layer in range(1,N):
            amp[:,i_layer] = np.dot(np.dot(np.linalg.inv((np.array([[cmath.exp(1j*k_z[i_layer]*z[i_layer-1]),0],[0,cmath.exp(-1j*k_z[i_layer]*z[i_layer-1])]])))\
            ,np.dot(np.linalg.inv((np.array([[1, 1],[(k_z[i_layer]/permit[i_layer])[0], (-k_z[i_layer]/permit[i_layer])[0]]])))\
            ,np.dot((np.array([[1, 1],[(k_z[i_layer-1]/permit[i_layer-1])[0], (-k_z[i_layer-1]/permit[i_layer-1])[0]]]))\
            ,(np.array([[cmath.exp(1j*(k_z[i_layer-1]*z[i_layer-1])[0]), 0],[0, cmath.exp(-1j*(k_z[i_layer-1]*z[i_layer-1])[0])]]))))), amp[:,i_layer - 1])
            
        amp[1][-1] = 0
        
        k_x = k_x.reshape(k_x.shape[1],k_x.shape[0])
        
        for i_point in range(0,N_z):
            temp = np.array([[cmath.exp(1j * k_z[z_layer[0][i_point] - 1] * z_obs[0][i_point]), cmath.exp(-1j * k_z[z_layer[0][i_point] - 1] * z_obs[0][i_point])]])
            temp2 = np.array([amp[:,(z_layer[0][i_point] - 1)]]).reshape(2,1)
            temp3 = np.array([np.exp(1j*beta*x[0])])           
            H_y[i_point,:] = np.dot(np.dot(temp,temp2),temp3)
                
            temp4 = np.array([[cmath.exp(1j * k_z[z_layer[0][i_point] - 1] * z_obs[0][i_point]), -cmath.exp(-1j * k_z[z_layer[0][i_point] - 1] * z_obs[0][i_point])]])
            temp5 = np.array([amp[:,(z_layer[0][i_point] - 1)]]).reshape(2,1)
            temp6 = np.array(k_z[(z_layer[0][i_point] - 1)]/permit[(z_layer[0][i_point] - 1)]/(epsilon_0**0.5)*(mu_0**0.5)/k_0[0][i_lambda]*[np.exp(1j*beta*x[0])])
            E_x[i_point,:] = np.dot(np.dot(temp4,temp5),temp6)
                
            temp7 = np.array([[-cmath.exp(1j * k_z[z_layer[0][i_point] - 1] * z_obs[0][i_point]), -cmath.exp(-1j * k_z[z_layer[0][i_point] - 1] * z_obs[0][i_point])]])
            temp8 = np.array([amp[:,(z_layer[0][i_point] - 1)]]).reshape(2,1)
            temp9 = np.array(beta/permit[(z_layer[0][i_point] - 1)]/(epsilon_0**0.5)*(mu_0**0.5)/k_0[0][i_lambda]*[np.exp(1j*beta*x[0])])
            E_z[i_point,:] = np.dot(np.dot(temp7,temp8),temp9)        
```


```python
plt.figure(figsize = (200,200))
fig, ax = plt.subplots()
q = ax.quiver(x_plot,z_plot,E_x.real,E_z.real,color = "Blue", scale = None , scale_units = 'inches')
plt.axhline(color='r', linestyle='-')
plt.title('Figure 1')
plt.xlabel(r"$ x/\lambda_0 $")
plt.ylabel(r"$ Z/\lambda_0  $")
plt.savefig('SPP.pdf')
plt.show()
```

For Figure 3:

```python
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
from scipy import special,linalg
import cmath
import math
```


```python
global n1,n2nc0,eps0,mu0,a,m,w,k0

n1 = 1.5
n2 = 1.0
c0 = 3e8
eps0 = 8.854187e-12
mu0 = 1.256637e-6
a = 5e-6
m = 1
```

```python
def tot_power(A,B,C,D,beta,p,q,jj):
    global n1,n2,c0,eps0,mu0,a,m,w,k0
    res_r = 800
    r = np.array([np.linspace(0.01*a,80*a,res_r)])
    dr = (r[0][len(r[0]) - 1]-r[0][0])/(res_r+1)
    power = 0
    for j in range(0,res_r):
        J1 = sp.special.jv(m, r[0][j]*p)
        K1 = sp.special.kv(m, r[0][j]*q)
        dJ1 = 0.5*(sp.special.jv(m-1,r[0][j]*p) - sp.special.jv(m+1,r[0][j]*p))
        dK1 = -0.5*(sp.special.kv(m-1,r[0][j]*q) + sp.special.kv(m+1,r[0][j]*q))
        point_r = r[0][j]
        point_phi = 0
        point_z = 0
        
        EH_array = field_calc(A,B,C,D,beta,p,q,J1,dJ1,K1,dK1,point_r,point_phi,point_z)
        
        E_r = EH_array[0]
        E_phi = EH_array[1]
        E_z = EH_array[2]
        H_r = EH_array[3]
        H_phi = EH_array[4]
        H_z = EH_array[5]
        
        power = power + 2*np.pi*r[0][j]*dr*(1/2)*(E_r*np.conj(H_phi)-E_phi*np.conj(H_r))
        
    return power
```

```python
def field_calc( A, B, C, D, beta, p, q, J1, dJ1, K1, dK1, point_r,point_phi,point_z ):
    global n1,n2,c0,eps0,mu0,a,m,w,k0

    E_r =  A*1j*beta/p[0][0]*dJ1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) - B[0][0]*mu0*w[0][0]*m/p[0][0]**2/point_r*J1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)
    
    E_phi=  -A*beta*m/p**2/point_r*J1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) - B*1j*mu0*w/p*dJ1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)

    E_z=    A*J1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z);

    H_r=    A*eps0*n1**2*w*m/p**2/point_r*J1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) + B*1j*beta/p*dJ1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)

    H_phi=  A*1j*eps0*n1**2*w/p*dJ1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) - B*beta*m/p**2/point_r*J1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)

    H_z=    B*J1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)
    
    if(type(point_r) == int):
        if(point_r > a):
            E_r = -C*1j*beta/q*dK1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) + D*mu0*w*m/q**2/point_r*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)
            E_phi = C*beta*m/q**2/point_r*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) + D*1j*mu0*w/q*dK1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)
            E_z= C*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)
            H_r = -C*eps0*n2**2*w*m/q**2/point_r*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) - D*1j*beta/q*dK1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)
            H_phi = -C*1j*eps0*n2**2*w/q*dK1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) + D*beta*m/q**2/point_r*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)
            H_z = D*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)
        else:
            E_r = E_r
            E_phi = E_phi
            E_z = E_z
            H_r = H_r
            H_phi = H_phi
            H_z = H_z
            
    else:
        dummy = (point_r > a).astype(int)
        dummy_inv = (point_r < a).astype(int)
        dummy = dummy*(-C*1j*beta/q*dK1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) + D*mu0*w*m/q**2/point_r*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z))
        dummy_inv = E_r * dummy_inv
        E_r = dummy + dummy_inv
        
        dummy = (point_r > a).astype(int)
        dummy_inv = (point_r < a).astype(int)
        dummy = dummy*(C*beta*m/q**2/point_r*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) + D*1j*mu0*w/q*dK1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z))
        dummy_inv = E_phi * dummy_inv
        E_phi = dummy + dummy_inv
        
        dummy = (point_r > a).astype(int)
        dummy_inv = (point_r < a).astype(int)
        dummy = dummy*(C*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z))
        dummy_inv = E_z * dummy_inv
        E_z = dummy + dummy_inv
        
        dummy = (point_r > a).astype(int)
        dummy_inv = (point_r < a).astype(int)
        dummy = dummy*(-C*eps0*n2**2*w*m/q**2/point_r*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) - D*1j*beta/q*dK1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z))
        dummy_inv = H_r * dummy_inv
        H_r = dummy + dummy_inv
        
        dummy = (point_r > a).astype(int)
        dummy_inv = (point_r < a).astype(int)
        dummy = dummy*(-C*1j*eps0*n2**2*w/q*dK1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) + D*beta*m/q**2/point_r*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z))
        dummy_inv = H_phi * dummy_inv
        H_phi = dummy + dummy_inv
        
        
        dummy = (point_r > a).astype(int)
        dummy_inv = (point_r < a).astype(int)
        dummy = dummy*(D*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z))
        dummy_inv = H_z * dummy_inv
        H_z = dummy + dummy_inv
        
        
        
    return [E_r,E_phi,E_z,H_r,H_phi,H_z,point_r,point_phi,point_z]
```

```python
alpha_ee = 1*(10e-9)**3
alpha_em = 0.01*(10e-9)**3
alpha_mm = 0.1*(10e-9)**3

res_r = 1
res_w = 1
res_beta = 20000
```

```python
if res_w == 1:
    mat_w = np.array([np.linspace(0.05e15, 0.1e15, 2)])
    mat_w = np.array([[mat_w[0][1]]])
else:
    mat_w = np.array([np.linspace(0.05e15, 0.1e15, res_w)])
    
mat_k0 = mat_w/c0

F1 = np.zeros((3,mat_w.shape[1]),dtype = complex)
F2 = np.zeros((3,mat_w.shape[1]),dtype = complex)
F3 = np.zeros((3,mat_w.shape[1]),dtype = complex)
F4 = np.zeros((3,mat_w.shape[1]),dtype = complex)
F5 = np.zeros((3,mat_w.shape[1]),dtype = complex)
F6 = np.zeros((3,mat_w.shape[1]),dtype = complex)
F7 = np.zeros((3,mat_w.shape[1]),dtype = complex)
F8 = np.zeros((3,mat_w.shape[1]),dtype = complex)
F_tot = np.zeros((3,mat_w.shape[1]), dtype= complex)
```

```python
for jj in range(0,mat_w.shape[1]):
    w = np.array([mat_w[jj]])
    lambda_matrix = np.array([2*np.pi*c0*w])
    k0 = w/c0
    mat_V = k0*a*(n1**2 - n2**2)**0.5
    
    beta = np.array([np.linspace(n2*k0[0][0]*1.0000001 , n1*k0[0][0]*0.999999, res_beta)])
    p = (n1**2*k0**2 -beta**2)**0.5
    q = (beta**2 - n2**2*k0**2)**0.5
    J1 = sp.special.jv(m,a*p)
    K1 = sp.special.kv(m,a*q)
    dJ1 = 0.5*(sp.special.jv(m-1,a*p) - sp.special.jv(m+1,a*p))
    dK1 = -0.5*(sp.special.kv(m-1,a*q) + sp.special.kv(m+1,a*q))
    
    error = (dJ1/p/J1 + dK1/q/K1)*(dJ1/p/J1 + n2**2/n1**2*dK1/q/K1) - m**2/a**2*(1/p**2 + 1/q**2)*(1/p**2 + n2**2/n1**2/q**2)
    
    error_min = np.amin(abs(error))
    error_min_index = np.where(abs(error) == error_min)
    out_beta = beta[(error_min_index[0][0])][(error_min_index[1][0])]
    
    beta = 0
    p = 0
    q = 0
    J1 = 0
    K1 = 0
    dJ1 = 0
    dK1 = 0
    
        
    beta = out_beta
    beta_normal = out_beta
    mat_b = (beta/k0 - n2)/(n1-n2)
    w = mat_w
    p = (n1**2*k0**2 -beta**2)**0.5
    q = (beta**2 - n2**2*k0**2)**0.5
    J1 = sp.special.jv(m,a*p)
    K1 = sp.special.kv(m,a*q)
    dJ1 = 0.5*(sp.special.jv(m-1,a*p) - sp.special.jv(m+1,a*p))
    dK1 = -0.5*(sp.special.kv(m-1,a*q) + sp.special.kv(m+1,a*q))
    
    A = 1
    B = -A*J1*beta*m*(1/p**2 + 1/q**2)/a/1j/w/mu0/(J1/K1*dK1/q+dJ1/p)
    C = A*J1/K1
    D = B*J1/K1
    
    mat_B = B
    mat_c = C
    mat_D = D
    
    power = tot_power( A, B, C, D, beta, p, q, jj );
    
    J1 = 0
    K1 = 0
    dJ1 = 0
    dK1 = 0
    
    r = 1.0*a+50e-9
    
    A = A/np.sqrt(power)
    B = B/np.sqrt(power)
    C = C/np.sqrt(power)
    D = D/np.sqrt(power)
    
    K1 = sp.special.kn(m,r*q)
    dK1 = -0.5*(sp.special.kn(m-1, r*q) + sp.special.kn(m+1,r*q))
    ddK1 = 0.25*(sp.special.kn(m-2,r*q) + 2*sp.special.kn(m,r*q) + sp.special.kn(m+2 ,r*q))
    
    F1[0][jj] = 1/4*(alpha_ee.real*( 2*abs(C)**2*q*K1*dK1+2*(abs(C)**2*beta**2/q+abs(D)**2*mu0**2*w**2/q)*dK1*ddK1+2*(abs(D)**2*mu0**2*w**2*m**2/q**4+abs(C)**2*beta**2*m**2/q**4)*(-K1**2/r**3+q*K1*dK1/r**2)+(-2*C*np.conj(D)+2*np.conj(C)*D)*1j*beta*mu0*w*m/q**3*(-K1*dK1/r**2+q*dK1**2/r+q*K1*ddK1/r))\
                + alpha_mm.real*( 2*abs(D)**2*q*K1*dK1+2*(abs(D)**2*beta**2/q+abs(C)**2*eps0**2*n2**4*w**2/q)*dK1*ddK1+2*(abs(C)**2*eps0**2*n2**4*w**2*m**2/q**4+abs(D)**2*beta**2*m**2/q**4)*(-K1**2/r**3+q*K1*dK1/r**2)+(-2*C*np.conj(D)+2*np.conj(C)*D)*1j*beta*eps0*n2**2*w*m/q**3*(-K1*dK1/r**2+q*dK1**2/r+q*K1*ddK1/r))\
                - 2*alpha_em.real*((-2*abs(C)**2*eps0*n2**2-2*abs(D)**2*mu0)*1j*beta*w*m/q**3*(-K1*dK1/r**2+q*dK1**2/r+q*K1*ddK1/r)+2*(np.conj(C)*D*beta**2/q**2-C*np.conj(D)*mu0*w**2*eps0*n2**2/q**2)*q*dK1*ddK1+2*(-C*np.conj(D)*mu0*w**2*eps0*n2**2*m**2/q**4+np.conj(C)*D*beta**2*m**2/q**4)*(-K1**2/r**3+q*K1*dK1/r**2)+2*np.conj(C)*D*q*K1*dK1).imag)
    F1[1][jj] = 0
    F1[2][jj] = 0
    
    F2[0][jj] = 0
    F2[1][jj] = (k0/eps0*(alpha_ee).imag+k0/mu0*(alpha_mm).imag-c0**2*k0**4/6/np.pi*((alpha_ee*np.conj(alpha_mm)).real+abs(alpha_em)**2))/c0*0.5*((-abs(C)**2*eps0*n2**2*w*m/q**2/r-abs(D)**2*mu0*w*m/q**2/r)*K1**2+2*C*np.conj(D)*1j*beta/q*K1*dK1)
    F2[2][jj] = (k0/eps0*(alpha_ee).imag+k0/mu0*(alpha_mm).imag-c0**2*k0**4/6/np.pi*((alpha_ee*np.conj(alpha_mm)).real+abs(alpha_em)**2))/c0*0.5*((abs(C)**2*eps0*n2**2+abs(D)**2*mu0)*beta*w/q**2*(dK1**2+m**2*K1**2/q**2/r**2)+2*(-C*np.conj(D)*beta**2+D*np.conj(C)*mu0*eps0*n2**2*w**2)*1j*m/q**3/r*K1*dK1)
    
    F3[0][jj] = 0
    F3[1][jj] = -(alpha_em).imag*(-1/2)*( 2*(abs(C)**2*eps0*n2**2+abs(D)**2*mu0)*beta*w/q**2*(q*dK1*ddK1+m**2/q/r**2*K1*dK1-m**2*K1**2/q**2/r**3)+2*(-C*np.conj(D)*beta**2+D*np.conj(C)*mu0*eps0*n2**2*w**2)*1j*m/q**3*(q*dK1**2/r+q*K1*ddK1/r-K1*dK1/r**2) )
    F3[2][jj] = -(alpha_em).imag/r*( -(abs(C)**2*eps0*n2**2+abs(D)**2*mu0)*w*m/q*K1*dK1+C*np.conj(D)*1j*beta*(1/q*K1*dK1+r*dK1**2+r*K1*ddK1) )
    
    F4[0][jj] = 0 
    F4[1][jj] = c0*k0/eps0*(alpha_ee).imag*(-eps0)/2/w/1j*( (-abs(C)**2*beta**2-abs(D)**2*mu0**2*w**2)*1j*m/q**3*(q*dK1**2/r+q*K1*ddK1/r-K1*dK1/r**2)-2*C*np.conj(D)*beta*mu0*w/q*dK1*ddK1+2*D*np.conj(C)*mu0*beta*w*m**2/q**4*(q*K1*dK1/r**2-K1**2/r**3) )
    F4[2][jj] = c0*k0/eps0*(alpha_ee).imag*eps0/2/w/1j/r*( abs(C)**2*1j*beta/q*(q*r*dK1**2+q*r*K1*ddK1+K1*dK1)+2*C*np.conj(D)*mu0*w*m/q*K1*dK1 )
    
    F5[0][jj] = 0
    F5[1][jj] = c0*k0/mu0*(alpha_mm).imag*(-mu0)/2/w/1j*( (-abs(C)**2*eps0**2*n2**4*w**2-abs(D)**2*beta**2)*1j*m/q**3*(q*dK1**2/r+q*K1*ddK1/r-K1*dK1/r**2)-2*C*np.conj(D)*eps0*n2**2*w*beta*m**2/q**4*(q*K1*dK1/r**2-K1**2/r**3)+2*D*np.conj(C)*eps0*n2**2*beta*w/q*K1*dK1 )
    F5[2][jj] = c0*k0/mu0*(alpha_mm).imag*mu0/2/w/1j/r*( abs(D)**2*1j*beta/q*(q*r*dK1**2+q*r*K1*ddK1+K1*dK1)-2*D*np.conj(C)*eps0*n2**2*w*m/q*K1*dK1 )
    
    F6[0][jj] = 0
    F6[1][jj] = w*(-2*w*(alpha_em).imag+c0*k0**4/3/np.pi/eps0*(alpha_ee*np.conj(alpha_em)).real)*eps0/2/w/1j*(abs(C)**2*1j*beta/q*K1*dK1+C*np.conj(D)*mu0*w*m/q**2/r*K1**2)
    F6[2][jj] = w*(-2*w*(alpha_em).imag+c0*k0**4/3/np.pi/eps0*(alpha_ee*np.conj(alpha_em)).real)*eps0/2/w/1j*( (-abs(C)**2*beta**2-abs(D)**2*mu0**2*w**2)*1j*m/q**3/r*K1*dK1-C*np.conj(D)*beta*mu0*w/q**2*dK1**2+D*np.conj(C)*beta*mu0*w*m**2/q**4/r**2*K1**2 )
    
    F7[0][jj] = 0
    F7[1][jj] = w*(-2*w*(alpha_em).imag+c0*k0**4/3/np.pi/mu0*(alpha_mm*np.conj(alpha_em)).real)*mu0/2/w/1j*(-D*np.conj(C)*eps0*n2**2*w*m/q**2/r*K1**2+abs(D)**2*1j*beta/q*K1*dK1)
    F7[2][jj] = w*(-2*w*(alpha_em).imag+c0*k0**4/3/np.pi/mu0*(alpha_mm*np.conj(alpha_em)).real)*mu0/2/w/1j*((-abs(C)**2*eps0**2*n2**4*w**2-abs(D)**2*beta**2)*1j*m/q**3/r*K1*dK1-C*np.conj(D)*eps0*n2**2*w*beta*m**2/q**4/r**2*K1**2+D*np.conj(C)*beta*eps0*n2**2*w/q**2*dK1**2 )
    
    F8[0][jj] = c0*k0**4/12/np.pi*(alpha_ee*np.conj(alpha_mm)).imag*(abs(D)**2*mu0-abs(C)**2*eps0*n2**2)*w/q*K1*dK1;
    F8[1][jj] = 0
    F8[2][jj] = 0
    
    F_tot[0][jj] = F1[0, jj]+F2[0, jj]+F3[0, jj]+F4[0, jj]+F5[0, jj]+F6[0, jj]+F7[0, jj]+F8[0, jj]
    F_tot[1][jj] = F1[1, jj]+F2[1, jj]+F3[1, jj]+F4[1, jj]+F5[1, jj]+F6[1, jj]+F7[1, jj]+F8[1, jj]
    F_tot[2][jj] = F1[2, jj]+F2[2, jj]+F3[2, jj]+F4[2, jj]+F5[2, jj]+F6[2, jj]+F7[2, jj]+F8[2, jj]
```

```python
res_long_z = 20
res_long_x = 20
res_trns_x = 20
res_trns_y = 20

x_long = np.linspace(-2*a, 2*a, res_long_x)
z_long = np.linspace(0*a, 4*a, res_long_z)
x_trns = np.linspace(-2*a, 2*a, res_trns_x)
y_trns = np.linspace(-2*a, 2*a, res_trns_y)

mat_long_z,mat_long_x = np.meshgrid(z_long, x_long)
mat_trns_x,mat_trns_y = np.meshgrid(x_trns, y_trns)

point_r = abs(mat_long_x);
point_phi = np.arctan2(0, mat_long_x);
point_z = mat_long_z;

J1 = sp.special.jv(m, point_r*p)
K1 = sp.special.kn(m, point_r*q)
dJ1 = 0.5*(sp.special.jv(m-1, point_r*p)- sp.special.jv(m+1, point_r*p))
dK1 = -0.5*(sp.special.kn(m-1, point_r*q)+sp.special.kn(m+1, point_r*q))

EH_array = field_calc(A,B,C,D,beta,p,q,J1,dJ1,K1,dK1,point_r,point_phi,point_z)

E_long_r = EH_array[0]
E_long_phi = EH_array[1]
E_long_z = EH_array[2]
H_long_r = EH_array[3]
H_long_phi = EH_array[4]
H_long_z = EH_array[5]
```

```python
E_long_x=    E_long_r*np.cos(point_phi)-E_long_phi*np.sin(point_phi)
E_long_y=    E_long_r*np.sin(point_phi)+E_long_phi*np.cos(point_phi)
H_long_x=    H_long_r*np.cos(point_phi)-H_long_phi*np.sin(point_phi)
H_long_y=    H_long_r*np.sin(point_phi)+H_long_phi*np.cos(point_phi)

Le_long_x=   eps0/4/w/1j*(E_long_y*np.conj(E_long_z)-E_long_z*np.conj(E_long_y))
Le_long_y=   eps0/4/w/1j*(E_long_z*np.conj(E_long_x)-E_long_x*np.conj(E_long_z))
Le_long_z=   eps0/4/w/1j*(E_long_x*np.conj(E_long_y)-E_long_y*np.conj(E_long_x))

Lm_long_x =   mu0/4/w/1j*(H_long_y*np.conj(H_long_z)-H_long_z*np.conj(H_long_y))
Lm_long_y =   mu0/4/w/1j*(H_long_z*np.conj(H_long_x)-H_long_x*np.conj(H_long_z))
Lm_long_z =   mu0/4/w/1j*(H_long_x*np.conj(H_long_y)-H_long_y*np.conj(H_long_x))

L_long_x=   Le_long_x+Lm_long_x
L_long_y=   Le_long_y+Lm_long_y
L_long_z=   Le_long_z+Lm_long_z

S_long_x=   1/2*(E_long_y*np.conj(H_long_z)-E_long_z*np.conj(H_long_y))
S_long_y=   1/2*(E_long_z*np.conj(H_long_x)-E_long_x*np.conj(H_long_z))
S_long_z=   1/2*(E_long_x*np.conj(H_long_y)-E_long_y*np.conj(H_long_x))
```

```python
point_r=    abs((mat_trns_x**2+mat_trns_y**2)**0.5)
point_phi=  np.arctan2(mat_trns_y, mat_trns_x)
point_z=    np.zeros((mat_trns_x.shape),dtype = complex)

J1 = sp.special.jv(m, point_r*p)
K1 = sp.special.kn(m, point_r*q)
dJ1 = 0.5*(sp.special.jv(m-1, point_r*p)- sp.special.jv(m+1, point_r*p))
dK1 = -0.5*(sp.special.kn(m-1, point_r*q)+sp.special.kn(m+1, point_r*q))

EH_array = field_calc(A,B,C,D,beta,p,q,J1,dJ1,K1,dK1,point_r,point_phi,point_z)

E_trns_r = EH_array[0]
E_trns_phi =EH_array[1]
E_trns_z = EH_array[2]
H_trns_r =EH_array[3]
H_trns_phi =EH_array[4]
H_trns_z =EH_array[5]
```


```python
E_trns_x=    E_trns_r*np.cos(point_phi)-E_trns_phi*np.sin(point_phi)
E_trns_y=    E_trns_r*np.sin(point_phi)+E_trns_phi*np.cos(point_phi)
H_trns_x=    H_trns_r*np.cos(point_phi)-H_trns_phi*np.sin(point_phi)
H_trns_y=    H_trns_r*np.sin(point_phi)+H_trns_phi*np.cos(point_phi)

Le_trns_x=   eps0/4/w/1j*(E_trns_y*np.conj(E_trns_z)-E_trns_z*np.conj(E_trns_y))
Le_trns_y=   eps0/4/w/1j*(E_trns_z*np.conj(E_trns_x)-E_trns_x*np.conj(E_trns_z))
Le_trns_z=   eps0/4/w/1j*(E_trns_x*np.conj(E_trns_y)-E_trns_y*np.conj(E_trns_x))

Lm_trns_x=   mu0/4/w/1j*(H_trns_y*np.conj(H_trns_z)-H_trns_z*np.conj(H_trns_y))
Lm_trns_y=   mu0/4/w/1j*(H_trns_z*np.conj(H_trns_x)-H_trns_x*np.conj(H_trns_z))
Lm_trns_z=   mu0/4/w/1j*(H_trns_x*np.conj(H_trns_y)-H_trns_y*np.conj(H_trns_x))

L_trns_x=   Le_trns_x+Lm_trns_x
L_trns_y=   Le_trns_y+Lm_trns_y
L_trns_z=   Le_trns_z+Lm_trns_z
```

```python
S_trns_x=   1/2*(E_trns_y*np.conj(H_trns_z)-E_trns_z*np.conj(H_trns_y))
S_trns_y=   1/2*(E_trns_z*np.conj(H_trns_x)-E_trns_x*np.conj(H_trns_z))
S_trns_z=   1/2*(E_trns_x*np.conj(H_trns_y)-E_trns_y*np.conj(H_trns_x))
```


```python
plt.figure(figsize = (200,200))
fig, ax = plt.subplots()
q = ax.quiver(mat_long_z/a, mat_long_x/a,(E_long_z).real,(E_long_x).real,color = "Blue")
plt.title('E')
plt.xlabel('z/a')
plt.ylabel('x/a')
plt.savefig('hegeneral.pdf')
plt.show()
```

For Figure 4,5 and 6:

```python
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
from scipy import special,linalg,interpolate
import cmath
import math
```
```python
global n1,n2nc0,eps0,mu0,a,m,w,k0

n1 = 1.5
n2 = 1.0
c0 = 3e8
eps0 = 8.854187e-12
mu0 = 1.256637e-6
#a = 5e-6
m = 1
```

```python
def tot_power(A,B,C,D,beta,p,q,jj):
    global n1,n2,c0,eps0,mu0,a,m,w,k0
    res_r = 800
    r = np.array([np.linspace(0.01*a,80*a,res_r)])
    dr = (r[0][len(r[0]) - 1]-r[0][0])/(res_r+1)
    power = 0
    for j in range(0,res_r):
        J1 = sp.special.jv(m, r[0][j]*p)
        K1 = sp.special.kv(m, r[0][j]*q)
        dJ1 = 0.5*(sp.special.jv(m-1,r[0][j]*p) - sp.special.jv(m+1,r[0][j]*p))
        dK1 = -0.5*(sp.special.kv(m-1,r[0][j]*q) + sp.special.kv(m+1,r[0][j]*q))
        point_r = r[0][j]
        point_phi = 0
        point_z = 0
        
        EH_array = field_calc(A,B,C,D,beta,p,q,J1,dJ1,K1,dK1,point_r,point_phi,point_z)
        
        E_r = EH_array[0]
        E_phi = EH_array[1]
        E_z = EH_array[2]
        H_r = EH_array[3]
        H_phi = EH_array[4]
        H_z = EH_array[5]
        
        power = power + 2*np.pi*r[0][j]*dr*(1/2)*(E_r*np.conj(H_phi)-E_phi*np.conj(H_r))
        
    return power
```

```python
def field_calc( A, B, C, D, beta, p, q, J1, dJ1, K1, dK1, point_r,point_phi,point_z ):
    global n1,n2,c0,eps0,mu0,a,m,w,k0

    E_r =  A*1j*beta/p[0][0]*dJ1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) - B[0][0]*mu0*w[0][0]*m/p[0][0]**2/point_r*J1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)
    
    E_phi=  -A*beta*m/p**2/point_r*J1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) - B*1j*mu0*w/p*dJ1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)

    E_z=    A*J1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z);

    H_r=    A*eps0*n1**2*w*m/p**2/point_r*J1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) + B*1j*beta/p*dJ1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)

    H_phi=  A*1j*eps0*n1**2*w/p*dJ1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) - B*beta*m/p**2/point_r*J1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)

    H_z=    B*J1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)
    
    if(type(point_r) == int):
        if(point_r > a):
            E_r = -C*1j*beta/q*dK1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) + D*mu0*w*m/q**2/point_r*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)
            E_phi = C*beta*m/q**2/point_r*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) + D*1j*mu0*w/q*dK1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)
            E_z= C*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)
            H_r = -C*eps0*n2**2*w*m/q**2/point_r*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) - D*1j*beta/q*dK1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)
            H_phi = -C*1j*eps0*n2**2*w/q*dK1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) + D*beta*m/q**2/point_r*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)
            H_z = D*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z)
        else:
            E_r = E_r
            E_phi = E_phi
            E_z = E_z
            H_r = H_r
            H_phi = H_phi
            H_z = H_z
            
    else:
        dummy = (point_r > a).astype(int)
        dummy_inv = (point_r < a).astype(int)
        dummy = dummy*(-C*1j*beta/q*dK1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) + D*mu0*w*m/q**2/point_r*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z))
        dummy_inv = E_r * dummy_inv
        E_r = dummy + dummy_inv
        
        dummy = (point_r > a).astype(int)
        dummy_inv = (point_r < a).astype(int)
        dummy = dummy*(C*beta*m/q**2/point_r*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) + D*1j*mu0*w/q*dK1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z))
        dummy_inv = E_phi * dummy_inv
        E_phi = dummy + dummy_inv
        
        dummy = (point_r > a).astype(int)
        dummy_inv = (point_r < a).astype(int)
        dummy = dummy*(C*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z))
        dummy_inv = E_z * dummy_inv
        E_z = dummy + dummy_inv
        
        dummy = (point_r > a).astype(int)
        dummy_inv = (point_r < a).astype(int)
        dummy = dummy*(-C*eps0*n2**2*w*m/q**2/point_r*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) - D*1j*beta/q*dK1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z))
        dummy_inv = H_r * dummy_inv
        H_r = dummy + dummy_inv
        
        dummy = (point_r > a).astype(int)
        dummy_inv = (point_r < a).astype(int)
        dummy = dummy*(-C*1j*eps0*n2**2*w/q*dK1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z) + D*beta*m/q**2/point_r*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z))
        dummy_inv = H_phi * dummy_inv
        H_phi = dummy + dummy_inv
        
        
        dummy = (point_r > a).astype(int)
        dummy_inv = (point_r < a).astype(int)
        dummy = dummy*(D*K1*np.exp(1j*m*point_phi)*np.exp(1j*beta*point_z))
        dummy_inv = H_z * dummy_inv
        H_z = dummy + dummy_inv
        
        
        
    return [E_r,E_phi,E_z,H_r,H_phi,H_z,point_r,point_phi,point_z]
```

```python
alpha_ee = (1.2+1j*0.1)*(10e-9)**3
alpha_em = 0.01*(10e-9)**3
alpha_mm = 0.0002*(10e-9)**3

res_a = 50
mat_a = np.array([np.linspace(0.7e-6,2e-6,res_a)])
res_w = 1
res_beta = 20000

if res_w == 1:
    mat_w = np.array([np.linspace(0.05e15, 2*np.pi*c0/5e-6, 2)])
    mat_w = np.array([[mat_w[0][1]]])
else:
    mat_w = np.array([np.linspace(0.05e15, 2*np.pi*c0/5e-6, res_w)])
```

```python
mat_k0 = mat_w/c0

F1 = np.zeros((3,mat_a.shape[1]),dtype = complex)
F2 = np.zeros((3,mat_a.shape[1]),dtype = complex)
F3 = np.zeros((3,mat_a.shape[1]),dtype = complex)
F4 = np.zeros((3,mat_a.shape[1]),dtype = complex)
F5 = np.zeros((3,mat_a.shape[1]),dtype = complex)
F6 = np.zeros((3,mat_a.shape[1]),dtype = complex)
F7 = np.zeros((3,mat_a.shape[1]),dtype = complex)
F8 = np.zeros((3,mat_a.shape[1]),dtype = complex)
F_tot = np.zeros((3,mat_a.shape[1]), dtype= complex)

error_min_dummy = np.array([[]])
error_min_index_dummy = np.array([[]],dtype = int)
out_beta = np.array([[]])
```

```python
for jj in range(0,mat_a.shape[1]):
    a = mat_a[0][jj]
    w = mat_w
    lambda_matrix = np.array([2*np.pi*c0*w])
    k0 = w/c0
    mat_V = k0*a*(n1**2 - n2**2)**0.5
    
    beta = np.array([np.linspace(n2*k0[0][0]*1.0000001 , n1*k0[0][0]*0.999999, res_beta)])
    p = (n1**2*k0**2 -beta**2)**0.5
    q = (beta**2 - n2**2*k0**2)**0.5
    J1 = sp.special.jv(m,a*p)
    K1 = sp.special.kv(m,a*q)
    dJ1 = 0.5*(sp.special.jv(m-1,a*p) - sp.special.jv(m+1,a*p))
    dK1 = -0.5*(sp.special.kv(m-1,a*q) + sp.special.kv(m+1,a*q))
    
    error = (dJ1/p/J1 + dK1/q/K1)*(dJ1/p/J1 + n2**2/n1**2*dK1/q/K1) - m**2/a**2*(1/p**2 + 1/q**2)*(1/p**2 + n2**2/n1**2/q**2)
    error_min = np.amin(abs(error))
    error_min_dummy = np.append(error_min_dummy, error_min)
    error_min_index =  np.where(abs(error) == error_min)
    error_min_index_dummy = np.append(error_min_index_dummy,error_min_index[1][0])
    out_beta = np.append(out_beta,beta[0][error_min_index_dummy[jj]])

    J1 = 0
    p =0
    q = 0
    beta = 0
    K1 = 0
    dJ1 = 0
    dK1 = 0
    
    beta = out_beta[jj]
    beta_normal = out_beta[jj] 
    mat_b = (beta/k0 - n2)/(n1-n2)
    
    p = (n1**2*k0**2 -beta**2)**0.5
    q = (beta**2 - n2**2*k0**2)**0.5
    J1 = sp.special.jv(m,a*p)
    K1 = sp.special.kv(m,a*q)
    dJ1 = 0.5*(sp.special.jv(m-1,a*p) - sp.special.jv(m+1,a*p))
    dK1 = -0.5*(sp.special.kv(m-1,a*q) + sp.special.kv(m+1,a*q))
    
    A = 1
    B = -A*J1*beta*m*(1/p**2 + 1/q**2)/a/1j/w/mu0/(np.dot(J1/K1,dK1/q) + dJ1/p)
    C = A*J1/K1
    D = B*J1/K1

    mat_b = B
    mat_c = C
    matD = D
    
    power = tot_power(A,B,C,D,beta,p,q,jj)
    
    J1 = 0
    K1 = 0
    dJ1 = 0
    dK1 = 0
    
    r = 1.0*a+100e-9
    
    A = A/np.sqrt(1e6*power)
    B = B/np.sqrt(1e6*power)
    C = C/np.sqrt(1e6*power)
    D = D/np.sqrt(1e6*power)
    
    K1 = sp.special.kn(m,r*q)
    dK1 = -0.5*(sp.special.kn(m-1, r*q) + sp.special.kn(m+1,r*q))
    ddK1 = 0.25*(sp.special.kn(m-2,r*q) + 2*sp.special.kn(m,r*q) + sp.special.kn(m+2 ,r*q))
    F1[0][jj] = 1/4*(alpha_ee.real*( 2*abs(C)**2*q*K1*dK1+2*(abs(C)**2*beta**2/q+abs(D)**2*mu0**2*w**2/q)*dK1*ddK1+2*(abs(D)**2*mu0**2*w**2*m**2/q**4+abs(C)**2*beta**2*m**2/q**4)*(-K1**2/r**3+q*K1*dK1/r**2)+(-2*C*np.conj(D)+2*np.conj(C)*D)*1j*beta*mu0*w*m/q**3*(-K1*dK1/r**2+q*dK1**2/r+q*K1*ddK1/r))\
                + alpha_mm.real*( 2*abs(D)**2*q*K1*dK1+2*(abs(D)**2*beta**2/q+abs(C)**2*eps0**2*n2**4*w**2/q)*dK1*ddK1+2*(abs(C)**2*eps0**2*n2**4*w**2*m**2/q**4+abs(D)**2*beta**2*m**2/q**4)*(-K1**2/r**3+q*K1*dK1/r**2)+(-2*C*np.conj(D)+2*np.conj(C)*D)*1j*beta*eps0*n2**2*w*m/q**3*(-K1*dK1/r**2+q*dK1**2/r+q*K1*ddK1/r))\
                - 2*alpha_em.real*((-2*abs(C)**2*eps0*n2**2-2*abs(D)**2*mu0)*1j*beta*w*m/q**3*(-K1*dK1/r**2+q*dK1**2/r+q*K1*ddK1/r)+2*(np.conj(C)*D*beta**2/q**2-C*np.conj(D)*mu0*w**2*eps0*n2**2/q**2)*q*dK1*ddK1+2*(-C*np.conj(D)*mu0*w**2*eps0*n2**2*m**2/q**4+np.conj(C)*D*beta**2*m**2/q**4)*(-K1**2/r**3+q*K1*dK1/r**2)+2*np.conj(C)*D*q*K1*dK1).imag)
    F1[1][jj] = 0
    F1[2][jj] = 0
    
    F2[0][jj] = 0
    F2[1][jj] = (k0/eps0*(alpha_ee).imag+k0/mu0*(alpha_mm).imag-c0**2*k0**4/6/np.pi*((alpha_ee*np.conj(alpha_mm)).real+abs(alpha_em)**2))/c0*0.5*((-abs(C)**2*eps0*n2**2*w*m/q**2/r-abs(D)**2*mu0*w*m/q**2/r)*K1**2+2*C*np.conj(D)*1j*beta/q*K1*dK1)
    F2[2][jj] = (k0/eps0*(alpha_ee).imag+k0/mu0*(alpha_mm).imag-c0**2*k0**4/6/np.pi*((alpha_ee*np.conj(alpha_mm)).real+abs(alpha_em)**2))/c0*0.5*((abs(C)**2*eps0*n2**2+abs(D)**2*mu0)*beta*w/q**2*(dK1**2+m**2*K1**2/q**2/r**2)+2*(-C*np.conj(D)*beta**2+D*np.conj(C)*mu0*eps0*n2**2*w**2)*1j*m/q**3/r*K1*dK1)
    
    F3[0][jj] = 0
    F3[1][jj] = -(alpha_em).imag*(-1/2)*( 2*(abs(C)**2*eps0*n2**2+abs(D)**2*mu0)*beta*w/q**2*(q*dK1*ddK1+m**2/q/r**2*K1*dK1-m**2*K1**2/q**2/r**3)+2*(-C*np.conj(D)*beta**2+D*np.conj(C)*mu0*eps0*n2**2*w**2)*1j*m/q**3*(q*dK1**2/r+q*K1*ddK1/r-K1*dK1/r**2) )
    F3[2][jj] = -(alpha_em).imag/r*( -(abs(C)**2*eps0*n2**2+abs(D)**2*mu0)*w*m/q*K1*dK1+C*np.conj(D)*1j*beta*(1/q*K1*dK1+r*dK1**2+r*K1*ddK1) )
    
    F4[0][jj] = 0 
    F4[1][jj] = c0*k0/eps0*(alpha_ee).imag*(-eps0)/2/w/1j*( (-abs(C)**2*beta**2-abs(D)**2*mu0**2*w**2)*1j*m/q**3*(q*dK1**2/r+q*K1*ddK1/r-K1*dK1/r**2)-2*C*np.conj(D)*beta*mu0*w/q*dK1*ddK1+2*D*np.conj(C)*mu0*beta*w*m**2/q**4*(q*K1*dK1/r**2-K1**2/r**3) )
    F4[2][jj] = c0*k0/eps0*(alpha_ee).imag*eps0/2/w/1j/r*( abs(C)**2*1j*beta/q*(q*r*dK1**2+q*r*K1*ddK1+K1*dK1)+2*C*np.conj(D)*mu0*w*m/q*K1*dK1 )
    
    F5[0][jj] = 0
    F5[1][jj] = c0*k0/mu0*(alpha_mm).imag*(-mu0)/2/w/1j*( (-abs(C)**2*eps0**2*n2**4*w**2-abs(D)**2*beta**2)*1j*m/q**3*(q*dK1**2/r+q*K1*ddK1/r-K1*dK1/r**2)-2*C*np.conj(D)*eps0*n2**2*w*beta*m**2/q**4*(q*K1*dK1/r**2-K1**2/r**3)+2*D*np.conj(C)*eps0*n2**2*beta*w/q*K1*dK1 )
    F5[2][jj] = c0*k0/mu0*(alpha_mm).imag*mu0/2/w/1j/r*( abs(D)**2*1j*beta/q*(q*r*dK1**2+q*r*K1*ddK1+K1*dK1)-2*D*np.conj(C)*eps0*n2**2*w*m/q*K1*dK1 )
    
    F6[0][jj] = 0
    F6[1][jj] = w*(-2*w*(alpha_em).imag+c0*k0**4/3/np.pi/eps0*(alpha_ee*np.conj(alpha_em)).real)*eps0/2/w/1j*(abs(C)**2*1j*beta/q*K1*dK1+C*np.conj(D)*mu0*w*m/q**2/r*K1**2)
    F6[2][jj] = w*(-2*w*(alpha_em).imag+c0*k0**4/3/np.pi/eps0*(alpha_ee*np.conj(alpha_em)).real)*eps0/2/w/1j*( (-abs(C)**2*beta**2-abs(D)**2*mu0**2*w**2)*1j*m/q**3/r*K1*dK1-C*np.conj(D)*beta*mu0*w/q**2*dK1**2+D*np.conj(C)*beta*mu0*w*m**2/q**4/r**2*K1**2 )
    
    F7[0][jj] = 0
    F7[1][jj] = w*(-2*w*(alpha_em).imag+c0*k0**4/3/np.pi/mu0*(alpha_mm*np.conj(alpha_em)).real)*mu0/2/w/1j*(-D*np.conj(C)*eps0*n2**2*w*m/q**2/r*K1**2+abs(D)**2*1j*beta/q*K1*dK1)
    F7[2][jj] = w*(-2*w*(alpha_em).imag+c0*k0**4/3/np.pi/mu0*(alpha_mm*np.conj(alpha_em)).real)*mu0/2/w/1j*((-abs(C)**2*eps0**2*n2**4*w**2-abs(D)**2*beta**2)*1j*m/q**3/r*K1*dK1-C*np.conj(D)*eps0*n2**2*w*beta*m**2/q**4/r**2*K1**2+D*np.conj(C)*beta*eps0*n2**2*w/q**2*dK1**2 )
    
    F8[0][jj] = c0*k0**4/12/np.pi*(alpha_ee*np.conj(alpha_mm)).imag*(abs(D)**2*mu0-abs(C)**2*eps0*n2**2)*w/q*K1*dK1;
    F8[1][jj] = 0
    F8[2][jj] = 0
    
    F_tot[0][jj] = F1[0, jj]+F2[0, jj]+F3[0, jj]+F4[0, jj]+F5[0, jj]+F6[0, jj]+F7[0, jj]+F8[0, jj]
    F_tot[1][jj] = F1[1, jj]+F2[1, jj]+F3[1, jj]+F4[1, jj]+F5[1, jj]+F6[1, jj]+F7[1, jj]+F8[1, jj]
    F_tot[2][jj] = F1[2, jj]+F2[2, jj]+F3[2, jj]+F4[2, jj]+F5[2, jj]+F6[2, jj]+F7[2, jj]+F8[2, jj]
```


```python
F0_plot = F_tot[0,:].reshape(1,50)
F1_plot = F_tot[1,:].reshape(1,50)
F2_plot = F_tot[2,:].reshape(1,50)

plt.figure(figsize = (200,200))
fig, ax = plt.subplots()
q = ax.scatter(1e6*mat_a,1e12*(F0_plot),color = "Blue")
plt.title('F_r')
plt.xlabel('$ a(\mu m)$')
plt.ylabel('$F_r(pN)$')
plt.savefig('force1.pdf')
plt.show()
```


```python
plt.figure(figsize = (200,200))
fig, ax = plt.subplots()
q = ax.scatter(1e6*mat_a,1e12*(F1_plot),color = "Blue")
plt.title('$F_\phi $')
plt.xlabel('$a(\mu m)$')
plt.ylabel('$F_\phi (pN)$')
plt.savefig('force2.pdf')
plt.show()
```


```python
plt.figure(figsize = (200,200))
fig, ax = plt.subplots()
q = ax.scatter(1e6*mat_a,1e12*(F2_plot),color = "Blue")
plt.title('F_z')
plt.xlabel('$a(\mu m)$')
plt.ylabel('F_z (pN)')
plt.savefig('force3.pdf')
plt.show()
```
