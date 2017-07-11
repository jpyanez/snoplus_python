import numpy as np

greydisc_radius = 137.7

def pmtContactTime( x0, y0, z0, u, v, w):
    # x0, y0, z0 are the starting point
    # u, v, w are the velocity components
    # ae, ce are ellipsoid characteristics
    ae = 203.7/2.
    ce = 150.3/2.
    a = w**2 / ce**2 + (u**2+v**2)/ae**2
    b = 2*((u*x0+v*y0)/ae**2 + w*z0/ce**2)
    c = (x0**2 + y0**2)/ae**2 + z0**2/ce**2 -1
    
    inner_sqrt = b**2 -4*a*c
    #print inner_sqrt
    if inner_sqrt < 0:
        return 0.
    
    sol1 = (-b+np.sqrt(inner_sqrt))/(2*a)
    sol2 = (-b-np.sqrt(inner_sqrt))/(2*a)
    
    if sol1 > 0 and sol2 < 0:
        return sol1
    if sol2 > 0 and sol1 < 0:
        return sol2
    if sol1 > 0 and sol2 > 0:
        return np.min([sol1, sol2])
    if sol1 < 0 and sol2 < 0:
        return 0.

    return
    
    
    

def pmtContactDepth(z0, w, t):
    return z0 + t*w
