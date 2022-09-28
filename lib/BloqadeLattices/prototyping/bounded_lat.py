import numpy as np
import matplotlib.pyplot as plt

#### Bravais lattice and unit cell vectors
theta = np.deg2rad(60)
## BLVs
a1 = np.array([1,0])
a2 = np.array([np.cos(theta),np.sin(theta)])
## Unit cell vectors
t1 = np.zeros(2)
t2 = np.array([0.0,1/np.sqrt(3)])




# generate tile vectors
zero = np.zeros(2)

(m1,m2) = (5,-2)
(n1,n2) = (2,5)

v1 = m1*a1 + m2*a2
v2 = n1*a1 + n2*a2
v3 = v1+v2
# calculat inverse for checks later on
V = np.array([v1,v2]).T
V_inv = np.linalg.inv(V)

# generate bounds  around tile
Xmin = min(m1,n1,0,m1+n1)
Xmax = max(m1,n1,0,m1+n1)
Ymin = min(m2,n2,0,m2+n2)
Ymax = max(m2,n2,0,m2+n2)
# generate grid based on bounds
X,Y = np.mgrid[Xmin:Xmax+1,Ymin:Ymax+1]
# generate actual BL points in 2-D space
# e.g. r = X*a1+Y*a2 
x = np.outer(X.ravel(),a1) + np.outer(Y.ravel(),a2)
# add in unit cells for full lattice
x = np.vstack([x+t1,x+t2])


# bounding for tile, calculate coeffs for each lattice site
c = np.around(x.dot(V_inv.T),14)
# print out truth table for bounds 0 <= c < 1 for both coeff of a1 and a2
print((c>=0.0)*(c<1.0))
# generate masks based on whether point is inside the tile or not
m = np.prod((c>=0.0)*(c<1.0),axis=1).astype(bool)
mm = np.logical_not(m)
# plot the points, color is based on inside or outside of tile
plt.scatter(x[m,0],x[m,1],color="blue")
plt.scatter(x[mm,0],x[mm,1],color="red")
# draw tile boundary
for a,b in [(zero,v1),(zero,v2),(v1,v3),(v2,v3)]:
    plt.plot([a[0],b[0]],[a[1],b[1]],color="black",linestyle=":")
# draw BL Vectors at origin
for a,b in [(zero,a1),(zero,a2)]:
    plt.plot([a[0],b[0]],[a[1],b[1]],color="orange")
# set Aspect ratio to prevent scaling along x/y axes
ax = plt.gca()
ax.set_aspect("equal")
plt.show() # show lattice