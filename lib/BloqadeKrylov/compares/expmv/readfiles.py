import numpy as np 
import scipy as sp

def readCSC(fname,dtype):
	m,n = np.fromfile(fname + ".size",np.int64)
	data = np.fromfile(fname+".nzval",dtype)
	indices = np.fromfile(fname+".rowval",np.int64)
	indptr = np.fromfile(fname+".colptr",np.int64)
	print(m," ",n)
	print(indices)
	print(indptr)
	print(data)
	return sp.sparse.csc_matrix((data,indices-1,indptr-1),shape = (m,n))


O = readCSC("MHt",np.complex128)
#print(np.linalg.norm(O.todense(),ord=1))
#print(sp.sparse.linalg.onenormest(O))
#print(sp.sparse.linalg.norm(O,ord=1))
v = np.arange(O.shape[0])
sp.sparse.linalg.expm_multiply(0.0029j*O,v,traceA=(0.0029j*O).trace())

