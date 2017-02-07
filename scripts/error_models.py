#!/usr/bin/python

import numpy as np
import sys

d = np.loadtxt(sys.argv[1], skiprows=3) / 100
a = d[:,7]
t = d[:,8]
g = d[:,9]
c = d[:,10]
x = d[:,12]
A = np.array([np.flatnonzero(a), np.ones(np.count_nonzero(a))])
T = np.array([np.flatnonzero(t), np.ones(np.count_nonzero(t))])
G = np.array([np.flatnonzero(g), np.ones(np.count_nonzero(g))])
C = np.array([np.flatnonzero(c), np.ones(np.count_nonzero(c))])
X = np.array([np.flatnonzero(x), np.ones(np.count_nonzero(x))])
log_a = np.log(a[np.flatnonzero(a)])
log_t = np.log(t[np.flatnonzero(t)])
log_g = np.log(g[np.flatnonzero(g)])
log_c = np.log(c[np.flatnonzero(c)])
log_x = np.log(x[np.flatnonzero(x)])
a_model = np.linalg.lstsq(A.T, log_a)[0]
t_model = np.linalg.lstsq(T.T, log_t)[0]
g_model = np.linalg.lstsq(G.T, log_g)[0]
c_model = np.linalg.lstsq(C.T, log_c)[0]
x_model = np.linalg.lstsq(X.T, log_x)[0]

print "(Intercept)\tA.nonzero" 
print str(a_model[1]) + "\t" + str(a_model[0])
print "(Intercept)\tT.nonzero" 
print str(t_model[1]) + "\t" + str(t_model[0])
print "(Intercept)\tG.nonzero"
print str(g_model[1]) + "\t" + str(g_model[0])
print "(Intercept)\tC.nonzero"
print str(c_model[1]) + "\t" + str(c_model[0])
print "(Intercept)\tX.nonzero"
print str(x_model[1]) + "\t" + str(x_model[0])

