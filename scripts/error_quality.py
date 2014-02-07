#!/usr/bin/python

import numpy as np
import sys

d = np.loadtxt(sys.argv[1], skiprows=0)
axis = d[:,0]
a = d[:,1]
t = d[:,2]
g = d[:,3]
c = d[:,4]
x = d[:,5]
a_model, a_resid = np.polyfit(np.flatnonzero(a),a[np.flatnonzero(a)],2,full=True)[:2]
r2 = 1 - a_resid / (a[np.flatnonzero(a)].size * a[np.flatnonzero(a)].var())
t_model, t_resid = np.polyfit(np.flatnonzero(t),a[np.flatnonzero(t)],2,full=True)[:2]
g_model, g_resid = np.polyfit(np.flatnonzero(g),a[np.flatnonzero(g)],2,full=True)[:2]
c_model, c_resid = np.polyfit(np.flatnonzero(c),a[np.flatnonzero(c)],2,full=True)[:2]
x_model, x_resid = np.polyfit(np.flatnonzero(x),a[np.flatnonzero(x)],2,full=True)[:2]


print "Intercept\tpoly(A.nonzero)1\tpoly(A.nonzero)2" 
print str(a_model[2]) + "\t" + str(a_model[1]) + "\t" + str(a_model[0])
print "Intercept\tpoly(T.nonzero)1\tpoly(T.nonzero)2"
print str(t_model[2]) + "\t" + str(t_model[1]) + "\t" + str(t_model[0])
print "Intercept\tpoly(G.nonzero)1\tpoly(G.nonzero)2"
print str(g_model[2]) + "\t" + str(g_model[1]) + "\t" + str(g_model[0])
print "Intercept\tpoly(C.nonzero)1\tpoly(C.nonzero)2"
print str(c_model[2]) + "\t" + str(c_model[1]) + "\t" + str(c_model[0])
print "Intercept\tpoly(X.nonzero)1\tpoly(X.nonzero)2"
print str(x_model[2]) + "\t" + str(x_model[1]) + "\t" + str(x_model[0])

