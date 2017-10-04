#!/usr/bin/env python

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
from scipy.stats import linregress

np.set_printoptions(precision=15)
np.set_printoptions(linewidth=132)

testcase = 1
finalTime = 0.5
elements = np.arange(25,201,25)
elements = np.array([10, 20, 40, 80, 160, 220])
elements = np.arange(25,301,25)
orders = range(1,9)
nn = 7
if(testcase == 0): nn = 2
if(testcase == 0): orders = range(1,9)
elements = 2*np.logspace(0,nn,num=nn+1,base=2)
#elements = 10*np.linspace(1,nn,num=nn)
if(len(sys.argv) == 1):
    inf = 'input.in'
    outf = 'output.dat'

    os.system('make')

    h = 1.0/elements
    error = np.zeros([len(orders), len(elements)])
    error1 = np.zeros([len(orders), len(elements)])
    error2 = np.zeros([len(orders), len(elements)])
    for i, order in enumerate(orders):
        for j, nele in enumerate(elements):
            fin = open(inf, 'w')
            fin.write('%d %d' % (order, nele))
            fin.close()
            n = (order+1) * nele

            os.system('./advection')

            x, u = np.loadtxt(outf, unpack = True, dtype=np.float128)
            ref = x[2*n:2*n+order+1]
            x = x[0:n]
            u0 = u[0:n]
            uf = u[n:2*n]
            weights = u[2*n:2*n+order+1]
            if j==0:
                print('integral')
                print(ref)
                print(weights)
                print(np.dot(ref**2.0,weights))


            xe = np.linspace(x[0],x[-1],300)
            ue = 0.0*xe
            for ix, xx in enumerate(xe):
                if -1.0+finalTime < xx < 1.0+finalTime:
                    ue[ix] = np.exp(-1.0/(1.0-(xx-finalTime)**2))
            if(testcase == 0): ue = np.sin(xe-finalTime)

#           plt.figure()
#           plt.plot(x,u0,'-bo',mec='blue',mfc='None')
#           plt.plot(x,uf,'-ro',mec='red',mfc='None')
#           plt.plot(xe,ue,'-k')
#           plt.show()

            ue = 0.0*x
            for ix, xx in enumerate(x):
                if -1.0+finalTime < xx < 1.0+finalTime:
                    ue[ix] = np.exp(-1.0/(1.0-(xx-finalTime)**2))
            if(testcase == 0): ue = np.sin(x-finalTime)

            error[i,j] = np.linalg.norm((ue[:] - uf[:])) / np.linalg.norm(ue[:])
            #error[i,j] = np.linalg.norm((ue[:] - uf[:]),2)
            uer = np.reshape(ue, (nele,order+1))
            ufr = np.reshape(uf, (nele,order+1))
            
            num = 0
            den = 0
            for iele in range(int(nele)):
                for iref in range(order+1):
                    err = uer[iele,iref] - ufr[iele,iref]
                    #num += weights[iref]*abs(err)**2
                    num += weights[iref]*abs(err)**2
                    den += weights[iref]*abs(uer[iele,iref])**2
            error[i,j] = np.sqrt(abs(num / den))
            err = (ue[:] - uf[:])
            error1[i,j] = np.linalg.norm(err,1) / np.linalg.norm(ue[:],1)
            error2[i,j] = np.linalg.norm(err,2) / np.linalg.norm(ue[:])
            #error1[i,j] = np.linalg.norm(err/ue[:],1) / nele
            #error2[i,j] = np.linalg.norm(err/ue[:],2) / nele

            print('order: %d \t nele: %d \t  quad error: %e \t l1e: %e \t  l2e: %e' % (order,nele,error[i,j], error1[i,j], error2[i,j]))
        print('**************')

    np.savez('save.npz', h=h, error=error, error1=error1, error2=error2)
else:
    data = np.load('save.npz')
    h = data['h']
    error = data['error']
    error1 = data['error1']
    error2 = data['error2']


pdfName = 'Figures.pdf'
pp=PdfPages(pdfName)
npts= -2
for i, order in enumerate(orders):
    #plt.loglog(h, error[i,:],'-o', label='p=%d'%order)
    logh = np.log10(h)
    loge2 = np.log10(error2[i,:])
    loge1 = np.log10(error1[i,:])
    logeq = np.log10(error[i,:])
    slope1, intercept, r_value, p_value, std_err = linregress(logh[npts:],loge1[npts:])
    slope2, intercept, r_value, p_value, std_err = linregress(logh[npts:],loge2[npts:])
    slopeq, intercept, r_value, p_value, std_err = linregress(logh[npts:],logeq[npts:])

    loge = np.log10(error[i,:])
    slope, intercept, r_value, p_value, std_err = linregress(logh[npts:],loge[npts:])
    plt.plot(logh, loge, '-o', label='p=%d,s=%f'%(order,slope))
    xfid = np.linspace(np.log10(h[0]),np.log10(h[-1]))
    plt.plot(xfid, xfid*slope+intercept,'-k')
    print('order=%d, sq=%f, s1=%f, s2=%f' %(order,slopeq,slope1,slope2))
plt.xlabel('log(h)')
plt.ylabel('log(err)')
plt.legend(loc='best',fontsize='small')
plt.show()
pp.savefig(bbx_inches='tight')
pp.close()

#plt.close()
#print(error)
#print(np.shape(error))
#for i, ele in enumerate(elements):
#    #plt.loglog(h, error[i,:],'-o', label='p=%d'%order)
#    e = error[:,i]
#    print(orders)
#    print(e)
#    plt.semilogy(orders, e, '-o', label='n=%d'%(ele))
#plt.xlabel('order p')
#plt.ylabel('log(err)')
#plt.legend(loc='best')
#plt.show()

