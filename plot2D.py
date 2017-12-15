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
from scipy.interpolate import lagrange

np.set_printoptions(precision=3)
np.set_printoptions(linewidth=132)

colors =['b', 'g', 'r', 'c', 'm', 'y', 'k']
markers=['o','v','^','<','>','1','2','3','4','8','s','p','P','*','h','H','+','x','X','D','d','|','_']

pdfName = 'Figures_2.pdf'
pp=PdfPages(pdfName)
testcase = 0
finalTime = 1.5
nn = 3
orders = np.arange(2,5,dtype=np.int32)
if(testcase == 0): nn = 4
if(testcase == 0): orders = np.arange(1,7,dtype=np.int32)
if(testcase == 1): nn = 6
if(testcase == 1): orders = np.arange(2,15,dtype=np.int32)

nn = 3
orders = np.arange(1,4,dtype=np.int32)
elements = 2*np.logspace(1,nn,num=nn+1,base=2,dtype=np.int32)
#elements = 10*np.linspace(1,nn,num=nn)
#elements = 32*np.logspace(0,1,1,base=2,dtype=np.int32)
    
orders = [1,2,3]
nn = 4
elements = 2*np.logspace(1,nn,num=nn+1,base=2,dtype=np.int32)
if(len(sys.argv) == 1):
    inf = 'input.in'
    outf = 'output.dat'

    os.system('make')

    h = 1.0/elements
    error = np.zeros([len(orders), len(elements)])
    slopes = np.zeros([len(orders), len(elements)])
    for i, order in enumerate(orders):
        for j, nvert in enumerate(elements):
            fin = open(inf, 'w')
            fin.write('%d %d %d' % (order, nvert, nvert))
            fin.close()

            nele = nvert-1

            nref = order + 1
            quad_npts = nref
            npts = nref
            npts = quad_npts
            n = (npts) * nele
            os.system('./advection')

            x0, u0 = np.loadtxt('init.x', unpack = True, dtype=np.float128,skiprows=0)
            xf, uf = np.loadtxt('final.x', unpack = True, dtype=np.float128,skiprows=0)

            x = x0

            quadR, quadW = np.loadtxt('quad.x', unpack = True, dtype=np.float128,skiprows=0)
            # Check quadrature
            if j==0:
                print('integral x**2', np.dot(quadR**2.0,quadW), '=2/3')
                print(quadR)
                print(quadW)

            xe = np.linspace(x[0],x[-1],300)
            xe = np.linspace(-2*3.1416,2*3.1416,300)
            ue = 0.0*xe
            for ix, xx in enumerate(xe):
                if -1.0+finalTime < xx < 1.0+finalTime:
                    ue[ix] = np.exp(-1.0/(1.0-(xx-finalTime)**2))
            if(testcase == 0): ue = np.sin(xe-finalTime)
            if(testcase == 3): ue = np.cos(xe)

            finalTime = 0.5
            ue = np.sin(xe-finalTime)

            plt.figure()
            plt.title('Order = %d, Elements = %d'%(order, nele))
            # Plotting numerical solution element by element
            for ix in range(nele):
                poly = lagrange(x[ix*npts:(ix+1)*npts],uf[ix*npts:(ix+1)*npts])
                xr = np.linspace(x[ix*npts],x[(ix+1)*npts-1],100,endpoint=True)
                ur = np.polyval(poly,xr)
                plt.plot(xr,ur,colors[ix%len(colors)],marker=None)#,mec=colors[ix%len(colors)],mfc='None')
                plt.plot(x[ix*npts:(ix+1)*npts],uf[ix*npts:(ix+1)*npts],'o',mec=colors[ix%len(colors)],mfc='None')
            # Exact Solution
            plt.plot(xe,ue,'-k')
            pp.savefig(bbx_inches='tight')
            plt.close()

            ue = 0.0*x
            ue = np.sin(x-finalTime)
            uer = np.reshape(ue, (nele,npts))
            ufr = np.reshape(uf, (nele,npts))
            num = 0
            den = 0
            for iele in range(int(nele)):
                for iref in range(npts):
                    err = uer[iele,iref] - ufr[iele,iref]
                    num += quadW[iref]*abs(err)**2
                    den += quadW[iref]*abs(uer[iele,iref])**2
            error[i,j] = np.sqrt(abs(num / (den+1e-15)))

            logh = np.log10(h)
            loge = np.log10(abs(error[i,:])+1e-14)
            slope = 0.0
            if j > 0:
                slope, intercept, r_value, p_value, std_err = linregress(logh[j-1:j+1],loge[j-1:j+1])
            slopes[i,j] = slope

            print('order: %d \t nele: %d \t  quad error: %e \t slope: %e' % (order,nele,error[i,j], slope))
        print('**************')
pp.close()

