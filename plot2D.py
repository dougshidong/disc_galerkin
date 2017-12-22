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
import scipy.integrate as inte

np.set_printoptions(precision=3)
np.set_printoptions(linewidth=132)

colors =['b', 'g', 'r', 'c', 'm', 'y', 'k']
markers=['o','v','^','<','>','1','2','3','4','8','s','p','P','*','h','H','+','x','X','D','d','|','_']

pdfName = 'Figures_2.pdf'
pp=PdfPages(pdfName)
nn = 1
orders = np.arange(2,5,dtype=np.int32)
#if(testcase == 0): nn = 5
#if(testcase == 0): orders = np.arange(1,9,dtype=np.int32)
#if(testcase == 1): nn = 5
#if(testcase == 1): orders = np.arange(2,5,dtype=np.int32)
#if(testcase == 3): nn = 5
#if(testcase == 3): orders = np.arange(2,5,dtype=np.int32)

wavespeed = [1.0, 1.0]
testcase = 1
finalTime = 0.1
ndim = 1
ngrids = 10
maxorders = 6
elements = (2*np.logspace(0,ngrids-1,num=ngrids,base=1.5,dtype=np.float64)+1)
elements = elements.astype(int)
elements = (2*np.logspace(2,ngrids-1,num=ngrids,base=1.5,dtype=np.float64)+1)
elements = elements.astype(int)
#if(ndim==2): elements = np.array([3,4,5,6,7,8,9,10,11,12,13,14,15,16])
#if(ndim==2): elements = np.array([3,4,5,6,7,8,9,10])
if(ndim==2): maxorders = 5
orders = np.arange(1,maxorders+1,dtype=np.int32)
print(elements)


if(testcase == -1):
    xmin = -1.0;    xmax = 1.0;    ymin = -1.0*np.pi;    ymax = 1.0*np.pi
if(testcase == 0):
    xmin = -1.0*np.pi;    xmax = 1.0*np.pi;    ymin = -1.0*np.pi;    ymax = 1.0*np.pi
if(testcase == 1):
    xmin = -3.0;    xmax = 3.0;    ymin = -3.0;    ymax = 3.0;
if(testcase == 3):
    xmin = -1.0;    xmax = 1.0;    ymin = -1.0;    ymax = 1.0
if(testcase == 4):
    xmin = -1.0*np.pi;    xmax = 1.0*np.pi;    ymin = -1.0*np.pi;    ymax = 1.0*np.pi

if(len(sys.argv) == 1):
    inf = 'input.in'
    outf = 'output.dat'

    os.system('make')

    if(ndim==1): h = (xmax-xmin)/(elements-1.0)
    if(ndim==2): h = (xmax-xmin)/(elements-1.0)
    error = np.zeros([len(orders), len(elements)])
    slopes = np.zeros([len(orders), len(elements)])
    for i, order in enumerate(orders):
        for j, nvert in enumerate(elements):
            fin = open(inf, 'w')
            fin.write('%d %d %d' % (order, nvert, nvert))
            fin.close()

            nele = nvert-1
            if(ndim == 2): nele = (nvert-1)**2

            nref = order + 1
            quad_npts = nref
            npts = nref
            if(ndim == 2): npts = nref**2
            dof = nele*npts

            os.system('./advection')

            if ndim == 1:
                x0, u0 = np.loadtxt('init.x', unpack = True, dtype=np.float128,skiprows=0)
                xf, uf = np.loadtxt('final.x', unpack = True, dtype=np.float128,skiprows=0)
                x = x0
            if ndim == 2:
                x0, y0, u0 = np.loadtxt('init.x', unpack = True, dtype=np.float128,skiprows=0)
                xf, yf, uf = np.loadtxt('final.x', unpack = True, dtype=np.float64,skiprows=0)
                x = x0
                y = y0


            quadR, quadW = np.loadtxt('quad.x', unpack = True, dtype=np.float128,skiprows=0)
            cubR, cubS, cubW = np.loadtxt('cub2D.x', unpack = True, dtype=np.float128,skiprows=0)
            # Check quadrature
#            if j==0:
#                print('integral x**(2n-1)', np.dot(quadR**4.0,quadW), '=%f'%(2.0/(5.0)))
#                print('integral x**2', np.dot(quadR**2.0,quadW), '=2/3')
#                print(quadR)
#                print(quadW)

            # Check cubature
            #if ndim==2 and j==0:
            #     print('integral x**(2n-1)', np.dot(cubR**4.0*cubS**4.0,cubW), '=%f'%(4.0/(9.0)))
            #     print(quadR)
            #     print(quadW)


            if ndim == 1:
                xexact = np.linspace(x[0],x[-1],100)
                uexact = 0.0*xexact
                if(testcase == -1): uexact = np.sin(2.15*(xexact)+0.23)
                if(testcase == 0): uexact = np.sin(xexact-wavespeed[0]*finalTime)+np.cos(xexact-wavespeed[0]*finalTime)
                if(testcase == 1): uexact = np.exp(-10*(xexact-wavespeed[0]*finalTime)**2)
                if(testcase == 3): uexact = np.cos(xexact)
                if(testcase == 4): uexact = np.cos(xexact)+ np.sin(xexact)


            if ndim == 1:
                plt.figure()
                plt.title('Order = %d, Elements = %d, DoF = %d'%(order, nele, dof))
                # Plotting numerical solution element by element
                for iele in range(nele):
                    poly = lagrange(x[iele*npts:(iele+1)*npts],uf[iele*npts:(iele+1)*npts])
                    xr = np.linspace(x[iele*npts],x[(iele+1)*npts-1],100,endpoint=True)
                    a = (iele+1e-16)/(nele*1.0)*(xmax-xmin)+xmin
                    b = (iele+1.0)/(nele*1.0)*(xmax-xmin)+xmin
                    xr = np.linspace(a, b, 100,endpoint=True)
                    ur = np.polyval(poly,xr)
                    plt.plot(xr,ur,colors[iele%len(colors)],marker=None)#,mec=colors[iele%len(colors)],mfc='None')
                    plt.plot(x[iele*npts:(iele+1)*npts],uf[iele*npts:(iele+1)*npts],'o',mec=colors[iele%len(colors)],mfc='None')
                # Exact Solution
                plt.plot(xexact,uexact,'-k')
                pp.savefig(bbx_inches='tight')
                plt.close()
            if ndim == 2:
                if(j==0):
                    plt.figure()
                    xexact = (np.random.rand(5000)*2-1)*max(xf)
                    yexact = (np.random.rand(5000)*2-1)*max(yf)
                    #uexact = np.sin(xexact-finalTime) + np.cos(yexact-finalTime)
                    if(testcase == 0): uexact = np.sin(xexact-wavespeed[0]*finalTime) + np.cos(xexact-wavespeed[0]*finalTime)  \
                                               +np.sin(yexact-wavespeed[1]*finalTime) + np.cos(yexact-wavespeed[1]*finalTime)
                    if(testcase == 1): uexact = np.exp(-5.0*(xexact-wavespeed[0]*finalTime)**2 -5.0*(yexact-wavespeed[1]*finalTime)**2)
                    if(testcase == 3): uexact = np.cos(xexact) - np.sin(yexact)
                    if(testcase == 4): uexact = np.cos(xexact) - np.sin(yexact)
                    plt.tripcolor(xexact,yexact,uexact)
                    plt.colorbar()
                    plt.xlim([min(xf),max(xf)])
                    plt.ylim([min(yf),max(yf)])
                    pp.savefig(bbx_inches='tight')
                    plt.close()

                plt.figure()
                plt.title('Order = %d, Elements = %d, DoF = %d'%(order, nele, dof))
                #f, ax = plt.subplots(1,2, sharex=True, sharey=True)
                plt.tripcolor(xf,yf,uf)
                plt.colorbar()
                plt.plot(xf,yf, 'ko', ms=1)
                plt.xlim([min(xf),max(xf)])
                plt.ylim([min(yf),max(yf)])
#               ax[1].tripcolor(xexact,yexact,uexact) # choose 20 contour levels, just to show how good its interpolation is
                #ax[1].tricontourf(xexact,yexact,uexact,20) # choose 20 contour levels, just to show how good its interpolation is
                #ax[1].plot(x,y, 'ko ')
                pp.savefig(bbx_inches='tight')
                plt.close()

            if ndim == 1:
                xexact = x
                uexact = 0.0*xexact
                if(testcase == -1): uexact = np.sin(2.15*(xexact)+0.23)
                if(testcase == 0): uexact = np.sin(xexact-wavespeed[0]*finalTime) + np.cos(xexact-wavespeed[0]*finalTime)
                if(testcase == 1): uexact = np.exp(-10*(xexact-wavespeed[0]*finalTime)**2)
                if(testcase == 3): uexact = np.cos(xexact)
                if(testcase == 4): uexact = np.cos(xexact)+np.sin(xexact)
            if ndim == 2: 
                xexact = x
                yexact = y
                #xe = x
                #ye = y
                #xe, ye = np.meshgrid(xe,ye)
                #xexact = np.reshape(xe,(100**2))
                #yexact = np.reshape(ye,(100**2))
                uexact = 0.0*xexact
                if(testcase == 0): uexact = np.sin(xexact-wavespeed[0]*finalTime) + np.cos(xexact-wavespeed[0]*finalTime)  \
                                           +np.sin(yexact-wavespeed[1]*finalTime) + np.cos(yexact-wavespeed[1]*finalTime)
                if(testcase == 1): uexact = np.exp(-5.0*(xexact-wavespeed[0]*finalTime)**2 -5.0*(yexact-wavespeed[1]*finalTime)**2)
                if(testcase == 3): uexact = np.cos(xexact) - np.sin(yexact)
                if(testcase == 4): uexact = np.cos(xexact) - np.sin(yexact)

            u_exact_reshaped = np.reshape(uexact, (nele,npts))
            u_final_reshaped = np.reshape(uf, (nele,npts))
            num = 0
            den = 0
            for iele in range(int(nele)):
                for iref in range(npts):
                    err = u_exact_reshaped[iele,iref] - u_final_reshaped[iele,iref]
                    if(ndim == 1):
                        num += quadW[iref]*abs(err)**2
                        den += quadW[iref]*abs(u_exact_reshaped[iele,iref])**2
                    if(ndim == 2):
                        num += cubW[iref]*abs(err)**2
                        den += cubW[iref]*abs(u_exact_reshaped[iele,iref])**2
            error[i,j] = np.sqrt(abs(num / (den+1e-16)))


#           relerr = 0.0
#           for iele in range(nele):
#               poly = lagrange(x[iele*npts:(iele+1)*npts],uf[iele*npts:(iele+1)*npts])
#               def integrand1(x):
#                   if(testcase == 0): u = np.sin(x-finalTime)
#                   if(testcase == 1): u = np.exp(-10*(x-finalTime)**2)
#                   if(testcase == 3): u = np.cos(x)
#                   if(testcase == 4): u = np.cos(x)
#                   return ((np.polyval(poly,x) - u)**2)
#               def integrand2(x):
#                   return ((np.sin(x-finalTime))**2)
#               a = (iele+1e-16)/(nele*1.0)*(xmax-xmin)+xmin
#               b = (iele+1.0)/(nele*1.0)*(xmax-xmin)+xmin
#               num = inte.quad(integrand1, a, b, epsrel=1e-15)
#               den = inte.quad(integrand2, a, b, epsrel=1e-15)
#               relerr = relerr + num[0] / (den[0])
#           error[i,j] = np.sqrt(relerr)

            logh = np.log10(h)
            loge = np.log10(abs(error[i,:])+1e-16)
            slope = 0.0
            if j > 0:
                slope, intercept, r_value, p_value, std_err = linregress(logh[j-1:j+1],loge[j-1:j+1])
            slopes[i,j] = slope
            slope2 = 0.0
            if j > 1:
                slope2, intercept, r_value, p_value, std_err = linregress(logh[j-2:j+1],loge[j-2:j+1])

            print('order: %d \t nele: %d \t  quad error: %e \t slope: %e \t slope2: %e' % (order,nele,error[i,j], slope, slope2))
        print('**************')
pp.close()

