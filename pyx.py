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

colors =['b', 'g', 'r', 'c', 'm', 'y', 'k' 
        ,'b', 'g', 'r', 'c', 'm', 'y', 'k' 
        ,'b', 'g', 'r', 'c', 'm', 'y', 'k' 
        ,'b', 'g', 'r', 'c', 'm', 'y', 'k' 
        ,'b', 'g', 'r', 'c', 'm', 'y', 'k' 
        ,'b', 'g', 'r', 'c', 'm', 'y', 'k' 
        ,'b', 'g', 'r', 'c', 'm', 'y', 'k' 
        ,'b', 'g', 'r', 'c', 'm', 'y', 'k' 
        ,'b', 'g', 'r', 'c', 'm', 'y', 'k' 
        ,'b', 'g', 'r', 'c', 'm', 'y', 'k' 
        ,'b', 'g', 'r', 'c', 'm', 'y', 'k' 
        ,'b', 'g', 'r', 'c', 'm', 'y', 'k' 
        ,'b', 'g', 'r', 'c', 'm', 'y', 'k' ]
markers=['o','v','^','<','>','1','2','3','4','8','s','p','P','*','h','H','+','x','X','D','d','|','_']

pdfName = 'Figures_2.pdf'
pp=PdfPages(pdfName)
testcase = 0
finalTime = 0.5
nn = 3
orders = np.arange(2,5,dtype=np.int32)
if(testcase == 0): nn = 4
if(testcase == 0): orders = np.arange(1,7,dtype=np.int32)
if(testcase == 1): nn = 6
if(testcase == 1): orders = np.arange(2,15,dtype=np.int32)
elements = 2*np.logspace(0,nn,num=nn+1,base=2,dtype=np.int32)
#elements = 10*np.linspace(1,nn,num=nn)
#elements = 32*np.logspace(0,1,1,base=2,dtype=np.int32)
if(len(sys.argv) == 1):
    inf = 'input.in'
    outf = 'output.dat'

    os.system('make')

    h = 1.0/elements
    slopes = np.zeros([len(orders), len(elements)])
    error = np.zeros([len(orders), len(elements)])
    error1 = np.zeros([len(orders), len(elements)])
    error2 = np.zeros([len(orders), len(elements)])
    for i, order in enumerate(orders):
        for j, nele in enumerate(elements):
            fin = open(inf, 'w')
            fin.write('%d %d' % (order, nele))
            fin.close()


            nref = order + 1
            quad_npts = nref

            npts = nref
            npts = quad_npts

            n = (npts) * nele

            os.system('./advection')

            xread, uread = np.loadtxt(outf, unpack = True, dtype=np.float128,skiprows=0)
            idata = 0
            x = xread[0:n]
            u0 = uread[0:n]
            idata += n
            uf = uread[idata:idata+n]
            idata += n
#           basis_fun = np.zeros([order+1,n])
#           for ix in range(order+1):
#               basis_fun[ix,:] = uread[idata:idata+n]
#               idata += n
#           plt.figure()
#           for ix in range(order+1):
#               plt.title('nbasis = %d, order = %d'%(order+1, max(order+1-14,2)))
#               #plt.title('nbasis = %d, order = %d'%(order+1, order))
#               plt.plot(x,basis_fun[ix,:],colors[ix], ms = 8, mec = colors[ix], marker=markers[ix], mfc='None',label=('ibasis=%d'%ix))
#               plt.legend()
#           pp.savefig(bbx_inches='tight')

            ref = xread[idata:idata+npts]
            weights = uread[idata:idata+npts]
            # Check quadrature
#           if j==0:
#               print('integral x**2', np.dot(ref**2.0,weights), '=2/3')
#               print(ref)
#               print(weights)

            xe = np.linspace(x[0],x[-1],300)
            ue = 0.0*xe
            for ix, xx in enumerate(xe):
                if -1.0+finalTime < xx < 1.0+finalTime:
                    ue[ix] = np.exp(-1.0/(1.0-(xx-finalTime)**2))
            if(testcase == 0): ue = np.sin(xe-finalTime)
            if(testcase == 3): ue = np.cos(xe)

            plt.figure()
            plt.title('Order = %d, Elements = %d'%(order, nele))
#           plt.plot(x,u0,'-bo',mec='blue',mfc='None')
#           for ix in range(nele):
#               plt.plot(x[ix*npts:(ix+1)*npts],u0[ix*npts:(ix+1)*npts],colors[ix],marker='o',mec=colors[ix],mfc='None')
            for ix in range(nele):
                poly = lagrange(x[ix*npts:(ix+1)*npts],uf[ix*npts:(ix+1)*npts])
                xr = np.linspace(x[ix*npts],x[(ix+1)*npts-1],100,endpoint=True)
                ur = np.polyval(poly,xr)
                plt.plot(xr,ur,colors[ix],marker=None)#,mec=colors[ix],mfc='None')
                plt.plot(x[ix*npts:(ix+1)*npts],uf[ix*npts:(ix+1)*npts],'o',mec=colors[ix],mfc='None')
                #plt.plot(x[ix*npts:(ix+1)*npts],uf[ix*npts:(ix+1)*npts],colors[ix],marker='o',mec=colors[ix],mfc='None')
            #plt.plot(x,uf,'-r^',mec='red',mfc='None')
            plt.plot(xe,ue,'-k')
            pp.savefig(bbx_inches='tight')
            plt.close()

            ue = 0.0*x
            for ix, xx in enumerate(x):
                if -1.0+finalTime < xx < 1.0+finalTime:
                    ue[ix] = np.exp(-1.0/(1.0-(xx-finalTime)**2))
            if(testcase == 0): ue = np.sin(x-finalTime)
            if(testcase == 3): ue = np.cos(x)

            error[i,j] = np.linalg.norm((ue[:] - uf[:])) / np.linalg.norm(ue[:]+1e-14)
            uer = np.reshape(ue, (nele,npts))
            ufr = np.reshape(uf, (nele,npts))

            num = 0
            den = 0
            for iele in range(int(nele)):
                for iref in range(npts):
                    err = uer[iele,iref] - ufr[iele,iref]
                    #num += weights[iref]*abs(err)**2
                    num += weights[iref]*abs(err)**2
                    den += weights[iref]*abs(uer[iele,iref])**2
            error[i,j] = np.sqrt(abs(num / (den+1e-15)))
            err = (ue[:] - uf[:])
            error1[i,j] = np.linalg.norm(err,1) / (np.linalg.norm(ue[:],1)+1e-14)
            error2[i,j] = np.linalg.norm(err,2) / (np.linalg.norm(ue[:])+1e-14)

            logh = np.log10(h)
            loge = np.log10(abs(error[i,:])+1e-14)
            slope = 0.0
            if j > 0:
                slope, intercept, r_value, p_value, std_err = linregress(logh[j-1:j+1],loge[j-1:j+1])
            slopes[i,j] = slope

            print('order: %d \t nele: %d \t  quad error: %e \t slope: %e' % (order,nele,error[i,j], slope))
        print('**************')

    np.savez('save.npz', h=h, error=error, error1=error1, error2=error2)
else:
    data = np.load('save.npz')
    h = data['h']
    error = data['error']
    error1 = data['error1']
    error2 = data['error2']
pp.close()


pdfName = 'Figures.pdf'
pp=PdfPages(pdfName)
slope_s= -2
slope_e= None
plt.figure()
for i, order in enumerate(orders):
    #plt.loglog(h, error[i,:],'-o', label='p=%d'%order)
    logh = np.log10(h)
    #logh = np.log10(elements)
    loge2 = np.log10(error2[i,:])
    loge1 = np.log10(error1[i,:])
    logeq = np.log10(error[i,:])
    slope1, intercept, r_value, p_value, std_err = linregress(logh[slope_s:slope_e],loge1[slope_s:slope_e])
    slope2, intercept, r_value, p_value, std_err = linregress(logh[slope_s:slope_e],loge2[slope_s:slope_e])
    slopeq, intercept, r_value, p_value, std_err = linregress(logh[slope_s:slope_e],logeq[slope_s:slope_e])

    loge = np.log10(error[i,:])
    #print(loge)
    #print(logh)
    slope, intercept, r_value, p_value, std_err = linregress(logh[slope_s:slope_e],loge[slope_s:slope_e])
    plt.plot(logh, loge, '-o', label='p=%d,s=%3.2f'%(order,slope))
    xfid = np.linspace(np.log10(h[0]),np.log10(h[-1]))
    #xfid = np.linspace(np.log10(elements[0]),np.log10(elements[-1]))
    plt.plot(xfid, xfid*slope+intercept,'-k')
    print('order=%d, sq=%f, s1=%f, s2=%f' %(order,slopeq,slope1,slope2))
plt.xlabel(r'$\log (h_k)$')
#plt.xlabel(r'log(number of elements)')
plt.ylabel(r'$\log(\|\| \mathbf{u}^* - \mathbf{u}_h \|\|_2)$')
#plt.ylabel(r'log(error)')
plt.title('Spectral Convergence of DG')
plt.legend(loc='best',fontsize='small')
pp.savefig(bbx_inches='tight')
pp.close()

def tab(a):
    """Returns a LaTeX bmatrix

    :a: numpy array
    :returns: LaTeX bmatrix as a string
    """
    if len(a.shape) > 2:
        raise ValueError('bmatrix can at most display two dimensions')
    lines = str(a).replace('[', '').replace(']', '').splitlines()
    rv = [r'\begin{bmatrix}']
    rv += ['  ' + ' & '.join(l.split()) + r' \\' for l in lines]
    #rv +=  [r'\end{bmatrix}']
    return '\n'.join(rv)
print(tab(np.transpose(error))+'\n')
print(tab(np.transpose(slopes))+'\n')


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

