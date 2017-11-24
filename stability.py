#!/usr/bin/env python

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from scipy.spatial import ConvexHull

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

pdfName = 'Figures.pdf'
pp=PdfPages(pdfName)
nele = 8
nk = 100
nrk = 3
orders = range(1,5)
def check_stab():
    if(len(sys.argv) == 1):
        inf = 'input.in'
        outf = 'out_stab.out'

        os.system('make')

        plt.figure()
        for i, order in enumerate(orders):
            fin = open(inf, 'w')
            fin.write('%d %d' % (order, nele))
            fin.close()

            nref = order + 1

            os.system('./advection > out_stab.out')

            real, imag = np.loadtxt(outf, unpack = True, dtype=np.float128,skiprows=0)

            istr = 0
            iend = nref*nk

            pts = np.empty([nref*nk,2],dtype=np.float64)
            print(istr,iend,np.shape(real))
            pts[:,0] = real[istr:iend]
            pts[:,1] = imag[istr:iend]
            hull = ConvexHull(pts)

            plt.plot(pts[hull.vertices, 0], pts[hull.vertices, 1],
                color=colors[i], marker='o', ms = 4, mec=colors[i],mfc=colors[i],label='p=%d'%order)

            plt.legend()
            plt.xlim([-20,0])
    pp.savefig(bbx_inches='tight')
    plt.close()
    pp.close()


def vonNeumann():
    inf = 'input.in'
    outf = 'out_stab.out'

    os.system('make')

    for i, order in enumerate(orders):
        fin = open(inf, 'w')
        fin.write('%d %d' % (order, nele))
        fin.close()

        nref = order + 1

        os.system('./advection > out_stab.out')

        real, imag = np.loadtxt(outf, unpack = True, dtype=np.float128,skiprows=0)

        plt.figure()
        plt.title('order = %d'%(order))
        for irk in range(nrk):
            istr = irk*nref*nk
            iend = (irk+1)*nref*nk

            pts = np.empty([nref*nk,2],dtype=np.float64)
            print(istr,iend,np.shape(real))
            pts[:,0] = real[istr:iend]
            pts[:,1] = imag[istr:iend]
            hull = ConvexHull(pts)

#           plt.plot(real[istr:iend],imag[istr:iend], \
#               linestyle='None', marker='o', ms = 4, mec=colors[irk],mfc=colors[irk],label='RK%d'%(irk+2))
#           plt.plot(real[istr:iend],imag[istr:iend], \
#               linestyle='None', marker='o', ms = 4, mec=colors[irk],mfc=colors[irk])
            plt.plot(pts[hull.vertices, 0], pts[hull.vertices, 1],
                color=colors[irk], marker='o', ms = 4, mec=colors[irk],mfc=colors[irk],
                label='RK%d'%(irk+2))
            plt.xlim([-3,0.5])
            x, y = np.meshgrid(np.arange(-6,6,0.01),np.arange(-3,3,0.01))
            z = x + y*1j
            R = 0.0
            def factorial(n):return reduce(lambda x,y:x*y,[1]+range(1,n+1))
            for i in range(irk+2+1):
                R = R + z**i / factorial(i)
            zlevel = abs(R)
            plt.contour(x,y,zlevel, [1], colors=colors[irk])
#           for simplex in hull.simplices:
#               plt.plot(pts[simplex, 0], pts[simplex, 1],
#                   color=colors[irk], marker='o', ms = 4, mec=colors[irk],mfc=colors[irk])


        plt.legend(loc=2)
        pp.savefig(bbx_inches='tight')
        plt.close()

    pp.close()

vonNeumann()
