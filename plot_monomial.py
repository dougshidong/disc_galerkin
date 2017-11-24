#!/usr/bin/env python

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

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


mono = lambda order, x : x**order

n = 100
x = np.arange(-1,1,1.0/n)

orders = range(2,3)

matplotlib.rcParams.update({'font.size': 16})
pdfName = 'Figures.pdf'
pp=PdfPages(pdfName)
for order in orders:
    nbasis = order + 1

    plt.figure(figsize=(16,12))
    plt.title('Representing a function with a polynomial')

    for i in range(nbasis):
        y = mono(i, x)
        plt.plot(x,y,'-o',label='$\psi_{%d}(x) = x^%d$'%(i+1,i))

    plt.plot(x,2*np.sin(x+1)-0.5, '-k', linewidth=4, label='$2\sin(x+1)-0.5$')
    plt.plot(x,1.18+0.975*x-0.782*x*x,'-o',label='$u_h(x) = \hat{u}_1\psi_1+\hat{u}_2\psi_2+\hat{u}_3\psi_3$')

    plt.legend(loc=4)
    pp.savefig(bbx_inches='tight')
    plt.close()

pp.close()
