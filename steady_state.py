from scipy.optimize import fsolve
from scipy import roots
from numpy import linspace, imag, real, sqrt, dot
import numpy as np
from matplotlib import pyplot
def main():
    a = 0.5
    zeta=0
    alpha=0.05
    D=linspace(0.1,0.9,80)
    # D = [0.4]
    threshold_offset = []
    xoffsets=[]

    xb=[]
    zb=[]
    xoffsets2=[]
    for d in D:
        Roots = roots([1.0/3.0, 0, d-d*d/(alpha+d), a-d*alpha*zeta/(alpha+d)])
        xbar=0
        zbar=0
        ybar=0
        count = 0
        for r in Roots:
            if imag(r) == 0:
                count += 1
                xbar =  real(r)
                zbar = (alpha*zeta+d*real(r))/(alpha+d)
                ybar =  xbar+a
        assert(count == 1)
        xb.append(xbar)
        zb.append(zbar)
        assert(zbar > xbar)


        b = d/(alpha+d)
        c = alpha*zeta/(alpha+d)
        kappa = sqrt(1+b*b)

        # eta = (xbar+b*(zbar-c))/(kappa*kappa)
        # Roots = roots([-1./3.*b**3/kappa**6,
        #                -2./3.*b**2/kappa**4*eta,
        #                b/kappa**2-2./3.*b/kappa**2*eta**2+d*(1-1/kappa**2-b/kappa**2),
        #                eta-1./3.*eta**3-ybar+d*(zbar*(1-1/kappa**2)+(b*xbar+c)/kappa**2-(xbar+b*(zbar-c))/kappa**2)])

        eta = (zbar+xbar-c)/(b+1)
        Roots = roots([-1./3./(b+1)**3,
                       -eta/(b+1)**2,
                       1./(b+1)-(b+1)*eta**2+d*(b-1.)/(b+1),
                       eta-1./3.*eta**3-ybar+d*c+d*(b-1)/(b+1.)*(zbar+xbar)])

        # eta = xbar*(1+b**2/kappa**2)+(c-zbar)*b/kappa**2
        # Roots = roots([1./3.*b**3/kappa**6,
        #                -2./3.*b**2/kappa**4*eta,
        #                -b/kappa**2+2./3.*b/kappa**2*eta**2+d*(1+1/kappa**2+b/kappa**2),
        #                eta-1./3.*eta**3-ybar+d*(zbar*(1+1/kappa**2)+(b*xbar+c)/kappa**2-eta)])

        realRoots=[real(x) for x in Roots]
        xoffsets.append(realRoots[0])
        Roots = roots([-1./3.,0,1.,-ybar])
        realRoots=[real(x)*2-zbar-xbar for x in Roots]
        xoffsets2.append(realRoots[2])



    pyplot.plot(D,xoffsets,'b')

    pyplot.plot(D,xoffsets2,'g')

    xoffset_c=[]
    Dvalues_c=[]
    for line in open('rst_0_0.05.txt').readlines():
        d,xof=[float(x) for x in line.split()]
        xoffset_c.append(xof)
        Dvalues_c.append(d)

    pyplot.plot(Dvalues_c,xoffset_c,"ro-")
    pyplot.show()


if __name__ == '__main__':
    main()
