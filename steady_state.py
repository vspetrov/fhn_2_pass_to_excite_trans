from scipy.optimize import fsolve
from scipy import roots
from numpy import linspace, imag, real, sqrt, dot
import numpy as np
from matplotlib import pyplot
def main():
    a = 0.5
    zeta=0.
    alpha=0.05
    D=linspace(0.1,0.9,80)
    # D = [0.4]
    threshold_offset = []
    xoffsets=[]
    xoffsets2=[]
    for d in D:
        Roots = roots([1.0/3.0, 0, d-d*d/(alpha+d), a-d*alpha*zeta/(alpha+d)])
        xbar=[]
        zbar=[]
        ybar=[]
        count = 0
        for r in Roots:
            if imag(r) == 0:
                count += 1
                xbar.append( real(r))
                zbar.append((alpha*zeta+d*real(r))/(alpha+d))
                ybar.append( xbar[-1]+a)

                c = zbar[-1]-xbar[-1]
                c = 0
                xm = -sqrt(1-d)
                xR = roots([1.0/3.0,0,-1,ybar[-1]+d*c])
                cc = 0
                v = []
                for xr in  xR:
                    if imag(xr) == 0:
                        # xoffsets.append(real(xr))
                        v.append(real(xr))
                tmp = []
                for vv in v:
                    V = 2*vv-xbar[-1]-zbar[-1]-c
                    if V > 0.1:
                        tmp.append(V)


                xoffsets.append(sorted(tmp)[0])

                # xoffsets2.append((ybar[-1]-xm+xm**3/3.0)/d-zbar[-1]-xbar[-1]+2*xm)
                th = abs(abs(ybar[-1])-abs(-2.0/3.0+d*(1+zbar[-1])))
                xoffsets2.append(th/d)
        assert(count == 1)
    pyplot.plot(D,xoffsets)
    pyplot.plot(D,xoffsets2)
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
