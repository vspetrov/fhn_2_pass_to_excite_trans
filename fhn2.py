import os, sys
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

import matplotlib.pyplot as plt
class FhnSystem:
    def __init__(self,aParam = 0):
        self.epsilon = 0.002
        self.a       = aParam
        self.initState(-aParam,0)
    def initState(self, x, y):
        self.x = x
        self.y = y
    def __str__(self):
        return "FHN (a,x,y) = ({0}, {1}, {2})".format(self.a, self.x, self.y)
    def __repr__(self):
        return self.__str__()
    def xFunction(self):
        return self.x-self.x*self.x*self.x/3.0-self.y
    def yFunction(self):
        return self.epsilon*(self.x-self.y+self.a)

class PassiveSystem:
    def __init__(self,state = -1, rate = 0.05):
        self.x = state
        self.state = state
        self.rate = rate
    def xFunction(self):
        return self.rate*(self.state - self.x)

class FHN_2:
    def __init__(self,aParam1 = 0, aParam2 = 0):
        self.e1 = FhnSystem(aParam1)
        self.e2 = PassiveSystem(aParam2)
        self.dt = 0.01
        self.D =  0.5
        self.visualize = False
        self.response = []
        self.X1=[]
        self.X2=[]
        self.Y1=[]
        self.T=[]
        self.time=0
        self.ts={}
    def get_response_ampl(self):
        return self.response
    def set_xoffset(self, offset):
        self.xoffset = offset
    def set_visualize(self, flag = True):
        self.visualize = flag
    def set_D(self, value):
        self.D = value
    def solve(self, Time = 100):
        counter = 0
        t = 0
        dt = self.dt
        X1 = self.X1
        X2 = self.X2
        Y1 = self.Y1
        T  = self.T
        x1_max = -100
        x2_max = -100
        y_max = -100
        while t<Time:
            self.time+=dt
            x1,y1,x2 = self.e1.x,self.e1.y,self.e2.x

            # if abs(t-100)<1e-10:
            #     x2 += self.xoffset
            #     # x1 = -0.4
            #     pass


            if  self.visualize and counter % 10 == 0:
                X1.append(self.e1.x)
                X2.append(self.e2.x)
                Y1.append(self.e1.y)
                T.append(self.time)

            if self.e2.x > x2_max:
                x2_max = self.e2.x

            if self.e1.x > x1_max:
                x1_max = self.e1.x
            if self.e1.y > y_max:
                y_max = self.e1.y

            rk1_0_x = dt*self.e1.xFunction()
            rk1_0_y = dt*self.e1.yFunction()
            rk2_0_x = dt*self.e2.xFunction()

            self.e1.x = x1+rk1_0_x/2.0
            self.e1.y = y1+rk1_0_y/2.0
            self.e2.x = x2+rk2_0_x/2.0
            rk1_1_x = dt*self.e1.xFunction()
            rk1_1_y = dt*self.e1.yFunction()
            rk2_1_x = dt*self.e2.xFunction()

            self.e1.x = x1+rk1_1_x/2.0
            self.e1.y = y1+rk1_1_y/2.0
            self.e2.x = x2+rk2_1_x/2.0
            rk1_2_x = dt*self.e1.xFunction()
            rk1_2_y = dt*self.e1.yFunction()
            rk2_2_x = dt*self.e2.xFunction()

            self.e1.x = x1+rk1_2_x
            self.e1.y = y1+rk1_2_y
            self.e2.x = x2+rk2_2_x
            rk1_3_x = dt*self.e1.xFunction()
            rk1_3_y = dt*self.e1.yFunction()
            rk2_3_x = dt*self.e2.xFunction()

            self.e1.x = x1 + (rk1_0_x+2*rk1_1_x+2*rk1_2_x+rk1_3_x)/6.0
            self.e1.y = y1 + (rk1_0_y+2*rk1_1_y+2*rk1_2_y+rk1_3_y)/6.0
            self.e2.x = x2 + (rk2_0_x+2*rk2_1_x+2*rk2_2_x+rk2_3_x)/6.0

            self.e1.x += dt*self.D*(x2-x1)
            self.e2.x += dt*self.D*(x1-x2)




            counter += 1
            t += dt

        self.ts['time'] = T
        self.ts['x']=X1
        self.ts['z']=X2
        self.ts['y']=Y1
        self.response = [x1_max,x2_max,y_max]

def build_curve():
    print "Starting"
    xoffset_min = 0.3
    xoffset_max = 1.5
    steps=100
    Dmax = 0.9
    Dmin = 0.1
    dsteps = 30
    result=[]
    x2_ss = 0
    x2_rate = 0.05
    for i in range(dsteps):
        xoffset = xoffset_min
        response_ampl = 0
        D = Dmin + (Dmax-Dmin)/float(dsteps-1)*i
        thresh = 1
        while response_ampl < thresh and xoffset < xoffset_max:
            f  = FHN_2(aParam1 = 0.5, aParam2 = x2_ss)
            f.e2.rate = x2_rate
            xoffset += (xoffset_max-xoffset_min)/steps
            f.set_xoffset(xoffset)
            f.set_D(D)
            f.solve(300)
            f.e2.x += xoffset
            f.solve(200)
            response_ampl = f.get_response_ampl()[0]


        if D < Dmax:
            print "FOUND: ","D=",D,"Xoffset=",xoffset
            result.append((str(D),str(xoffset)))
        else:
            print "NOT FOUND:", "D=",D
            result.append((str(D),"0"))

        print i
        sys.stdout.flush()

    print "opening file"
    rst_file = open("rst_" + str(x2_ss)+"_"+str(x2_rate)+".txt","w");
    print "file opened"
    print os.getcwd()
    rst_file.write("\n".join([" ".join(x) for x in result]))
    rst_file.close()
def single_run():
    a = 0.5
    zeta = 0
    alpha=0.05
    f = FHN_2(aParam1 = a, aParam2 = 0)
    f.e2.rate=alpha
    f.e2.state = zeta
    f.set_xoffset(0.6)
    f.set_D(0.9)
    f.solve(1000)
    f.set_visualize(True)
    f.e2.x += f.xoffset
    f.solve(1500)

    fig=plt.figure()
    ax = Axes3D(fig)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    ax.set_zticks([])
    ax.view_init(90,-90)
    X = np.arange(-2, 2, 0.1)
    Y = np.arange(-2, 2, 0.1)
    X, Y = np.meshgrid(X, Y)
    d=f.D
    Z = X-X*X*X/3.0+d*(Y-X)
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=True,alpha=0.6)
    Z = X+a
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=True,alpha=0.3)
    Z = (alpha*zeta+d*X)/(alpha+d)

    surf = ax.plot_surface(X, Z, Y, rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=True,alpha=0.8)


    x = -1.2
    z = (a+x**3/3.0+d*x)/d
    y = x+a
    p1 = [x, z, y]
    x = -0.9
    z = (a+x**3/3.0+d*x)/d
    y = x+a
    p2 = [x, z, y]
    ax.plot([p1[0],p2[0]],[p1[1],p2[1]],[p1[2],p2[2]],'b',lw=3)

    x = np.linspace(-2,2,50)
    z = (alpha*zeta+d*x)/(alpha+d)
    y = (x-x**3/3+d*(z-x))
    ax.plot(x,z,y,'m',lw=3)

    # xm = -np.sqrt(1-d)
    # p1 = [xm, -2, xm-xm**3/3.0+d*(-2-xm)]
    # p2 = [xm, 2, xm-xm**3/3.0+d*(2-xm)]
    # ax.plot([p1[0],p2[0]],[p1[1],p2[1]],[p1[2],p2[2]],'r')
    # p1 = [-2,-2,-2+8/3.0]
    # p2 = [2,2,2-8/3.0]
    # ax.plot([p1[0],p2[0]],[p1[1],p2[1]],[p1[2],p2[2]],'m')
    # fig.colorbar(surf, shrink=0.5, aspect=5)

    # plt.plot(f.ts['x'],f.ts['z'],f.ts['y'],'g',lw=2)
    rstx=[]
    rsty=[]
    rstz=[]
    xts = f.ts['x']
    yts = f.ts['y']
    zts = f.ts['z']
    xprev=0
    yprev=0
    yref_prev=0
    zprev=0
    for i in range(len(xts)):

        if i == 0:
            x_prev = xts[i]
            z_prev = zts[i]
            y_prev = yts[i]
            yref_prev = x_prev-x_prev**3/3+d*(z_prev-x_prev)
            rstx.append([x_prev])
            rsty.append([y_prev])
            rstz.append([z_prev])
            continue
        x = xts[i]
        z = zts[i]
        y = yts[i]
        y_ref = x-x**3/3+d*(z-x)

        if np.sign(y-y_ref)*np.sign(y_prev-yref_prev) > 0:
            rstx[-1].append(x)
            rsty[-1].append(y)
            rstz[-1].append(z)
        else:
            rstx.append([x])
            rsty.append([y])
            rstz.append([z])

        x_prev = x
        y_prev = y
        z_prev = z
        yref_prev = y_ref

    for i in range(len(rstx)):
        rx = rstx[i]
        ry = rsty[i]
        rz = rstz[i]
        x = rx[len(rx)/2]
        y = ry[len(ry)/2]
        z = rz[len(rz)/2]
        y_ref = x-x**3/3+d*(z-x)
        if y > y_ref:
            plt.plot(rx,rz,ry,'r',lw=2)
        else:
            plt.plot(rx,rz,ry,'g',lw=2)
    plt.plot([f.ts['x'][-1]],[f.ts['z'][-1]],[f.ts['y'][-1]],'r.',lw=5,markersize=20)
    plt.plot([f.ts['x'][-0]],[f.ts['z'][0]],[f.ts['y'][0]],'y.',lw=5,markersize=20)
    plt.show()

def main():
    single_run()
    # build_curve()

if __name__ == '__main__':
    main()
