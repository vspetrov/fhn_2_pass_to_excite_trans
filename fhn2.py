import os, sys
import matplotlib.pyplot as plt
class FhnSystem:
    def __init__(self,aParam = 0):
        self.epsilon = 0.1
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
    def set_xoffset(self, offset):
        self.xoffset = offset
    def solve(self, Time = 100):
        counter = 0
        t = 0
        dt = self.dt
        X1 = []
        X2 = []
        T  = []
        while t<Time:
            x1,y1,x2 = self.e1.x,self.e1.y,self.e2.x

            if abs(t-100)<1e-10:
                x2 = self.xoffset
                # x1 = -0.4
                pass

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

            if counter % 10 == 0:

                X1.append(self.e1.x)
                X2.append(self.e2.x)
                T.append(t)
            counter += 1
            t += dt

        plt.plot(T,X1,'-',lw=2)
        plt.plot(T,X2,'-',lw=2,color='red')
        plt.show()

def main():
    print "Starting"
    f  = FHN_2(aParam1 = 0.5, aParam2 = 0)
    f.set_xoffset(0.5)
    f.solve(200)


if __name__ == '__main__':
    main()
