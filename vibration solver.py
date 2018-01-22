#customary function trial
from scipy.integrate import ode
import matplotlib.pyplot as mpl
def function(t , y):
    return[y[1], -y[0] - 0.5*y[1  ] + 2 ]

r=ode(function, jac=None).set_integrator('dopri5', nsteps=100)
r.set_initial_value([10,0], 0)
t1=100
dt=0.1
i=0

disp=[]
vel=[]

while r.successful() and r.t<t1:
    b=r.integrate(r.t+dt)
    disp.append(b[0])
    vel.append(b[1])
    i=i+1
mpl.plot(disp)
mpl.plot(vel,'r')
mpl.show()

