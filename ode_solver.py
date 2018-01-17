#ODE trial
from scipy.integrate import ode

y0, t0 = [1.0j, 2.0], 0

def f(t, y, arg1):
    return [1j*arg1*y[0] + y[1], -arg1*y[1]**2]

#def jac(t,y,arg1):
#    return[ [1j*arg1, 1], [0, -arg1*2*y[1]]]

r=ode(f).set_integrator('zvode', method='bdf')
r.set_initial_value(y0,t0).set_f_params(2.0)
t1=10
dt=0.01

while r.successful() and r.t <t1:
    print(r.t+dt, r.integrate(r.t+dt))
