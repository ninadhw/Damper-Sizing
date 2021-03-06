#Integrated Vibration Code with damping moment, inertia and spring force
from scipy.integrate import ode
import matplotlib.pyplot as mpl
import numpy as np

#curve fitting fr stiffning damper curve fitting
#speed = np.array([-30,-25,-20,-15,-10,-4,0,4,10,15,20,25,30])    
#dampingforce=np.array([-500,-240,-120,-70,-30,-20,0,20,30,70,110,160,300])
speed = np.array([-25,-20,-15,-10,-4,0,4,10,15,20,25])          #mm/s units
dampingforce=np.array([-90,-70,-50,-30,-20,0,20,30,50,70,90])   #kg units
z = np.polyfit(speed, 10*dampingforce , 6)                      #newton and mm/s units
poly=np.poly1d(z)


#wind forcing function
Angle=np.array([-60,-50,-45,-40,-30,-20,-10,-5,0,5,10,20,30,40,45,50,60])
torque=np.array([11.80,10.51,10.43,9.78,8.89,13.74,11.16,10.51,-1.62,-10.03,-10.67,-10.31,-8.89,-10.03,-10.67,-11.40,-12.37])
Angle=np.deg2rad(Angle)
torque=torque*2000

wi=np.polyfit(Angle, torque, 16)
poly_wind=np.poly1d(z)



#damper geometry and kinematics calculations (user inputs)
GL=np.array([0,0])                          #Ground level coordinate  (base of VP) (x,y coordinate)
PL=np.array([0,1500])                       #Pivot level (centre of rotation) (x,y coordinate)
DMP=np.array([0,630])                       #Damper mounting point (fom the Pivot Level) (x,y coordinate)
DLA=140                                     #Damper lever arm in mm
LOS=76                                      #Lever arm offset
f=0.63                                      #Natural frequency in Hz
Th=5                                        #Maximum amplitude of torsional vibration in degree
k_tt= 6666.67                               #Newton and meter units

I1=39000                                    #torque tube
I2=8110336                                  #rail
I3=241661000                                #module
I=(I1+I2+I3)/1000000

Ipl=249.8

#####calculating machine parameters
OA = 90 - np.rad2deg(np.arctan2(DLA,LOS))   #Offset angle as the lever arm is mounted at the bottom of TT
OA = np.deg2rad(OA)                         #Offset angle below the positive x axis
DML=PL-DMP                                  #Damper mounting level (from pivot reference)
Th=np.deg2rad(Th)                           #convert to radian
l=np.sqrt(DLA**2+LOS**2)                    #Actual lever arm length from centre of rotation
#OMEGA = (2 * f * np.pi)

def motion_equation(t , y):
    AD = y[0]
    AS = y[1]

    #damper 1
    DA = y[0] + OA
    CPx= PL[0] + (l*np.cos(DA))
    CPy= PL[1] - (l*np.sin(DA))
    damper_l= np.sqrt( (CPx-DML[0])**2 + (CPy-DML[1])**2)
    Lin_Vel = -1* l * AS *  ((CPx - DML[0])*np.sin(DA) + (CPy - DML[1])*np.cos(DA))/(damper_l)
    f=poly(Lin_Vel)
    fx_cap = (CPx-DML[0])/damper_l              #unit vectors in force direction
    fy_cap = (CPy-DML[1])/damper_l
    rx_cap = (CPx - PL[0])/l                    #unit vector in radius direction (this is for cross product
    ry_cap = (CPy - PL[1])/l

    #damper 2
    DA2 = y[0] - OA
    CPx2= PL[0] + (l*np.cos(DA2))
    CPy2= PL[1] - (l*np.sin(DA2))
    damper_l2= np.sqrt( (CPx2-DML[0])**2 + (CPy2-DML[1])**2)
    Lin_Vel2 = 1* l * AS *  ((CPx2 - DML[0])*np.sin(DA2) + (CPy2 - DML[1])*np.cos(DA2))/(damper_l2)
    f2=poly(Lin_Vel2)
    fx_cap2 = (CPx2 -DML[0])/damper_l2
    fy_cap2 = (CPy2 -DML[1])/damper_l2
    rx_cap2 = (CPx2 - PL[0])/l
    ry_cap2 = (CPy2 - PL[1])/l
    
    D_M = (f * l * (fy_cap*rx_cap - fx_cap*ry_cap))/1000 - (f2 * l * (fy_cap2*rx_cap2 - fx_cap2*ry_cap2))/1000 #damper moment
    K_M = k_tt * AD                             #spring force moment
    W_tor=poly_wind(AD)
    #print('angle:',AD,'torque',W_tor*2000,'spring moment', K_M)

    return[y[1], (-K_M   - W_tor*2000 + D_M )/(Ipl)]

r=ode(motion_equation, jac=None).set_integrator('dopri5', nsteps=100000)
r.set_initial_value([np.deg2rad(15), 0], 0)
t1=100
dt=0.01

i=0
disp=[]
vel=[]

while r.successful() and r.t<t1:
    b=r.integrate(r.t+dt)
    disp.append(b[0])
    vel.append(b[1])
    i=i+1
tzone=np.linspace(0,t1,(t1/dt))
mpl.plot(tzone,disp,'r')
mpl.plot(tzone,vel,'b--')
mpl.show()
