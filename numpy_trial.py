import numpy as np
import matplotlib.pyplot as mpl

#curve fitting fr stiffning damper curve fitting
x=np.array([4,10,15,20,25,30])
y_extension=np.array([20,30,70,110,160,300])
y_compression=np.array([20,30,70,120,240,500])

p_extension = np.polyfit(x,y_extension,3)
p_compression = np.polyfit(x,y_compression,3)

#dampergeometry and kinematics calculations

GL=np.array([0,0])                          #Ground level coordinate  (base of VP) (x,y coordinate)
PL=np.array([0,1500])                       #Pivot level (centre of rotation) (x,y coordinate)

DMP=np.array([0,630])                       #Damper mounting point (fom the Pivot Level) (x,y coordinate)
DML=PL-DMP                                  #Damper mounting level (from ground reference
DLA=140                                     #Damper lever arm in mm
LOS=76                                      #Lever arm offset
f=0.63                                      #Natural frequency in Hz
Th=5                                       #Maximum amplitude of torsional vibration in degree

l=np.sqrt(DLA**2+LOS**2)                    #Actual lever arm length from centre of rotation
print (l)

OA = 90 - np.rad2deg(np.arctan2(DLA,LOS))   #Offset angle as the lever arm is mounted at the bottom of TT
print(OA)
OA = np.deg2rad(OA)                         #Offset angle below the positive x axis

t=np.linspace(0,(3/(0.63)),1000,True)
OMEGA=2*(np.pi)*f
THETA=OMEGA*t
AD= Th*np.sin(THETA)                        #Angular displacement: degree change through time
AS= Th*OMEGA*np.cos(THETA)                  #Angular Speed: a*omega*cos (theta)
AA= -Th*OMEGA**2*np.sin(THETA)              #Trial's sake angular accelaration

mpl.plot(t,AD)
mpl.plot(t,AS)
mpl.plot(t,AA)
mpl.show()

AD_prime= np.deg2rad(AD) + OA*np.ones((len(AD))) #Adding the offset angle for location of coupling point
print(AD_prime)

CPx=PL[0]+(l*np.cos(AD_prime))                  #coupling point x coordinate
CPy=PL[1]- (l*np.sin(AD_prime))                 #coupling point y coordinate
mpl.plot(CPx,CPy)
mpl.show()

#next few lines only for debugging
AD_deg=np.rad2deg(AD_prime)                     #convert to degrees
mpl.plot(AD_deg)
mpl.show()
print(np.mean(AD_deg))

alpha=np.arctan(((CPy-DML[1]*np.ones(len(CPy)))/(CPx-DML[0]*np.ones(len(CPx))))) #angle of damper axis with horizontal
print(np.rad2deg(AD_prime))
print(np.rad2deg(alpha))

IA=(180- (90*np.ones(len(alpha))-np.rad2deg(alpha)) - (90*np.ones(len(AD_prime))-np.rad2deg(AD_prime)))         #Internal angle in degrees between damper and lever arm
IA=np.deg2rad(IA)

damper_l= np.sqrt( (CPx-DML[0]*np.ones(len(CPx)))**2 + (CPy-DML[1]*np.ones(len(CPy)))**2)
Lin_Vel = -AS * l * ( (-np.sin (AD_prime)/np.cos (alpha)) + (np.cos(AD_prime)*np.sin(alpha)/(np.cos(alpha)*np.cos(alpha))*(((- l**2 + l* np.sqrt(DML[1]**2 + DML[0]**2)) * np.ones(len(CPx)))/(damper_l**2) ) ) )/1000#print(Lin_Vel)

mpl.figure(1)
mpl.plot(Lin_Vel)
mpl.plot(1000*np.rad2deg(alpha)-1000*np.mean(alpha)*np.ones(len(alpha)))
mpl.plot(1000*AD_deg - 1000*np.mean(AD_deg)*np.ones(len(AD_deg)))
mpl.show()

mpl.figure(2)
mpl.plot(np.rad2deg(alpha))
mpl.show()

mpl.figure(3)
mpl.plot(damper_l)
mpl.show()

l_damp= np.sqrt(((l*np.cos(AD_prime))**2+(630 - (l*np.sin(AD_prime)))**2))
mpl.figure(4)
mpl.plot(l_damp)
mpl.show()

t_interval = t[1]-t[0]
d_vel=np.empty(len(damper_l)-1)

for i in range(0, (len(damper_l)-1)):
    d_vel[i]= (damper_l[i+1]-damper_l[i])/t_interval
    #print('damperl',damper_l[i],'damperl+1',damper_l[i+1],'tinterval',t_interval,'damper velocity',d_vel[i])
    #print (d_vel[i])

print("done...")
print(d_vel)
mpl.figure(5)
mpl.plot(d_vel)
mpl.show()

i=0
f=np.empty(len(d_vel))

for x in d_vel:
    if x<0:
        x=np.abs(x)
        f[i]= (0.05045257* x**3)+(-1.49989757 * x**2)+(17.17375005 * x)+(-30.90373659) #cubic curve fit
        #f[i]= np.abs(f[i])
        #if f[i]<0:
        #    print('ping',i)
    else:
        f[i]= (0.02137117* x**3)+(-0.60748134* x**2)+(9.23398964* x)+(-11.76072502) #cubic curve fit
        #f[i]= np.abs(f[i])
        #if f[i]<0:
        #    print('pong',i)
    #if f[i]<0:
        #f[i]=0
    i=i+1

mpl.figure(6)
mpl.plot(f)
mpl.show()
