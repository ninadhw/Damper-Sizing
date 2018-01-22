import numpy as np
import matplotlib.pyplot as mpl

#curve fitting fr stiffning damper curve fitting
speed = np.array([-30,-25,-20,-15,-10,-4,0,4,10,15,20,25,30])    
dampingforce=np.array([-500,-240,-120,-70,-30,-20,0,20,30,70,110,160,300])
z = np.polyfit(speed, dampingforce , 6)
poly=np.poly1d(z)

#dampergeometry and kinematics calculations
GL=np.array([0,0])                          #Ground level coordinate  (base of VP) (x,y coordinate)
PL=np.array([0,1500])                       #Pivot level (centre of rotation) (x,y coordinate)

DMP=np.array([0,630])                       #Damper mounting point (fom the Pivot Level) (x,y coordinate)
DML=PL-DMP                                  #Damper mounting level (from pivot reference)
DLA=140                                     #Damper lever arm in mm
LOS=76                                      #Lever arm offset
f=0.63                                      #Natural frequency in Hz
Th=5                                        #Maximum amplitude of torsional vibration in degree
Th=np.deg2rad(Th)

l=np.sqrt(DLA**2+LOS**2)                    #Actual lever arm length from centre of rotation
#print (l)

OA = 90 - np.rad2deg(np.arctan2(DLA,LOS))   #Offset angle as the lever arm is mounted at the bottom of TT
OA = np.deg2rad(OA)                         #Offset angle below the positive x axis

t=np.linspace(0,(1/(0.63)),100,True)
OMEGA=2*(np.pi)*f
THETA=OMEGA*t
AD= Th*np.sin(THETA)                        #Angular displacement: degree change through time
AS= Th*OMEGA*np.cos(THETA)                  #Angular Speed: a*omega*cos (theta)
AA= -Th*OMEGA**2*np.sin(THETA)              #Trial's sake angular accelaration

AD_prime= (AD) + OA*np.ones((len(AD))) #Adding the offset angle for location of coupling point

AD_prime2= (AD) - OA*np.ones((len(AD)))

CPx=PL[0] + (l*np.cos(AD_prime))                  #coupling point x coordinate
CPy=PL[1] - (l*np.sin(AD_prime))                 #coupling point y coordinate

CPx2=PL[0] - (l*np.cos(AD_prime2))                  #coupling point x coordinate
CPy2=PL[1] + (l*np.sin(AD_prime2))                 #coupling point y coordinate


#next few lines only for debugging
#AD_deg=np.rad2deg(AD_prime)                     #convert to degrees
#alpha=np.arctan(((CPy-DML[1]*np.ones(len(CPy)))/(CPx-DML[0]*np.ones(len(CPx))))) #angle of damper axis with horizontal
#IA=(180- (90*np.ones(len(alpha))-np.rad2deg(alpha)) - (90*np.ones(len(AD_prime))-np.rad2deg(AD_prime)))         #Internal angle in degrees between damper and lever arm
#IA=np.deg2rad(IA)


#the damper length is calculated using 2 independant methods
damper_l= np.sqrt( (CPx-DML[0]*np.ones(len(CPx)))**2 + (CPy-DML[1]*np.ones(len(CPy)))**2)
damper_l2= np.sqrt( (CPx2-DML[0]*np.ones(len(CPx)))**2 + (CPy2-DML[1]*np.ones(len(CPy)))**2)
l_damp= np.sqrt(((l*np.cos(AD_prime))**2+(630 - (l*np.sin(AD_prime)))**2))

#damper velocity calculated using 2 independamt methods
Lin_Vel = -1* l * AS *  ((CPx - DML[0])*np.sin(AD_prime) + (CPy - DML[1])*np.cos(AD_prime))/(damper_l)
Lin_Vel2 =1* l * AS *  ((CPx2 - DML[0])*np.sin(AD_prime2) + (CPy2 - DML[1])*np.cos(AD_prime2))/(damper_l2)

d_vel=np.empty(len(damper_l))
for i in range(1, (len(damper_l)-1)):
    d_vel[i]= (damper_l[i+1]-damper_l[i-1])/(t[i+1]-t[i-1])

d_vel[(len(damper_l)-1)]=0

b=np.empty(len(d_vel))
for i in range(0, len(d_vel)-1):
    b[i]=(d_vel[i]/Lin_Vel[i])


i=0
f= np.empty(len(d_vel))
f2=np.empty(len(d_vel))
for x in d_vel:
    f[i]=poly(Lin_Vel[i])
    f2[i]=poly(Lin_Vel2[i])
    i=i+1

#plotting damper force response based on linear velocity
#mpl.figure(6)
#mpl.plot(f)
#mpl.show()

fx_cap = (CPx-DML[0]*np.ones(len(CPx)))/damper_l
fy_cap = (CPy-DML[1]*np.ones(len(CPy)))/damper_l

rx_cap = (CPx - PL[0])/l
ry_cap = (CPy - PL[1])/l

fx2_cap = (CPx2-DML[0]*np.ones(len(CPx)))/damper_l2
fy2_cap = (CPy2-DML[1]*np.ones(len(CPy)))/damper_l2

rx2_cap = (CPx2 - PL[0])/l
ry2_cap = (CPy2 - PL[1])/l



M = f * l * (fy_cap*rx_cap - fx_cap*ry_cap)
M2= f2 * l * (fy2_cap*rx2_cap - fx2_cap*ry2_cap)

Mt=M+M2

#plotting damper moment response against time
mpl.plot(t,damper_l,'b--')
mpl.plot(t,damper_l2,'r--')
mpl.show()
#mpl.plot(CPx2,CPy2,'r')
#mpl.plot(CPx,CPy)
#mpl.show()
