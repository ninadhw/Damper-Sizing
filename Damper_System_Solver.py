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

t=np.linspace(0,(3/(0.63)),10000,True)
OMEGA=2*(np.pi)*f
THETA=OMEGA*t
AD= Th*np.sin(THETA)                        #Angular displacement: degree change through time
AS= Th*OMEGA*np.cos(THETA)                  #Angular Speed: a*omega*cos (theta)
AA= -Th*OMEGA**2*np.sin(THETA)              #Trial's sake angular accelaration

AD_prime= (AD) + OA*np.ones((len(AD))) #Adding the offset angle for location of coupling point

CPx=PL[0] + (l*np.cos(AD_prime))                  #coupling point x coordinate
CPy=PL[1] - (l*np.sin(AD_prime))                 #coupling point y coordinate

#next few lines only for debugging
#AD_deg=np.rad2deg(AD_prime)                     #convert to degrees
#alpha=np.arctan(((CPy-DML[1]*np.ones(len(CPy)))/(CPx-DML[0]*np.ones(len(CPx))))) #angle of damper axis with horizontal
#IA=(180- (90*np.ones(len(alpha))-np.rad2deg(alpha)) - (90*np.ones(len(AD_prime))-np.rad2deg(AD_prime)))         #Internal angle in degrees between damper and lever arm
#IA=np.deg2rad(IA)


#the damper length is calculated using 2 independant methods
damper_l= np.sqrt( (CPx-DML[0]*np.ones(len(CPx)))**2 + (CPy-DML[1]*np.ones(len(CPy)))**2)
l_damp= np.sqrt(((l*np.cos(AD_prime))**2+(630 - (l*np.sin(AD_prime)))**2))

#damper velocity calculated using 2 independamt methods
Lin_Vel = -1* l * AS *  ((CPx - DML[0])*np.sin(AD_prime) + (CPy - DML[1])*np.cos(AD_prime))/(damper_l)

d_vel=np.empty(len(damper_l))
for i in range(1, (len(damper_l)-1)):
    d_vel[i]= (damper_l[i+1]-damper_l[i-1])/(t[i+1]-t[i-1])

d_vel[(len(damper_l)-1)]=0

b=np.empty(len(d_vel))
for i in range(0, len(d_vel)-1):
    b[i]=(d_vel[i]/Lin_Vel[i])


i=0
f=np.empty(len(d_vel))
for x in d_vel:
    f[i]=poly(Lin_Vel[i])
    i=i+1

#plotting damper force response based on linear velocity
mpl.figure(6)
mpl.plot(f)
mpl.show()
