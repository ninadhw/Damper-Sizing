#calculating forcing function due to wind:
import numpy as np
import matplotlib.pyplot as mpl

#suction function (coeffecient)
angle=np.array([-5,-10,-20,-30,-40,-45,-50,-60])      #angle in degrees
Cv_suction = np.array([-0.62, -0.66, -1, -1.1, -1.24, -1.32, -1.41, -1.53])
cop_suction= np.array([0.3,0.3,0.35,0.4,0.4,0.4,0.4,0.4])
pc_suction = np.polyfit(angle,Cv_suction,2)
pcop_suction = np.polyfit(angle, cop_suction,3)

#pressure function (coeffecient)
angle=np.array([0,5,10,20,30,40,45,50,60])      #angle in degrees
Cv_pressure = np.array([0.14,0.65,0.69,1.1,1.1,1.21,1.29,1.3,1.46])
cop_pressure = np.array([0.04,0.3,0.3,0.3,0.4,0.4,0.4,0.4,0.4])
pc_pressure = np.polyfit(angle,Cv_pressure,2)
pcop_pressure = np.polyfit(angle, cop_pressure,3)

f=0.63                                      #Natural frequency in Hz
Th=5                                       #Maximum amplitude of torsional vibration in degree
WindVelocity = 47                           #m/s
ModuleLength = 2                            #m (2 for 1P, 3 for 3L)
Table_Length = 60                           #m (table length)

###########################################################################################

t=np.linspace(0,(2/(0.63)),200,True)
OMEGA=2*(np.pi)*f
THETA=OMEGA*t
AD= Th*np.sin(THETA)                        #Angular displacement: degree change through time
AS= Th*OMEGA*np.cos(THETA)                  #Angular Speed: a*omega*cos (theta)
AA= -Th*OMEGA**2*np.sin(THETA)              #Trial's sake angular accelaration

###########################################################################################
#based on the angle of the table, the force caused by the wind
i=0
Load_Torque=[]
Cop_array=[]
CV_array=[]
for x in AD:
    if x<0:
        #case 1: suction
        Cv=pc_suction[0]*x**2+pc_suction[1]*x+pc_suction[2]
        #x=np.abs(x)
        CoP=pcop_suction[0]*x**3+pcop_suction[1]*x**2+pcop_suction[2]*x+pcop_suction[3]
        Load_Torque.append(0.613 * WindVelocity**2 * Cv * (0.5-CoP) * ModuleLength * Table_Length/2 * 2 )
        Cop_array.append(CoP)
        CV_array.append(Cv)
        #################################################################################                   
    else:
        #case 2: pressure
        Cv=pc_pressure[0]*x**2+pc_pressure[1]*x+pc_pressure[2]
        CoP=pcop_pressure[0]*x**3+pcop_pressure[1]*x**2+pcop_pressure[2]*x+pcop_pressure[3]
        Load_Torque.append(0.613 * WindVelocity**2 * Cv * (0.5-CoP) * ModuleLength * Table_Length/2 * 2 )
        Cop_array.append(CoP)
        CV_array.append(Cv)
        #################################################################################
mpl.plot(t,AD)
mpl.plot(t,Cop_array,'g')
mpl.plot(t,CV_array,'r')
#mpl.plot(t,Load_Torque)
mpl.show()

print('COP \t CV \t Angle')
for i in range(0, len(AD)-1):
    print(Cop_array[i],'\t',CV_array[i],'\t', AD[i])
