import numpy as np


Angle=np.array([-60,-50,-45,-40,-30,-20,-10,-5,0,5,10,20,30,40,45,50,60])

torque=np.array([11.80,10.51,10.43,9.78,8.89,13.74,11.16,10.51,-1.62,-10.03,-10.67,-10.31,-8.89,-10.03,-10.67,-11.40,-12.37])

wi=np.polyfit(Angle, 1000*torque, 16)
poly_wind=np.poly1d(wi)

t=poly_wind(0 )
