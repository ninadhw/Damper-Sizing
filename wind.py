import numpy as np
import matplotlib.pyplot as mpl

Angle=np.array([-60,-50,-45,-40,-30,-20,-10,-5,0,5,10,20,30,40,45,50,60])

torque=np.array([11.80,10.51,10.43,9.78,8.89,13.74,11.16,10.51,-1.62,-10.03,-10.67,-10.31,-8.89,-10.03,-10.67,-11.40,-12.37])


Angle=np.array([-1.04720,
-0.87266,
-0.78540,
-0.69813,
-0.52360,
-0.34907,
-0.17453,
-0.08727,
0.00000,
0.08727,
0.17453,
0.34907,
0.52360,
0.69813,
0.78540,
0.87266,
1.04720
])

torque=np.array([23608.02,
21020.84,
20859.15,
19565.55,
17786.87,
27488.80,
22314.43,
21020.84,
-3233.98,
-20050.65,
-21344.24,
-20616.60,
-17786.87,
-20050.65,
-21344.24,
-22799.53,
-24739.92
])





wi=np.polyfit( Angle, -1* torque, 9)
poly_wind=np.poly1d(wi)

angle=np.linspace(-1.05,1.05,100)
torque=[]

for x in angle:
    t=poly_wind(x)
    torque.append(t)

mpl.plot(angle,torque)
mpl.plot(angle, angle*6666.7,'r')
mpl.plot (angle, angle*66666.7,'g')
mpl.show()
t=poly_wind(0)
print(t)
