README:

1. Damper is modeled as a linear viscous damper with cubic polynomial curve representing the velocity-force relation. Find details in excel sheet "MAHINDRA SUSTEN -DAMPER FORCES.xls".

2. The geometry of the tracker table creates a non-linear relationship between the angle of the tracker table and position of the damper.

3. The geometry of the tracker table creates a non-linear relationship between the angular velocity of the table and linear velocity of the damper.
	3a. the linear velocity of the damper causes force (resistance to motion) in the damper. The relation in between the linear velocity and force cause is explained in point 1.
	3b. the non-linear relationship between the angular velocity of the tracker table and the linear velocity of the damper is a function of i. position of the table (theta) ii. angular velocity of table (theta dot). The equation is derived on paper from basic kinematics of 4 bar system.

4. The inertia of the system is calculated theoretcaly as the mass-moment of inertia of solar modules multiplied by the number of solar modules (about the axis of rotation of tracker). Find details in "Inertia and Stiffness data.xls".

5. The stiffness of the spring is calculated using classical mechanics formula for torsional rigidity of non-circular hollow sections. Find detail value in "Inertia and Stiffness data.xls"

6. "dopri5" solver is a scipy solver for solving first order linear (linearized) differential equations using 4th order Runge-Kutta method.

7. The analysis is done on half of tracker table assuming rigid fixity in the middle. The mathematical model is for a lumped spring-mass-damper system.

8. File "Damper_System_Solver".py" is a basic model for the geometry of the tracker table and the damper model. It gives the length of damper assuming a sinusoidal displacement to the tracker table. It also gives the force caused by this displacement and the moment caused by this moment for the two counter-placed dampers arranged accordinig to the prescribed geometry.

9. File "Integrated_Code_wind_trial.py" is an integrated code which has the spring, inertia and damper system working together. Also used is a crude wind forcing function. **The forcing function is in doubt currently.**

10. Run all codes in any python IDE environment which has "matplotlib","numpy" and "scipy" modules installed. All modules and the python environment are free and open-source.