-- Numerical Simulation of Rigid Body Rotation --

This project implements a numerical simulation of rigid body rotation dynamics using quaternions and the Runge-Kutta 4th order method. The code solves Euler's rotation equations coupled with quaternion kinematics to model the motion of a rotating body under external forces and moments.

-- Features --

1) Quaternion-based rotation representation
2) RK4 integration for both angular velocity and orientation
3) Energy and angular momentum calculations
4) 3D visualization of the rotating body's axes
5) GIF animation output of the simulation

-- Output --

The simulation produces:
1) A 3D animation showing the evolution of the body's principal axes
2) Console output of angular momentum and energy values at each frame
3) Energy conservation tracking (kinetic + potential)

-- Customization --

You can modify:
1) The body's inertial properties (I_values)
2) Initial conditions (q0, w0)
3) External forces and torques (F, M_ext)
4) Simulation parameters (t_max, dt)


-- Energy Conservation --

The code tracks and prints:
Kinetic energy (T), Potential energy (U), Total energy (E = T + U)
Energy conservation can be used to verify the numerical accuracy of the simulation.

-- Visualization --

The animation shows:
1) Red arrow: body's x-axis
2) Green arrow: body's y-axis
3) Blue arrow: body's z-axis
The axes rotate according to the computed quaternion orientation.
