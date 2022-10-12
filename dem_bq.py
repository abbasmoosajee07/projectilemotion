"""
===============================================================================
    # University of Birmingham
    # Labs and Data Analysis 2
    # Portfolio 2: particle Tracking
    # Written by:   Abbas Moosajee
    # Written on:   10/04/2022
    # Contact:      AHM080@student.bham.ac.uk
    # A DEM model using the Euler's Numerical Solver
===============================================================================
"""
import numpy as np


def dem_model(particle, initial_conditions, cor):

    # %% User Defined Data

    # Environmental Factors
    g = 9.81        # Gravitational acceleration in m s^-2
    f_rho = 1.225   # Fluid Density in kg m^-3
    # Time step in s
    dt = 0.01
    # Maximum time span of calculations in s
    t_max = 10000

    # particle Properties
    D = particle["Diameter"]
    p_m = particle["Mass"]  # Mass in kg
    CD = 0.47       # Drag coefficient of a ball is taken as 0.47
    RE =  cor[0]             # Restitution co-efficient, non-dimensional

    # Maximum number of iterations for loop
    N = int(abs(t_max/dt))
    tt, r_x, r_y, V = (np.zeros(N) for _ in range(4))

    # Initial conditions
    # Initial position of particle in m
    r_xy = initial_conditions["Coordinates"]
    # Velocity of particle in m s^-1
    V[0] = (initial_conditions["Velocity"])
    # Angle of particle relative to horizontal in degrees
    THETA_D = initial_conditions["Launch Angle"]
    # %%
    Vol = 4/3 * (np.pi*(D/2)**3)  # Volume of Sphere in m^3
    P_a = (np.pi)*(D/2)**2  # Projected Area of Sphere in m^2
    fb = f_rho*g*Vol       # Buoyancy force  (N)
    w = p_m*g             # particle Weight (N)

    # Initial Calculations-Velocity, Acceleration, Drag Force and Energy of Paticle
    tt[0] = 0
    theta = np.deg2rad(THETA_D)
    v_x = V * np.cos(theta)
    v_y = V * np.sin(theta)
    r_x[0] = r_xy[0]
    r_y[0] = r_xy[1]
    FD_x = 0.5*CD*f_rho*P_a*v_x*v_x
    FD_y = 0.5*CD*f_rho*P_a*v_y*v_y
    a_x = (FD_x)/p_m
    a_y = (w-(FD_y+fb))/p_m
    GPE = p_m*g*r_y
    KE = p_m*0.5*v_y*v_y
    TE = GPE + KE
    # %% Main Loop for projectile motion
    for ts in range(1, N):
        tt[ts] = tt[ts-1] + dt

        # New X velocity, updated using x acceleration
        v_x[ts] = v_x[ts-1] + (dt * -a_x[ts-1])
        # New Y velocity, updated using y acceleration
        # v_y[ts] = v_y[ts-1] + (dt * -a_y[ts-1])

        # New X position, updated using x velocity
        r_x[ts] = r_x[ts-1] + (dt * v_x[ts-1])
        # New Y position, updated using y velocity
        r_y[ts] = r_y[ts-1] + (dt * v_y[ts-1])

        # Gravitational Potential Energy in J
        GPE[ts] = p_m*g*r_y[ts]
        # Kinetic Energy in J
        KE[ts] = p_m*0.5*v_y[ts-1]**2
        # Total Energy of particle, sum of GPE and KE
        TE[ts] = GPE[ts]+KE[ts]

        # Drag force, used to find acceleration in x and y components
        FD_x[ts] = 0.5*CD*f_rho*P_a*v_x[ts]**2
        FD_y[ts] = 0.5*CD*f_rho*P_a*v_y[ts]**2

        a_x[ts] = (FD_x[ts])/p_m
        a_y[ts] = (w-(FD_y[ts]+fb))/p_m

        # If condition that determines when projectile reaches the ground
        if r_y[ts] < 0:
            # New y velocity, as particle bounces back, and redefines y position
            v_y[ts] = -RE * v_y[ts-1]
            r_y[ts] = 0   # The  y position, as particle bounces back is 0
        else:
            # Find new y velocity, if no impact occurs
            v_y[ts] = v_y[ts-1]+(-a_y[ts-1])*dt
        if TE[ts] < 5e-3:
            break
        # Determnes when two consecutive velocities are same, and thus terminal
        if v_y[ts] is v_y[ts-1]:
            break

    # Truncates all arrays to the number of iterations actually calculated
    tt = tt[0:ts+1]
    r_x = r_x[0:ts+1]
    r_y = r_y[0:ts+1]
    v_x = v_x[0:ts+1]
    v_y = v_y[0:ts+1]
    a_x = a_x[0:ts+1]
    a_y = a_y[0:ts+1]
    FD_x = FD_x[0:ts+1]
    FD_y = FD_y[0:ts+1]
    GPE = GPE[0:ts+1]
    KE = KE[0:ts+1]
    TE = TE[0:ts+1]
    Data_Table = [tt, r_x, r_y, v_x, v_y]
    return Data_Table
