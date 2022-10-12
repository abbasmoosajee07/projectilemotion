"""
===============================================================================
    # University of Birmingham
    # Labs and Data Analysis 2
    # Portfolio 2: Particle Tracking
    # Written by:   Abbas Moosajee
    # Written on:   10/04/2022
    # Contact:      AHM080@student.bham.ac.uk
    # An Empirical Model of a bouncing ball
===============================================================================
"""
import numpy as np


def empirical_model(initial_conditions, RE1):

    # Initial position of particle in m
    r_xy = initial_conditions["Coordinates"]
    # Velocity of Particle in m s^-1
    v = initial_conditions["Velocity"]
    # Angle of particle relative to horizontal in degrees
    theta_d = initial_conditions["Launch Angle"]
    # Restitution co-efficient, non-dimensional
    RE = RE1[0]
    N = 1000
    g = 9.81
    counter = 0

    # Calculates initial parameters at t=0
    t_c = 0
    theta = np.deg2rad(theta_d)
    vx = v * np.cos(theta)
    vy = v * np.sin(theta)
    rx = r_xy[0]
    ry = r_xy[1]

    # Save the simulated times, velocities and positions in growable lists

    tt = [t_c]
    v_x = [vx]
    v_y = [vy]
    r_x = [rx]
    r_y = [ry]

    while counter < N:
        counter = counter+1

        # Determines max height based on velocity
        ryn = ((vy**2)/(2*g))+ry 
        # Calculates total time based on y velocity
        tl = (vy)/g
        tn = t_c+tl   
        vyn = 0  # (tl*-g)            # y velocity when ball reaches ground
        vxn = vx
        rxn = rx + (vxn*tl)

        # Redefining new variables to current
        t_c = tn
        vx = vxn
        vy = vyn
        rx = rxn
        ry = ryn

        # Save velocities / positions in our lists
        tt.append(t_c)
        v_x.append(vx)
        v_y.append(vy)
        r_x.append(rx)
        r_y.append(ry)

        # determines time to fall to ground from max height
        tl = np.sqrt((2*ryn)/g)
        ryn = (0)  # new y position ground
        vyn = -RE * (tl*-g)  # velocity after bounce
        tn = t_c+tl
        vxn = vx
        rxn = rx + (vxn*tl)
        if abs(vyn) < 0.05:
            break
        # Redefining new variables to current
        t_c = tn
        vx = vxn
        vy = vyn
        rx = rxn
        ry = ryn

        # Save velocities / positions in our lists
        tt.append(t_c)
        v_x.append(vx)
        v_y.append(vy)
        r_x.append(rx)
        r_y.append(ry)

    # Transform the lists into NumPy arrays
    tt = np.array(tt)
    r_x = np.array(r_x)
    r_y = np.asarray(r_y)
    v_x = np.array(v_x)
    v_y = np.array(v_y)

    Data_Table = [tt, r_x, r_y, v_x, v_y]
    return Data_Table
# %%
