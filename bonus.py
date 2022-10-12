"""
===============================================================================
    # University of Birmingham
    # Labs and Data Analysis 2
    # Portfolio 2: Particle Tracking
    # Written by:   Abbas Moosajee
    # Written on:   10/04/2022
    # Contact:      AHM080@student.bham.ac.uk
    # Tracking a falling ball
===============================================================================
"""
import time
import numpy as np
import matplotlib.pyplot as plt
import seaborn; seaborn.set_style("whitegrid")

from tracker_func  import Particle_Tracker
from empirical_bq  import empirical_model
from dem_bq        import dem_model
start = time.time()

lt=1

Red=[(150, 100, 10),(185, 255, 255)] # Red Boundaries
Cricket_Ball={"Mass":160E-03, "Diameter":6.86E-02, "Color":Red}

R2,I2,RE2=Particle_Tracker('RDropTest_125.mp4',lt,Cricket_Ball,7,1.25,220)

Blue=[(94,100,100),(184, 255, 255)]  # Blue Boundaries
Rubber_Ball={"Mass":80E-03, "Diameter":6.86E-02, "Color":Blue}

R3,I3,RE3=Particle_Tracker('BDropTest_125.mp4',lt,Rubber_Ball,7,1.25,50)

Green=[(29, 86,  6),(64, 255, 255)]
Tennis_Ball={"Mass":56.93E-03, "Diameter":6.86E-02, "Color":Green}

R1,I1,RE1=Particle_Tracker('DropTest_W125_R2.mp4',lt,Tennis_Ball,2,1.25,88.7)


# %%
T1=empirical_model(I1,RE1)
T2=dem_model(Tennis_Ball,I1,RE1)
P1=R1
tp=P1[0]
ParticleData = plt.figure(figsize=[18, 18])
plt.subplot(3, 3, 1)
plt.plot(P1[1], P1[2], 'b-')
plt.xlabel('Distance(m)')
plt.ylabel('Height(m)')
plt.subplot(3, 3, 2)
plt.plot(tp, P1[1], 'b-')
plt.xlabel('Time(s)')
plt.ylabel('Distance(m)')
plt.subplot(3, 3, 3)
plt.plot(tp, P1[2], 'b-')
plt.xlabel('Time(s)')
plt.ylabel('Height(m)')

plt.subplot(3, 2, 3)
plt.plot(tp[1:], P1[4], 'c.')
plt.xlabel('Time(s)')
plt.ylabel('Vertical Velocity(m s^-1)')
plt.subplot(3, 2, 4)
plt.plot(tp[2:], P1[6], 'k.')
plt.xlabel('Time(s)')
plt.ylabel('Vertical Acceleration(m s^-2)')

plt.subplot(3, 1, 3)
plt.plot(tp[1:], P1[7], 'g-')
plt.plot(tp[1:], P1[8], 'r-')
plt.plot(tp[1:], P1[9], 'k-')
plt.xlabel('Time(s)')
plt.ylabel('Energy(J)')
plt.legend(['GPE', 'KE', 'TE'])

# %%
ParticleDataEB = plt.figure(figsize=[18, 9])
plt.subplot(1, 3, 1)
plt.plot(P1[1], P1[2], 'b-')
plt.plot(T1[1], T1[2], 'r-')
plt.plot(T2[1], T2[2], 'g-')

plt.xlabel('Distance(m)')
plt.ylabel('Height(m)')
plt.subplot(1, 3, 2)
plt.plot(tp, P1[2], 'b-')
plt.plot(T1[0], T1[2], 'r-')
plt.plot(T2[0], T2[2], 'g-')
plt.xlabel('Time(s)')
plt.ylabel('Height(m)') 
plt.legend(['Actual', 'Empirical', 'DEM'])

plt.subplot(1, 3, 3)
plt.plot(tp[1:], P1[4], 'b-')
plt.plot(T1[0], T1[4], 'ro')
plt.plot(T2[0], T2[4], 'g.')
plt.xlabel('Time(s)')
plt.ylabel('Velocity(m/s)')
end = time.time()
print(f"Execution took: {np.round((end - start),2)} s")
