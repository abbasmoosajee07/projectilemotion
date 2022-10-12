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
from pandas import pandas as pd
from tracker_func  import Particle_Tracker
start = time.time()
lt=0#live track
#%%
Green=[(29, 86,  6),(64, 255, 255)]
Tennis_Ball={"Mass":56.93E-03, "Diameter":6.86E-02, "Color":Green}
# Video Naming Convention- DropTest-Surface-Height-Run
# Data,Intial,RE=Particle_Tracker(Video,Particle,Iterationd,Height,ymin)
# Surface- Wooden
R1,I1,RE1=Particle_Tracker('DropTest_W125_R1.mp4',lt,Tennis_Ball,2,1.25,88.7)
R2,I2,RE2=Particle_Tracker('DropTest_W125_R2.mp4',lt,Tennis_Ball,2,1.25,88.7)
R3,I3,RE3=Particle_Tracker('DropTest_W75_R1.mp4',lt,Tennis_Ball,2,0.75,250)
R4,I4,RE4=Particle_Tracker('DropTest_W75_R2.mp4',lt,Tennis_Ball,2,0.75,250)
RE_WV125=(RE1[0]+RE2[0])/2
RE_WH125=(RE1[1]+RE2[1])/2
RE_WV75 =(RE3[0]+RE4[0])/2
RE_WH75 =(RE3[1]+RE4[1])/2

# Surface- Cardboard
R5,I5,RE5=Particle_Tracker('DropTest_Cb125_R1.mp4',lt,Tennis_Ball,2,1.25,173)
R6,I6,RE6=Particle_Tracker('DropTest_Cb125_R2.mp4',lt,Tennis_Ball,2,1.25,173)
R7,I7,RE7=Particle_Tracker('DropTest_Cb75_R1.mp4',lt,Tennis_Ball,2,0.75,300)
R8,I8,RE8=Particle_Tracker('DropTest_Cb75_R2.mp4',lt,Tennis_Ball,2,0.75,300)
RE_CbV125=(RE5[0]+RE6[0])/2
RE_CbH125=(RE5[1]+RE6[1])/2
RE_CbV75 =(RE7[0]+RE8[0])/2
RE_CbH75 =(RE7[1]+RE8[1])/2

# Surface- Carpet
R9,I9,RE9=Particle_Tracker('DropTest_Cp125_R1.mp4',lt,Tennis_Ball,2,1.25,240)
R10,I10,RE10=Particle_Tracker('DropTest_Cp125_R2.mp4',lt,Tennis_Ball,2,1.25,240)
R11,I12,RE11=Particle_Tracker('DropTest_Cp75_R1.mp4',lt,Tennis_Ball,2,0.75,240)
R12,I12,RE12=Particle_Tracker('DropTest_Cp75_R1.mp4',lt,Tennis_Ball,2,0.75,240)
RE_CpV125=(RE9[0]+RE10[0])/2
RE_CpH125=(RE9[1]+RE10[1])/2
RE_CpV75 =(RE11[0]+RE12[0])/2
RE_CpH75 =(RE11[1]+RE12[1])/2
# %%Restitution Table
RE_C=np.round([[RE_CbV75,RE_CbH75,RE_CbV125,RE_CbH125],
              [RE_CpV75,RE_CpH75,RE_CpV125,RE_CpH125],
              [RE_WV75,RE_WH75,RE_WV125,RE_WH125]],2)
column_names = pd.DataFrame([["75 cm", "Velocity"], 
                             ["75 cm", "Height"], 
                             ["125 cm", "Velocity"], 
                             ["125 cm", "Height"]], 
                             columns=["",""])

columns = pd.MultiIndex.from_frame(column_names)
index = ["Cardboard","Carpet","Wood" ]
df = pd.DataFrame(RE_C, columns=columns, index=index)
print(df)
# %%
MEH_H=(RE_CbH125+RE_WH125)/2
MEH_L=(RE_CbH75+RE_WH75)/2
MEH=MEH_H-MEH_L

MES_H=(RE_WH75+RE_WH125)/2
MES_L=(RE_CbH75+RE_CbH125)/2
MES=MES_H-MES_L
IE=((RE_WH125+RE_CbH75)/2)-((RE_CbH125+RE_WH75)/2)

RestitutionEffect = plt.figure(figsize=[7, 7])
plt.subplot(1, 1, 1)
plt.plot((0,1), (MEH_L,MEH_H), 'b-')
plt.plot((0,1), (MES_L,MES_H), 'r-')
plt.xlabel('Level')
plt.ylabel('Effect')
plt.legend(['Height Effect', 'Surface Effect'])
# %%Graphical Representation
P1=R1
tp=P1[0]
ParticleData = plt.figure(figsize=[18, 18])
plt.subplot(1, 3, 1)
plt.plot(P1[1], P1[2], 'b-')
plt.xlabel('Distance(m)')
plt.ylabel('Height(m)')
plt.subplot(1, 3, 2)
plt.scatter(tp[:-1], P1[2][:-1], c=P1[4])
plt.xlabel('Time(s)')
plt.ylabel('Height(m)')
plt.subplot(1, 3, 3)
plt.plot(tp[1:], P1[4], 'c.-')
plt.xlabel('Time(s)')
plt.ylabel('Vertical Velocity(m s^-1)')
end = time.time()
print(f"Execution took: {end - start} s")