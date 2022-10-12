"""
===============================================================================
    # University of Birmingham
    # Labs and Data Analysis 2
    # Portfolio 2: Particle Tracking
    # Written by:   Abbas Moosajee
    # Written on:   10/04/2022
    # Contact:      AHM080@student.bham.ac.uk
    # Tracks balls in videos, to see factors effecting Restitution Coefficient
===============================================================================
"""
# Setting up the workspace import the necessary packages
import numpy as np
import matplotlib.pyplot as plt
from cv2 import cv2
import seaborn; seaborn.set_style("whitegrid")
cv2.destroyAllWindows()
# %% construct the argument parse and parse the arguments
track = (1)
vs = cv2.VideoCapture("BDropTest_125.mp4.mp4")
g = 9.81       # Gravitational Acceleration
height = 125E-02    # Ball drop height in m
df = 150E-02   # Focal length of camera
Db = 6.86E-02  # Diameter of a tennis ball
m = 56.69E-03  # Mass of Tennis Ball
itr=7
ymin=240
count=0
# define the lower and upper boundaries of the ball in the HSV color space
# Lower = (29, 86,  6); Upper = (64, 255, 255)  # Green Boundaries
# Lower = (150, 100, 10); Upper = (185, 255, 255); # Red Boundaries
Lower = (94,100,100); Upper = (184, 255, 255); # Blue Boundaries
buffer = int(vs.get(cv2.CAP_PROP_FRAME_COUNT))
x,y,center  = (None for _ in range(3))
pts,yt,xt,ct= ([] for _ in range(4))

# %% keep looping
while True:
    oframe = vs.read()  # grab the current frame
    oframe = oframe[1]  # handle the frame from VideoCapture
    fps = vs.get(cv2.CAP_PROP_FPS)
    dt = 1/fps
    # if we are viewing a video and we did not grab a frame,
    # then we have reached the end of the video
    if oframe is None:
        break
    count=count+1
    # resize the frame, blur it, and convert it to the HSV color space
    frame = oframe[::-1, :, :]  # flips image so bottom left is origin
    frameR = cv2.resize(frame,dsize=(500,1000))
    blurred = cv2.GaussianBlur(frameR, (11, 11), 0)
    hsv = cv2.cvtColor(blurred, cv2.COLOR_BGR2HSV)

    # constructs a mask for the ball color, then perform
    # a series of dilations and erosions to remove irregularities
    maskR = cv2.inRange(hsv, Lower, Upper)
    maskE = cv2.erode(maskR, None, itr)
    maskD = cv2.dilate(maskE, None, itr)

    # find contours in the mask and initialize the circle from center
    contours, hierarchy = cv2.findContours(maskD.copy(),
                        mode=cv2.RETR_TREE, method=cv2.CHAIN_APPROX_NONE)

    if len(contours) > 0:  # only proceed if at least one contour was found
        # find the largest contour in the mask, then use it to compute
        # the minimum enclosing circle and centroid
        c = max(contours, key=cv2.contourArea)
        ((x, y), radius) = cv2.minEnclosingCircle(c)
        M = cv2.moments(c)
        center = (int(M["m10"] / M["m00"]), int(M["m01"] / M["m00"]))
        if track == 1:  # allows live tracking on video
            if radius > 3:  # only proceed if the radius meets a minimum size
                # draw the circle and centroid on the frame,
                # then update the list of tracked points
                cv2.circle(frameR, (int(x), int(y)),
                           int(radius), (0, 255, 255), 2)
                cv2.circle(frameR, center, 5, Upper, -1)
                  # update the points queue
            for i in range(1, len(pts)):  # loop of trackd points that creates the trail
                # if either of the tracked points are None, ignore them
                if pts[i - 1] is None or pts[i] is None:
                    continue
                # otherwise, compute the thickness draws connecting lines
                thickness = int(np.sqrt(buffer / float(i + 1)) * 2.5)
                cv2.line(frameR, pts[i - 1], pts[i], (0, 0, 255), thickness)
            # show the frame to screen
            cv2.imshow("Frame", frameR[::-1, :])
            if cv2.waitKey(1) & 0xFF == ord('q'):
                vs.release()
                break
    pts.append(center)
    ct.append(count)
    yt.append(y)
    xt.append(x)
cv2.destroyAllWindows()  # close all windows

# %%
inl=np.argwhere(np.array(yt)==None)
if len(inl)==0: inl=[1]
# inl=np.append(inl,1)
yt=np.array(yt[int(inl[-1]+1):])
xt=np.array(xt[int(inl[-1]+1):])
ct=np.array(ct[int(inl[-1]+1):])
print(yt)

r = height/(yt[0]-yt[-1])
thetad = np.rad2deg(np.arctan(yt[0]/xt[0]))
ypm = (yt-yt[-1])*r
xpm = xt*r
tp  = ct*dt
vxp = np.asarray([xpm[ti+1]-xpm[ti] for ti in range(len(xpm)-1)])/dt
vyp = np.asarray([ypm[ti+1]-ypm[ti] for ti in range(len(ypm)-1)])/dt
axp = np.asarray([vxp[ti+1]-vxp[ti] for ti in range(len(vxp)-1)])/dt
ayp = np.asarray([vyp[ti+1]-vyp[ti] for ti in range(len(vyp)-1)])/dt
vp  = np.sqrt((vxp**2)+(vyp**2))

GPE = m*g*ypm[1:]
KE  = (0.5)*m*(vyp**2)
TE  = GPE+KE

# %%Restitution Coefficient Calculations
Bts = np.argwhere(yt<(ymin))
# REa = (vyp[Bts[0:2]]/vyp[Bts[0:2]-1])
# REv = sum(REa)/len(REa)
REv = max(vyp)/min(vyp)
REy = np.sqrt(abs(max(ypm[int(Bts[0]):int(Bts[1])])/max(ypm[0:int(Bts[0])])))
REs=[abs(REv),abs(REy)]
r_xyi=[xpm[0],ypm[0]]
Data_Table=[tp,xpm,ypm,vxp,vyp,axp,ayp,GPE,KE,TE] # Store Data
Initial_Conditions={"Coordinates":r_xyi,
                    "Velocity":vp[0],
                    "Launch Angle":thetad}
# %%
ParticleData = plt.figure(figsize=[18, 18])
plt.subplot(2, 3, 1)
plt.plot(xpm, ypm,'bo-')
plt.xlabel('Distance(m)')
plt.ylabel('Height(m)')
plt.subplot(2, 3, 2)
plt.plot(tp, xpm, 'b-')
plt.xlabel('Time(s)')
plt.ylabel('Distance(m)')
plt.subplot(2, 3, 3)
plt.scatter(tp[1:], ypm[1:], c=vyp)
plt.xlabel('Time(s)')
plt.ylabel('Height(m)')

plt.subplot(2, 2, 3)
plt.plot(tp[1:], vyp, 'c.')
plt.xlabel('Time(s)')
plt.ylabel('Vertical Velocity(m s^-1)')
plt.subplot(2, 2, 4)
plt.plot(tp[2:], ayp, 'k.')
plt.xlabel('Time(s)')
plt.ylabel('Vertical Acceleration(m s^-2)')
print(REs)

