"""
===============================================================================
    # University of Birmingham
    # Labs and Data Analysis 2
    # Portfolio 2:  Particle Tracking
    # Written by:   Abbas Moosajee
    # Written on:   10/04/2022
    # Contact:      AHM080@student.bham.ac.uk
    # Function:     Track a falling ball with OpenCV
===============================================================================
"""

import numpy as np
from cv2 import cv2


def Particle_Tracker(video, track, Particle, itr, height, ymin):
    count = 0
    vs = cv2.VideoCapture(video)

    color = Particle["Color"]
    m = Particle["Mass"]  # Mass in kg
    Lower = color[0]      # Lower Color Boundary
    Upper = color[1]      # Upper Color Boundary
    g = 9.81              # Gravitational Acceleration
    buffer = int(vs.get(cv2.CAP_PROP_FRAME_COUNT))
    x, y, center = (None for _ in range(3))  # Defines position of particle
    pts, yt, xt, ct = ([] for _ in range(4))   # Creates growable lists

    # %% keep looping
    while True:
        oframe = vs.read()  # Grabs current frame
        oframe = oframe[1]
        fps = vs.get(cv2.CAP_PROP_FPS)  # Frame rate of video
        dt = 1/fps                     # Time step
        # Viewing a video and determines,when video reaches end
        if oframe is None:
            break
        count = count+1
        frame = oframe[::-1, :, :]  # Flips image so bottom left is origin
        # Resizes the frame, blurs it, and convert it to HSV color gamut
        frameR = cv2.resize(frame, dsize=(450, 650))
        blurred = cv2.GaussianBlur(frameR, (11, 11), 0)
        hsv = cv2.cvtColor(blurred, cv2.COLOR_BGR2HSV)

        # Constructs a mask for the ball color, then performs a series
        # of dilations and erosions to clean up the mask
        maskR = cv2.inRange(hsv, Lower, Upper)
        maskE = cv2.erode(maskR, None, itr)
        maskD = cv2.dilate(maskE, None, itr)

        # find contours in the mask and initialize the circle from center
        contours, hierarchy = cv2.findContours(maskD.copy(),
                         mode=cv2.RETR_TREE, method=cv2.CHAIN_APPROX_NONE)

        if len(contours) > 0:  # only proceed if at least one contour is found
            # in the maskD, and then use it to create an enclosing circle
            c = max(contours, key=cv2.contourArea)
            ((x, y), radius) = cv2.minEnclosingCircle(c)
            M = cv2.moments(c)
            center = (int(M["m10"] / M["m00"]), int(M["m01"] / M["m00"]))

            if track == 1:  # allows live tracking on video
                if radius > 3:  # only proceed if the radius meets a minimum size
                    # Draws circle on the frame
                    cv2.circle(frameR, (int(x), int(y)),
                               int(radius), (0, 255, 255), 2)
                    cv2.circle(frameR, center, 5, Upper, -1)

                for i in range(1, len(pts)):
                    if pts[i - 1] is None or pts[i] is None:
                        continue
                    # Computes the thickness draws connecting lines
                    thickness = int(np.sqrt(buffer / float(i + 1)) * 2.5)
                    cv2.line(frameR, pts[i - 1], pts[i],
                             (0, 0, 255), thickness)
                # show the frame to screen
                cv2.imshow("Frame", frameR[::-1, :])
                if cv2.waitKey(1) & 0xFF == ord('q'):
                    vs.release()
                    break

        # Adds value from obtained from curent frame to list
        pts.append(center)
        ct.append(count)
        yt.append(y)
        xt.append(x)
    cv2.destroyAllWindows()  # close all windows

    # %%
    # Identifies the positions in which the y position arrays have None values
    inl = np.argwhere(np.array(yt) == None)
    if len(inl) == 0:
        inl = [1]
    yt = np.array(yt[int(inl[-1]+1):])
    xt = np.array(xt[int(inl[-1]+1):])
    ct = np.array(ct[int(inl[-1]+1):])

    # Uses drop height to determine the relative conversion from pixels to m
    r = height/(yt[0]-yt[-1])
    # Calculates angle relative t horizontal using Pythagoras Theorem
    thetad = np.rad2deg(np.arctan(yt[0]/xt[0]))
    # Converts Pixel coordinates to values in m
    ypm = (yt-(yt[-1]))*r
    xpm = xt*r
    # Converts frame counter to time values by multiplying by dt
    tp = (ct-ct[0])*dt
    # Calculates velocity of particle between two frames, by dividing the
    # change in its dispalcement with the dt
    vxp = np.asarray([xpm[ti+1]-xpm[ti] for ti in range(len(xpm)-1)])/dt
    vyp = np.asarray([ypm[ti+1]-ypm[ti] for ti in range(len(ypm)-1)])/dt
    vp = np.sqrt((vxp**2)+(vyp**2))
    # Calculates acceleration by dividing the velocity of particle between two
    # points with the dt
    axp = np.asarray([vxp[ti+1]-vxp[ti] for ti in range(len(vxp)-1)])/dt
    ayp = np.asarray([vyp[ti+1]-vyp[ti] for ti in range(len(vyp)-1)])/dt

    # Calculate Gravitational Potential Energy, Kinetic Energy, and Total Energy
    GPE = m*g*ypm[1:]
    KE = (0.5)*m*(vyp**2)
    TE = GPE+KE

    # %% Restitution Coefficient Calculations
    # Determines the index in y arrays where particle has made contact with 
    # ground and changed direction. Note: ymin must be manually defined
    Bts = np.argwhere(yt < (ymin))

    # As the Velocity before impact is most negative and after most positive,
    # the COR for the velocity method can be found by dividing the max velocity
    # over the min, albeit only for the first bounce
    REv = max(vyp)/min(vyp)
    
    # Uses the indexed position Bts, to create two height arrays one of before
    # bounce and one of after bounce. Use the max heights of each to find COR
    REy = np.sqrt(
        abs(max(ypm[int(Bts[0]):int(Bts[1])])/max(ypm[0:int(Bts[0])])))

    REs = [abs(REv), abs(REy)]
    r_xyi = [xpm[0], ypm[0]]
    Data_Table = [tp, xpm, ypm, vxp, vyp, axp, ayp, GPE, KE, TE]
    Initial_Conditions = {"Coordinates": r_xyi,
                          "Velocity": vp[0],
                          "Launch Angle": thetad}

    print('Coefficient of Restitution based on Velocity:', np.round(REs[0], 2),
          'Coefficient of Restitution based on Height:', np.round(REs[1], 2))
    return(Data_Table, Initial_Conditions, REs) 
