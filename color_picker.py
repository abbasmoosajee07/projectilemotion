from cv2 import cv2
import numpy as np
# Pic=cv2.imread('Green_Ball.jpg')
# prevCircle = None
# frame = cv2.resize(Pic,dsize=(600,700))
# Lower = (29, 86, 6); Upper = (64, 255, 255); # Green Boundaries
# # Lower = (150, 10, 10); Upper = (185, 255, 255); # Red Boundaries
# # Lower = (94, 100, 100); Upper = (184, 255, 255); # Blue Boundaries

# blurred = cv2.GaussianBlur(frame, (11, 11), 0)
# hsv = cv2.cvtColor(blurred, cv2.COLOR_BGR2HSV)

# # constructs a mask for the ball color, then perform
# # a series of dilations and erosions to remove irregularities
# maskR = cv2.inRange(hsv, Lower, Upper)
# maskE = cv2.erode(maskR, None, iterations=7)
# maskD = cv2.dilate(maskE, None, iterations=7)

# # find contours in the mask and initialize the circle from center 
# cnts, hierarchy = cv2.findContours(maskD.copy(),
#                     mode=cv2.RETR_TREE, method=cv2.CHAIN_APPROX_NONE)

# cv2.imshow("Circles", maskD)
# key = cv2.waitKey(7000) & 0xFF
# cv2.destroyAllWindows()

# import the necessary packages

# load the image, clone it for output, and then convert it to grayscale
image = cv2.imread('Green_Ball.jpg')
output = image.copy()
gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
# detect circles in the image
circles = cv2.HoughCircles(gray, cv2.HOUGH_GRADIENT, 1.2, 100)# ensure at least some circles were found
if circles is not None:
	# convert the (x, y) coordinates and radius of the circles to integers
	circles = np.round(circles[0, :]).astype("int")
	# loop over the (x, y) coordinates and radius of the circles
	for (x, y, r) in circles:
		# draw the circle in the output image, then draw a rectangle
		# corresponding to the center of the circle
		cv2.circle(output, (x, y), r, (0, 255, 0), 4)
		cv2.rectangle(output, (x - 5, y - 5), (x + 5, y + 5), (0, 128, 255), -1)
	# show the output image
	cv2.imshow("output", image)
	cv2.waitKey(0)