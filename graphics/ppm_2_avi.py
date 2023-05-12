
#https://askubuntu.com/questions/971119/convert-a-sequence-of-ppm-images-to-avi-video


import cv2
import glob

img1 = cv2.imread('bat/FalseColor/CS585Bats-FalseColor_frame000000750.ppm')
height, width, layers = img1.shape

video = cv2.VideoWriter('video.avi', -1, 1, (width, height))

filenames = glob.glob('bat/FalseColor/*.ppm')
	for filename in filenames:
		img = cv2.imread(filename)
		video.write(img)

cv2.destroyAllWindows()
video.release()
