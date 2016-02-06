import numpy as np

class Window:
	def __init__(self, lowleft, upright, wall, frame=40, sill=280):
		self.lowleft = lowleft
		self.upright = upright
		self.wall = wall
		self.frame = frame
		self.sill = 280

#	def partitions(self, size=(32,32)):
		#for i in range(5):
		#	yield (0,0,0)
