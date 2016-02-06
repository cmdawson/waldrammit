import numpy as np
from utils import moller_trombore

EPSILON = 1.0E-8

class Ray:
	def __init__(self, start, end):
		self.color = [1,0,0,1]
		self.start = np.array(start)
		self.end = np.array(end)

		self.direction = self.end - self.start
		self.direction = self.direction / np.sqrt(np.dot(self.direction,self.direction))

	def plot(self, ax):
		ax.plot(*zip(self.start,self.end), c=self.color, lw=0.5)
