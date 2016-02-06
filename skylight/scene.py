import numpy as np
import pylab as pl
import mpl_toolkits.mplot3d as a3
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Circle
from shapely.geometry import Point, Polygon
import pickle

from mesh import Mesh

__all__ = ['Scene', 'Mesh']

GREY = [0.5,0.5,0.5,1.0]
LIGHTGREY = [0.6,0.6,0.6,1.0]
BEIGE = [0.955, 0.955, 0.85, 1.0]

EPSILON=1.0E-6

class Scene:
	"""Scene definition for ROL illuminance calculations"""

	def __init__(self):
		"""Initialize a scene for working plane illuminance calculations
		"""
		self.objects = {}
		self.abstractions = []
		self.hidden = set()

	def __getitem__(self, label):
		if label in self.objects:
			return self.objects[label]
		return None

	def add_object(self, obj, label=None):
		"""Add a face (external wall or whatever"""
		obj.scene = self
		self.objects[label] = obj

	def add_abstraction(self, abstract):
		self.abstractions.append(abstract)

	def is_hidden(self, label):
		return label in self.hidden

	def hide(self, label):
		if label in self.objects:
			self.hidden.add(label)

	def unhide(self, label):
		if label in self.hidden:
			self.hidden.remove(label)

	def save(self, filename):
		pickle.dump(self,open(filename,'wb'))

	@staticmethod
	def load(filename):
		return pickle.load(open(filename,'rb'))

	def plot(self, viewpoint=None, normals=False):
		"""Plot the scene using matplotlib"""
		ax = a3.Axes3D(pl.figure())

		# Render 
		for lbl,obj in self.objects.iteritems():
			obj.plot(ax,normals)

		for ab in self.abstractions:
			ab.plot(ax)

		# azim=-90 gives you right behind the origin
		#ax.view_init(elev=7.5, azim=-90)
		#ax.dist=7

		# TODO: calculate minimum we can get to fit everything in and keep aspect
		W = 4500
	
		ax.set_zlim(0,2*W)
		ax.set_xlim(-W,W) 
		ax.set_ylim(-W,W)

		ax.set_xlabel('X')
		ax.set_ylabel('Y')
		ax.set_zlabel('Z')

		#ax.set_axis_off()

		pl.show()
		#pl.savefig('demo.png', dpi=300)

	@staticmethod
	def plot_perspective(polygons):
		ax = pl.subplot('111')
		i = 1

		xmin,xmax = 1E6,-1E6
		ymin,ymax = 1E6,-1E6
    
		for pgon in polygons:
			if pgon[0].geom_type == 'GeometryCollection':
				print pgon[0].geoms
				continue

			verts = [v for v in pgon[0].exterior.coords[:-1]] + [(0,0)]
			codes = [Path.MOVETO]+([Path.LINETO]*(len(verts)-2))+[Path.CLOSEPOLY]

			path = Path(verts, codes)    
			patch = PathPatch(path, fc=pgon[1], ec='black')
			ax.add_patch(patch)    

			bb = pgon[0].bounds
			xmin = bb[0] if bb[0] < xmin else xmin
			xmax = bb[2] if bb[2] > xmax else xmax
			ymin = bb[1] if bb[1] < ymin else ymin
			ymax = bb[3] if bb[3] > ymax else ymax 

			i = i+1

		ax.set_xlim(int(1.1*xmin-1), int(1.1*xmax+1))
		ax.set_ylim(int(1.1*ymin-1), int(1.1*ymax+1))

		pl.show()
		# save figure?

