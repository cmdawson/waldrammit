import mpl_toolkits.mplot3d as a3
import numpy as np

GREY = [0.4, 0.4, 0.4, 1.0]

class Mesh:
	def __init__(self, vertices, edges, faces, color=None, edgecolor=None):
		"""A generic mesh."""
		self.scene = None

		self.vertices = vertices
		self.edges = edges
		self.faces = faces

		self.midpoint = np.zeros(3)
		for v in self.vertices:
			self.midpoint = self.midpoint + v
		self.midpoint = self.midpoint / len(self.vertices)

		if edgecolor is None:
			self.edgecolor = np.array([0.0,0.0,0.0,1.0])
		else:
			self.edgecolor = edgecolor

		self._normals = {}
		for f in self.faces:
			n = np.cross(vertices[f[1]]-vertices[f[0]], vertices[f[2]]-vertices[f[0]]) 
			self._normals[f] = n / np.sqrt(np.dot(n,n))

		if color is None:
			self.color = np.array([1.0,1.0,1.0,1.0])
		else:
			self.color = color

	def transform(self, R=None, b=None, coords='local'):
		"""Transform in local coordinate system"""
		if R is None:
			R = np.eye(3,3)
		if b is None:
			b = np.zeros(3)

		if coords.lower() == 'global':
			for i in range(len(self.vertices)):
				self.vertices[i] = np.dot(self.vertices[i],R) + b
			self.midpoint = np.dot(self.midpoint,R) + b
		else:
			for i in range(len(self.vertices)):
				self.vertices[i] = np.dot(self.vertices[i]-self.midpoint,R) + self.midpoint + b
			self.midpoint = self.midpoint + b

		for k,v in self._normals.iteritems():
			self._normals[k] = np.dot(v, R)

	def normal(self, face):
		return self._normals[face]
		
	def plot(self, ax, normals=False):
		for face in self.faces:  
			v = np.array([self.vertices[i] for i in face])

			shape = a3.art3d.Poly3DCollection([v])
			shape.set_color(self.color)
			shape.set_edgecolor(self.edgecolor)

			ax.add_collection3d(shape)

			if normals:
				n0 = sum(v) / len(v)
				nv = 1000*self.normal(face)
				ax.plot([n0[0],n0[0]+nv[0]],[n0[1],n0[1]+nv[1]],[n0[2],n0[2]+nv[2]],'black',lw=1)

