import numpy as np
from shapely.geometry import Point, Polygon
from utils import *
from ray import Ray
from collections import namedtuple
from tesselator import Tesselator

class LightMapper:
	def __init__(self, room, height=850):
		self.room = room
		self.window = room.window
		self.scene = room.scene
		self.height = height

		self.environment = []
		for obj in self.scene.objects.values():
			if obj == self.room:
				continue
	
			for face in obj.faces:
				self.environment.append([obj.vertices[i] for i in face])

		for sill in self.room.sills:
			for face in sill.faces:
				self.environment.append([sill.vertices[i] for i in face])

	def _calcOmega(self, point, windlet):
		"""Calculate the solid angle subtended by the sky visible through a windlet
		"""
		corners = [intersects_any((point,c-point),self.environment) for c in windlet]

		if np.all(corners):
			return 0.0
		elif not np.any(corners):
			return solid_angle(windlet, point)

		# Otherwise sample a 16x16 grid
		N, m = 16, 0
		dx = (windlet[1] - windlet[0]) * 1.0/N
		dy = (windlet[2] - windlet[1]) * 1.0/N 
		for i in range(N):
			for j in range(N):
				wv = windlet[0] + (0.5+i)*dx + (0.5*j)*dy
				if not intersects_any((point,wv-point),self.environment):
					m = m + 1

		return (solid_angle(windlet,point) * m) / N
		

	def luminance(self, point):
		"""Calculate the luminance received from the visible sky at point"""
		window = self.window
		room = self.room
		w0 = window.midpoint

		window = [self.window.vertices[j] for j in self.window.faces[0]]
		omega = 0.0

		for wp in partition(window, size=(128,48)):
			omega += self._calcOmega(point, wp)

		return omega
		
	def map(self, resolution, xlim=None, ylim=None):
		"""Build a lightmap of the room at the given resolution"""
		lightmap = []
		for point in self.grid(resolution):
			try:
				lightmap.append((point[0],point[1],point[2],self.skylight(point)))
			except Exception as e:
				raise RuntimeError('Exception %s at point (%d,%d,%d)' \
					% (str(e),point[0],point[1],point[2]))

		return lightmap


	def grid(self, resolution):
		"""Generate regular grid of points on the working plane"""
		vertices = self.room.floor.vertices

		zdir = np.cross(vertices[1]-vertices[0], vertices[2]-vertices[0])
		zdir = zdir / np.sqrt(np.dot(zdir,zdir))

		p0 = vertices[0] + self.height*zdir

		vx = vertices[1] - vertices[0]
		vy = vertices[2] - vertices[1]

		nx = int(np.sqrt(np.dot(vx,vx))/resolution)
		ny = int(np.sqrt(np.dot(vy,vy))/resolution)

		dx = vx * 1.0 / nx
		dy = vy * 1.0 / ny

		for i in range(nx):
			for j in range(ny):
				yield p0 + (0.5+i)*dx + (0.5+j)*dy


	def perspective(self, point, pplane=None):
		"""Perspective from the given point"""
		window = self.window
		room = self.room

		w0 = window.midpoint
		if pplane is None:
			pplane = projection_plane(point, w0)

		#ray = Ray(point, w0 + 4*(w0-point))
		#self.scene.add_abstraction(ray)

		polygons = []	

		wgon = project2d(point, window, w0, pplane, union=False)

		# TODO: sort by distance from viewpoint?
		for label,mesh in self.scene.objects.iteritems():
			if self.scene.is_hidden(label) or mesh == room:
				continue

			mgons = project2d(point, mesh, w0, pplane, union=False)

			if mgons is None:
				#print label, "(None)"
				continue
			elif len(mesh.faces) == 1:
				#print label, "(single face)"
				#print wgon, mgons
				vgon = polygon_intersect(wgon[0],mgons[0])
				polygons.append((vgon,[0.75,0.1,0.1,1]))

		for sill in reversed(room.sills):
			sp = project2d(point, sill, w0, pplane, hidden=True)
			if not sp is None:
				polygons.append((sp,[1,1,1,1]))

		polygons.append((wgon[0],[0,0,1,0.1]))

		return polygons

	def skylight(self, point):
		"""Calculate the solid angle subtended by the sky visible from a given point
		Almost identical to 
		"""
		window = self.window
		room = self.room

		w0 = window.midpoint
		pplane = projection_plane(point, w0)

		polygons = self.perspective(point,pplane=pplane)
		wgon = polygons[-1][0]
		blobs = []

		for p in polygons[:-1]:
			if not wgon.intersects(p[0]):
				continue

			#foo = polygon_intersect(wgon,p[0])
			#print foo.area
			#print foo

			blobs.append(polygon_intersect(wgon,p[0]))
			
		bigblob = polygon_union(blobs)

		#return [(bigblob, [0,0,1,1])]

		T = Tesselator()

		triangles = T.tesselate(bigblob, clear=True)

		omega = solid_angle(window.vertices, point)

		#print "Omega_0", omega
		for tri in triangles:
			tri3d = project3d(tri, w0, pplane)
			omega -= solid_angle(tri3d, point)

		return 0.5*omega/np.pi

		#print bigblob.area
		#print wgon.area

		# Only take the (x,y) as returned by tesselate!
		#return [(Polygon([v[:2] for v in t]), [0,0,1,1]) for t in triangles]

		#print "NumBlobs:", len(blobs)
		#return [(blob,[0,0,1,1]) for blob in blobs]

	def test(self):
		pass
