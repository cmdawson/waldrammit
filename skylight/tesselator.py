from OpenGL.GLU import *
from OpenGL.GL import GL_TRIANGLE_FAN, GL_TRIANGLE_STRIP, GL_TRIANGLES
import numpy as np

__all__ = ['Tesselator']

class Tesselator:
	"""Tesselate Shapely polygon using OpenGL tesselation utility."""
	def __init__(self):
		self.winding = GLU_TESS_WINDING_ODD
		self.geometry = None
		self.vertices = list()

		self.controller = gluNewTess()

		gluTessCallback(self.controller, GLU_TESS_COMBINE, self._glu_combine)
		gluTessCallback(self.controller, GLU_TESS_VERTEX, self._glu_vertex)
		gluTessCallback(self.controller, GLU_TESS_BEGIN, self._glu_begin)
		gluTessCallback(self.controller, GLU_TESS_END, self._glu_end)
		gluTessProperty(self.controller, GLU_TESS_WINDING_RULE, self.winding)
				
	def clear(self):
		self.triangles = list()  


	def tesselate(self, polygon, clear=True):
		"""Generate list of triangles forming a tesselation of a polygon"""
		if clear:
			self.clear()

		if polygon.geom_type == 'MultiPolygon':
			for subgon in polygon:
				self.tesselate(subgon, clear=False)
			return self.triangles
		elif polygon.geom_type != 'Polygon':
			raise NotImplementedError('Unknown Shapely geometry type: %s' % polygon.geom_type)

		gluTessBeginPolygon(self.controller, self)

		gluTessBeginContour(self.controller)
		for coord in polygon.exterior.coords[:-1]:
			vertex = coord+(0.,)
			gluTessVertex(self.controller, vertex, vertex) 
		gluTessEndContour(self.controller)

		for interior in polygon.interiors:
			if interior is None:
				continue
			
			gluTessBeginContour(self.controller)
			icoords = [v for v in interior.coords[:-1]]
			icoords.reverse()
			for coord in icoords:
				vertex = coord+(0.,)
				gluTessVertex(self.controller, vertex, vertex) 
			gluTessEndContour(self.controller)

		gluTessEndPolygon(self.controller)

		return self.triangles

	def _glu_begin(self, geometry, data=None ):
		self.geometry =  geometry
		self.vertices = list()
    
	def _glu_vertex(self, vertex, data=None):
		self.vertices.append(vertex)
    
	def _glu_combine( self, newpos, vertices, weights, data=None):
		newvertex = np.array((0.0,0.0,0.0))
		for vertex,weight in zip(vertices,weights):
			if not vertex is None:
				newvertex = newvertex + weight * np.array(vertex)

		return tuple(newvertex)

	def _glu_end(self, *args, **namedargs):
		verts = self.vertices
        
		if self.geometry == GL_TRIANGLE_FAN:
			v1 = verts[1]
			for v in verts[2:]:
				self.triangles.append((verts[0],v1,v))
				v1 = v            
		elif self.geometry == GL_TRIANGLE_STRIP:
			for ii in range(len(verts)-2):
				if ii % 2: 
					self.triangles.append((verts[ii], verts[ii+1], verts[ii+2]))
				else:
					self.triangles.append((verts[ii+1], verts[ii], verts[ii+2]))    
		elif self.geometry == GL_TRIANGLES:
			if len(verts) % 3 != 0:
				raise NotImplementedError(r'Expected 3xn vertices for n GL_TRIANGLES')
			ii = 0
			while ii < len(verts):
				self.triangles.append(tuple(verts[ii:ii+3]))
				ii += 3
		else:
			raise NotImplementedError(r'Unknown geometry %d' % self.geometry)
