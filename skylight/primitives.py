from mesh import Mesh
import numpy as np
from utils import rotation_matrix

LIGHTBLUE = [0.1, 0.1, 0.4, 1.0]
LIGHTGREY = [0.6,0.6,0.6,1.0]

EX = np.array([1.0,0.0,0.0])
EY = np.array([0.0,1.0,0.0])
EZ = np.array([0.0,0.0,1.0])

def build_rectangle(xlim, ylim, z=0, zdir='z'):
	"""Build a simple face in the XY plane given the"""
	vertices = np.zeros((4,3))
	faces = [(0,1,2,3)]
	edges = [(i,(i+1)%4) for i in range(4)]

	xi,yi,zi = (0,1,2)
	if zdir.upper() == 'Y':
		zi,yi = (1,2)
	elif zdir.upper() == 'X':
		zi,xi = (0,2)

	for i in range(4):
		a,b = i%2, i/2
		vertices[i][xi] = xlim[(a+b)%2]
		vertices[i][yi] = ylim[b]
		vertices[i][zi] = z

	return Mesh(vertices, edges, faces)

def build_box(width, length, height):
	"""Build a box mesh"""
	vertices = np.zeros((8,3))
	edges =[]
	faces =[]

	for i in range(4):
		a,b = i%2, i/2
		vertices[i][0] = vertices[4+i][0] = width * ((a+b)%2)
		vertices[i][2] = vertices[4+i][2] = height * b
		vertices[i][1] = 0.0
		vertices[4+i][1] = length

		edges.append((i,(i+1)%4))
		edges.append((4+i,(4+i+1)%8))
		edges.append((i,i+4))

		faces.append((i,i+4,4+((i+1)%4),(i+1)%4))

	faces.append((0,1,2,3))
	faces.append((4,7,6,5))

	return Mesh(vertices, edges, faces)

def extrude_polygon(polygon, length, interior=None):
	"""Take a 2d polygon in the YZ plane and extrude it a given length"""
	plen = len(polygon)

	vertices = np.zeros((2*plen,3))
	edges = []
	faces = []

	for i in range(plen):
		vertices[i][0] = 0.0
		vertices[i][1] = vertices[i+plen][1] = polygon[i][0]
		vertices[i][2] = vertices[i+plen][2] = polygon[i][1]
		vertices[i+plen][0] = length

		edges.append((i,(i+1)%plen))
		edges.append((plen+i,(plen+i+1)%(2*plen)))
		edges.append((i,i+plen))

		faces.append((i,i+plen,plen+((i+1)%plen),(i+1)%plen))

	faces.append(tuple(range(plen)))
	faces.append(tuple(range(plen,2*plen)))

	return Mesh(vertices, edges, faces)




