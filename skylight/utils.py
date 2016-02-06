import numpy as np
from shapely.geometry import Point, Polygon
from shapely.ops import cascaded_union

EX = np.array([1.0,0.0,0.0])
EZ = np.array([0.0,0.0,1.0])
EPSILON = 1.0E-8
MACHINE_EPS = np.finfo(float).eps

#PSCALE = 10**-int(np.log10(MACHINE_EPS)+3)
PSCALE = 10**6

def scale_polygon(polygon, how='up'):
	return _scale_polygon(polygon, how)

def _scale_polygon(polygon, how='up'):
	if polygon.type == 'MultiPolygon':
		return cascaded_union([_scale_polygon(sp,how) for sp in polygon])
    
	fscale = lambda x: int(0.5 + PSCALE*x) if how=='up' else 1.0*x/PSCALE
	shell_ = [tuple(fscale(e) for e in p) for p in polygon.exterior.coords]
	holes_ = [[tuple(fscale(e) for e in p) for p in interior.coords] \
					for interior in polygon.interiors]

	return Polygon(shell=shell_, holes=holes_)

def polygon_union(polygons):
	"""Takes the union of shapely polygons attempting to avoid floating point issues"""
	if len(polygons) == 1:
		return polygons[0]
		
	U = cascaded_union([_scale_polygon(poly,'up') for poly in polygons])
	return _scale_polygon(U,'down')


def polygon_diff(A, B):
	"""Takes the difference of shapely polygons attempting to avoid floating point issues"""
	D =	_scale_polygon(A,'up') - _scale_polygon(B,'up')
	if D.is_empty:
		return D
	else:
		return _scale_polygon(D,'down')

def polygon_intersect(A, B):
	"""Takes the union of shapely polygons attempting to avoid floating point issues"""
	#print "intersect(", A, ",", B, ")"
	I = _scale_polygon(A,'up').intersection(_scale_polygon(B,'up'))
	if I.is_empty:
		return I
	else:
		return _scale_polygon(I,'down')

def rotation_matrix(axis, theta):
    """Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians."""
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2)
    b, c, d = -axis * np.sin(theta/2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)], [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
		[2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def partition(rect, size=(8,8)):
	"""Generate regular grid on the given rectangle"""
	p0 = rect[0]
	dx = (rect[1] - rect[0]) * 1.0/size[0]
	dy = (rect[2] - rect[1]) * 1.0/size[1]

	for i in range(size[0]):
		for j in range(size[1]):
			pij = p0 + i*dx + j*dy
			yield [pij, pij+dx, pij+dx+dy, pij+dy]

def solid_angle(polygon, origin=np.zeros(3)):
	"""Calculate the solid angle subtended by a convex polygon as viewed from origin.
	No check is made to ensure the polygon actually is convex"""
	omega = 0.0
	verts = [v for v in polygon]
	while len(verts) > 2:
		omega = omega + solid_triangle(verts[:3],origin)
		del verts[1]
	return omega


def solid_triangle(triangle, origin=np.zeros(3)):
	"""Calculate the solid angle subtended by a triangle as seen from origin. Uses the
	tetrahedron solid-angle formula."""
	norm2 = lambda x: np.sqrt(np.dot(x,x))

	va = triangle[0] - origin
	vb = triangle[1] - origin
	vc = triangle[2] - origin

	a,b,c = norm2(va), norm2(vb), norm2(vc)

	alpha = np.dot(np.cross(va,vb),vc) 
	beta = a*b*c + np.dot(vc,va)*b + np.dot(va,vb)*c + np.dot(vb,vc)*a

	return 2.0*np.arctan2(np.abs(alpha), beta) 


def moller_trombore(ray, triangle):
	"""Moller-Trombore algorithm to determine if a ray intersects a triangle.
	
	Args:
		ray (2x3 array) : ray described by origin and (unit) direction
		triangle (3x3 array) : Triangle described by three vertices
	"""
	e1 = triangle[1] - triangle[0]
	e2 = triangle[2] - triangle[0]

	# direction of ray
	rd = ray[1]  - ray[0]
	p = np.cross(rd, e2)
	det  = np.dot(e1, p)

	# if determinant is near zero, ray lies in plane of triangle
	if  np.abs(det) < EPSILON:
		return False

	idet = 1.0 / det

	# calculate distance from vertex 0 to ray origin
	t = ray[0] - triangle[0]

	# calculate u parameter and test bound
	u = np.dot(t,p) * idet

	# we're definately outside
	if  u < 0.0 or u > 1.0:
		return False

	# otherwise test the v parameter
	q = np.cross(t,e1)

	# and calculate v parameter
	v  = np.dot(rd,q) * idet

	if v  < 0.0 or (u+v)>1.0:
		return False

	t = np.dot(e2,q) * idet

	if  t > EPSILON:
		return True

	return False


def intersects_any(vector, faces):
	"""Returns True if the (directed) vector intersects any of the given faces at any point.
	"""
	for face in faces:
		if len(face) == 3:
			return moller_trombore(vector, face)
		elif len(face) == 4:
			return moller_trombore(vector, face[:3]) or moller_trombore(vector, (face[0],face[2],face[3]))
		else:
			raise NotImplementedError("Face must be of degree < 5")

def clip_face(clipplane, face):
	"""Clip a face according to the given clipping plane
	Args:
		clipplane (point,normal): Clipping plane defined by a face and normal
		face: array of vertices
	"""
	c0 = clipplane[0]
	cN = clipplane[1]
	clip = [np.dot(v-c0,cN) <= 0.0 for v in face]
	result = []

	if np.all(clip):
		return result
	elif not np.any(clip):
		return face

	n = len(face)

	for j in range(n):
		if clip[j]:
			if not clip[(j+1)%n]:
				result.append(line_plane_intersection((face[(j+1)%n],face[j]),c0,cN)) 
		else:	
			result.append(face[j])
			if clip[(j+1)%n]:
				result.append(line_plane_intersection((face[j],face[(j+1)%n]),c0,cN)) 

	return result

def project2d(viewpoint, mesh, origin, axes, union=True, hidden=False, clip=True):
	"""Project visible faces of given mesh as viewed from a given point onto a plane and return
	the result as polygon(s) in the plane's 2D coordinate system""" 
	p0 = origin
	pN = axes[2]

	projection = None
	polygons = list()

	# define the clipping plane
	cN = origin - viewpoint
	cN = cN / np.sqrt(np.dot(cN,cN))
	c0 = viewpoint + 0.02*cN

	for face in mesh.faces:
		fv = np.array([mesh.vertices[i] for i in face])
		fmid = sum(fv) / len(fv)
		fn = mesh.normal(face)

		if np.dot(fn,fmid-viewpoint) >= 0.0 and not hidden:
			continue # (face is not visible)

		cf = clip_face((c0,cN),fv)
		if len(cf) == 0:
			continue

		# projected
		proj = [line_plane_intersection((viewpoint,v), p0, pN)-p0 for v in cf]
		polygons.append(Polygon([(np.dot(v,axes[0]), np.dot(v,axes[1])) for v in proj]))

	if len(polygons) > 0:
		if union:
			return polygon_union(polygons)
		else:
			return polygons

	return projection


def project3d(vertices, origin, axes):
	"""Given vertices in a 2d coordinate system defined by (origin,axes) transform them back
	into 3d coordinates"""
	return [origin + v[1]*axes[0] + v[0]*axes[1] for v in vertices]


def line_plane_intersection(line, point=None, normal=None, triangle=None):
	"""Return the point of intersection between a plane defined by either point within the plane
	and a normal or alternatively as triangle within the plane."""
	if point is None:
		point = triangle[0]
		normal = np.cross(triangle[1]-p0, triangle[2]-p0)

	l = line[1] - line[0]
	nl = np.dot(normal,l)
	if abs(nl) < EPSILON:
		return None

	return line[0] + (np.dot(normal,point-line[0])/nl)*l 

def projection_plane(point, origin):
	ez = point - origin
	ez = ez / np.sqrt(np.dot(ez,ez))
	ey = line_plane_intersection((origin+EZ,point), origin, ez) - origin
	ey = ey / np.sqrt(np.dot(ey,ey))
	ex = np.cross(ez,ey)
	return (-ex,ey,ez)


def intersects_triangle(ray, triangle):
	"""Moller-Trombore algorithm to determine if a ray intersects a triangle
	
	Args:
	 ray (2x3 array) : ray described by origin and direction
	 triangle (3x3 array) : Triangle described by three vertices
	"""
	e1 = triangle[1]  - triangle[0]
	e2 = triangle[2]  - triangle[0]

	# direction of ray
	rd = ray[1]  - ray[0]
	rd = rd / np.linalg.norm(rd)
	p = np.cross(rd, e2)

	det  = np.dot(e1, p)

	# if determinant is near zero, ray lies in plane of triangle
	if  np.abs(det) < EPSILON:
		return False

	idet = 1.0 / det

	# calculate distance from vertex 0 to ray origin
	t = ray[0] - triangle[0]

	# calculate u parameter and test bound
	u = np.dot(t,p) * idet

	# we're definately outside
	if  u < 0.0 or u > 1.0:
		return False

	# otherwise test the v parameter
	q = np.cross(t,e1)

	# and calculate v parameter
	v  = np.dot(rd,q) * idet

	if v < 0.0 or (u+v)>1.0:
		return False

	t = np.dot(e2,q) * idet

	if  t > EPSILON:
		return True

