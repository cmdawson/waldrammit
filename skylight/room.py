from utils import *
from primitives import *
from collections import namedtuple

class Room:
	"""Room analyzed for right-of-light calculations"""

	def __init__(self, floor, window):
		"""Rectangular floor is passed as (width,length). Window given by two corners
		in 2D and a wall index where:
			0: wall with y = 0
			1: wall with x = width
			2: wall with y = length
			3: wall with x = 0
		"""
		self.scene = None

		width,length = floor[0], floor[1]
		self.floor = build_rectangle(xlim=(0,width), ylim=(0,length))
		self.floor.color = LIGHTGREY

		ll,ur = window.lowleft, window.upright
		fw = window.frame

		self.window = build_rectangle(xlim=(ll[0]+fw,ur[0]-fw), ylim=(ll[1]+fw,ur[1]-fw), zdir='Y')
		self.window.color = LIGHTBLUE[0:3] + [0.15]
		self.window.edgecolor = LIGHTBLUE

		self.sills = []
		silldepth = window.sill

		if silldepth > 1.0E-6:
			self.sills.append(build_rectangle(xlim=(ll[0],ur[0]), ylim=(0,silldepth)))
			self.sills.append(build_rectangle(xlim=(ll[0],ur[0]), ylim=(0,silldepth)))	
			ya,yb = -0.5*(ur[1]-ll[1]), 0.5*(ur[1]-ll[1])
			self.sills.append(build_rectangle(xlim=(ya,yb), ylim=(0,silldepth)))
			self.sills.append(build_rectangle(xlim=(ya,yb), ylim=(0,silldepth)))

			self.sills[0].transform(b=ll[1]*EZ)
			self.sills[1].transform(rotation_matrix(EY,np.pi), ur[1]*EZ)
			self.sills[2].transform(rotation_matrix(EY,np.pi/2), 0.5*(ur[1]+ll[1])*EZ + ur[0]*EX)
			self.sills[3].transform(rotation_matrix(EY,-np.pi/2), 0.5*(ur[1]+ll[1])*EZ + ll[0]*EX) 

		# Position the window on the correct wall
		A = np.eye(3,3)
		b = (0,0,0)

		if window.wall == 0:
			A = rotation_matrix(EZ, np.pi)
			b = (width,0,0)
		elif window.wall == 1:
			A = rotation_matrix(EZ, 0.5*np.pi)
			b = (width,length,0)
		elif window.wall == 2:
			b = (0,length,0)
		elif window.wall == 3:
			A = rotation_matrix(EZ, -0.5*np.pi)
		else:
			raise RuntimeError("Unknown wall index")

		self.window.transform(A,b,'global')
		for sill in self.sills:
			sill.transform(A,b,'global')
	
	def grid(size=(8,8), z=850):
		
		p0 = rect[0]
		dx = (rect[1] - rect[0]) * 1.0/size[0]
		dy = (rect[2] - rect[1]) * 1.0/size[1]

		for i in range(size[0]):
			for j in range(size[1]):
				pij = p0 + i*dx + j*dy
				yield [pij, pij+dx, pij+dx+dy, pij+dy]

	
	def transform(self, R=None, b=None):
		"""Transform in local coordinate system"""
		if R is None:
			R = np.eye(3,3)
		if b is None:
			b = np.zeros(3)

		self.floor.transform(R,b)
		self.window.transform(R,b)
		for sill in self.sills:
			sill.transform(R,b)


	def plot(self, ax, normals=False):
		self.floor.plot(ax, normals)
		self.window.plot(ax, normals)
		for sill in self.sills:
			sill.plot(ax, normals)


	# Interior walls if rendering in blender or something.
	# opaque off-white patch for the interior wall
	# (bit pointless with matplotlib)
	#fw = self.framewidth
	#maxz = 1.2*self.window.vertices[2][2]
	#p0,p1 = self.plane[0], self.plane[1]

	#wpath = Path( \
		#[(p0[0],0), (p1[0],0), (p1[0],maxz), (p0[0],maxz), (np.nan,np.nan)] \
		 #+ [(v[0],v[2]) for v in self.window.vertices[::-1]] \
		 #+ [(np.nan,np.nan)], \
		#[\
			#Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY, \
			#Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY  \
		#])
	#wpatch = PathPatch(wpath, fc=BEIGE, ec=GREY) #, hatch='/')
	#ax.add_patch(wpatch)
	#art3d.pathpatch_2d_to_3d(wpatch, z=self.plane[2][1], zdir='y')


	#def _sanity_check(self):
		# calculate the solid angle subtended by the window (i) directly
		# and (ii) summing up the solid angles of the individual positions
		#spoints = [p for p in self.ppoints()]
		#pidx = 26

		#for pt in spoints[26:55]:
			#momega = solid_angle(self.window.vertices, pt)
			#somega = 0.0

			#for wp in self.window_partitions:
				#rect = [self.window_points[i] for i in wp]
				#somega = somega + solid_angle(rect, pt)

			#print pidx, momega, r'(%.2f%%)' % (50*momega/np.pi), somega
