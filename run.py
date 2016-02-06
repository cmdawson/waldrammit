#!/usr/bin/python
import pandas as pd, numpy as np, sys
from skylight import *

# room
floor2 = (2800,2200)
window21 = Window(lowleft=(1450,830), upright=(2180,2140), wall=2, frame=40, sill=280)
#window22 = Window(lowleft=(200,1550), upright=(550,2140), wall=2, frame=40, sill=280)
bathroom = Room(floor2, window21)
bathroom.transform(b=(200,-3240,2600))

sr29 = build_rectangle((0,6630), (0,4900), z=1210, zdir='y')

or29 = build_rectangle((0,3320), (4900,6300), z=1810, zdir='y')
or29_2 = Mesh( \
		vertices = [ \
			np.array((0,0,0)), np.array((3320,0,0)), \
			np.array((3320,950,950)), np.array((0,950,950)) \
		], \
		edges = [(0,1),(1,2),(2,3),(3,0)], \
		faces = [(0,1,2,3)])
or29_2.transform(b=np.array([0,1810,6300]))

or27 = build_rectangle((0,3320), (4900,7300), z=5620, zdir='y')

ch29 = build_rectangle((4550,5400), (4900,6000), z=1210, zdir='y')
ch29_2 = build_rectangle((4900,6000), (1210,1710), z=4550, zdir='x')

rear_wall = build_rectangle((0,8000), (8000,-1500), z=0, zdir='x')

bk_road = build_rectangle((0,8000), (50000,-10000), z=37000, zdir='x')

S = Scene()

S.add_object(bathroom,'bathroom')

existing = False
resolution = 35

# External environment
S.add_object(sr29,'sr_29')

if existing:
	S.add_object(ch29,'ch_29')
	S.add_object(ch29_2,'ch_29_2')
	S.add_object(or27,'or_27')
else:
	S.add_object(or29,'or_29')
	S.add_object(or29_2,'or_29_2')

S.add_object(rear_wall,'rear_wall')
S.add_object(bk_road,'bk_road')

lightmapper = LightMapper(bathroom, height=850)

S.plot(normals=True)

#mapf = pd.DataFrame(lightmapper.map(resolution), columns=['A','B','C','Omega'])
#mapf.to_csv('bathroom_%d_%s.csv' % (resolution,'e' if existing else 'p'))


