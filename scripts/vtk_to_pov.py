import sys

import numpy
import vtk
from vtk.util.numpy_support import vtk_to_numpy

reader = vtk.vtkPolyDataReader()
reader.SetFileName(sys.argv[1])
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()

output=reader.GetOutput()
positions=output.GetPoints().GetData()
n=positions.GetNumberOfTuples()
assert(positions.GetNumberOfComponents()==3)

bead_types=output.GetPointData().GetArray("bead_type")


sys.stdout.write("""
#include "colors.inc"					
#include "stones.inc"					
background {color White}				    
camera{ location < 50, -50, 22 >
        look_at  < 50, 50, 22 >	}

light_source{ < 30, 30, 30 > color White	
            }									

light_source{ < -30, 0, 0 > color White	
            }									

light_source{ < 0, -30, 0 > color White	
            }									

light_source{ < 0, 0, -30 > color White	
            }									

//box{ < 0,0,0 > < 100,100,48 >
//     texture{ pigment{ color rgbf < 0.9,0.9,0.9,0.9 > } } }
""")

colours=sys.argv[2].split(",")

#finish="""finish { ambient .2 diffuse .4 roughness .001 }"""
finish=""

for i in range(0,n):
    pos=positions.GetTuple(i)
    bt=int(bead_types.GetTuple1(i))

    colour=colours[bt]
    if colour!="None":
        sys.stdout.write(f"sphere {{ < {pos[0]}, {pos[1]}, {pos[2]} >, 0.5 texture {{ pigment {{color {colour} }} {finish} }} }}\n")

