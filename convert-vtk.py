from paraview.simple import *
import sys

r = LegacyVTKReader( FileNames=[sys.argv[1]] )
SaveData(sys.argv[2], proxy=r)