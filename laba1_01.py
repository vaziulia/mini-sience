import gmsh
import math
import os
import sys

gmsh.initialize()

gmsh.model.add("torus")
tor2 = gmsh.model.occ.addTorus(0, 0, 0, 3, 1, 8)
tor21 = gmsh.model.occ.addTorus(0, 0, 0, 3, 0.5, 2)
fluid2 = gmsh.model.occ.cut([(3, 8)], [(3, 2)])

gmsh.model.occ.synchronize()

gmsh.model.mesh.generate(3)

gmsh.write("torus.msh")
gmsh.write("torus.geo_unrolled")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
