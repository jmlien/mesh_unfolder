# mesh_unfolder
Edge unfolding a triangular mesh

The program unfolder allows you to generate a polyhedra net for a given 3D mesh.

by 
Jyh-Ming Lien, jmlien@cs.gmu.edu, George Mason University
Zhonghua Xi, xizhonghua@gmail.com, George Mason University
Hao Yue, yhao3@gmu.edu, George Mason University

Project webpage: http://masc.cs.gmu.edu/wiki/Origami



This archive contains:
- README.txt
- mesh unfolder
- several 3d models in obj or off format

Compile:
This archive uses CMake(https://cmake.org/) compiling software and works across linux/Mac OS X/Windows platform.
linux and Mac OS X users can run gen.sh for a quick setup in build/release or build/debug.
Open GL and glut is required.

Usage:
./unfolder [options] 3d_mesh.obj

Options:
Heuristic Methods
  -h heuristic | use heuristic method
      s        | STEEPEST_EDGE
      f        | FLAT_TREE (default)
     uf        | UNFLAT_TREE
      p        | MINIMUM_PERIMETER
     mp        | MAXIMUM_PERIMETER
      r        | RANDOM
     ga        | Genetic Algorithm
      c        | Clustering based method
     lp        | Linear Programming branch and bound based method

Clustering
  -c filename  | specify the clustering config file. defualt = 'unfolding.cluster'
  -k           | override the number of components.
  -i           | maximum iterations.
  -l filename  | using external labels.
  -w           | using weighted distance between faces.
  -f           | do not use cache file.

Unfolding
  -ga filename | specify the ga config file. default = 'unfolding.ga'
  -lp filename | specify the LP config file. default = 'unfolding.lp'
  -s seed      | specify random seed
  -r times     | retry times, default is 100
  -weights fn  | using the specify weights to unfold the mesh.
  -q           | quite mode.
  -bf          | specify the base face.
  -rb          | random base face.
  -pc          | using pixel checker for overlapping estimation.
  -lc          | less cuts / don't cut flat edges.
  -g           | disable GUI. dump outputs.
  -tab         | add tabs in the net.
  -nbb         | do not find best base face.

Net Optimization
  -objective # | specify the objective for optimization
     hull_area
     cut_length

Dumping SVG
  -scale       | scale factor applied. both *.svg and *.ori will be affected.
  -nl          | do not dump labels in SVG file
  -lfs         | label font size scale [default=1.0]
  -ecx         | extra cut x [default=0.7]
  -ecy         | extra cut x [default=0.1]

Here is a list of keys to control the openGL visualization tool.

Displaying:
  b: show background
  c: toggle random colors
  e: show edge types
  n: show numbers
  o: show overlapping
  w: show wire frame [Model View Mode Only]
  u: toggle showing unfolding
  C: toggle collision detection
  R: show animation in reverse order
  S: only show sharp edges [Model View Mode Only]
  T: show spanning tree
z/Z: rotation model around z-axis

Control:
  r: fully folded
  t: fully unfolded
 sp: toggle animation
Esc: exit

Dumping:
  D: dump unfolding (ori/svg)
  M: dump current (folded) model (obj)
  U: dump flat (unfolded) model (obj)
  W: dump components (wrl)
