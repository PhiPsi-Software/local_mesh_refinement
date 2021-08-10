## 1 About this program

local_mesh_refinement is a program written in Fortran to perform local refinement of XFEM-enriched 3D hexahedral elements. Given the original mesh (*.elem and *.node files) and a file (*.ennd) that describes enriched nodes, this program generates the refined mesh (output.elem and output.ndoe).

The main program is local_mesh_refinement.f90.

Two Matlab scripts, matlab_plot_original_mesh.m and matlab_plot_refined_mesh.m, are provided to plot the original mesh and refined mesh, respectively.

## 2 Input files description

### 2.1 *.node file

This file contains the coordinate system of each node, example lines are given as follows:

        0.00000000E+00   0.50000000E+00   0.00000000E+00
        0.00000000E+00   0.00000000E+00   0.00000000E+00
        0.00000000E+00   0.45454545E+00   0.00000000E+00
        0.00000000E+00   0.40909091E+00   0.00000000E+00
        0.00000000E+00   0.36363636E+00   0.00000000E+00

### 2.2 *.elem file

This file contains nodes number (columns 1 to 8) and material type number (column 9) of each element, example lines are given as follows:

      2      33      41      12     121     123     161     160       1
     33      34      42      41     123     124     171     161       1
     34      35      43      42     124     125     181     171       1
     35      36      44      43     125     126     191     181       1
     36      37      45      44     126     127     201     191       1

### 2.3 *.ennd file

This file contains enriched nodes. If a node is an enriched node, then its corresponding line equals "1". Example lines are given as follows:

        0
        0
        1
        1
        0

## 3 Please cite the following papers

+ Shi F., Liu J. A fully coupled hydromechanical XFEM model for the simulation of 3D non-planar fluid-driven fracture propagation. Computers and Geotechnics, 2021, 132, 103971.

+ Shi F. XFEM-based numerical modeling of well performance considering proppant transport, embedment, crushing and rock creep in shale gas reservoirs. Journal of Petroleum Science and Engineering, 2021, 201, 108523.

+ Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An XFEM-based numerical model to calculate conductivity of propped fracture considering proppant transport, embedment and crushing. Journal of Petroleum Science and Engineering, 2018, 167: 615-626.

+ Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An XFEM-based method with reduction technique for modeling hydraulic fracture propagation in formations containing frictional natural fractures. Engineering Fracture Mechanics, 2017, 173: 64-90.

+ Shi F., Wang X.L., Liu C., Liu H., Wu H.A. A coupled extended finite element approach for modeling hydraulic fracturing in consideration of proppant. Journal of Natural Gas Science and Engineering, 2016, 33: 885-897.