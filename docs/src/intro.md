# JUDYN

2021/12/14

Computational efficiency of the JUDYN code  has been drastically improved through different
 major changes.

1) The user-defined types Vec3, RV3 and Mat3 of the MyAlgebra.jl  library are directly defined
in terms of Vector{Float64} and Array{Float64,2} types instead of using the 
StaticArrays package. 
As a consequence, It is no longer recommended to define Vec3 and RV3 quantities 
in the forms Vec3([...]) and RV3([...]) since there is no control of the vector size. 
They should be input in the form Vec3(a,b,c) and RV3(x,y,z).
2) The computation of the rotation and tangent operators and associated operations has been revised.
3) The exact type of function arguments (e.g. Vector{Float64}, Vector{Int}, Array{Float64,2}, etc.) 
is always specified  in the defintition of the function calling sequence.
4) The "save_freq" parameter  has been introduced to define the stepping interval at which the solution 
is archived in the .h5 output file.
5) Archiving on the .MAT file has been removed. 
6) an "init_judyn_postpro.jl" file has been  created to define  the post-processing environment.

Two examples have been updated so far  to take the modifications introduced : "Top" and "Haug".

The "Haug" beam with 20 nonlinear elements has been used as a benchmark to measure the improvement
 in computational efficiency. On this specific example, the elapsed time for the solution phase solve()
 is reduced from 144 to 32 s. (Julia 1.6.4 Linux version under WSL2 on Windows, Open BLAS library, 
Intel(R) Core(TM) i7 CPU  960  @ 3.20GHz).

In the postprocessing files, a Bool::save_fig parameter has been introduced to control the graphic output 
to .pdf files.

2021/11/08 

upate 01/08/2021

This new version of JUDYN  contains many new features:

introduction of natural strain measure in superelement.

nonlinear beam element

static solution

and may other corrections

update 27/09

This new version of JUDYN has been strongly re-organized.

- For better readability, all the modules (MyAlgebra, Nodes, SetElements, BuildElements, SuperElements, Solve, Frames, InputFunctions)  have been split into small files.

- A "super-beam" nonlinear element constructed  from a highly-flexible linear kernel has been introduced.

- The set of benchmark problems is being enlarged, modified and improved. In particular, on 27/09,  the "NL_beam" and "super-beam" solutions of the "Yigit" problem have been introduced  and a post-processing file "postpro_Yigit_compare.jl"  has been written to compare the solutions.

- New functionalities have been introduced (e.g. saving of the problem topology allowing delayed post-processing)

- New examples and modifications to the other existing ones will be introduced soon.
