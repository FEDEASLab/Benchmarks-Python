# [1] M. A. Crisfield, “Accelerating and damping the modified Newton-Raphson method,” 
#     Computers & Structures, vol. 18, no. 3, pp. 395–407, Jan. 1984, 
#     doi: 10.1016/0045-7949(84)90059-2.
#

model basic -ndm 1 -ndf 2

node  1 0

element Rosenbrock 1 1

pattern Plain 1 Constant  {
  load 1  2.421  3
}


set n 1
numberer Plain
system BandGen
algorithm ModifiedNewton
integrator LoadControl [expr 1.0/double($n)]
analysis Static
test Residual 1e-8 3 4
analyze $n
printA
printB

