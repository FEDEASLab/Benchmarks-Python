set element PrismFrame

model Basic -ndm 3 -ndf 6

node 1 0.0 0.0 0.0 
node 2 0.2 0.0 0.0 
node 3 0.4 0.0 0.0 
node 4 0.6 0.0 0.0 
node 5 0.8 0.0 0.0 
node 6 1.0 0.0 0.0 

fix 1 1 1 1 1 1 1 
fix 6 0 0 1 1 1 0 

section FrameElastic 1 -E 1.0 -G 1.0 -A 2.0 -J 2.0 -Iy 2.0 -Iz 2.0 -Ay 2.0 -Az 2.0 

geomTransf Corotational01 1 0 0 1 

element $element 1  1 2  -section 1 -transform 1 -shear 1
element $element 2  2 3  -section 1 -transform 1 -shear 1
element $element 3  3 4  -section 1 -transform 1 -shear 1
element $element 4  4 5  -section 1 -transform 1 -shear 1
element $element 5  5 6  -section 1 -transform 1 -shear 1


test EnergyIncr 1e-12 10 0
pattern Plain 1 Linear {
  load 6 0 0 0 0 0 [expr -4.0*acos(-1)];
}

integrator LoadControl 0.1
constraints Plain
numberer RCM
system BandGen
algorithm Newton
analysis Static 
verify value [analyze 10] 0
verify value [nodeDisp 6 1] -1 1e-4

