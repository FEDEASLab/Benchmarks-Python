pragma openseespy
model  -ndm 3 -ndf 6
section ElasticFrame 1 -E 10000.0 -G 10000.0 -A 1 -Ay 1 -Az 1 -Iy 0.01 -Iz 0.01 -J 0.01
geomTransf Corotational01 1 0 0 1 

node 0 0.0 0 0 
node 1 1.0 0 0 
node 2 2.0 0 0 
node 3 3.0 0 0 
node 4 4.0 0 0 
node 5 5.0 0 0 
node 6 6.0 0 0 
node 7 7.0 0 0 
node 8 8.0 0 0 
node 9 9.0 0 0 
node 10 10.0 0 0 
element ForceFrame 1 0 1 -section 1 -transform 1 -shear 1
element ForceFrame 2 1 2 -section 1 -transform 1 -shear 1
element ForceFrame 3 2 3 -section 1 -transform 1 -shear 1
element ForceFrame 4 3 4 -section 1 -transform 1 -shear 1
element ForceFrame 5 4 5 -section 1 -transform 1 -shear 1
element ForceFrame 6 5 6 -section 1 -transform 1 -shear 1
element ForceFrame 7 6 7 -section 1 -transform 1 -shear 1
element ForceFrame 8 7 8 -section 1 -transform 1 -shear 1
element ForceFrame 9 8 9 -section 1 -transform 1 -shear 1
element ForceFrame 10 9 10 -section 1 -transform 1 -shear 1

fix 0 1 1 1 1 1 1 
fix 10 0 0 0 0 0 0 

pattern Plain 1 Linear 
nodalLoad 10 0 0 25 0 0 314.1592653589793 -pattern 1 

system Umfpack 
integrator LoadControl 0.0025 
test NormUnbalance 1e-10 55 0 
algorithm Newton 
analysis Static 
nodeCoord 10 1 
nodeDisp 10 1 
nodeDisp 10 2 
analyze 1 
getTime  
nodeDisp 10 
