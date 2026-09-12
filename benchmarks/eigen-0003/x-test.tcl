
set A   216
set Iz 2592
set Iy 5832
set J  6085.12013626099
set E  3600
set G  1500
set L  96.0

model Basic -ndm 3 -ndf 6

node 1  0 0 0 ; # Joint=1
node 2 $L 0 0 ; # Joint=2

fix 1  1 1 1 1 1 1 
fix 2  1 0 0 1 0 0; # -dof 1, 4


section FrameElastic 1 -A $A -Ay $A -Az $A -Iz $Iz -Iy $Iy -J $J -E $E -G $G

# Elements
geomTransf Linear 1 0.0 0.0 1.0 
element ForceFrame 1 1 2 -section 1 -n 8 -gauss_type Legendre -transform 1 -mass 4.968e-05 -shear 0; # Frame 1


constraints Transformation 
system BandGeneral 
numberer RCM 
number


set ev [eigen 2 -solver -genBandArpack]

verify error [lindex $ev 0] 13268.512 1e-6
verify error [lindex $ev 1] 29854.153 1e-6

