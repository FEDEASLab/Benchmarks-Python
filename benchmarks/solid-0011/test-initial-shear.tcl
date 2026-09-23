model BasicBuilder -ndm 3 -ndf 3
node 1 0 0 0 
node 2 1 0 0 
node 3 0 1 0 
node 4 0 0 1 
fixZ 0.0 1 1 1 

nDMaterial ElasticIsotropic 100 1000.0 0.0 
nDMaterial InitStrain 300 100 0.0 0.0 0.1 0.0 0.0 0.0 
element FourNodeTetrahedron 1 1 2 3 4 300 


timeSeries Constant 1 
pattern Plain 1 1 
constraints Transformation 
numberer RCM 
system FullGeneral 
test NormDispIncr 1e-05 100 0 
algorithm Newton 
integrator LoadControl 1.0 
analysis Static 

verify value [analyze 1] 0

verify value [nodeDisp 4 3] -0.1

