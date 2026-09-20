# TODO:  SparseSPD still fails this
set E 29e3
set I 1400
set A   20
set L  100

foreach system {Umfpack BandSPD FullGen ProfileSPD SparseGen BandGen} {
  model basic -ndm 2 -ndf 3
  node 1  0 0
  node 2 $L 0
  fix 1 1 1 1
  fix 2 1 1 1
  geomTransf Linear 1 
  element PrismFrame 1 1 2 $A $E $I 1 

  pattern Plain 1 "Constant" {
    sp 2 2 1.0 
  }

  puts "  $system"
  constraints Transformation 
  system $system 
  analysis Static 
  verify value [analyze 1] 0
  reactions

  verify value [nodeReaction 2 2] [expr 12*$E*$I/($L*$L*$L)]
  verify value [nodeReaction 2 3] [expr -6*$E*$I/($L*$L)]

  wipe
}
