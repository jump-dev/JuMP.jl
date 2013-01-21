for s in [100,1000,5000,10000]
  println("100 100 $s")
  system("julia pmedian.jl 100 100 $s")
  system("julia pmedian.jl 100 100 $s")
  system("julia pmedian.jl 100 100 $s")
end
