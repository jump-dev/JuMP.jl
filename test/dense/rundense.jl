for s in {[10000,1000],[1000,1000],[1000,100],[100,10]}
  println("$(s[1]) $(s[2])")
  system("julia dense.jl $(s[1]) $(s[2])")
  system("julia dense.jl $(s[1]) $(s[2])")
  system("julia dense.jl $(s[1]) $(s[2])")
end


