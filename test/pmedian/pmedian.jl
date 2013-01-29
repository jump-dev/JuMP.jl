# pmedian.jl
# Solves the P-Median problem

using MathProg

function doTest(numFacility::Int,numCustomer::Int,numLocation::Int)
	
	
	
	#numFacility = 100
	#numCustomer = 100
	#numLocation = 1000
    srand(10)	
	#customerLocations = [randi(numLocation) for a = 1:numCustomer ]
	customerLocations = [rand(1:numLocation) for a = 1:numCustomer ]
	
	tic()
	m = Model("max")
	
	# Facility locations
	#s = [ Variable(m,0,1,0,"s$i") for i=1:numLocation ]
	s = [ Variable(m,0,1,0,@sprintf("s%d",i)) for i=1:numLocation ]
	#s = addVars(m, 0, 1, JUMP_CONTINUOUS, numLocation, "s")

	# Aux. variable: x_a,i = 1 iff the closest facility to a is at i
	#x = [ Variable(m,0,1,0,"x$a,$i") for i = 1:numLocation, a = 1:numCustomer]
	#x = [ Variable(m,0,1,0,@sprintf("x%d,%d",a,i)) for i = 1:numLocation, a = 1:numCustomer]
	x = addVars(m,0,1,JUMP_CONTINUOUS, (numLocation, numCustomer), "x")
	
	# Objective: min distance
	#vars::Array{Variable,1} = reshape([x[i,a] for a = 1:numCustomer, i = 1:numLocation], (numCustomer*numLocation,))
	#coef::Array{Float64,1} = reshape([abs(customerLocations[a]-i) for a = 1:numCustomer, i = 1:numLocation], (numCustomer*numLocation,))
	#setObjective(m, AffExpr(vars, coef, 0.))
    @setObjective(m, sum{abs(customerLocations[a]-i)*x[i,a], a = 1:numCustomer, i = 1:numLocation} )
	
	# Constraints
	
	for a in 1:numCustomer
	  # Subject to linking x with s
	  for i in 1:numLocation
		#addConstraint(m, 1.0*x[i,a] + (-1.0*s[i]) <= 0 )
		#addConstraint(m, AffExpr([x[i,a],s[i]],[1.,-1.],0.) <= 0)
		@addConstraint(m, x[i,a] - s[i] <= 0)
	  end
	  # Subject to one of x must be 1
	  @addConstraint(m, sum{x[i,a],i=1:numLocation} == 1 )
	end
	
	# Subject to must allocate all facilities
	@addConstraint(m, sum{s[i],i=1:numLocation} == numFacility )	
	toc()

	tic()
	#writeMPS(m,"pmedian_jl.mps")
	writeMPS(m,"/dev/null")
	toc()
end

@assert length(ARGS) == 3

numFacility = int(ARGS[1])
numCustomer = int(ARGS[2])
numLocation = int(ARGS[3])

doTest(numFacility,numCustomer,numLocation)

#m += lpSum( x[(a,i)] * abs(customerLocations[a]-i) for a in range(numCustomer) for i in range(numLocation) )
	#setObjective(m, @sumExpr(abs(customerLocations[a]-i) * x[i,a] for a = 1:numCustomer, i = 1:numLocation) )
	#x = [ [ Variable(m,0,1,0,"x$a,$i") for a = 1:numCustomer ] for i = 1:numLocation ]
