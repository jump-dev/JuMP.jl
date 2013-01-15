# pmedian.jl
# Solves the P-Median problem

require("../src/julp.jl")
using Julp

using ArgParse

function doTest()
	
	#s = ArgParseSettings()
    #@add_arg_table s begin
    #    "arg1"
    #    "arg2"
    #    "arg3"
    #end
    #parsed_args = parse_args(s)    
	#numFacility::Int = int(parsed_args["arg1"])
	#numCustomer::Int = int(parsed_args["arg2"])
	#numLocation::Int = int(parsed_args["arg3"])
	
	numFacility = 100
	numCustomer = 100
	numLocation = 1000
	
	customerLocations = [randi(numLocation) for a = 1:numCustomer ]
	
	tic()
	m = Model("max")
	
	# Facility locations
	s = [ Variable(m,0,1,0,"s$i") for i=1:numLocation ]

	# Aux. variable: x_a,i = 1 iff the closest facility to a is at i
	x = [ Variable(m,0,1,0,"x$a,$i") for i = 1:numLocation, a = 1:numCustomer]
	
	# Objective: min distance
	vars::Array{Variable,1} = reshape([x[i,a] for a = 1:numCustomer, i = 1:numLocation], (numCustomer*numLocation,))
	coef::Array{Float64,1} = reshape([abs(customerLocations[a]-i) for a = 1:numCustomer, i = 1:numLocation], (numCustomer*numLocation,))
	setObjective(m, AffExpr(vars, coef, 0.))
	
	# Constraints
	
	for a in 1:numCustomer
	  # Subject to linking x with s
	  for i in 1:numLocation
		#addConstraint(m, 1.0*x[i,a] + (-1.0*s[i]) <= 0 )
		addConstraint(m, AffExpr([x[i,a],s[i]],[1.,-1.],0.) <= 0)
	  end
	  # Subject to one of x must be 1
	  addConstraint(m, @sumExpr([1.0*x[i,a] for i = 1:numLocation]) == 1 )
	end
	
	# Subject to must allocate all facilities
	addConstraint(m, @sumExpr([1.0*s[i] for i = 1:numLocation]) == numFacility )	
	toc()

	tic()
	writeMPS(m,"pmedian_jl.mps")
	toc()
end

doTest()


#m += lpSum( x[(a,i)] * abs(customerLocations[a]-i) for a in range(numCustomer) for i in range(numLocation) )
	#setObjective(m, @sumExpr(abs(customerLocations[a]-i) * x[i,a] for a = 1:numCustomer, i = 1:numLocation) )
	#x = [ [ Variable(m,0,1,0,"x$a,$i") for a = 1:numCustomer ] for i = 1:numLocation ]
