#############################################################################
# Author:  Louis Luangkesorn <lugerpitt@gmail.com>
# Date:  February 26, 2015
#
# This model determines a set of workforce levels that will most economically
# meet demands and inventory requirements over time. The formulation is
# motivated by the experiences of a large producer in the United States.  The
# data are for three products and 13 numperiods
#
# Problem taken from the Appendix C of the expanded version of Fourer, Gay,
# and Kernighan, A Modeling Language for Mathematical Programming
#############################################################################

function PrintSolution(status, CREWS, HIRE, LAYOFF)
    println("RESULTS:")
    if status == :Optimal
      println("Crews")
      for t = 0:length(CREWS)-1
        print(" $(getValue(CREWS[t])) ")
      end
      println()
      println("Hire")
      for t = 1:length(HIRE)
        print(" $(getValue(HIRE[t])) ")
      end
      println()
      println("Layoff")
      for t = 1:length(LAYOFF)
        print(" $(getValue(LAYOFF[t])) ")
      end
      println()
    else
        println("  No solution")
    end
    println("")
end


using JuMP
#using DataFrames
#using GLPKMathProgInterface
#using Gurobi

####  PRODUCTION SETS AND PARAMETERS  ###

prd = ["18REG"  "24REG" "24PRO"]
# Members of the product group
numprd = length(prd)
pt=	[1.194,	1.509,	1.509]
# Crew-hours to produce 1000 units
pc=	[2304,	2920,	2910]
# Nominal production cost per 1000, used
# to compute inventory and shortage costs


###  TIME PERIOD SETS AND PARAMETERS  ###

firstperiod =  int(1)
			# Index of first production period to be modeled
lastperiod  = int(13)
			# Index of last production period to be modeled
numperiods = firstperiod:lastperiod
# 'planning horizon' := first..last;

###  EMPLOYMENT PARAMETERS  ###

cs = 18
			# Workers per crew
sl =  8
			# Regular-time hours per shift
rtr = 16.00
			# Wage per hour for regular-time labor
otr = 43.85
			# Wage per hour for overtime labor
iw =  8
			# Crews employed at start of first period
dpp =	 [19.5,	19,	20,	19,	19.5,	19,	19,	20,	19,	20,	20,	18,	18]
			# Regular working days in a production period
ol =	 [96,	96,	96,	96,	96,	96,	96,	96,	96,	96,	96,	96,	96]
			# Maximum crew-hours of overtime in a period
cmin =	[0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0]
			# Lower limit on average employment in a period
cmax =	[8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8]
			# Upper limit on average employment in a period
hc =	 [7500,	7500,	7500,	7500,	15000,	15000,	15000,	15000,	15000,	15000,	7500,	7500,	7500]
			# Penalty cost of hiring a crew
lc =	 [7500,	7500,	7500,	7500,	15000,	15000,	15000,	15000,	15000,	15000,	7500,	7500,	7500]
			# Penalty cost of laying off a crew

###  DEMAND PARAMETERS  ###

d18REG=[63.8,	76,	88.4,	913.8,	115,	133.8,	79.6,	111,	121.6,	470,	78.4,	99.4,	140.4,	63.8]
d24REG = [1212,	306.2,	319,	208.4,	298,	328.2,	959.6,	257.6,	335.6,	118,	284.8,	970,	343.8,	1212]
d24PRO = [0,	0,	0,	0,	0,	0,	0,	0,	0,	1102,	0,	0,	0,	0]
dem = Array[d18REG, d24REG, d24PRO]
			# Requirements (in 1000s)
			# to be met from current production and inventory

pro = Array[
  [0,	0,	0,	1,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0],
  [1,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	1,	1],
  [0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0]]
			# true if product will be the subject
			# of a special promotion in the period
###  INVENTORY AND SHORTAGE PARAMETERS  ###

rir =  0.75
			# Proportion of non-promoted demand
			# that must be in inventory the previous period
pir =  0.80
			# Proportion of promoted demand
			# that must be in inventory the previous period
life  =  2
			# Upper limit on number of periods that
			# any product may sit in inventory
cri	=[0.015,	0.015,	0.015]
			# Inventory cost per 1000 units is
			# cri times nominal production cost
crs	=[1.1,	1.1,	1.1]
			# Shortage cost per 1000 units is
			# crs times nominal production cost
iinv	=[82,	792.2,	0]
			# Inventory at start of first period; age unknown
iil = [[max(0,iinv[p] - sum([dem[p][v] for v=firstperiod:t]))
        for t = numperiods]
        for p=1:numprd]
			# Initial inventory still available for allocation
			# at end of period t
function checkpro(product, timeperiod, production, promotionalrate, regularrate)
  if production[product][timeperiod + 1]==1
    value=promotionalrate
  else
    value =regularrate
  end
  value
end
minv = [[dem[p][t+1] * checkpro(p,t, pro, pir, rir) for t=numperiods] for p=1:numprd]
#param minv 'minimum inventory' {p in prd, t in time}
#	      := dem[p,t+1] * (if pro[p,t+1] then pir else rir);
			# Lower limit on inventory at end of period t

prod = Model()

###  VARIABLES  ###
@defVar(prod, Crews[0:lastperiod] >= 0)
			# Average number of crews employed in each period

@defVar(prod, Hire[numperiods] >= 0)
#var Hire{time} >= 0;    # Crews hired from previous to current period

@defVar(prod, Layoff[numperiods]>= 0)
#var Layoff{time} >= 0;  # Crews laid off from previous to current period

@defVar(prod, Rprd[1:numprd, numperiods] >=0)
			# Production using regular-time labor, in 1000s

@defVar(prod, Oprd[1:numprd, numperiods]>=0)
			# Production using overtime labor, in 1000s

@defVar(prod, Inv[1:numprd, numperiods, 1:life] >=0)
			# Inv[p,t,a] is the amount of product p that is
			# a numperiods old -- produced in period (t+1)-a --
			# and still in storage at the end of period t
@defVar(prod, Short[1:numprd, numperiods]>=0)
			# Accumulated unsatisfied demand at the end of period t

###  CONSTRAINTS  ###

@addConstraint(prod, xyconstr[t=numperiods],
               sum{pt[p] * Rprd[p,t], p=1:numprd} <= sl * dpp[t] * Crews[t])
# rlim 'regular-time limit' {t in time}:
			# Hours needed to accomplish all regular-time
			# production in a period must not exceed
			# hours available on all shifts

@addConstraint(prod, xyconstr[t=numperiods],
               sum{pt[p] * Oprd[p,t], p:1:numprd}  <= ol[t])
#olim 'overtime limit' {t in time}:
			# Hours needed to accomplish all overtime
			# production in a period must not exceed
			# the specified overtime limit

@addConstraint(prod, Crews[firstperiod-1]==iw)
#empl0 'initial crew level':  Crews[first-1] = iw;
			# Use given initial workforce

@addConstraint(prod,xyconstr[t=numperiods],
               Crews[t]== Crews[t-1] + Hire[t]-Layoff[t])
#empl 'crew levels' {t in time}:  Crews[t] = Crews[t-1] + Hire[t] - Layoff[t];
			# Workforce changes by hiring or layoffs

@addConstraint(prod, xyconstr[t=numperiods], cmin[t] <= Crews[t])
@addConstraint(prod, xyconstr[t=numperiods], Crews[t] <= cmax[t])
#emplbnd 'crew limits' {t in time}:  cmin[t] <= Crews[t] <= cmax[t];
			# Workforce must remain within specified bounds

@addConstraint(prod, xyconstr[p=1:numprd],
               Rprd[p, firstperiod] + Oprd[p,firstperiod] + Short[p, firstperiod]-Inv[p, firstperiod, 1]
                 == max(0,dem[p][firstperiod]-iinv[p]))
#dreq1 'first demand requirement' {p in prd}:

# NOTE: JuMP xyconstr[] requires that indices be integer at compile time,
# so firstperiod +1 could not be an index within xycontr or triconstr
for t=(firstperiod+1:lastperiod)
  @addConstraint(prod, xyconstr[p=1:numprd],
                 Rprd[p,t] + Oprd[p,t] + Short[p,t] - Short[p,t-1] +
                   sum{Inv[p, t-1, a] - Inv[p,t,a], a=1:life} ==
                   max(0, dem[p][t]-iil[p][t-1]))
end
#dreq 'demand requirements' {p in prd, t in first+1..last}:
			# Production plus increase in shortage plus
			# decrease in inventory must equal demand

@addConstraint(prod, xyconstr[p=1:numprd, t=numperiods],
               sum{Inv[p,t,a] + iil[p][t], a=1:life} >= minv[p][t])
#ireq 'inventory requirements' {p in prd, t in time}:
			# Inventory in storage at end of period t
			# must meet specified minimum

@addConstraint(prod, triconstr[p=1:numprd, v=1:(life-1), a=v+1:life], Inv[p, firstperiod+v-1, a] ==0)

#izero 'impossible inventories' {p in prd, v in 1..life-1, a in v+1..life}:
    #Inv[p,first+v-1,a] = 0;
			# In the vth period (starting from first)
			# no inventory may be more than v numperiods old
			# (initial inventories are handled separately)

@addConstraint(prod, xyconstr[p=1:numprd, t=numperiods],
               Inv[p,t,1] <= Rprd[p,t]+Oprd[p,t])
#ilim1 'new-inventory limits' {p in prd, t in time}:
    #Inv[p,t,1] <= Rprd[p,t] + Oprd[p,t];
			# New inventory cannot exceed
			# production in the most recent period

secondperiod = firstperiod + 1
@addConstraint(prod, triconstr[p=1:numprd, t=2:lastperiod, a=2:life],
                     Inv[p,t,a] <= Inv[p,t-1,a-1])
#ilim 'inventory limits' {p in prd, t in first+1..last, a in 2..life}:
    #Inv[p,t,a] <= Inv[p,t-1,a-1];
			# Inventory left from period (t+1)-p
			# can only decrease as time goes on

###  OBJECTIVE  ###

@setObjective(prod, Min,
                sum{rtr * sl * dpp[t] * cs * Crews[t], t=numperiods} +
                sum{hc[t] * Hire[t], t=numperiods} +
                sum{lc[t] * Layoff[t], t=numperiods} +
                sum{sum{ otr * cs * pt[p] * Oprd[p,t], t=numperiods}, p=1:numprd} +
                sum{sum{sum{cri[p] * pc[p] * Inv[p,t,a], t=numperiods}, p=1:numprd}, a=1:life} +
                sum{sum{ crs[p] * pc[p] * Short[p,t], t=numperiods},p=1:numprd})

#minimize cost:

#    sum {t in time} rtr * sl * dpp[t] * cs * Crews[t] +
#    sum {t in time} hc[t] * Hire[t] +
#    sum {t in time} lc[t] * Layoff[t] +
#    sum {t in time, p in prd} otr * cs * pt[p] * Oprd[p,t] +
#    sum {t in time, p in prd, a in 1..life} cri[p] * pc[p] * Inv[p,t,a] +
#    sum {t in time, p in prd} crs[p] * pc[p] * Short[p,t];

			# Full regular wages for all crews employed, plus
			# penalties for hiring and layoffs, plus
			# wages for any overtime worked, plus
			# inventory and shortage costs

			# (All other production costs are assumed
			# to depend on initial inventory and on demands,
			# and so are not included explicitly.)


production = solve(prod)

PrintSolution(production, Crews, Hire, Layoff)

