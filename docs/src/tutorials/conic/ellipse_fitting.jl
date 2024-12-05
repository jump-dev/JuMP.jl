# # Example: least-squares fitting of circles and ellipses

# Ellipse fitting is a common task in data analysis and computer vision and is of 
# key importance in many application areas. In this tutorial we show how to fit an
# ellipse to a set of points using a conic optimization approach.

# ## Required packages

# This tutorial uses the following packages:

using JuMP
import LinearAlgebra
import Plots
import SCS

# ## Problem formulation
