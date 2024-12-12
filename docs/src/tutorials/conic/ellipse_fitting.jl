# # Example: least-squares fitting of circles and ellipses

# Ellipse fitting is a common task in data analysis and computer vision and is of 
# key importance in many application areas. In this tutorial we show how to fit an
# ellipse to a set of points using a conic optimization approach.

# ## Required packages

# This tutorial uses the following packages:

using JuMP
import SCS
import LinearAlgebra
import Images
import Wavelets
using LinearAlgebra
using DSP
using Clustering
using Plots
using LinearOperatorCollection
using RegularizedLeastSquares

# ## Parametrization of an ellipse

# An ellipse is a set of the form:
# ```math
# \begin{equation}
# \mathcal{E} = \{ \xi : (\xi - c)^\top D (\xi - c) = r^2 \}
# \end{equation}
# ```
# where ``c \in \mathbb{R}^2`` is the center of the ellipse, 
# ``D \in \mathbb{R}^{2 \times 2} \succ 0`` is a symmetric positive definite matrix and 
# ``r > 0``.

# We can setup a coordinate system ``(x, y) \in \mathbb{X} \times \mathbb{Y}`` with 
# ``x, y \geq 0``. We use definition (1) to write an ellipse as the root of a quadratic 
# form in homogeneous coordinates:
# ```math
# \begin{equation} 
#     \begin{bmatrix} \xi \\ 1 \end{bmatrix}^T \begin{bmatrix} Q & d \\ d^T & e \end{bmatrix} \begin{bmatrix} \xi \\ 1 \end{bmatrix} = 0
# \end{equation}
# ```
# Where:
# ```math
# \begin{align}
#     Q &= D\\
#     d &= -D c\\
#     e &= c^T D c - r^2
# \end{align}
# ```
# The residual distance ``r_0`` of a random point ``\xi_0 = (x_0, y_0)`` to the ellipse is 
# then given by:
# ```math
# \begin{equation}
#     r_0 \triangleq \begin{bmatrix} \xi_0 \\ 1 \end{bmatrix}^T \begin{bmatrix} Q & d \\ d^T & e \end{bmatrix} \begin{bmatrix} \xi_0 \\ 1 \end{bmatrix}
# \end{equation}
# ```
# The value of ``r_0`` is positive if the point is outside the ellipse, zero if it is on
# the ellipse and negative if it is inside the ellipse. We also see we only need 6 
# parameters to uniquely define an ellipse.

# ## Helper functions

# We define some helper functions to help us visualize the results.
# <details>
# <summary>Functions</summary>
function plot_dwt(x, sz = (500, 500))
    return plt = heatmap(
        x;
        color = :grays,
        aspect_ratio = 1,
        cbar = false,
        xlims = (0, size(x, 2)),
        ylims = (0, size(x, 1)),
        size = sz,
        dpi = 300,
    )
end

function normalize(x)
    x = (x .- minimum(x)) ./ (maximum(x) - minimum(x))
    return x
end
# </details>

# ## Reading the test image

# To test our ellipse-fitting algorithm we need a test image with elliptical features.
# For our test image we will use an image of the cartwheel galaxy, 
# [captured by the James Webb Space Telescope](https://webbtelescope.org/contents/media/images/2022/039/01G8JXN0K2VBQP112RNSQWTCTH).
# Galaxies [come in many shapes and sizes](https://science.nasa.gov/universe/galaxies/types/), elliptical being one of them.

# This is just a toy problem with little scientific value, but you can imagine how the 
# rotation and position of elliptical galaxies can be useful information to astronomers.
img = load("docs/src/assets/cartwheel_galaxy.png");

# We convert the image to gray scale so that we can work with a single channel.

img_gray = Gray.(img)
mosaicview(img, img_gray; nrow = 1)

# Instead of operating on the entire image, we select a region of interest (roi) which 
# is a subset of ``\mathbb{X} \times \mathbb{Y}``.
sz = 256
X_c = 600
Y_c = 140

X = X_c:X_c+sz-1
Y = Y_c:Y_c+sz-1
roi = (X, Y)

img_roi = img[roi...]
img_gray_roi = img_gray[roi...]
mosaicview(img_roi, img_gray_roi; nrow = 1)

# ## Extracting image features

# We cannot directly fit ellipses to the image, so we need to extract features that
# enable us to find the elliptical galaxies. 

# The first step is to find a sparse representation of the image. We will use the 
# discrete wavelet transform (DWT) in combination with the Iterative Shrinking and 
# Thresholding (ISTA) algorithm to denoise the image and find a sparse representation.
# This will remove redundant information and make it much easier to detect the edges of
# galaxies.

# Finding a sparse representation amounts to solving the following optimization problem:
# ```math
# \begin{equation}
#     \min_{x} \frac{1}{2} \left\| y - \Phi x \right\|_2^2 + \lambda \left\| x \right\|_1
# \end{equation}
# ```
# where ``y`` is the noisy image, ``\Phi`` is the sparsifying basis, ``x`` is the sparse
# representation of our image, and ``\lambda`` is the regularization parameter which 
# we set to 0.1.

# To work with our image we must first convert it to Float64.
x = convert(Array{Float64}, img_gray_roi)

# We then use ISTA in combination with our wavelet sparsifying basis ``\Psi`` obtained 
# from the [family of Daubechies wavelets](https://en.wikipedia.org/wiki/Daubechies_wavelet). 
# We use the db4 wavelet which has 4 vanishing moments. We set the number of iterations 
# to 15.
reg = L1Regularization(0.1);
Φ = WaveletOp(Float64; shape = size(x), wt = wavelet(WT.db4));
solver = createLinearSolver(OptISTA, Φ; reg = reg, iterations = 15);

# The sampled image in wavelet domain is given by:
b = Φ * vec(x);

# We can now solve the optimization problem to find the sparse representation of the image.
x_approx = solve!(solver, b)
x_approx = reshape(x_approx, size(x));
x_final = normalize(x_approx)

mosaicview(x, Gray.(x_final); nrow = 1)

# We then use a binarization algorithm to map each grayscale pixel 
# ``(x_i, y_i)`` to a binary value so ``x_i, y_i \to \{0, 1\}``.
x_bin = binarize(x_final, Otsu(); nbins = 128)
x_bin = convert(Array{Bool}, x_bin)

plt = plot_dwt(img_roi, (500, 500))
heatmap!(x_bin; color = :grays, alpha = 0.45)

# ## Edge detection and clustering

# Now that we have our binary image, we can use edge detection to find the edges of the
# galaxies. We will use the [Sobel operator](https://en.wikipedia.org/wiki/Sobel_operator) for this task.
# <details>
# <summary>Sobel operator</summary>
function edge_detector(
    f_smooth::Matrix{Float64},
    d1::Float64 = 0.1,
    d2::Float64 = 0.1,
)
    rows, cols = size(f_smooth)

    gradient_magnitude = zeros(Float64, rows, cols)
    laplacian_magnitude = zeros(Float64, rows, cols)
    edge_matrix = zeros(Bool, rows, cols)  # Binary matrix for edges

    sobel_x = [-1 0 1; -2 0 2; -1 0 1]
    sobel_y = [-1 -2 -1; 0 0 0; 1 2 1]
    sobel_xx = [-1 2 -1; 2 -4 2; -1 2 -1]
    sobel_yy = [-1 2 -1; 2 -4 2; -1 2 -1]

    gradient_x = conv(f_smooth, sobel_x)
    gradient_y = conv(f_smooth, sobel_y)

    gradient_magnitude = sqrt.(gradient_x .^ 2 + gradient_y .^ 2)

    gradient_xx = conv(f_smooth, sobel_xx)
    gradient_yy = conv(f_smooth, sobel_yy)
    laplacian_magnitude = sqrt.(gradient_xx .^ 2 + gradient_yy .^ 2)

    for i in 1:rows
        for j in 1:cols
            if gradient_magnitude[i, j] > d1 && laplacian_magnitude[i, j] < d2
                edge_matrix[i, j] = true
            end
        end
    end

    return edge_matrix
end
# </details>

# We apply the Sobel operator to the binary image:
edges = edge_detector(convert(Matrix{Float64}, x_bin), 1e-1, 1e2)
edges = thinning(edges; algo = GuoAlgo())

# And finally we cluster the edges using dbscan so we can fit ellipses to individual 
# galaxies. We can control the minimum size of galaxies by changing the minimum cluster 
# size.
points = findall(edges)
points = getfield.(points, :I)
points = hcat([p[1] for p in points], [p[2] for p in points])'
result = dbscan(
    convert(Matrix{Float64}, points),
    3.0;
    min_neighbors = 2,
    min_cluster_size = 20,
)

# The result of the clustering is a list of clusters to which we will assign a unique 
# color. Each cluster is a list of points that belong to the same galaxy.
clusters = result.clusters
N_clusters = length(clusters)

# <details>
# <summary>Plotting code</summary>
colors = distinguishable_colors(N_clusters + 1)[2:end]
plt = plot_dwt(img_roi, (600, 600))

for (i, cluster) in enumerate(clusters)
    p_cluster = points[:, cluster.core_indices]
    scatter!(
        p_cluster[2, :],
        p_cluster[1, :];
        color = colors[i],
        label = false,
        markerstrokewidth = 0,
        markersize = 1.5,
    )
end

plot!(
    plt;
    axis = false,
    legend = :topleft,
    legendcolumns = 1,
    legendfontsize = 12,
)
# </details>

# ## Fitting ellipses

# Now that we have all the ingredients we can finally start fitting ellipses. We will
# use a conic optimization approach to do so since it is a very natural way to represent
# ellipses.

# First, we define the residual distance definition (6) of a point to an ellipse in JuMP:
function ellipse(Ξ::Array{Tuple{Int,Int},1}, ϵ = 1e-4)
    model = Model(SCS.Optimizer)
    N = length(Ξ)

    @variable(model, Q[1:2, 1:2], PSD)
    @variable(model, d[1:2])
    @variable(model, e)

    @constraint(model, Q - ϵ * I in PSDCone())

    @expression(
        model,
        r[i = 1:N],
        [Ξ[i][1], Ξ[i][2], 1]' * [Q d; d' e] * [Ξ[i][1], Ξ[i][2], 1]
    )

    return model
end

# ## Objective 1: Minimize the total squared distance

# For our first objective we will minimize the total squared distance of all points to
# the ellipse. Hence we will use the sum of the squared distances as our objective 
# function, also known as the ``L^2`` norm.
# ```math
# \begin{equation}
#     \min_{Q, d, e} P_\text{res}(\mathcal{E}) = \sum_{i \in N} d^2_{\text{res}}(\xi_i, \mathcal{E}) = ||d_{\text{res}}||^2_2
# \end{equation}
# ```

# This problem is equivalent to:
# ```math
# \begin{align}
#     \min_{Q, d, e, ζ} &ζ \\
# \text{s.t.} \quad ζ \geq d^2_{\text{res}}(\xi_i, \mathcal{E}) &\quad \forall i \in N
# \end{align}
# ```
# And hence can be modelled as a second-order cone program (SOCP) using 
# [`MOI.RotatedSecondOrderCone`](@ref) as follows:

ellipses_C1 = Vector{Dict{Symbol,Any}}()

for (i, cluster) in enumerate(clusters)
    p_cluster = points[:, cluster.core_indices]
    Ξ = [(point[1], point[2]) for point in eachcol(p_cluster)]

    model = ellipse(Ξ)
    @variable(model, ζ >= 0)
    @constraint(
        model,
        [1 / 2; ζ; model[:r]] in
        MOI.RotatedSecondOrderCone(2 + length(model[:r]))
    )
    obj = @objective(model, Min, ζ)
    set_objective_function(model, obj)
    optimize!(model)

    Q, d, e = value.(model[:Q]), value.(model[:d]), value.(model[:e])
    push!(ellipses_C1, Dict(:Q => Q, :d => d, :e => e))
end

# <details>
# <summary>Plotting code</summary>
W, H = size(img_roi)
x_range = 0:1:W
y_range = 0:1:H
X, Y = [x for x in x_range], [y for y in y_range]
Z = zeros(length(X), length(Y))

function ellipse_eq(x, y, Q, d, e)
    for i in 1:length(x)
        for j in 1:length(y)
            ξ = [x[i], y[j]]
            Z[i, j] = [ξ; 1.0]' * [Q d; d' e] * [ξ; 1.0]
        end
    end
    return Z
end

for ellipse in ellipses_C1
    Q, d, e = ellipse[:Q], ellipse[:d], ellipse[:e]
    Z_sq = ellipse_eq(X, Y, Q, d, e)
    contour!(
        plt,
        x_range,
        y_range,
        Z_sq;
        levels = [0.0],
        linewidth = 2,
        color = :blue,
        cbar = false,
    )
end

display(plt)
# </details>

# ## Objective 2: Minimize the maximum residual distance

# For our second objective we will minimize the maximum residual distance of all points 
# to the ellipse.
# ```math	
# \begin{equation}
#     \min \max_{\xi_i \in \mathcal{F}} d_\text{res}(\xi_i, \mathcal{E}) = 
# \min ||d_\text{res}||_\infty
# \end{equation}
# ```

# This objective can be implemented in JuMP using [`MOI.NormInfinityCone`](@ref) as 
# follows:
ellipses_C2 = Vector{Dict{Symbol,Any}}()

for (i, cluster) in enumerate(clusters)
    p_cluster = points[:, cluster.core_indices]
    Ξ = [(point[1], point[2]) for point in eachcol(p_cluster)]

    model = ellipse(Ξ)
    N = length(Ξ)

    @variable(model, ζ)
    @constraint(
        model,
        [ζ; model[:r]] in MOI.NormInfinityCone(1 + length(model[:r]))
    )

    obj = @objective(model, Min, ζ)
    set_objective_function(model, obj)
    optimize!(model)

    Q, d, e = value.(model[:Q]), value.(model[:d]), value.(model[:e])
    push!(ellipses_C2, Dict(:Q => Q, :d => d, :e => e))
end

# <details>
# <summary>Plotting code</summary>
for ellipse in ellipses_C2
    Q, d, e = ellipse[:Q], ellipse[:d], ellipse[:e]
    Z_sq = ellipse_eq(X, Y, Q, d, e)
    contour!(
        plt,
        x_range,
        y_range,
        Z_sq;
        levels = [0.0],
        linewidth = 2,
        color = :red,
        cbar = false,
    )
end
scatter!([0], [0]; color = :blue, label = "Squared (Objective 1)")
scatter!([0], [0]; color = :red, label = "Min-Max (Objective 2)")
display(plt)
# </details>