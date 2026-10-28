# This script rebuilds the JuMP logos. Install Luxor using
#
#   import Pkg; Pkg.pkg"add Luxor@3"

using Luxor

"""
    _logo_no_text(; color::String)

Creates the three Julia dots with the constraints above.
"""
function _logo_no_text(; color::String)
    setcolor(color)
    setline(8)
    setlinecap("round")
    line(Point(20, 68), Point(120, 28), :stroke)
    line(Point(4, 120), Point(120, 4), :stroke)
    pts = box(Point(86, 84), 40, 40; vertices = true)
    @layer begin
        setcolor(Luxor.julia_red)
        circle(pts[1], 17, :fill)
        setcolor(Luxor.julia_green)
        circle(pts[3], 17, :fill)
        setcolor(Luxor.julia_purple)
        circle(pts[4], 17, :fill)
    end
    return
end

"""
    _jump_text_outlines(; color::String)

Creates the text part of the JuMP logo with the top-aligned "u".
"""
function _jump_text_outlines(; color::String)
    setcolor(color)
    @layer begin
        # "J"
        newsubpath()
        move(Point(34.018, 6.789))
        line(Point(46.768, 6.789))
        line(Point(46.768, 49.619))
        curve(
            Point(46.768, 53.479),
            Point(47.058, 57.329),
            Point(46.118, 61.119),
        )
        curve(
            Point(45.022, 65.714),
            Point(42.273, 69.748),
            Point(38.398, 72.449),
        )
        curve(
            Point(35.706, 74.325),
            Point(32.589, 75.501),
            Point(29.328, 75.869),
        )
        curve(Point(18.823, 77.054), Point(8.995, 70.07), Point(6.658, 59.759))
        line(Point(19.038, 56.599))
        curve(
            Point(19.632, 59.365),
            Point(21.576, 61.655),
            Point(24.208, 62.689),
        )
        curve(
            Point(25.38, 63.124),
            Point(26.637, 63.278),
            Point(27.878, 63.139),
        )
        curve(
            Point(29.023, 63.018),
            Point(30.122, 62.62),
            Point(31.078, 61.979),
        )
        curve(
            Point(33.848, 60.049),
            Point(34.018, 56.639),
            Point(34.018, 53.639),
        )
        closepath()
        # "U"
        move(Point(60.478, 6.789))
        line(Point(72.998, 6.789))
        line(Point(72.998, 30.599))
        curve(
            Point(72.998, 35.209),
            Point(73.308, 38.439),
            Point(73.948, 40.229),
        )
        curve(
            Point(74.504, 41.912),
            Point(75.575, 43.378),
            Point(77.008, 44.419),
        )
        curve(
            Point(78.538, 45.459),
            Point(80.359, 45.988),
            Point(82.208, 45.929),
        )
        curve(
            Point(84.063, 45.982),
            Point(85.891, 45.469),
            Point(87.448, 44.459),
        )
        curve(
            Point(88.948, 43.389),
            Point(90.07, 41.869),
            Point(90.648, 40.119),
        )
        curve(
            Point(91.155, 38.699),
            Point(91.408, 35.659),
            Point(91.408, 30.999),
        )
        line(Point(91.408, 6.819))
        line(Point(103.838, 6.819))
        line(Point(103.838, 27.739))
        curve(
            Point(103.838, 36.349),
            Point(103.148, 42.259),
            Point(101.838, 45.429),
        )
        curve(
            Point(100.362, 49.077),
            Point(97.799, 52.185),
            Point(94.498, 54.329),
        )
        curve(
            Point(91.258, 56.403),
            Point(87.145, 57.439),
            Point(82.158, 57.439),
        )
        curve(
            Point(76.738, 57.439),
            Point(72.362, 56.233),
            Point(69.028, 53.819),
        )
        curve(
            Point(65.63, 51.311),
            Point(63.149, 47.755),
            Point(61.968, 43.699),
        )
        curve(
            Point(61.028, 40.699),
            Point(60.518, 35.259),
            Point(60.518, 27.379),
        )
        closepath()
        # "M"
        move(Point(175.278, 73.459))
        line(Point(187.868, 73.459))
        line(Point(187.868, 6.779))
        line(Point(173.468, 6.779))
        line(Point(152.578, 35.449))
        line(Point(131.688, 6.779))
        line(Point(117.188, 6.779))
        line(Point(117.188, 73.459))
        line(Point(129.778, 73.459))
        line(Point(129.778, 24.969))
        line(Point(151.048, 54.029))
        line(Point(153.528, 54.029))
        line(Point(175.278, 25.059))
        closepath()
        # "P"
        move(Point(199.348, 6.789))
        line(Point(212.808, 6.789))
        curve(
            Point(220.141, 6.789),
            Point(225.395, 7.456),
            Point(228.568, 8.789),
        )
        curve(
            Point(231.747, 10.093),
            Point(234.423, 12.387),
            Point(236.198, 15.329),
        )
        curve(
            Point(238.123, 18.596),
            Point(239.09, 22.339),
            Point(238.988, 26.129),
        )
        curve(
            Point(239.137, 30.302),
            Point(237.848, 34.402),
            Point(235.338, 37.739),
        )
        curve(
            Point(232.823, 40.92),
            Point(229.346, 43.204),
            Point(225.428, 44.249),
        )
        curve(
            Point(222.988, 44.916),
            Point(218.531, 45.249),
            Point(212.058, 45.249),
        )
        line(Point(212.058, 73.509))
        line(Point(199.348, 73.509))
        closepath()
        newsubpath()
        move(Point(212.058, 32.869))
        line(Point(216.058, 32.869))
        curve(
            Point(218.29, 32.974),
            Point(220.524, 32.745),
            Point(222.688, 32.189),
        )
        curve(
            Point(223.879, 31.776),
            Point(224.904, 30.986),
            Point(225.608, 29.939),
        )
        curve(
            Point(226.346, 28.816),
            Point(226.72, 27.493),
            Point(226.678, 26.149),
        )
        curve(
            Point(226.815, 23.864),
            Point(225.653, 21.687),
            Point(223.678, 20.529),
        )
        curve(
            Point(222.258, 19.659),
            Point(219.568, 19.249),
            Point(215.628, 19.249),
        )
        line(Point(212.058, 19.249))
        closepath()
    end
    fillpath()
    return
end

"""
    logo_with_text(;
        filename::String,
        color::String,
        background::String = "",
        verbose::Bool = false,
    )

Create the JuMP logo with the text on the right-hand side.

If `verbose`, print the margin lines.
"""
function logo_with_text(;
    filename::String,
    color::String,
    background::String = "",
    verbose::Bool = false,
)
    width, height, margin = 415, 150, 10
    Drawing(width, height, filename)
    if !isempty(background)
        @layer begin
            translate(O)
            translate(width / 2, height / 2)
            setcolor(background)
            box(O, width - 1, height - 1, 10, :fill)
        end
    end
    @layer begin
        translate(O)
        if verbose
            line(Point(0, margin), Point(width, margin), :stroke)
            line(Point(margin, 0), Point(margin, height), :stroke)
            w = width - margin
            line(Point(w, 0), Point(w, height), :stroke)
            h = height - margin
            line(Point(0, h), Point(width, h), :stroke)
        end
        translate(10, 13)
        _logo_no_text(; color = color)
        translate(130, 40)
        scale(1.1)
        _jump_text_outlines(; color = color)
    end
    finish()
    add_license(filename)
    return
end

"""
    logo_square(; filename::String, color::String, verbose::Bool = false)

Create the JuMP logo with the text vertically beneath the pictogram.

If `verbose`, print the margin lines.
"""
function logo_square(; filename::String, color::String, verbose::Bool = false)
    N, margin, w = 200, 10, 110
    Drawing(N, N, filename)
    @layer begin
        translate(O)
        if verbose
            line(Point(0, margin), Point(N, margin), :stroke)
            line(Point(margin, 0), Point(margin, N), :stroke)
            line(Point(N - margin, 0), Point(N - margin, N), :stroke)
            line(Point(0, N - margin), Point(N, N - margin), :stroke)
            lhs = (N - w) / 2
            line(Point(lhs, 0), Point(lhs, N), :stroke)
            line(Point(N - lhs, 0), Point(N - lhs, N), :stroke)
        end
        translate((N - w) / 2, margin)
        scale(0.9)
        _logo_no_text(; color = color)
        translate(-50, 135)
        scale(0.9)
        _jump_text_outlines(; color = color)
    end
    finish()
    add_license(filename)
    return
end

"""
    logo_no_text(; filename::String, color::String)

Create the JuMP logo without the "JuMP" text.
"""
function logo_no_text(; filename::String, color::String)
    Drawing(150, 150, filename)
    @layer begin
        translate(O)
        translate(10, 10)
        _logo_no_text(; color = color)
    end
    finish()
    add_license(filename)
    return
end

function add_license(filename)
    header = """
    <svg xmlns:svg="http://www.w3.org/2000/svg"
      xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
      xmlns:cc="http://creativecommons.org/ns#"
      xmlns:dc="http://purl.org/dc/elements/1.1/"
    """
    footer = """
    <metadata
        id="metadata21">
        <rdf:RDF>
            <cc:Work
                rdf:about="">
            <dc:creator>
                <cc:Agent>
                <dc:title>JuMP project</dc:title>
                </cc:Agent>
            </dc:creator>
            <dc:title>Logo for JuMP</dc:title>
            <cc:license
                rdf:resource="http://creativecommons.org/licenses/by/4.0/" />
            <dc:description>Logo for JuMP, an algebraic modeling language and a collection of supporting packages for mathematical optimization embedded in the Julia programming language</dc:description>
            <dc:date>August 2018</dc:date>
            </cc:Work>
            <cc:License
                rdf:about="http://creativecommons.org/licenses/by/4.0/">
            <cc:permits
                rdf:resource="http://creativecommons.org/ns#Reproduction" />
            <cc:permits
                rdf:resource="http://creativecommons.org/ns#Distribution" />
            <cc:requires
                rdf:resource="http://creativecommons.org/ns#Notice" />
            <cc:requires
                rdf:resource="http://creativecommons.org/ns#Attribution" />
            <cc:permits
                rdf:resource="http://creativecommons.org/ns#DerivativeWorks" />
            </cc:License>
        </rdf:RDF>
    </metadata>
    </svg>
    """
    file = read(filename, String)
    file = replace(file, "<svg " => header)
    file = replace(file, "</svg>" => footer)
    write(filename, file)
    return
end

color(prefix) = isempty(prefix) ? "#333" : "#ffe"
dark_color(prefix) = isempty(prefix) ? "#ffe" : "#333"

for prefix in ("", "-dark")
    path = joinpath(@__DIR__, "logo$(prefix)")
    logo_with_text(; filename = "$(path)-with-text.svg", color = color(prefix))
    logo_with_text(;
        filename = "$(path)-with-text-background.svg",
        color = color(prefix),
        background = dark_color(prefix),
    )
    logo_square(; filename = "$(path).svg", color = color(prefix))
    logo_no_text(; filename = "$(path)-without-text.svg", color = color(prefix))
end
