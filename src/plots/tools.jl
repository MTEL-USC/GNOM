# Only keep obs part of the cost function
# and interpolate model onto their location
# Function to extract indices of those model boxes (or obs.) that go in comparison
function _locations(xmodel, obs, ux, grd=grd)
    M = interpolationmatrix(grd, obs)
    xmodelatobs = M * xmodel
    xdepthatobs = M * depthvec(grd)
    iobswet = findall(iswet(grd, obs))
    xwetobs = uconvert.(ux, obs.value[iobswet] * upreferred(ux))
    return xmodelatobs, xdepthatobs, iobswet, xwetobs
end




#####
# A recipe for edgy heatmaps

# my recipe for unitful content transects
@recipe(MakieTransect, x1D, ct) do scene
    Theme()
end
midpoints(x) = [(x[i] + x[i+1])/2 for i in eachindex(x)[1:end-1]]
function Makie.plot!(p::MakieTransect)
    x1D = p[:x1D][]
    ct = p[:ct][]
    ex, ey, z = verticalsection2(x1D, grd, ct=ct)
    ex = ustrip.(ex)
    ey = ustrip.(ey)
    z = ustrip.(z)
    #heatmap!(p, ex, ey, permutedims(z); p.attributes...)
    # replacing heatmap with contours to test

    x, y = midpoints(ex), midpoints(ey)
    z_painted = inpaint(z)
    z_itp = LinearInterpolation((y, x), z_painted, extrapolation_bc=Interpolations.Line())
    x_extended = sort([x;ex])
    y_extended = sort([y;ey])
    z_extended = [z_itp(y,x) for y in y_extended, x in x_extended]
    z_extended[2:2:end, 2:2:end] .= z
    contourf!(p, x_extended, y_extended, permutedims(z_extended); p.attributes...)
    #contour!(p, x, y, permutedims(z); color=:black, linestyle=:dashed, linewidth=0.5, levels=p.attributes.levels[][1:5:end])
    p
end
@recipe(MakieScatterTransect, t) do scene
    Theme()
end
function Makie.plot!(p::MakieScatterTransect)
    t = p[:t][]
    x, y, v = scattertransect(t)
    x = ustrip.(x)
    y = ustrip.(y)
    v = ustrip.(v)
    Makie.scatter!(p, x, y; color=v, p.attributes...)
    Makie.scatter!(p, x, y; color=v, strokewidth=0, p.attributes...)
    p
end

# function for polygons of land
mypolys(cl) = [y for s in shp for x in polyvectors(s) for y in polysplit(x, cl)]
function mapit!(ax, cl::Number, polys; kwargs...)
    [poly!(ax, x; kwargs...) for x in polys]
    xlims!(ax, cl .+ (-180,180))
    ylims!(ax, (-90,90))
    ax
end

# function for plotting all transects
function maptransects!(ax, cts, ts, tcol, wlon)
    for it in eachindex(ts)
        ct = sort(OceanographyCruises.shiftlon(CruiseTrack(ts[it]), baselon=wlon))
        ctlon, ctlat = [s.lon for s in ct.stations], [s.lat for s in ct.stations]
        tmp = Makie.lines!(ax, ctlon, ctlat, color=(tcol[it], 0.5), linewidth=4)
        tmp = scatter!(ax, ctlon, ctlat, color=tcol[it], markersize=2)
        tmp = scatter!(ax, ctlon, ctlat, color=tcol[it], markersize=2, strokewidth=0)
        tmp = scatter!(ax, [ctlon[1]], [ctlat[1]], color=tcol[it])
        #text!(ax, ts[it].cruise, position=(ctlon[1]-10, ctlat[1]), textsize=10, font="DejaVu Sans", align=(:right,:center))
        #Label(fig, bbox = ax.scene.px_area, "a", textsize=20, halign=:left, valign=:bottom, padding=(10,10,5,10), font=labelfont, color=labelcol)
        push!(cts, tmp)
    end
end

function R2ε(R)
    R = ustrip.(R)
    #return [(isnan(r) ? 0.0 : (r / R_CHUR - 1)) |> per10000 for r in R]
    return [(r / R_CHUR - 1) |> per10000 for r in R]
end

function myscatter!(ax, args...; kwargs...)
    out = scatter!(ax, args...; kwargs...)
    scatter!(ax, args...; kwargs..., strokewidth=0)
    out
end


# Function to plot ribbon
#function xribbon!(ax, μ, σ, y; color=:red, strokecolor=:red, linewidth = 5, ylims=(5000,0), kwargs...)
#    xl, xr = μ - σ, μ + σ + 0.001abs.(μ) # Adding a small amount for cases where σ=0
#    X = [xl; reverse(xr)]
#    Y = [y; reverse(y)]
#    isempty(X) && return nothing
#    σribbon = poly!(ax, Point2f0.(zip(X,Y)); color=(color, 0.3), strokewidth=0)
#    μline = lines!(ax, μ, y; color=strokecolor, linewidth=linewidth, kwargs...)
#    ylims!(ax, ylims)
#    return nothing
#end
function xribbon!(ax, μ::Vector{T}, σ, y; color=:red, strokecolor=color, linewidth = 10, ylims=(maximum(bottomdepthvec(grd)),0), αribbon=0.3, kwargs...) where T
    u = Unitful.unit(T)
    μ, σ = ustrip.(μ), ustrip.(σ)
    ikeep = .!isnan.(μ)
    μ, σ, y = μ[ikeep], σ[ikeep], y[ikeep]
    xl, xr = μ - σ, μ + σ + 0.001abs.(μ) # Adding a small amount for cases where σ=0
    X = [xl; reverse(xr)]
    Y = [y; reverse(y)]
    isempty(X) && return u
    σribbon = poly!(ax, Point2f0.(zip(X,Y)); color=(color, αribbon), strokewidth=0)
    μline = lines!(ax, μ, y; color=strokecolor, linewidth=linewidth, kwargs...)
    ylims!(ax, ylims) # as tuple, rever
    return u
end

# Function for errobars
function xerrorbars!(ax, μ, σ, y; color=:red, linewidth=10)
    isempty(μ) && return nothing
    lines!(ax, μ, y; color=color, linewidth=linewidth)
    errorbars!(ax, μ, y, σ, σ, color=color, linewidth=linewidth, whiskerwidth=5, direction=:x)
    return nothing
end


# Make vector of model value where obs are and reassign them to the model?
# TODO, maybe comput the mismatch using the nearest neighbor
function mismatch_along_transect(x::AbstractVector{T}, grd::OceanGrid, t::Transect) where {T}
    unitx = Unitful.unit(T)
    itpx = OceanGrids.interpolate(ustrip.(x), grd)
    pros = [mismatch_profile(itpx, unitx, p) for p in t.profiles]
    Transect(tracer=replace(t.tracer, "Observed"=>"model-obs"), cruise=t.cruise, profiles=pros)
end
function mismatch_profile(itpx, unitx, p)
    values = itpx(p.station.lat, p.station.lon, p.depths) * unitx .- p.values
    ikeep = findall(!isnan, values)
    DepthProfile(depths=p.depths[ikeep], station=p.station, values=values[ikeep])
end

function mismatch_map(model, obs, grd)
    ux = Unitful.unit(eltype(model))
    xmodelatobs, xdepthatobs, iobswet, xwetobs = _locations(model, obs, ux)
    x, y = ustrip.(xwetobs), ustrip.(xmodelatobs)
    values = y - x
    return DataFrame(value=values, lat=obs.lat[iobswet], lon=obs.lon[iobswet], depth=obs.depth[iobswet])
end

function generic_ZA!(ax, x2D::Array{T,2}, lat, depth; kwargs...) where T
    u = Unitful.unit(T)
    x2D = ustrip.(x2D)
    #hm = heatmap!(ax, lat, depth, x2D; kwargs...)
    kwargsf = if haskey(kwargs, :filllevels)
        (;kwargs..., levels=kwargs[:filllevels])
    else
        kwargs
    end
    hm = contourf!(ax, lat, depth, x2D; kwargsf...)
    heatmap!(ax, lat, depth, x2D; kwargs...)
    ylims!(ax, get(kwargs, :ylims, (6000,0)))
    contour!(ax, lat, depth, x2D; kwargs..., color=:black)
    return u, hm
end
function generic_ZA!(ax, v::Vector{T}, grd; mask=1, kwargs...) where T
    u = Unitful.unit(T)
    v2D = zonalaverage(ustrip.(v), grd, mask)
    #depth = [0; cumsum(ustrip.(grd.δdepth))]
    #lat = [ustrip(grd.lat[1]) - ustrip.(grd.δlat[1])/2; ustrip(grd.lat) + ustrip.(grd.δlat)/2]
    depth = ustrip.(grd.depth)
    lat = ustrip.(grd.lat)
    return generic_ZA!(ax, v2D, lat, depth; kwargs...)
end



#=============================#
#      Tools for nice maps    #
#=============================#

# Download using RemoteFiles (from Shapefile's tests)
url_physical = "https://github.com/nvkelso/natural-earth-vector/raw/v4.1.0/110m_physical"
natural_earth_dir = joinpath(DataDeps.determine_save_path(""), "natural_earth")
# Note that this is not a DataDep although I use DataDeps to get a good save path
# that should work as long as DataDeps has been used before, which it should have been.

@RemoteFileSet natural_earth "Natural Earth 110m physical" begin
    # polygon
    ne_land_shp = @RemoteFile joinpath(url_physical, "ne_110m_land.shp") dir = natural_earth_dir
    ne_land_shx = @RemoteFile joinpath(url_physical, "ne_110m_land.shx") dir = natural_earth_dir
    ne_land_dbf = @RemoteFile joinpath(url_physical, "ne_110m_land.dbf") dir = natural_earth_dir
    # linestring
    ne_coastline_shp = @RemoteFile joinpath(url_physical, "ne_110m_coastline.shp") dir = natural_earth_dir
    ne_coastline_shx = @RemoteFile joinpath(url_physical, "ne_110m_coastline.shx") dir = natural_earth_dir
    ne_coastline_dbf = @RemoteFile joinpath(url_physical, "ne_110m_coastline.dbf") dir = natural_earth_dir
end
RemoteFiles.download(natural_earth)
shp = Shapefile.Table(joinpath(natural_earth_dir, "ne_110m_land.shp")) |> Shapefile.shapes
# polyvectors function from Julius Krumbiegel
function polyvectors(shp)
    npoints = length(shp.points)
    offsets = shp.parts .+ 1 # 1 based
    map(zip(offsets, vcat(offsets[2:end], npoints+1))) do (from, to)
        [Point(p.x, p.y) for p in view(shp.points, from:to-1)]
    end
end
# some extra functions for polygons from my notebook
function polysplit(p, central_longitude)
	p2 = shift(p)
	if isleft(p, central_longitude - 180)
		[p2]
	elseif isright(p2, central_longitude + 180)
		[p]
	else
		[p, p2]
	end
end
shift(p) = [Point(pt[1]+360, pt[2]) for pt in p] # shifts p 360° to the east
isright(p, lon) = all(pt[1] ≥ lon for pt in p) # check ≥ or >
isleft(p, lon) = all(pt[1] ≤ lon for pt in p) # check ≤ or <



#===================================#
#      Other Presets and tools      #
#=================------============#

# Shared clims and colormaps
εcmap = make_colorscheme([colorant"cyan", colorant"darkblue", colorant"white", colorant"darkred", colorant"yellow"], 2^8)
εclims = (-30, 10)
εlevels = range(extrema(εclims)..., step=0.5)

Ndcmap = cgrad(:viridis)
Ndclims = (0, 40)
Ndlevels = range(extrema(Ndclims)..., step=1)

# shared clims and colormaps for mismatch plots
εcmap2 = :grayC
δεcmap = cgrad(:PiYG_11, rev=true)
δεclims = (-5, 5)
δεlevels = range(extrema(δεclims)..., step=1)

Ndcmap2 = :grayC
δNdcmap = cgrad(:PuOr_10, rev=true)
δNdclims = (-10, 10)
δNdlevels = range(extrema(δNdclims)..., step=1)

# α colormap
αcmap = cgrad(:solar)

# σ cmap ? (TODO comment what σ is)
logσcmap = cgrad(:magma)
logσclims = (-10, -6)
σcmap = cgrad(:inferno)
σclims = (0, 1) # TODO check this value

# Shared lon y z lims and ticks
ylims = (6000,0)
lons, lats = ustrip.(grd.lon), ustrip.(grd.lat)
clon = 204 # central longitdue to avoid split ocean basins
wlon = clon - 180 # west-most lon for lon centered on clon
centerlon(lon, wlon=wlon) = mod(lon - wlon, 360) + wlon
clons = @. centerlon.(lons)
ilon = sortperm(clons)       # sorting permutation for longitudes
sclons = view(clons, ilon) # sorted centered grd longitudes


# Transect plots stuff
# sort transects by distance so that longests are at bottom
Nd_t_sort = [11, 10, 8, 9, 7, 1, 2, 4, 5, 3, 6]
ε_t_sort = [10, 9, 7, 8, 1, 2, 4, 5, 6]
t_ext = 1 # transect extension in km (so it goes from -t_ext to total distance + t_ext)

cl0 = 0 # centered central longitude for dust-region map

clonATL = -40 # centered central longitude for ATL zoom in
wlonATL = clonATL - 180 # west-most lon for lon centered on clon
clonsATL = @. centerlon.(lons, wlonATL)
ilonATL = sortperm(clonsATL)       # sorting permutation for longitudes
sclonsATL = view(clonsATL, ilonATL) # sorted centered grd longitudes

# trying to add nice continents
# Common lat lon ticks and labels
lonticks30 = -180:30:500
lonticks45 = -180:45:500
lonticks60 = -180:60:500
function lonlabelfun(lon)
    lon = mod(round(Int, lon), 360)
    if lon == 0
        "0°"
    elseif lon ≤ 180
        "$(lon)°E"
    else
        "$(360-lon)°W"
    end
end
lonlabelfun(lons::AbstractVector) = lonlabelfun.(lons)
latticks30 = -90:30:90
latticks45 = -90:45:90
function latlabelfun(lat)
    lat = round(Int, lat)
    if lat == 0
        "EQ"
    elseif lat == 90
        "NP"
    elseif lat == -90
        "SP"
    elseif lat > 0
        "$(lat)°N"
    else
        "$(-lat)°S"
    end
end
latlabelfun(lats::AbstractVector) = latlabelfun.(lats)
function mylatlons!(ax, lats, lons)
    ax.xticks = lons
    ax.yticks = lats
    ax.xtickformat = lonlabelfun
    ax.ytickformat = latlabelfun
end
function myxlats!(ax, lats)
    ax.xticks = lats
    ax.xtickformat = latlabelfun
end


# parameters
function plot_param!(ax, p, s, powers_of_ten=0; color=:black, density_color=:lightgray, plotval=true)
    param_unit = units(p, s)
    param_value = getproperty(p, s)
    param_dist = prior(p, s) #/ (ustrip(param_unit, 1.0upreferred(param_unit)))
    xlims = Makie.default_range(param_dist, 0.005)
    xleft = min(xlims.left, param_value)
    xright = max(xlims.right, param_value)
    xs = range(xleft, xright, length=1001)
    ys = pdf.(param_dist, xs)
    band!(ax, xs ./ exp10(powers_of_ten), zeros(length(xs)), ys, color=density_color)
    plotval && vlines!(ax, param_value / exp10(powers_of_ten); color)
    ylims!(ax; low=0.0)
    pow_str = (powers_of_ten==0) ? "" : "×10" * Unitful.superscript(powers_of_ten)
    u_str = (param_unit==NoUnits) ? "" : param_unit
    paren_str = (u_str=="" && pow_str=="") ? "" : (u_str=="" || pow_str=="") ? "($pow_str$u_str)" : "($pow_str $u_str)"
    ax.xlabel = "$s $paren_str"
    #hideydecorations!(ax)
    (u_str == "‱") && xlims!(ax, s == :σ_ε ? (0,20) : (-35,15))
    #(u_str == "pM") && xlims!(ax, (0,300))
    return maximum(ys)
end

function plot_param_dirty!(ax, p, s, initialparms, finalparams, run_num, powers_of_ten=0; color=:black, density_color=:lightgray)
    ymax = plot_param!(ax, p, s, powers_of_ten; color, density_color, plotval=false)
    runs = keys(initialparams)
    refrun = Symbol("run", run_num)

    for run in runs
        initialparam_value = getproperty(getproperty(initialparms, run), s)
        finalparam_value = getproperty(getproperty(finalparams, run), s)
        lines!(ax, kron([finalparam_value, initialparam_value], [1, 1]) / exp10(powers_of_ten), range(0, 1.1ymax, length=4);
            color=(run==refrun ? RGBA(0,0,1,0.3) : RGBA(0,0,0,0.1)),
            linewidth=(run==refrun ? 3 : 1))
    end

    ylims!(ax, (0, 1.1ymax))
end






# ribbon / profiles
αribbon = 0.3


# some options for the plot
labelcol = :black
labelfont = "Dejavu Sans"
panellabels = map(x -> string("(",x,")"), 'a':'l')
#nan_color = :red # really? (It does not show)
nan_color = :transparent ; #:gray20
land_color = :gray40
seafloor_color = :gray20
continent_color = :gray50
water_color = RGBf(0.8, 0.8, 1.0)
markersize = 7
markersize_maps = 4

outer_padding = 30

# density colors
density_colors = ColorSchemes.Pastel1_8[[1,2,3,4,5,7,8]]

# diagnosis plots
maskscmap = ColorSchemes.okabe_ito[[8,6,1,5,2,4,7,3]]
masks2cmap = [ColorSchemes.okabe_ito[[6,5]]; :white]
Ωcmap = ColorSchemes.okabe_ito[[6,5,3]]
Ωcmap2 = ColorSchemes.tableau_10[[1, 3, 6, 7, 5]]
Ωcmap3 = [ColorSchemes.tableau_10[[6, 3, 1, 7, 2, 4]]; :lightgray]


Ndpartitionlvls = 0:10:80
Ndpartitionclims = extrema(Ndpartitionlvls)




# table string and LaTeX stuff
# and print a LaTeX table for paper
function symbol2latex(s)
    s = latexstring(s)
    s = replace(s, "σ" => "\\sigma")
    s = replace(s, "τ" => "\\tau")
    s = replace(s, "α" => "\\alpha")
    s = replace(s, "β" => "\\beta")
    s = replace(s, "ϕ" => "\\phi")
    s = replace(s, "ε" => "\\varepsilon")
    s = replace(s, "′" => "'") # replace spaces with small spaces
    s = replace(s, r"_([A-Za-z]+)_dust" => s"_\\textrm{\1}")
    s = replace(s, r"_([A-Za-z]+)" => s"_\\textrm{\1}")
    s = replace(s, r"_([0-9]+)" => s"_\1")
    s = replace(s, r"_∞" => s"_{\\infty}")
    s = replace(s, r"z₀" => "z_0")
    s = replace(s, r"DSi" => "[\\DSi]")
    s = replace(s, r"₀" => "")
end
function latexify(U)
    str = string(U)
    str = replace(str, "⁻" => s"^-")
    str = replace(str, "‱" => s"\textpertenthousand")
    str = replace(str, "¹" => "1")
    str = replace(str, "²" => "2")
    str = replace(str, r"\^-(?<exp>\d+)" => s"$^{-\g<exp>}$") # add brackets around exponents
    str = replace(str, r"\s" => s"\\,")
    return str
end
latexify(::Nothing) = "~"
latexify(D::Distribution) = string("(", latexify(support(D).lb), ",", latexify(support(D).ub), ")")
latexify(x::Float64) = isfinite(x) ? string(round(Int,x)) : replace(string(x), r"Inf" => s"\\infty")
latexeNd(s) = replace(s, "εNd"=>"\\eNd{}")
latexbool(optimized) = optimized ? "\\checkmark" : ""
function numformat(s)
    s = replace(s, r"e\+0([0-9])" => s"\\times 10^{\1}")
    s = replace(s, r"e\+([0-9]+)" => s"\\times 10^{\1}")
end

retools = false # flag to avoid reloading tools