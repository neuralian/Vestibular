using MakieLayout
using Makie
using Colors

scene = Scene(resolution = (1200, 900), camera=campixel!)

function scattergon(ax, x, y, n;
                    markersize = 1, color = :red,
                    strokewidth=.1, strokecolor=:black)
    # scatterplot points x as n-gons
    # size re. x-axis, aspect ratio adjusted to compensate for dataaspectratio
    # (workaround for bug(?) in Makielayout)
    # NB use movegon to reposition these markers (for animation)

    axlim = decompose(Point2f0, ax.limits[])
    xlim = axlim[4][1] - axlim[1][1]
    ylim = axlim[4][2] - axlim[1][2]
    aspect = 2.0*ylim/xlim

    ngonx = [0.5*markersize*cos(2.0*π*i/n) for i in 1:n]
    ngony = [0.5*markersize*aspect*sin(2.0*π*i/n) for i in 1:n]

    ngon = [[ Point2f0(x[1]+ ngonx[j], y[1] + ngony[j]) for j in 1:n]]

    for i in 2:length(x)
        push!(ngon,[ Point2f0(x[i]+ ngonx[j], y[i] + ngony[j]) for j in 1:n] )
    end

    h = poly!(ax, ngon, color = color,
                strokewidth=strokewidth, strokecolor=strokecolor)

    return h[1]    # observable marker vertex coordinates
end

function movegon(ngon, i, x,y)
    # move polygon created by scattergon() to (x,y)
    # nb ngon returned from scattergon() is an observable
    #    1D array whose entries are
    #    nx2 arrays defining n-gon markers. Index i specifies which of these
    #    to move to (x,y).  This version moves only 1 marker.

    vertex = ngon[][i]
    n = size(vertex,1)
    oldpos = (sum(vertex, dims=1)/n)[1]  # marker location Point2f0
    newpos = Point2f0(x,y)
    for j in 1:n
        vertex[j] = vertex[j] - oldpos + newpos
    end
    ngon[] = ngon[]   # triggers re-draw (ngon is observable)
end



layout = GridLayout(
    scene, 2, 2,
    colsizes = [Relative(0.75), Relative(0.25)],
    rowsizes = [Relative(0.5), Relative(0.5)],
    alignmode = Outside(30, 30, 30, 30))

layout[1,1] = ax = LAxis(scene)

lines!(ax, [0.0, 81.0], [0.0, 1.0], scale_plot = false)

# scatter!(ax, rand(4), rand(4), markersize = .25)

h = scattergon(ax, [10 20], [.1 .2], 6, markersize = 5, color = :red,
        strokewidth = 1, strokecolor = :green)


display(scene)
sleep(1)
movegon(h, 1, 50., .8)
display(scene)
