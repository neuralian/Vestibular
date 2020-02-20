# Interactive animation of vestibular hair cell transduction channel
# gating as a function of kinocillium deflection
# Mike Paulin University of Otago 2019
#


using Makie
using MakieLayout
using AbstractPlotting
using Colors
using Distributions


# Biophysical parameters
kᵦ = 1.38e-23  # Boltzmann constant J/K or m^2 kg ^-2 K^-1
T  = 300.      # temperature K
z  = 40.e-15   # Gating force 40 fN (Howard, Roberts & Hudspeth 1988)
d  = 3.5e-9    # Gate swing distance 3.5nm
pᵣ = 0.15      # resting/spontaneous open state probability
Nch = 48       # number of gating channels
nm = 1e-9  # nanometers

# simulation parameters

# plot parameters
kcx0 = 60.   # location of kinocilium in scene (pixels from BL)
kcy0 = 0.6   # ""
pRange = 1e-7  # range of probabilities to plot (pRange, 1-pRange)
hairScale = 0.05 # scale deflection from plot to gate state animation
maxTime = 1000   # duration of time series plots (haircell state and spikes)

min_deflect = -500.0
max_deflect = 1000.0
deflect_range = collect(min_deflect:max_deflect)

# time series plot time base
t = collect(0:maxTime)


# solve p₀(x₀)= 1/2 (deflection when open state prob = 1/2)
x₀ =  kᵦ*T*log( (1-pᵣ)/pᵣ)/z

# solve p₀(xRange)= pRange to find plot range
xRange =  kᵦ*T*log( (1-pRange)/pRange)/z

# open state probability as a function of bundle deflection
p_open(x) = 1.0./(1.0 .+ exp.(-z*(x.-x₀)/(kᵦ*T)))

CleanAxis(scene) = LAxis(scene,
                          titlevisible = false,
                          xticksvisible = false,
                          xticklabelsvisible = false,
                          xlabelvisible = false,
                          yticksvisible = false,
                          yticklabelsvisible = false,
                          ylabelvisible = false,
                          )

# create a scene
scene = Scene(resolution = (1000,800), camera=campixel!)
nPlots = 2
nRows = nPlots+1
scene_layout = GridLayout(scene, nRows, 2,
                    colsizes = [Relative(0.5), Relative(0.5)],
                    rowsizes = [Relative(0.5), Relative(0.25), Relative(0.25)],
                    alignmode = Outside(30, 30, 30, 30))

animation_layout = GridLayout(2,1, rowsizes = [Relative(.05), Relative(.95)])
animation_layout[2, 1] = hc_animation_axis = LAxis(scene)

# Hair cell animation pane
scene_layout[1,2]  = animation_layout
hc_animation_axis.xlabel  = "Deflection /nm"
hc_animation_axis.ylabel = "Open Probability"
hc_animation_axis.xgridvisible = false
hc_animation_axis.ygridvisible = false
# plot open state probability as a function of kinocilium deflection
lines!(hc_animation_axis, deflect_range,  p_open(deflect_range*nm),
           linewidth =4,
           color = :darkcyan,
           leg = false,
           limits = FRect(-500., -0.1, 1500., 1.2)
      )
# hc_animation_axis.limits[] = FRect(-500., -0.1, 1500., 1.2)

# hair cell state (depolarization atm)
scene_layout[2, 1:2] = hc_state_axis = CleanAxis(scene)
hc_state = fill(0.0, maxTime+1)
hc_state_plot = lines!(hc_state_axis, t, hc_state, color = :darkcyan)
hc_state_axis.limits[] = FRect(0., -515., 1001., 1530.)

display(scene)

scene_layout[3, 1:2] = afferent_spike_axis = CleanAxis(scene)

# # x-axis for animation pane
# nPts = 100.
# xScale = 1e-9    # x-axis in nm
# x = (x₀ .+ collect((-nPts/2.):(nPts/2.))/nPts*lengthUnit)/xScale

# slider to control kinocilium deflection
kinocilium_slider  = LSlider(scene,
                             range = LinRange(min_deflect, max_deflect, 100))
animation_layout[1,1] = kinocilium_slider

display(scene)

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
    aspect = 1.5*ylim/xlim

    ngonx = [0.5*markersize*cos(2.0*π*i/n) for i in 1:n]
    ngony = [0.5*markersize*aspect*sin(2.0*π*i/n) for i in 1:n]

    ngon = [[ Point2f0(x[1]+ ngonx[j], y[1] + ngony[j]) for j in 1:n]]

    for i in 2:length(x)
        push!(ngon,[ Point2f0(x[i]+ ngonx[j], y[i] + ngony[j]) for j in 1:n] )
    end

    h = poly!(ax, ngon, color = color,
                strokewidth=strokewidth, strokecolor=strokecolor, zorder = 1)

    return h    # observable marker vertex coordinates
end

function movegon(ngon, i, x,y)
    # move polygon created by scattergon() to (x,y)
    # nb ngon returned from scattergon() is an observable
    #    1D array whose entries are
    #    nx2 arrays defining n-gon markers. Index i specifies which of these
    #    to move to (x,y).  This version moves only 1 marker.

    vertex = ngon[1][][i]
    n = size(vertex,1)
    oldpos = (sum(vertex, dims=1)/n)[1]  # marker location Point2f0
    newpos = Point2f0(x,y)
    for j in 1:n
        vertex[j] = vertex[j] - oldpos + newpos
    end
    ngon[1][] = ngon[1][]   # triggers re-draw (ngon is observable)
end

function drawHairCell(panel, x0,y0, state)

  dx = 75.
  dy = .055

  # kinocilium, drawn in scene at (x0, y0)
  # outline (stays in place)
  scattergon(panel, [x0],[y0], 6, markersize = 64,
          color = :white, strokecolor =:black,  strokewidth=.75)
  kinocilium_handle =   scattergon(panel, [x0],[y0], 6, markersize = 64,
            color = RGBA(.75,.25,.5,.5), strokecolor =:black,  strokewidth=.75)

  x = zeros(48)
  y = zeros(48)

  # (x,y) coordinates of stereocilia.
  # 5 columns of 6,   # centre + 2 on each side
  for i in 1:6
    x[i] = x0-i*dx; y[i] = y0;
    x[6+i] = x0-i*dx + dx/2.0; y[6+i]=y0+dy
    x[12+i] = x0-i*dx + dx/2.0; y[12+i]=y0-dy
    x[18+i] = x0 - i*dx; y[18+i] = y0 + 2.0*dy
    x[24+i] = x0 - i*dx; y[24+i] = y0 - 2.0*dy
  end
  # two columns of 5 (one each side)
  for i in 1:5
    x[30+i] = x0 - i*dx - dx/2.0; y[30+i] = y0 + 3.0*dy
    x[35+i] = x0 - i*dx - dx/2.0; y[35+i] = y0 - 3.0*dy
  end
  # ...and two columns of 4
  for i in 1:4
    x[40+i] = x0 - (i+1)*dx; y[40+i] = y0 + 4.0*dy
    x[44+i] = x0 - (i+1)*dx; y[44+i] = y0 - 4.0*dy
  end

  # colours
  c = [state[i] ? :gold1 : :dodgerblue1 for i in 1:48]
  channel_handle = scattergon(panel, x,y, 16,
                              markersize = 52, color = c,
                              strokewidth = .75, strokecolor = :black)

  # return observable handles to bundle and kinocilium
  return (channel_handle,  kinocilium_handle)
end

# draw hair cell
(channel_handle, kinocilium_handle) =
                drawHairCell(hc_animation_axis,kcx0, kcy0, rand(48).<pᵣ)

# draw state tracker
tracker_handle = scattergon(hc_animation_axis, [0],[p_open(0.)], 16,
                              markersize = 42, color = RGBA(.75,.25,.5),
                              strokecolor =:black,  strokewidth=.75)

display(scene)

# animate gate states
# gates flicker open (yellow) and closed (blue)
# WARNING: These graphical objects (including callbacks) persist
#          unless scene is closed before re-running the script
@async while isopen(scene) # run this block as parallel thread
                       # while scene (window) is open

  # random (Normal) Brownian perturbation to deflection, RMS 2nm
  # nb deflection is an Observable whose (observed) value is deflection[]
  # Similarly randn(1) is a 1-element array of random numbers
  #    and randn(1)[] (or randn(1)[1]) is a random number
  Δk = kinocilium_slider.value[] +2.0*randn(1)[]

  p = p_open(Δk*nm)
  gateState = rand(48).<p
  channel_handle[:color] = [gateState[i] ? :gold1 : :dodgerblue1 for i in 1:48]

  movegon(kinocilium_handle, 1, kcx0+Δk*hairScale, kcy0)
  movegon(tracker_handle, 1, Δk, p)

  dScale = .5
  push!(deleteat!(hc_state,1), Δk)
  hc_state_plot[2] = hc_state
  # haircell_handle[1][] = haircell_handle[1][]
  # display(scene)
  sleep(.005)

  yield() # allow code below this block to run
          # while continuing to run this block
end

RecordEvents(scene, "output")
