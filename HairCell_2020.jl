# Interactive animation of vestibular hair cell transduction channel
# gating as a function of kinocillium deflection
# Mike Paulin University of Otago 2019
#


using Makie
using MakieLayout
using StatsMakie
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

receptor_alpha = Float32(1e-1)
receptor_channel_conductance = Float32(5.e-3)



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

#**************************************************************
#
#   CONSTRUCT GUI
#
#**************************************************************

CleanAxis(scene) = LAxis(scene,
                          titlevisible = false,
                          xticksvisible = false,
                          xticklabelsvisible = false,
                          xlabelvisible = true,
                          yticksvisible = false,
                          yticklabelsvisible = false,
                          ylabelvisible = false,
                          )

# Scene with 2 panels:
#    upper panel for GUI and animation control
#    lower panel for time series plots
scene = Scene(resolution = (1000,800), camera=campixel!)
scene_layout = GridLayout(scene, 2, 1,
                    rowsizes = [Relative(0.5), Relative(0.5)],
                    alignmode = Outside(30, 30, 30, 30))

# Control panel has 2 columns
#    left for hair cell animation
#    right for Exwald model interface
controlpanel_layout = GridLayout(1,2, colsizes = [Relative(.5), Relative(.5)])

# Hair cell animation pane
haircell_animation_layout= GridLayout(2,1, rowsizes = [Relative(.025), Relative(.975)])
haircell_animation_layout[2, 1] = hc_animation_axis = LAxis(scene)
controlpanel_layout[1,1] = haircell_animation_layout
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
# slider to control kinocilium deflection
kinocilium_slider  = LSlider(scene,      range = LinRange(min_deflect, max_deflect, 100))
haircell_animation_layout[1,1] = kinocilium_slider

# GUI for Exwald model (Goes in controlpanel[1,2])
exwald_layout = GridLayout(3,1,
                  rowsizes = [Relative(.05), Relative(.05), Relative(.9)])

# sliders for Exwald neuron parameters
exwald_layout[1,1] = tau_slider  =
                     LSlider(scene,  range = LinRange(-2.0, 2.0, 100))
exwald_layout[2,1] = lambda_slider =
                     LSlider(scene, range = LinRange(0.0, 5.0, 100))
exwald_layout[3,1] = exwald_axis = LAxis(scene)
exwald_axis.xlabel = "ISI Distribution"
exwald_axis.xgridvisible = false
exwald_axis.ygridvisible = false



controlpanel_layout[1,2] = exwald_layout




# insert control panel in scene
scene_layout[1,1] = controlpanel_layout


# time series plot pane
timeseries_layout = GridLayout(3,2,
                    colsizes = [Relative(.2), Relative(.8)],
                    rowsizes = [Relative(.33), Relative(.34), Relative(.33)])

# receptor current
timeseries_layout[1,2] = receptor_current_axis = CleanAxis(scene)
receptor_current_axis.xlabel = "receptor current"
receptor_current_trace = fill(0.0f0, maxTime+1)
lines!(receptor_current_axis, t, zeros(length(t)), color = :darkred)
receptor_current_plothandle =
         lines!(receptor_current_axis, t, receptor_current_trace, color = :darkcyan)
receptor_current_axis.limits[] = FRect(0., -15., 1001., 16.)

# receptor_potential
timeseries_layout[2,2] = receptor_potential_axis = CleanAxis(scene)
receptor_potential_axis.xlabel = "receptor potential"
RRP = -40.0f0   # resting receptor potential
receptor_potential_trace = fill(RRP, maxTime+1)
lines!(receptor_potential_axis, t, RRP*ones(length(t)), color = :darkred)
receptor_potential_plothandle =
         lines!(receptor_potential_axis, t, receptor_potential_trace,
                color = :darkcyan)
receptor_potential_axis.limits[] = FRect(0., -42.0, 1001., 12.0)

# afferent spike train axis
timeseries_layout[3,2] = spike_axis = CleanAxis(scene)
spike_axis.xlabel = "afferent spike train"
spike_trace = fill(0.0f0, maxTime+1)
spike_plothandle =
         lines!(spike_axis, t, spike_trace,
                color = :darkcyan)
spike_axis.limits[] = FRect(0., 0., 1001., 1.25)

# distributions
timeseries_layout[1, 1] = channel_distn_axis = CleanAxis(scene)
channel_distn_axis.xlabel = "H=0"
channel_distn_histogram = plot!(channel_distn_axis, collect(1:10), zeros(10))
channel_distn_axis.limits[] = FRect(-1., 0., 51., 1.)
timeseries_layout[2, 1] = receptor_distn_axis = CleanAxis(scene)
receptor_distn_axis.xlabel = "H=0"
timeseries_layout[3, 1] = ISI_distn_axis = CleanAxis(scene)
ISI_distn_axis.xlabel = "H=0"

# insert time series pane in scene
scene_layout[2,1] = timeseries_layout

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
    aspect = 1.4*ylim/xlim

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
    @inbounds for j in 1:n
        vertex[j] = vertex[j] - oldpos + newpos
    end
    ngon[1][] = ngon[1][]   # triggers re-draw (ngon is observable)
end

function drawHairCell(panel, x0,y0, state)

  dx = 75.
  dy = .05

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

function shiftinsert(X::Array{Float32,1}, x::Float32 )
  # shift X left and insert x at the end
   n = size(X,1)
   @inbounds for i in 2:n
     X[i-1] =  X[i]
   end
   X[n]=x
   return X
end

function afferentNeuron(u, mu, lambda, tau)
  # Integrate and fire with


end

# animate gate states
# gates flicker open (yellow) and closed (blue)
# WARNING: These graphical objects (including callbacks) persist
#          unless scene is closed before re-running the script
framecount = 0
X = fill(Point2f0(0.,0.), 40)  # buffer for distribution data
Iaff = 0.0   # integrate and fire neuron state
@async while isopen(scene) # run this block as parallel thread
                       # while scene (window) is open

  global framecount = framecount + 1
  global receptor_alpha
  global receptor_channel_conductance
  global RRP
  global Iaff

  threshold = 250.

  # random (Normal) Brownian perturbation to deflection, RMS 2nm
  # nb deflection is an Observable whose (observed) value is deflection[]
  # Similarly randn(1) is a 1-element array of random numbers
  #    and randn(1)[] (or randn(1)[1]) is a random number

  Δk = kinocilium_slider.value[] + 5.0f0*randn(Float32,1)[]

  p = p_open(Δk*nm)          # open probability

  gateOpen = rand(48).<p     # gate states

  receptor_conductance = Float32(sum(gateOpen))*receptor_channel_conductance

  receptor_potential = receptor_potential_trace[end]

  receptor_current = receptor_potential*receptor_conductance

  receptor_potential = RRP + (1.0f0-receptor_alpha)*(receptor_potential - RRP)-
                              receptor_alpha*receptor_current

  # Exwald neuron
  spike = (rand(1)[]< 0.01*p) ? Float32(1.0) : Float32(0.0)


  # # channel state distribution
  # if framecount > 100
  #     D = histogram(receptor_current_trace, nbins=20)
  #     n = length(D.weights)
  #     b = collect(D.edges[1])
  #     w = sum(D.weights)
  #     pdf = D.weights./w
  #     inotzero = findall(x-> x>1.0e-6 && x<(1.0-1.0e-6), pdf)
  #     channel_entropy = -sum(pdf[inotzero].*log.(pdf[inotzero],2))
  #      # println("Channel entropy: ", channel_entropy)
  #
  #     for i in 1:n
  #         X[i] = Point2f0((b[i]+b[i+1])/2., pdf[i])
  #     end
  #     channel_distn_histogram[1][] = X[1:n]
  #     framecount = 0
  # end

  # update display
  movegon(kinocilium_handle, 1, kcx0+Δk*hairScale, kcy0)
  movegon(tracker_handle, 1, Δk, p)
  channel_handle[:color] =
        @inbounds [gateOpen[i] ? :gold1 : :dodgerblue1 for i in 1:48]
  shiftinsert(receptor_current_trace, receptor_current )
  receptor_current_plothandle[2] = receptor_current_trace
  shiftinsert(receptor_potential_trace, receptor_potential )
  receptor_potential_plothandle[2] = receptor_potential_trace
  shiftinsert(spike_trace, spike )
  spike_plothandle[2] = spike_trace



  yield() # allow code below this block to run
          # while continuing to run this block
end

RecordEvents(scene, "output")
