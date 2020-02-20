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
kcx0 = 600.   # location of kinocilium in scene (pixels from BL)
kcy0 = 600.   # ""
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
                    colsizes = [Relative(0.3), Relative(0.7)],
                    rowsizes = [Relative(0.6), Relative(0.2), Relative(0.2)],
                    alignmode = Outside(30, 30, 30, 30))


# Hair cell animation pane
scene_layout[1,2]    = hc_animation_axis = LAxis(scene)
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
scene_layout[1,1] = kinocilium_slider



display(scene)

function drawHairCell(panel, x0,y0, state)

  dx = 20.
  dy = 16.

  # kinocilium, drawn in scene at (x0, y0)
  scatter!(panel, [x0],[y0],
    marker=:hexagon,
    markersize = 18,
    color =  :white,
    strokecolor =:black,
    strokewidth=.75)

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
  scatter!(panel, x,y,
        marker=:circle,
        color = c,
        markersize = 16,
        strokewidth = .5,
        strokecolor=:black)


  return panel[end]  # return handle to hair cell bundle

end

# draw hair cell (resting state)

haircell_handle = drawHairCell(scene,kcx0, kcy0, rand(48).<pᵣ)

# kinocilium_deflection = kinocilium_slider[end][:value]
# vbox(s1, parent=control_panel)

# draw kinocillium deflection indicator
scatter!(scene, [kcx0+kinocilium_slider.value[]*hairScale, kcy0+0.5],
                [kcx0+0.5,kcy0+p₀(kinocilium_slider.value[])],
                marker = [:hexagon,:hexagon],
                color = [RGB(.5,0.,.5), RGBA(.5,0.,.5,.25)],
                markersize = [18, 16],
                strokewidth = 0.5,
                strokecolor = :black)
kinocilium_handle = scene[end]  # Array{Point{2,Float32},1} coordinates

display(scene)
#
#
#
# #depolarizationPane
# S = hbox(deflectionPane, scene, s1,
#  sizes = [.3, .55, .05], parent = Scene(resolution = (1000, 800)));
#
#
#
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
  haircell_handle[:color] = [gateState[i] ? :gold1 : :dodgerblue1 for i in 1:48]

  kinocilium_handle[1][] = [Point2f0(kcx0+Δk*hairScale, kcy0+0.5),
                            Point2f0(kcx0+Δk,  kcy0+p)]

  dScale = .5
  push!(deleteat!(hc_state,1), Δk)
  hc_state_plot[2] = hc_state
  sleep(.005)


  yield() # allow code below this block to run
          # while continuing to run this block
end

RecordEvents(scene, "output")
