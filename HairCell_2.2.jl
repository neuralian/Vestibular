# Interactive animation of vestibular hair cell transduction channel
# gating as a function of kinocillium deflection
# with Exwald afferent model
# Mike Paulin University of Otago 2020
#


using Makie
using MakieLayout
using StatsMakie
using AbstractPlotting: px
using Colors
using Distributions
using Printf


# Hair cell biophysical parameters

const kᵦ = 1.38e-23  # Boltzmann constant J/K or m^2 kg ^-2 K^-1
const T  = 300.      # temperature K
const z  = 40.e-15   # Gating force 40 fN (Howard, Roberts & Hudspeth 1988)
const d  = 3.5e-9    # Gate swing distance 3.5nm
const pᵣ = 0.15      # resting/spontaneous open state probability
const nano = Float32(1e-9)     #  conversion factor

const haircell_resting_potential = Float32(-0.06) # -60mV (Corey and Hudspeth 1979)
const haircell_input_resistance  = 2.5e8  # 250Mohm (Corey and Hudspeth 1979)
const haircell_capacitance = 3.0e-11 # 30pF (Roberts, Howard and Hudspeth, 1988)
const haircell_transduction_reversal_potential = Float32(-.002) # -2mv (C&H 1979)
const haircell_single_channel_conductance = Float32(500.e-12) # Geleoc &c 1997

# simulation parameters
const dt = 1.0e-4    # 100 microsecond steps

# plot parameters
const kcx0 = 60.   # location of kinocilium in animation axis
const kcy0 = 0.6   # ""
const pRange = 1e-7  # range of probabilities to plot (pRange, 1-pRange)
const hairScale = 0.05 # scale deflection fr]om plot to gate state animation
maxTime = 0.1   # duration of time series plots /s
BGcolor = RGB(.995, .995, .975)
SCcolor = RGB(.95, .95, .95)
open_color = :gold1
closed_color = :cornflowerblue
#RGB(0.5, 0.75, .9)

# exwald pdf compute range needs to extend over all possible interval lengths
# (so convolution of exp and wald can be normalized)
# but we don't want to see the whole range
exwaldpdf_xrange = collect(0.0:0.025:100.0)
exwaldpdf_displayrange = collect(0.0:.025:50.0)

const min_deflect = -500.0
const max_deflect = 1000.0
const adaptation_motor_time_constant = 0.050;
const deflect_range = collect(min_deflect:max_deflect)


# time series plot time base
const t = collect(0:dt:maxTime)
numt = length(t)

# solve p₀(x₀)= 1/2 (deflection when open state prob = 1/2)
const x₀ =  kᵦ*T*log( (1-pᵣ)/pᵣ)/z

# solve p₀(xRange)= pRange to find plot range
const xRange =  kᵦ*T*log( (1-pRange)/pRange)/z

"""
  # Stereocilia channel open probability as a function of kinocilium deflection
"""
p_open(x) = 1.0./(1.0 .+ exp.(-z*(x.-x₀)/(kᵦ*T)))

"""
  # Sample from Exwald density
"""
function exwaldSample(μ, λ, τ)

  return rand(InverseGaussian(μ,λ)) + rand(Exponential(τ))

end

"""
  # Exwald probability density function
"""
function exwaldpdf(μ, λ, τ, x)

    n = length(x)
    dx = x[2]-x[1]  # assuming equal spacing
    e = pdf.(Exponential(τ), x)
    w = pdf.(InverseGaussian(μ,λ), x)
    c = zeros(n)
    for i in 1:n
      for j in 1:i-1
        c[i] = c[i] + e[i-j]*w[j]
      end
    end
    return c./(sum(c)*dx)
end



"""
 # Hair cell type
"""
struct HairCell

  # state variables
   potential::Array{Float32,1}  # receptor potential
   p_open::Array{Float32,1}     # gate open probability
   gateOpen::Array{Bool,1}      # gate states (true=open)
   current::Array{Float32,1}    # receptor channel current
   Δk::Array{Float32,1}         # kinocilium deflection
   Δk_::Array{Float32,1}        # previous kinocilium deflection
   x::Array{Float32,1}          # gating spring extension

   # parameters
   resting_potential::Float32
   resistance::Float32          # passive input resistance
   capacitance::Float32
   nCh::Int64                   # number of transduction channels
   conductance::Float32         # conductance per channel
   reversal_potential::Float32 # for channel currents
   τ_adaptationmotor::Float32

 end

"""
  # Hair cell constructor (overloaded)
"""
function HairCell(resting_potential)

    p = p_open(0.0)
    return HairCell([resting_potential],
                    [p],
                    rand(48).<p,
                    [0.0],
                    [0.0],
                    [0.0],
                    [0.0],

                    haircell_resting_potential,
                    haircell_input_resistance,
                    haircell_capacitance,
                    48,
                    haircell_single_channel_conductance,
                    haircell_transduction_reversal_potential,
                    adaptation_motor_time_constant)
end

"""
  #  Hair cell behaviour
"""
function haircell_stateupdate!(haircell::HairCell, kinocilium_deflection)

  global dt  # not necessary if dt is const, but documents where dt comes from
  global nano  # ditto

  # change in kinocilium defection including 2nm RMS noise
  haircell.Δk_[] = haircell.Δk[]  # previous deflection
  haircell.Δk[] = (kinocilium_deflection + 2.0f0*randn(Float32,1)[])*nano

  haircell.x[] += (haircell.Δk[] - haircell.Δk_[]) -
                dt*haircell.x[]/haircell.τ_adaptationmotor


  haircell.p_open[] = p_open(haircell.x[])  # gate open probability

  haircell.gateOpen[:] = rand(48) .< haircell.p_open[]     # gate states

  conductance = Float32(sum(haircell.gateOpen))*haircell.conductance

  # leak current
  leak_current = (haircell.potential[] - haircell.resting_potential)/
                                                   haircell.resistance
  # channel current
  haircell.current[] = (haircell.potential[]-haircell.reversal_potential)*
                                                                 conductance

  # receptor potential
  haircell.potential[] -= dt*(leak_current +haircell.current[])/haircell.capacitance
end

# Construct a hair cell
haircell = HairCell(0.0)

"""
 # Shannon Entropy (Channel Capacity) (in bits) estimated from a sample signal
 # treating a quantized empirical distribution (relative frequency histogram)
 # as the probability distribution of states.
"""
function sample_entropy(sample)


  # empirical pdf
  D = histogram(sample, nbins=12)
  n = length(D.weights)
  b = collect(D.edges[1])
  w = sum(D.weights)
  pdf = D.weights./w

  # bins with non-zero frequencies
  inonzero = findall(x-> x>1.0e-6 && x<(1.0-1.0e-6), pdf)

  entropy = -sum(pdf[inonzero].*log.(2, pdf[inonzero]))

  return entropy
end
#**************************************************************
#
#   CONSTRUCT GUI
#
#**************************************************************

# Scene with 2 panels:
#    upper panel for GUI and animation control
#    lower panel for time series plots
scene = Scene(resolution = (1000,800), camera=campixel!)
scene_layout = GridLayout(scene, 2, 1,
                    rowsizes = [Relative(0.5), Relative(0.5)],
                    alignmode = Outside(30, 30, 30, 30))
scene.backgroundcolor[] = SCcolor

# Control panel has 2 columns
#    left for hair cell animation
#    right for Exwald model interface
controlpanel_layout = GridLayout(1,2,
        colsizes = [Relative(.55), Relative(.45)])

# Hair cell animation pane
haircell_animation_layout= GridLayout(3,3,
                colsizes = [Relative(.05), Relative(.90), Relative(.05)],
                rowsizes = [Relative(.025), Relative(.025), Relative(.95)])
haircell_animation_layout[3, 1:3] = hc_animation_axis = LAxis(scene)
hc_animation_axis.backgroundcolor[] = BGcolor
hc_animation_axis.xlabel  = "Kinocilium deflection Δk   (nm)"
hc_animation_axis.ylabel = "Open Probability"
hc_animation_axis.xgridvisible = false
hc_animation_axis.ygridvisible = false
# plot open state probability as a function of kinocilium deflection
lines!(hc_animation_axis, deflect_range,  p_open(deflect_range*nano),
           linewidth =3,
           color = RGB(.4, .4, .4),
           leg = false
      )
# slider to control kinocilium deflection
kinocilium_slider  = LSlider(scene,
                    range = LinRange(min_deflect, max_deflect, 100)
                  )
haircell_animation_layout[2,2] = kinocilium_slider
haircell_animation_layout[1, 2] = LText(scene,
        "HAIR CELL", textsize = 14)
haircell_animation_layout[2, 1] = LText(scene, "Δk", textsize = 20)
controlpanel_layout[1,1] = haircell_animation_layout
# GUI for Exwald model
#   located in right sub-panel of control panel
#   2 sliders (log tau, log lambda) in upper sub-layout
#   Exwald pdf in lower sub-layout
exwald_layout = GridLayout(4,3,
                  colsizes = [Relative(.025), Relative(.875), Relative(.1)],
                  rowsizes = [Relative(.025), Relative(.025),
                              Relative(.025), Relative(.9)])

# EXWALD MODEL CONTROL PANEL (TOP RIGHT)
exwald_layout[1,1:3] = LText(scene,
         "AFFERENT INTER-SPIKE INTERVAL DISTRIBUTION", textsize = 14)

# τ slider
exwald_layout[2,1] = LText(scene, "τ", textsize = 20)
exwald_layout[2,2] = tau_slider = LSlider(scene, range=LinRange(-2.0, 2.0, 101))
tau_slider.value[] = 0.0
exwald_layout[2,3] = LText(scene, lift(x->@sprintf("%.2f", 10.0^x),
                    tau_slider.value), textsize = 14)

# λ slider
exwald_layout[3,1] = LText(scene, "λ", textsize = 20)
exwald_layout[3,2] = lam_slider = LSlider(scene,
                                  range=LinRange(2.0, 4.0, 101),
                                  startvalue = 3.5)
exwald_layout[3,3] = LText(scene, lift(x->@sprintf("%.0f", 10.0^x),
                           lam_slider.value), textsize = 14)


# plot of exwald distribution with parameters from sliders
exwald_layout[4,1:3]= exwald_plot_axis = LAxis(scene,
                      yticksvisible = false, yticklabelsvisible = false)
# exwald_plot_axis.limits[]  =FRect(0.0, 0.0, 60.0, 1.0)
exwald_plot_axis.backgroundcolor[] = BGcolor
exwald_plot_axis.xlabel = "Interval (ms)"
exwald_plot_axis.ylabel = "probability density"


exwaldpdfplot_xrange = collect(0.0:.025:250.0)

#  two sliders control exwald parameters & update plot interactively
exwald_pdf_plothandle = plot!(exwald_plot_axis,
    exwaldpdf_displayrange,
    lift((μ_control, λ_control, τ_control) ->
      exwaldpdf( 1.9/p_open(μ_control*nano),
      10.0^λ_control,
      10.0^τ_control, exwaldpdf_xrange)[1:length(exwaldpdf_displayrange)],
      kinocilium_slider.value,
      lam_slider.value,
      tau_slider.value ), linewidth = 3, color = RGB(.4, .4, .4) )

# insert exwald panel into control panel
controlpanel_layout[1,2] = exwald_layout

# insert control panel into scene
scene_layout[1,1] = controlpanel_layout

# time series plot pane
scene_layout[2,1] = timeseries_layout = GridLayout(3,1,
                    rowsizes = [Relative(.33), Relative(.34), Relative(.33)])

# receptor current
timeseries_layout[1,1] = receptor_current_axis = LAxis(scene,
                          titlevisible = false,
                          xticksvisible = false,
                          xticklabelsvisible = false,
                          xlabelvisible = true,
                          yticksvisible = true,
                          yticklabelsvisible = true,
                          ylabelvisible = false
                          )
receptor_current_axis.backgroundcolor[] = BGcolor
receptor_current_axis.xlabel = "receptor current  (nA)"
receptor_current_trace = fill(-0.5f0, numt)
lines!(receptor_current_axis, t, zeros(length(t)), color = :darkred)
receptor_current_plothandle =
          lines!(receptor_current_axis, t,
                 receptor_current_trace, color = :darkcyan )
receptor_current_axis.yticks[] = AutoLinearTicks{Int64}(2)
# tightlimits!(receptor_current_axis)
# receptor_potential
timeseries_layout[2,1] = receptor_potential_axis = LAxis(scene,
                          titlevisible = false,
                          xticksvisible = false,
                          xticklabelsvisible = false,
                          xlabelvisible = true,
                          yticksvisible = true,
                          yticklabelsvisible = true,
                          ylabelvisible = false,
                          )
receptor_potential_axis.backgroundcolor[] = BGcolor
receptor_potential_axis.xlabel = "receptor potential (mV)"
receptor_potential_trace = fill(haircell_resting_potential, numt)
lines!(receptor_potential_axis,
       t, 1000.0*haircell_resting_potential*ones(numt), color = :darkred)
       #
       # display(scene)
       # sleep(5000)
receptor_potential_plothandle =
         lines!(receptor_potential_axis, t, receptor_potential_trace,
                color = :darkcyan)
receptor_potential_axis.yticks[] = AutoLinearTicks{Int64}(2)



# afferent spike train axis
timeseries_layout[3,1] = spike_axis = LAxis(scene,
                          titlevisible = false,
                          xticksvisible = true,
                          xticklabelsvisible = true,
                          xlabelvisible = true,
                          yticksvisible = false,
                          yticklabelsvisible = false,
                          ylabelvisible = false,
                          )
spike_axis.backgroundcolor[] = BGcolor
# spike_axis.xlabel = @sprintf("afferent spike train (%.2fs)", maxTime)
spike_axis.xlabel = "afferent spike train (seconds)"
spike_trace = fill(0.0f0, numt)
spike_plothandle =
         lines!(spike_axis, t, spike_trace, color = :darkcyan)
spike_axis.xticks[] = AutoLinearTicks{Int64}(2)


spike_axis.limits[] = FRect(0., 0., maxTime, 1.25)
receptor_current_axis.limits[] =   FRect(0., -0.5, maxTime, 0.6)
receptor_potential_axis.limits[] = FRect(0., -65.0, maxTime, 65.0)
hc_animation_axis.limits[] = FRect(-500., -0.1, 2500., 1.2)


display(scene)

# """
#    Kluge for problem that scatter() gets the aspect ratio wrong in sublayouts
#    and draws distorted markers. This draws scatterplots uisng n-gon markers.
#    NB The initial code worked out the aspect ratio from the axis limits,
#    which worked in initial tests but for reasons yet to be determined is
#    not correct.  The "aspect" factor needs an additional kluge factor. This
#    seems to depend on the scene layout but at present I adjust it manually.
# """
# function scattergon(ax, x, y, n;
#                     markersize = 1, color = :red,
#                     strokewidth=.1, strokecolor=:black)
#     # scatterplot points x as n-gons
#     # size re. x-axis, aspect ratio adjusted to compensate for dataaspectratio
#     # (workaround for bug(?) in Makielayout)
#     # NB use movegon to reposition these markers (for animation)
#
#     axlim = decompose(Point2f0, ax.limits[])
#     xlim = axlim[4][1] - axlim[1][1]
#     ylim = axlim[4][2] - axlim[1][2]
#     aspect = 1.5*ylim/xlim
#
#     ngonx = [0.5*markersize*cos(2.0*π*i/n) for i in 1:n]
#     ngony = [0.5*markersize*aspect*sin(2.0*π*i/n) for i in 1:n]
#
#     ngon = [[ Point2f0(x[1]+ ngonx[j], y[1] + ngony[j]) for j in 1:n]]
#
#     for i in 2:length(x)
#         push!(ngon,[ Point2f0(x[i]+ ngonx[j], y[i] + ngony[j]) for j in 1:n] )
#     end
#
#     h = poly!(ax, ngon, color = color,
#                 strokewidth=strokewidth, strokecolor=strokecolor, zorder = 1)
#
#     return h    # observable marker vertex coordinates
# end

# """
#   Move polygon created by scattergon() to (x,y)
#   NB an ngon returned from scattergon() is an observable
#     1D array whose entries are nx2 arrays defining n-gon markers.
#     Index i specifies which of these to move to (x,y).
#     This version moves only 1 marker.
# """
# function movegon(ngon, i, x,y)
#
#     vertex = ngon[1][][i]
#     n = size(vertex,1)
#     oldpos = (sum(vertex, dims=1)/n)[1]  # marker location Point2f0
#     newpos = Point2f0(x,y)
#     @inbounds for j in 1:n
#         vertex[j] = vertex[j] - oldpos + newpos
#     end
#     ngon[1][] = ngon[1][]   # triggers re-draw (ngon is observable)
# end

"""
    Function to draw sterocilia (from above) as an array of discs &
    kinocilium as a hexagon
"""
function drawHairCell(panel, x0,y0, state)

  dx = 75.
  dy = .065

  # kinocilium, drawn in scene at (x0, y0)
  # outline (stays in place)
  scatter!(panel, [x0],[y0], markersize = 14px, marker = :hexagon,
          color = :white, strokecolor =:black,  strokewidth=.75)
  kinocilium_handle =   scatter!(panel, [x0],[y0], markersize = 14px,
          marker = :hexagon,  color = RGBA(.75,.25,.5,.5),
          strokecolor =:black,  strokewidth=.75)

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
  c = [state[i] ? open_color : closed_color for i in 1:48]

  channel_handle = scatter!(panel, x,y, markersize = 12px,
    color = c, strokecolor = :black, strokewidth=1)

  # return observable handles to bundle and kinocilium
  return (channel_handle,  kinocilium_handle)
end

# draw hair cell
(channel_handle, kinocilium_handle) =
                drawHairCell(hc_animation_axis,kcx0, kcy0, rand(48).<pᵣ)

# draw state tracker
tracker_handle = scatter!(hc_animation_axis, [0],[p_open(0.)],
                              markersize = 8px, color = RGBA(.75,.25,.5),
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


#


# animate gate states
# gates flicker open (yellow) and closed (blue)
# WARNING: These graphical objects (including callbacks) persist
#          unless scene is closed before re-running the script
framecount = 0
sweepCount = 0
sum_current_entropy = 0.0
sum_potential_entropy = 0.0
X = fill(Point2f0(0.,0.), 40)  # buffer for distribution data

# main simulation loop
interval = 0.0
@async while isopen(scene) # run this block as parallel thread
                       # while scene (window) is open
  #
  global framecount = framecount + 1
  global interval


  haircell_stateupdate!(haircell, kinocilium_slider.value[])
  threshold = 250.

  # Exwald neuron
  interval -= dt  # time to next spike
  spike = (interval <=0) ? true : false
  if spike
    μ = (1.0/3000.0)/(haircell.potential[] - haircell.resting_potential)
    # nb -3.0 in exponent because the sliders are calibrated in ms
    #     but simulation time unit is second
    interval += exwaldSample(μ,
                            10.0^(lam_slider.value[]-3.0),
                            10.0^(tau_slider.value[]-3.0) )
  end

  # update display
  kinocilium_handle[1][] = [Point2f0(kcx0+haircell.Δk[]*hairScale/nano, kcy0)]
  tracker_handle[1][] = [Point2f0(haircell.x[]/nano, haircell.p_open[])]
  #movegon(tracker_handle, 1, haircell.x[]/nano, haircell.p_open[])

  # gate marker is gold if open, blue if closed
  channel_handle[:color] =
         [haircell.gateOpen[i] ? open_color : closed_color for i in 1:48]

  # shift display buffers left, insert new values on right
  shiftinsert(receptor_current_trace, haircell.current[]/nano )
  shiftinsert(receptor_potential_trace, Float32(1000.0*haircell.potential[]) )
  shiftinsert(spike_trace, Float32(spike) )

  # update signal traces on display
  # nb the plot handles are observables, which trigger a redraw when they change
  receptor_current_plothandle[2]   = receptor_current_trace
  receptor_potential_plothandle[2] = receptor_potential_trace
  spike_plothandle[2]              = spike_trace

  yield() # allow other processes to run
end

# redraw the display
RecordEvents(scene, "output")

# nb the main simulation loop and Recordevents() are parallel processes
# that run while the scene window is open.  Recordevents runs when the main
# loop yields, updates the scene if an observable object has changed (which
# is always), and then yields back to the main loop (or another waiting process,
# if there was one). This kind af capablity is why I use Makie for graphics -
# it allows real-time animation of simulations with user (mouse/keyboard)
# interaction.  The downside is painfully slow time-to-first plot ... it takes
# a while to get started.
