using Makie
using MakieLayout
using Distributions

scene = Scene(resolution = (600,600), camera=campixel!)
scene_layout = GridLayout(scene, 2, 1,
                    rowsizes = [Relative(0.25), Relative(0.75)],
                    alignmode = Outside(30, 30, 30, 30))

wald_layout = GridLayout(3,1,
         rowsizes = [Relative(.05), Relative(.05), Relative(.9)])

wald_layout[1,1] = tau_slider = LSlider(scene, range=LinRange(-2.0, 2.0, 101))
tau_slider.value[] = 0.0
wald_layout[2,1] = lam_slider = LSlider(scene, range=LinRange(2.0, 4.0, 101))
wald_layout[3,1]= wald_plot_axis = LAxis(scene,
                                        yticksvisible = false,
                                        yticklabelsvisible = false)
wald_plot_axis.xlabel = "Interval (ms)"
wald_plot_axis.ylabel = "probability density"

scene_layout[2,1] = wald_layout

function waldpdf(μ, λ, x)

    return pdf.(InverseGaussian(10.0^μ,10.0^λ), x)
end


wald_x = collect(0.0:.5:100.0)

exwald_pdf_plothandle = plot!(wald_plot_axis, wald_x,
                  lift((x,y) ->  waldpdf( x, y, wald_x),
                  tau_slider.value, lam_slider.value ) )

RecordEvents(scene, "output")
