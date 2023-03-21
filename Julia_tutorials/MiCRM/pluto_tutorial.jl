### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 55cdaa2d-aaf8-4fcf-946f-4f21a8b2ac93
using Pkg, Random

# ╔═╡ 618fbb5e-642d-4459-a86e-8391badd2898
# Activate Environment
Pkg.activate(".")

# ╔═╡ e8eec266-342b-402d-aa59-6c4f4cf56794
Pkg.add("Distributions")

# ╔═╡ 474974b4-834e-4c22-8e8f-d1e2db53d856
Pkg.add("PlutoUI")

# ╔═╡ 766fa534-e57e-4dd1-a711-022bbbb3e5e4
Pkg.add("ModelingToolkit")

# ╔═╡ f12769a5-86ee-4cee-a26d-6b69f4aed8ac
using MiCRM

# ╔═╡ 1399e02b-46c1-4f35-9c8c-2be205491f4c
using Distributions

# ╔═╡ 87cd2475-cfa5-4113-8494-10d799540cb9
using DifferentialEquations, Plots, PlutoUI; plotly()

# ╔═╡ 755758bd-32ba-4751-b2e2-5eaed5ebbebd
md"# MiCRM Model exploration and parameter simulation"

# ╔═╡ ad4420f9-6ab3-45c8-a99e-71c2b88ddeda
md"## Installing, loading and activating environment!!!"

# ╔═╡ 9add2123-041c-4fcf-b3b0-b5063f6e1e49
#Random.seed!(1)

# ╔═╡ a7524bd1-0374-4412-aada-0b8a9eaa556c
#Pkg.instantiate()

# ╔═╡ 17a6bcbc-f894-485e-a815-36d0f10f9ccb
# Check working directory
pwd()

# ╔═╡ 337a95c2-b6d9-4112-b894-6104d08cc867


# ╔═╡ 4066f57f-3b37-4e78-acc6-d00e6d6ed2a2
md"## MiCRM Walkthrough"

# ╔═╡ a8e77d1b-84ad-401c-a48c-880cadf8d855
md"### Adjust the parameter settings to interact with Plots"

# ╔═╡ 0dd53345-c083-40c5-abcf-1979168fdd12
md"##### Any change will automatically generate a solution"

# ╔═╡ fe152d3d-54cb-45b9-9bd0-d5b1313c3632
md"Time Span $(@bind tspan Slider(1.0:5000.0, 1000.0, true))"

# ╔═╡ 47e0f8e0-ac81-482f-8adb-0eee013f783b
md"N = $(@bind N Slider(1:20, 10, true))"

# ╔═╡ 8206b882-be11-4a1a-b74b-2e28cecd7bad
md"M = $(@bind M Slider(1:50, 10, true))"

# ╔═╡ 8557fdad-32da-489f-b8ea-eea51a7c8069
md"leakage = $(@bind leakage Slider(0:0.01:1, 0.5, true))"

# ╔═╡ 4aecd61f-4bb5-43b7-a427-80ef501144c3
md"'Generalists' $(@bind gs_ratio Slider(1.0:0.5:10.0, 1.0, true)) 'Specialists'"

# ╔═╡ 1af130d3-b15e-46d9-baff-32292fffff75
md"Below is defined the distribution of metabolic by-products:"

# ╔═╡ 52cb1686-0946-4576-985a-7eddf0bb7a02
md"'Uniform' $(@bind by_product_ratio Slider(1.0:0.5:10.0, 1.0, true)) 'Structured'"

# ╔═╡ aba26bbc-c8cf-4d66-9334-47e64f8c3b0c
mm(N,M,kw) = repeat([0.2], N)

# ╔═╡ dd2adbfb-8037-407a-984e-0add743af01a
md"### Below we can inspect the specific values generated for the model"

# ╔═╡ 85fc8f14-4f12-46ef-aaa5-a91d4fa9c7b2


# ╔═╡ 0ce53e19-00c6-426f-9bd0-0d68eb2e3d30
md"### ODEProblem is defined below:"

# ╔═╡ 28b78450-b8ea-4077-b3aa-eae4c2e72ff3
md"##### Defining Consumers Biomass and Resource Concentrations"

# ╔═╡ ac7243ee-2029-4ec2-bd24-f4da56bd397b
x0 = ones(N+M)

# ╔═╡ 758a53f1-6aa6-4721-85e6-0029266c3e4b
typeof(x0)

# ╔═╡ 751908c5-a61d-42a5-be48-07e93ebc5bef


# ╔═╡ 03c88c00-a959-44e0-967e-f52606d6a1a9
md"### Extra parameters, such as `mortality`, modular `leakage` and `uptake`, are generated below:"

# ╔═╡ 159de021-3a9f-44ec-a6ca-7f8be7928433
md"**Mortality:**"

# ╔═╡ 360a4869-4b26-4d8f-b6cb-963b774c1f96
#mortality_rate = rand(Uniform(0.3, 0.5), N)'

# ╔═╡ 1b49575f-7c4c-4ce3-a428-5d4c61f90350
#mm(N,M,kw) = copy(mortality_rate)

# ╔═╡ 10456478-3acd-4632-8f8c-bcb6dca086a0
# Defining mortality parameter
#mm(N,M,kw) = copy(rand(Uniform(0.0, 0.5), N)')

# ╔═╡ 77f32659-99f9-44e9-8a58-9e53e731bec8


# ╔═╡ 37805ecf-c17e-480f-a16f-3b408c954c66
md"**Uptake:**"

# ╔═╡ 6fe70ddf-26aa-40b7-9db8-8eac583a1635
mod_up(N,M,kw) = MiCRM.Parameters.modular_uptake(N,M, N_modules = 1, s_ratio = gs_ratio)

# ╔═╡ 56bb9573-920b-493e-9fa5-5f7a3f8c83fb


# ╔═╡ 23e811fb-49b6-4590-9006-320981c32fda
md"**leakage:**"

# ╔═╡ 6c578beb-3f85-4d68-86bc-a4cec36ec923
mod_leak(N,M,kw) = MiCRM.Parameters.modular_leakage(M, N_modules = 1, s_ratio = by_product_ratio, λ = leakage)

# ╔═╡ 8bbe837b-335f-49df-95bb-ad7f7d25d662
p = MiCRM.Parameters.generate_params(N, M, f_m = mm,
											f_u = mod_up,
											f_l = mod_leak,
											λ = leakage)

# ╔═╡ 6df65794-66f1-46bb-85a0-1ef42aff1aaf
prob = ODEProblem(MiCRM.Simulations.dx!, x0, tspan, p)

# ╔═╡ 304832ac-e1c6-42a1-8c05-9382e66671a6
sol = solve(prob,AutoTsit5(Rosenbrock23()))

# ╔═╡ eb821c59-66ea-439e-b0bd-ee33d15d0272
plot(sol,
	widen = true,
	ylims = (0.0, 15),
	xlims = (0.0, tspan),
	linewidth = 2)

# ╔═╡ 46f739e4-1941-4d4d-b28d-9b9c49fd4dc8


# ╔═╡ 6171b71a-1bb7-486f-8fc8-13f4176d8272
md"## Wandering with the model..."

# ╔═╡ Cell order:
# ╟─755758bd-32ba-4751-b2e2-5eaed5ebbebd
# ╠═ad4420f9-6ab3-45c8-a99e-71c2b88ddeda
# ╠═55cdaa2d-aaf8-4fcf-946f-4f21a8b2ac93
# ╠═e8eec266-342b-402d-aa59-6c4f4cf56794
# ╠═474974b4-834e-4c22-8e8f-d1e2db53d856
# ╠═9add2123-041c-4fcf-b3b0-b5063f6e1e49
# ╠═a7524bd1-0374-4412-aada-0b8a9eaa556c
# ╠═f12769a5-86ee-4cee-a26d-6b69f4aed8ac
# ╠═1399e02b-46c1-4f35-9c8c-2be205491f4c
# ╠═87cd2475-cfa5-4113-8494-10d799540cb9
# ╠═17a6bcbc-f894-485e-a815-36d0f10f9ccb
# ╠═618fbb5e-642d-4459-a86e-8391badd2898
# ╟─337a95c2-b6d9-4112-b894-6104d08cc867
# ╟─4066f57f-3b37-4e78-acc6-d00e6d6ed2a2
# ╟─a8e77d1b-84ad-401c-a48c-880cadf8d855
# ╟─0dd53345-c083-40c5-abcf-1979168fdd12
# ╟─fe152d3d-54cb-45b9-9bd0-d5b1313c3632
# ╟─47e0f8e0-ac81-482f-8adb-0eee013f783b
# ╟─8206b882-be11-4a1a-b74b-2e28cecd7bad
# ╟─8557fdad-32da-489f-b8ea-eea51a7c8069
# ╟─4aecd61f-4bb5-43b7-a427-80ef501144c3
# ╟─1af130d3-b15e-46d9-baff-32292fffff75
# ╟─52cb1686-0946-4576-985a-7eddf0bb7a02
# ╠═aba26bbc-c8cf-4d66-9334-47e64f8c3b0c
# ╠═eb821c59-66ea-439e-b0bd-ee33d15d0272
# ╟─dd2adbfb-8037-407a-984e-0add743af01a
# ╠═8bbe837b-335f-49df-95bb-ad7f7d25d662
# ╟─85fc8f14-4f12-46ef-aaa5-a91d4fa9c7b2
# ╟─0ce53e19-00c6-426f-9bd0-0d68eb2e3d30
# ╟─28b78450-b8ea-4077-b3aa-eae4c2e72ff3
# ╟─ac7243ee-2029-4ec2-bd24-f4da56bd397b
# ╠═758a53f1-6aa6-4721-85e6-0029266c3e4b
# ╠═6df65794-66f1-46bb-85a0-1ef42aff1aaf
# ╠═304832ac-e1c6-42a1-8c05-9382e66671a6
# ╟─751908c5-a61d-42a5-be48-07e93ebc5bef
# ╟─03c88c00-a959-44e0-967e-f52606d6a1a9
# ╟─159de021-3a9f-44ec-a6ca-7f8be7928433
# ╠═360a4869-4b26-4d8f-b6cb-963b774c1f96
# ╠═1b49575f-7c4c-4ce3-a428-5d4c61f90350
# ╠═10456478-3acd-4632-8f8c-bcb6dca086a0
# ╟─77f32659-99f9-44e9-8a58-9e53e731bec8
# ╟─37805ecf-c17e-480f-a16f-3b408c954c66
# ╠═6fe70ddf-26aa-40b7-9db8-8eac583a1635
# ╟─56bb9573-920b-493e-9fa5-5f7a3f8c83fb
# ╟─23e811fb-49b6-4590-9006-320981c32fda
# ╠═6c578beb-3f85-4d68-86bc-a4cec36ec923
# ╟─46f739e4-1941-4d4d-b28d-9b9c49fd4dc8
# ╟─6171b71a-1bb7-486f-8fc8-13f4176d8272
# ╠═766fa534-e57e-4dd1-a711-022bbbb3e5e4
