# obtain a reasonable initial guess for the optimization route in fitting adsorption 
#  models to adsorption isotherm data
function _guess(
    df::DataFrame,
    pressure_col_name::Symbol,
    loading_col_name::Symbol,
    model::Symbol
)
    n = df[!, loading_col_name]
    p = df[!, pressure_col_name]
    if model == :langmuir || model == :henry
        # Henry coefficient H guess as the secant line between origin and first data point 
        #  of the adsorption isotherm
        if isapprox(n[1], 0.0)
            H0 = n[2] / p[2]
        else
            H0 = n[1] / p[1]
        end
        if model == :henry
            return Dict("H0" => H0)
        else
            # saturation loading M is largest value of adsorption observed plus 10%
            M0 = n[end] * 1.1
            K0 = H0 / M0
            return Dict("M0" => M0, "K0" => K0)
        end
    else
        error("Model not available. Currently only `:langmuir` and `:henry` are available")
    end
end

"""
    params = fit_adsorption_isotherm(df, pressure_col_name, loading_col_name, model)

Takes in a DataFrame `df` containing adsorption isotherm data and fits an analytical model
to the data to identify its parameters of best fit, returned as a dictionary.
Available models are `:henry` and `:langmuir`

The Henry model takes the following form:
N = HP
The identified Henry coefficient is `params["H"]`.

The Langmuir model takes the following form:
N = (MKP)/(1+KP)
where N is the total adsorption, M is the maximum monolayer coverage, K is the Langmuir constant. and P is the pressure of the gas.

# Arguments

  - `df::DataFrame`: The DataFrame containing the pressure and adsorption data for the isotherm
  - `pressure_col_name::Symbol`: The header of the pressure column. Can be found with `names(df)`
  - `loading_col_name::Symbol`: The header of the loading/adsorption column. Can be found with `names(df)`
  - `model::Symbol`: The model chosen to fit to the adsorption isotherm data

# Returns

  - `params::Dict{AbstractString, Float64}`: A Dictionary with the parameters corresponding to each model along with the MSE of the fit. `:langmuir` contains "M" and "K". `:henry` contains "H".
"""
function fit_adsorption_isotherm(
    df::DataFrame,
    pressure_col_name::Symbol,
    loading_col_name::Symbol,
    model::Symbol,
    options::Optim.Options=Optim.Options()
)
    _df = sort(df, [pressure_col_name])
    n = _df[!, loading_col_name]
    p = _df[!, pressure_col_name]
    θ0 = _guess(_df, pressure_col_name, loading_col_name, model)

    if model == :langmuir
        function objective_function_langmuir(θ)
            return sum([
                (n[i] - θ[1] * θ[2] * p[i] / (1 + θ[2] * p[i]))^2 for i in eachindex(n)
            ])
        end
        res = optimize(
            objective_function_langmuir,
            [θ0["M0"], θ0["K0"]],
            NelderMead(),
            options
        )
        if !Optim.converged(res)
            error("Optimization algorithm failed!")
        end
        M, K = res.minimizer
        mse = res.minimum / length(n)
        return Dict("M" => M, "K" => K, "MSE" => mse)
    elseif model == :henry
        objective_function_henry(θ) = sum([(n[i] - θ[1] * p[i])^2 for i in eachindex(n)])
        res = optimize(objective_function_henry, [θ0["H0"]], LBFGS(), options)
        if !Optim.converged(res)
            error("Optimization algorithm failed!")
        end
        H = res.minimizer
        mse = res.minimum / length(n)
        return Dict("H" => H[1], "MSE" => mse)
    end
end
