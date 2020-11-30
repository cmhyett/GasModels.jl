function add_lp_relaxation(
    f::Function,
    partition::Vector{Float64},
    # error_tolerance::Float64,
    m::JuMP.AbstractModel;
    name = "1",
)
    relaxation, _ = construct_lp_relaxation(
        f,
        partition,
        length_tolerance = 1e-5,
        derivative_tolerance = 1e-5,
        error_tolerance = 1e-5,
    )
    num_variables = get_num_variables(relaxation)
    lb, ub = get_variable_bounds(relaxation)
    names = get_variable_names(relaxation)
    x = JuMP.@variable(
        m,
        [i = 1:num_variables],
        lower_bound = lb[i],
        upper_bound = ub[i],
        base_name = names[i] * "_" * name
    )
    A, b = get_eq_constraint_matrices(relaxation)
    JuMP.@constraint(m, A * x .== b)
    return x, relaxation.x_index, relaxation.y_index
end


function add_milp_relaxation(
    f::Function,
    partition::Vector{Float64},
    # error_tolerance::Float64,
    m::JuMP.AbstractModel;
    name = "1",
)
    relaxation, _ = construct_milp_relaxation(
        f,
        partition,
        length_tolerance = 1e-5,
        derivative_tolerance = 1e-5,
        error_tolerance = 1e-5,
    )
    num_variables = get_num_variables(relaxation)
    lb, ub = get_variable_bounds(relaxation)
    names = get_variable_names(relaxation)
    binary = relaxation.binary
    x = JuMP.@variable(
        m,
        [i = 1:num_variables],
        lower_bound = lb[i],
        upper_bound = ub[i],
        binary = Bool(binary[i]),
        base_name = names[i] * "_" * name
    )
    A, b = get_eq_constraint_matrices(relaxation)
    JuMP.@constraint(m, A * x .== b)
    A, b = get_leq_constraint_matrices(relaxation)
    JuMP.@constraint(m, A * x .<= b)
    return x, relaxation.x_index, relaxation.y_index
end
