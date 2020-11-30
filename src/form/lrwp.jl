# Define LRWP implementations of Gas Models

######################################################################################################
## Constraints
######################################################################################################
# new functions:
    # variables
        # fmodf_pipe_lifted

"Constraint: Weymouth equation--not applicable for LRWP models"
function constraint_pipe_weymouth(gm::AbstractLRWPModel, n::Int, k, i, j, f_min, f_max, w, pd_min, pd_max)
    pi = var(gm, n, :psqr, i)
    pj = var(gm, n, :psqr, j)
    f = var(gm, n, :f_pipe, k)

    # variable_pipe_fmod_f()
    # if get(var(gm, n, :fmodf_pipe_lifted),k) == false
        var(gm, n)[:fmodf_pipe_lifted] = Dict()
    # end
    var(gm, n, :fmodf_pipe_lifted)[k] = JuMP.@variable(gm.model, base_name="fmodf_pipe_lifted_$(k)")
    fmodf_lifted = var(gm, n, :fmodf_pipe_lifted, k)

    _add_constraint!(gm, n, :weymouth1, k, JuMP.@constraint(gm.model, w * (pi - pj) == fmodf_lifted))


    #relaxation for fmodf
        fmodf = x->x*(abs(x))
        if(f_min<=0 && f_max>=0)
            partition = [f_min,0,f_max]
        else
            partition = [f_min, f_max]
        end
        relaxation_data = add_lp_relaxation(fmodf, partition,gm.model,name="pipe_fmodf_lp_relaxation")
        vars_in_relax = relaxation_data[1];
        f_index = relaxation_data[2];
        fmodf_index = relaxation_data[3];
        _add_constraint!(gm, n, :weymouth2, k, JuMP.@constraint(gm.model, f == vars_in_relax[f_index]))
        _add_constraint!(gm, n, :weymouth3, k, JuMP.@constraint(gm.model, f_modf_lifted == vars_in_relax[fmodf_index]))

end

"Auxillary function: relaxation for pressure square at junctions incdident on resistors"
function constraint_resistor_junction_pressure_squared_relaxation(gm::AbstractLRWPModel, n::Int, i::Int, p_min, p_max)
    p = var(gm, n, :p, i)
    psqr = var(gm, n, :psqr, i)
    if get(var(gm, n), :pp_lifted) == false
        var(gm, n)[:pp_lifted] = Dict()
    end
    var(gm, n, :pp_lifted)[i] = JuMP.@variable(gm.model, base_name="pp_lifted_$(i)")
    pp_lifted = var(gm, n, :pp_lifted, i)

    _add_constraint!(gm, n, :sqrd_pressure_relaxation_1, i, JuMP.@constraint(gm.model, psqr == pp_lifted))

    # relaxation for p*p
    pp = x->x*x
    partition = [p_min, p_max]

    relaxation_data = add_lp_relaxation(pp, partition,gm.model,name="pp_lp_relaxation")
    vars_in_relax = relaxation_data[1];
    p_index = relaxation_data[2];
    pp_index = relaxation_data[3];
    _add_constraint!(gm, n, :sqrd_pressure_relaxation_2, k, JuMP.@constraint(gm.model, p == vars_in_relax[p_index]))
    _add_constraint!(gm, n, :sqrd_pressure_relaxation_3, k, JuMP.@constraint(gm.model, pp_lifted == vars_in_relax[pp_index]))

end

"Constraint: Darcy-Weisbach equation--not applicable for LRWP models"
function constraint_resistor_darcy_weisbach(gm::AbstractLRWPModel, n::Int, k, i, j, f_min, f_max, w, pd_min, pd_max)
    p_i = var(gm, n, :p, i)
    p_j = var(gm, n, :p, j)
    f = var(gm, n, :f_resistor, k)
    # variable_resistor_y_linear()
    if get(var(gm, n), :y_resistor_linear) == false
        var(gm, n)[:y_resistor_linear] = Dict()
    end
    var(gm, n, :y_resistor_linear)[k] = JuMP.@variable(gm.model, base_name="y_resistor_linear_$(k)")
    y_linear = var(gm, n, :y_resistor_linear, k)

    # variable_resistor_fmod_f()
    if get(var(gm, n), :fmodf_resistor_lifted) == false
        var(gm, n)[:fmodf_resistor_lifted] = Dict()
    end
    var(gm, n, :fmodf_resistor_lifted)[k] = JuMP.@variable(gm.model, base_name="fmodf__resistor_lifted_$(k)")
    fmodf_lifted = var(gm, n, :fmodf__resistor_lifted, k)

    _add_constraint!(gm, n, :darcy_weisbach_1, k, JuMP.@constraint(gm.model, (1.0/w)* (p_i - p_j) == f_modf_lifted))

    #relaxation for fmodf
        fmodf = x->x*(abs(x))
        if(f_min<=0 && f_max>=0)
            partition = [f_min,0,f_max]
        else
            partition = [f_min, f_max]
        end

        relaxation_data = add_lp_relaxation(fmodf, partition,gm.model,name="resistor_fmodf_lp_relaxation")
        vars_in_relax = relaxation_data[1];
        f_index = relaxation_data[2];
        fmodf_index = relaxation_data[3];
        _add_constraint!(gm, n, :darcy_weisbach_2, k, JuMP.@constraint(gm.model, f == vars_in_relax[f_index]))
        _add_constraint!(gm, n, :darcy_weisbach_3, k, JuMP.@constraint(gm.model, f_modf_lifted == vars_in_relax[fmodf_index]))

        constraint_resistor_junction_pressure_squared_relaxation(gm,n,i,ref(gm, n, :junction, i)["p_min"],ref(gm, n, :junction, i)["p_max"])
        constraint_resistor_junction_pressure_squared_relaxation(gm,n,j,ref(gm, n, :junction, j)["p_min"],ref(gm, n, :junction, j)["p_max"])
end


"Constraint: Define pressures across a resistor"
function constraint_resistor_pressure(gm::AbstractLRWPModel, n::Int, k::Int, i::Int, j::Int, pd_min::Float64, pd_max::Float64)
end


"Constraint: Constraints which define pressure drop across a loss resistor"
function constraint_loss_resistor_pressure(gm::AbstractLRWPModel, n::Int, k::Int, i::Int, j::Int, pd::Float64)
    p_i = var(gm, n, :p, i)
    p_j = var(gm, n, :p, j)
    f = var(gm, n, :f_loss_resistor, k)
    # variable_loss_resistor_y_linear()
    if get(var(gm, n), :y_loss_resistor_linear) == false
        var(gm, n)[:y_loss_resistor_linear] = Dict()
    end
    var(gm, n, :y_loss_resistor_linear)[k] = JuMP.@variable(gm.model, base_name="y_loss_resistor_linear_$(k)")
    y_linear = var(gm, n, :y_loss_resistor_linear, k)

    constraint_resistor_junction_pressure_squared_relaxation(gm,n,i,ref(gm, n, :junction, i)["p_min"],ref(gm, n, :junction, i)["p_max"])
    constraint_resistor_junction_pressure_squared_relaxation(gm,n,j,ref(gm, n, :junction, j)["p_min"],ref(gm, n, :junction, j)["p_max"])

    _add_constraint!(gm, n, :loss_resistor_1, k, JuMP.@constraint(gm.model, p_i -p_j == (2*y_linear -1)*pd))

    f_min = ref(gm, n, :loss_resistor, k)["f_min"]
    f_max = ref(gm, n, :loss_resistor, k)["f_max"]
    _add_constraint!(gm, n, :loss_resistor_2, k, JuMP.@constraint(gm.model, (1-y_linear)*f_min <= f))
    _add_constraint!(gm, n, :loss_resistor_3, k, JuMP.@constraint(gm.model, f<= y_linear*f_max))
end


"Constraint: Compressor ratio constraints on pressure differentials--not applicable for LRWP models"
function constraint_compressor_ratios(gm::AbstractLRWPModel, n::Int, k, i, j, min_ratio, max_ratio, i_pmin, i_pmax, j_pmin, j_pmax, type)
    pi = var(gm, n, :psqr, i)
    pj = var(gm, n, :psqr, j)
    f = var(gm, n, :f_compressor, k)

    # compression in both directions
    if type == 0
        if (min_ratio <= 1.0 && max_ratio >= 1)
            pk = JuMP.@variable(
                gm.model,
                upper_bound = max(JuMP.upper_bound(pi), JuMP.upper_bound(pj)),
                lower_bound = min(JuMP.lower_bound(pi), JuMP.lower_bound(pj))
            )
            pik = JuMP.@variable(
                gm.model,
                upper_bound = JuMP.upper_bound(pi) - JuMP.lower_bound(pk),
                lower_bound = JuMP.lower_bound(pi) - JuMP.upper_bound(pk)
            )
            pjk = JuMP.@variable(
                gm.model,
                upper_bound = JuMP.upper_bound(pj) - JuMP.lower_bound(pk),
                lower_bound = JuMP.lower_bound(pj) - JuMP.upper_bound(pk)
            )
            fpik = JuMP.@variable(gm.model)
            fpjk = JuMP.@variable(gm.model)

            _add_constraint!(gm, n, :compressor_ratios1, k, JuMP.@constraint(gm.model, pk - max_ratio^2 * pi <= 0))
            _add_constraint!(gm, n, :compressor_ratios2, k, JuMP.@constraint(gm.model, min_ratio^2 * pi - pk <= 0))
            _add_constraint!(gm, n, :compressor_ratios3, k, JuMP.@constraint(gm.model, fpik <= 0))
            _add_constraint!(gm, n, :compressor_ratios4, k, JuMP.@constraint(gm.model, pk - max_ratio^2 * pj <= 0))
            _add_constraint!(gm, n, :compressor_ratios5, k, JuMP.@constraint(gm.model, min_ratio^2 * pj - pk <= 0))
            _add_constraint!(gm, n, :compressor_ratios6, k, JuMP.@constraint(gm.model, -fpjk <= 0))
            _add_constraint!(gm, n, :compressor_ratios7, k, JuMP.@constraint(gm.model, pjk == pj - pk))
            _add_constraint!(gm, n, :compressor_ratios8, k, JuMP.@constraint(gm.model, pik == pi - pk))
            _IM.relaxation_product(gm.model, f, pik, fpik)
            _IM.relaxation_product(gm.model, f, pjk, fpjk)

            # There is a disjunction, so we have to use a binary variable for this one
        else
            y = get_compressor_y(gm, n, k)
            _add_constraint!(gm, n, :on_off_compressor_ratios1, k, JuMP.@constraint(gm.model, pj - max_ratio^2 * pi <= (1 - y) * (j_pmax^2)))
            _add_constraint!(gm, n, :on_off_compressor_ratios2, k, JuMP.@constraint(gm.model, min_ratio^2 * pi - pj <= (1 - y) * (min_ratio^2 * i_pmax^2)))
            _add_constraint!(gm, n, :on_off_compressor_ratios3, k, JuMP.@constraint(gm.model, pi - max_ratio^2 * pj <= y * (i_pmax^2)))
            _add_constraint!(gm, n, :on_off_compressor_ratios4, k, JuMP.@constraint(gm.model, min_ratio^2 * pj - pi <= y * (min_ratio^2 * j_pmax^2)))
        end
        # compression when flow is from i to j.  No flow in reverse, so nothing to model in that direction
    elseif type == 1
        _add_constraint!(gm, n, :compressor_ratios1, k, JuMP.@constraint(gm.model, pj - max_ratio^2 * pi <= 0))
        _add_constraint!(gm, n, :compressor_ratios2, k, JuMP.@constraint(gm.model, min_ratio^2 * pi - pj <= 0))
        # compression when flow is from i to j.  no compression when flow is from j to i. min_ratio = 1
    else # type 2
        if min_ratio == 1
            pij = JuMP.@variable(
                gm.model,
                upper_bound = JuMP.upper_bound(pi) - JuMP.lower_bound(pj),
                lower_bound = JuMP.lower_bound(pi) - JuMP.upper_bound(pj)
            )
            fpij = JuMP.@variable(gm.model)

            _add_constraint!(gm, n, :compressor_ratios1, k, JuMP.@constraint(gm.model, pj - max_ratio^2 * pi <= 0))
            _add_constraint!(gm, n, :compressor_ratios2, k, JuMP.@constraint(gm.model, min_ratio^2 * pi - pj <= 0))
            _add_constraint!(gm, n, :compressor_ratios3, k, JuMP.@constraint(gm.model, fpij <= 0))
            _add_constraint!(gm, n, :compressor_ratios4, k, JuMP.@constraint(gm.model, pij == pi - pj))
            _IM.relaxation_product(gm.model, f, pij, fpij)
            # compression when flow is from i to j.  no compression when flow is from j to i. min_ratio != 1. This is a disjunctive model
        else
            y = get_compressor_y(gm, n, k)
            _add_constraint!(gm, n, :on_off_compressor_ratios1, k, JuMP.@constraint(gm.model, pj - max_ratio^2 * pi <= (1 - y) * (j_pmax^2)))
            _add_constraint!(gm, n, :on_off_compressor_ratios2, k, JuMP.@constraint(gm.model, min_ratio^2 * pi - pj <= (1 - y) * (i_pmax^2)))
            _add_constraint!(gm, n, :on_off_compressor_ratios3, k, JuMP.@constraint(gm.model, pi - pj <= y * (i_pmax^2)))
            _add_constraint!(gm, n, :on_off_compressor_ratios4, k, JuMP.@constraint(gm.model, pj - pi <= y * (j_pmax^2)))
        end
    end
end


"constraints on pressure drop across control valves--not applicable for LRWP models"
function constraint_on_off_regulator_pressure(gm::AbstractLRWPModel, n::Int, k, i, j, min_ratio, max_ratio, f_min, i_pmin, i_pmax, j_pmin, j_pmax)
    pi = var(gm, n, :psqr, i)
    pj = var(gm, n, :psqr, j)
    v = var(gm, n, :v_regulator, k)
    f = var(gm, n, :f_regulator, k)

    if max_ratio == 1
        pij = JuMP.@variable(
            gm.model,
            upper_bound = JuMP.upper_bound(pi) - JuMP.lower_bound(pj),
            lower_bound = JuMP.lower_bound(pi) - JuMP.upper_bound(pj)
        )
        fpij = JuMP.@variable(gm.model)

        _add_constraint!(gm, n, :regulator_pressure_drop1, k, JuMP.@constraint(gm.model, pj - max_ratio^2 * pi <= (1 - v) * j_pmax^2))
        _add_constraint!(gm, n, :regulator_pressure_drop2, k, JuMP.@constraint(gm.model, min_ratio^2 * pi - pj <= (1 - v) * (min_ratio * i_pmax^2)))
        _add_constraint!(gm, n, :regulator_pressure_drop3, k, JuMP.@constraint(gm.model, fpij >= 0))
        _add_constraint!(gm, n, :regulator_pressure_drop4, k, JuMP.@constraint(gm.model, pij == pi - pj))
        _IM.relaxation_product(gm.model, f, pij, fpij)
    elseif f_min >= 0
        _add_constraint!(gm, n, :regulator_pressure_drop1, k, JuMP.@constraint(gm.model, pj - max_ratio^2 * pi <= (1 - v) * j_pmax^2))
        _add_constraint!(gm, n, :regulator_pressure_drop2, k, JuMP.@constraint(gm.model, min_ratio^2 * pi - pj <= (1 - v) * (min_ratio * i_pmax^2)))
    else
        # There condition implies a disjunction when flow is reversed
        y = JuMP.@variable(gm.model, binary = true)
        _add_constraint!(gm, n, :regulator_pressure_drop1, k, JuMP.@constraint(gm.model, pj - (max_ratio^2 * pi) <= (2 - y - v) * j_pmax^2))
        _add_constraint!(gm, n, :regulator_pressure_drop2, k, JuMP.@constraint(gm.model, (min_ratio^2 * pi) - pj <= (2 - y - v) * i_pmax^2))
        _add_constraint!(gm, n, :regulator_pressure_drop3, k, JuMP.@constraint(gm.model, pj - pi <= (1 + y - v) * j_pmax^2))
        _add_constraint!(gm, n, :regulator_pressure_drop4, k, JuMP.@constraint(gm.model, pi - pj <= (1 + y - v) * i_pmax^2))
    end
end


"Constraint: Weymouth equation"
function constraint_pipe_weymouth_ne(gm::AbstractLRWPModel, n::Int, k, i, j, w, f_min, f_max, pd_min, pd_max)
    #TODO Linear convex hull equations in wp.jl
end


"Constraint: compressor ratios on a new compressor"
function constraint_compressor_ratios_ne(gm::AbstractLRWPModel, n::Int, k, i, j, min_ratio, max_ratio, i_pmin, i_pmax, j_pmin, j_pmax, type)
    pi = var(gm, n, :psqr, i)
    pj = var(gm, n, :psqr, j)
    zc = var(gm, n, :zc, k)
    f = var(gm, n, :f_ne_compressor, k)

    MR = max(i_pmax^2, j_pmax^2) / min(i_pmin^2, j_pmin^2)

    # compression in both directions
    if type == 0
        if (min_ratio <= 1.0 && max_ratio >= 1)
            k_pmax = max(i_pmax, j_pmax)
            k_pmin = min(i_pmin, j_pmin)
            pk = JuMP.@variable(gm.model, upper_bound = k_pmax, lower_bound = k_pmin)
            pik = JuMP.@variable(
                gm.model,
                upper_bound = JuMP.upper_bound(pi) - JuMP.lower_bound(pk),
                lower_bound = JuMP.lower_bound(pi) - JuMP.upper_bound(pk)
            )
            pjk = JuMP.@variable(
                gm.model,
                upper_bound = JuMP.upper_bound(pj) - JuMP.lower_bound(pk),
                lower_bound = JuMP.lower_bound(pj) - JuMP.upper_bound(pk)
            )
            fpik = JuMP.@variable(gm.model)
            fpjk = JuMP.@variable(gm.model)

            _add_constraint!(gm, n, :compressor_ratios_ne1, k, JuMP.@constraint(gm.model, pk - max_ratio^2 * pi <= (1 - zc) * k_pmax^2))
            _add_constraint!(gm, n, :compressor_ratios_ne2, k, JuMP.@constraint(gm.model, min_ratio^2 * pi - pk <= (1 - zc) * i_pmax^2 * min_ratio^2))
            _add_constraint!(gm, n, :compressor_ratios_ne3, k, JuMP.@constraint(gm.model, fpik <= 0)) # zc = 0 implies f = 0 so always true then
            _add_constraint!(gm, n, :compressor_ratios_ne4, k, JuMP.@constraint(gm.model, pk - max_ratio^2 * pj <= (1 - zc) * k_pmax^2))
            _add_constraint!(gm, n, :compressor_ratios_ne5, k, JuMP.@constraint(gm.model, min_ratio^2 * pj - pk <= (1 - zc) * j_pmax^2 * min_ratio^2))
            _add_constraint!(gm, n, :compressor_ratios_ne6, k, JuMP.@constraint(gm.model, -fpjk <= 0))  # zc = 0 implies f = 0 so always true then
            _add_constraint!(gm, n, :compressor_ratios_ne7, k, JuMP.@constraint(gm.model, pjk == pj - pk))
            _add_constraint!(gm, n, :compressor_ratios_ne8, k, JuMP.@constraint(gm.model, pik == pi - pk))
            _IM.relaxation_product(gm.model, f, pik, fpik)
            _IM.relaxation_product(gm.model, f, pjk, fpjk)
            # There is a disjunction, so we have to use a binary variable for this one
        else
            y = get_ne_compressor_y_lrwp(gm, n, k)
            _add_constraint!(gm, n, :on_off_compressor_ratios_ne1, k, JuMP.@constraint(gm.model, pj - (max_ratio^2 * pi) <= (2 - y - zc) * j_pmax^2))
            _add_constraint!(gm, n, :on_off_compressor_ratios_ne2, k, JuMP.@constraint(gm.model, (min_ratio^2 * pi) - pj <= (2 - y - zc) * (min_ratio^2 * i_pmax^2)))
            _add_constraint!(gm, n, :on_off_compressor_ratios_ne3, k, JuMP.@constraint(gm.model, pi - (max_ratio^2 * pj) <= (1 + y - zc) * i_pmax^2))
            _add_constraint!(gm, n, :on_off_compressor_ratios_ne4, k, JuMP.@constraint(gm.model, (min_ratio^2 * pj) - pi <= (1 + y - zc) * (min_ratio^2 * j_pmax^2)))
        end
        # compression when flow is from i to j.  No flow in reverse, so nothing to model in that direction
    elseif type == 1
        _add_constraint!(gm, n, :on_off_compressor_ratios_ne1, k, JuMP.@constraint(gm.model, pj - (max_ratio^2 * pi) <= (1 - zc) * j_pmax^2))
        _add_constraint!(gm, n, :on_off_compressor_ratios_ne2, k, JuMP.@constraint(gm.model, (min_ratio^2 * pi) - pj <= (1 - zc) * i_pmax^2 * min_ratio^2))
        # compression when flow is from i to j.  no compression when flow is from j to i. min_ratio = 1
    else # type 2
        if min_ratio == 1
            pij = JuMP.@variable(
                gm.model,
                upper_bound = JuMP.upper_bound(pi) - JuMP.lower_bound(pj),
                lower_bound = JuMP.lower_bound(pi) - JuMP.upper_bound(pj)
            )
            fpij = JuMP.@variable(gm.model)

            _add_constraint!(gm, n, :compressor_ratios_ne1, k, JuMP.@constraint(gm.model, pj - max_ratio^2 * pi <= (1 - zc) * j_pmax^2))
            _add_constraint!(gm, n, :compressor_ratios_ne2, k, JuMP.@constraint(gm.model, min_ratio^2 * pi - pj <= (1 - zc) * (min_ratio * i_pmax^2)))
            # z_c = 0 implies f = 0 (constraint_compressor_ne), so 0 <= 1. So constraint is off
            # z_c = 1 implies f * (pi - pj) <= 0, which is the constraint we want when the edge is actve
            _add_constraint!(gm, n, :compressor_ratios_ne3, k, JuMP.@constraint(gm.model, fpij <= (1 - zc)))
            _add_constraint!(gm, n, :compressor_ratios_ne4, k, JuMP.@constraint(gm.model, pij == pi - pj))
            _IM.relaxation_product(gm.model, f, pij, fpij)
            # compression when flow is from i to j.  no compression when flow is from j to i. min_ratio != 1. This is a disjunctive model
        else
            y = get_ne_compressor_y_lrwp(gm, n, k)
            _add_constraint!(gm, n, :on_off_compressor_ratios_ne1, k, JuMP.@constraint(gm.model, pj - (max_ratio^2 * pi) <= (2 - y - zc) * j_pmax^2))
            _add_constraint!(gm, n, :on_off_compressor_ratios_ne2, k, JuMP.@constraint(gm.model, (min_ratio^2 * pi) - pj <= (2 - y - zc) * (min_ratio^2 * i_pmax^2)))
            _add_constraint!(gm, n, :on_off_compressor_ratios3, k, JuMP.@constraint(gm.model, pi - pj <= (1 + y - zc) * (i_pmax^2)))
            _add_constraint!(gm, n, :on_off_compressor_ratios4, k, JuMP.@constraint(gm.model, pj - pi <= (1 + y - zc) * (j_pmax^2)))
        end
    end
end

"Constraint: constrains the ratio to be ``p_i \\cdot \\alpha = p_j``"
function constraint_compressor_ratio_value(gm::AbstractLRWPModel, n::Int, k, i, j, type, i_pmax, j_pmax, max_ratio)
    pi = var(gm, n, :psqr, i)
    pj = var(gm, n, :psqr, j)
    r = var(gm, n, :rsqr, k)

    if type == 0
        y = get_compressor_y(gm, n, k)
        rpi = JuMP.@variable(gm.model)
        rpj = JuMP.@variable(gm.model)

        _IM.relaxation_product(gm.model, pi, r, rpi)
        _IM.relaxation_product(gm.model, pj, r, rpj)
        _add_constraint!(gm, n, :compressor_ratio_value1, k, JuMP.@constraint(gm.model, rpi <= pj + (1 - y) * i_pmax^2 * max_ratio^2))
        _add_constraint!(gm, n, :compressor_ratio_value2, k, JuMP.@constraint(gm.model, rpi >= pj - (1 - y) * j_pmax^2))
        _add_constraint!(gm, n, :compressor_ratio_value3, k, JuMP.@constraint(gm.model, rpj <= pi + y * j_pmax^2 * max_ratio^2))
        _add_constraint!(gm, n, :compressor_ratio_value4, k, JuMP.@constraint(gm.model, rpj >= pi - y * i_pmax^2))
    else
        _IM.relaxation_product(gm.model, pi, r, pj)
    end
end


"Constraint: constrains the ratio to be ``p_i \\cdot \\alpha = p_j``"
function constraint_compressor_ratio_value_ne(gm::AbstractLRWPModel, n::Int, k, i, j, type, i_pmax, j_pmax, max_ratio)
    pi = var(gm, n, :psqr, i)
    pj = var(gm, n, :psqr, j)
    r = var(gm, n, :rsqr_ne, k)
    z = var(gm, n, :zc, k)

    if type == 0
        y = get_ne_compressor_y_wp(gm, n, k)
        rpi = JuMP.@variable(gm.model)
        rpj = JuMP.@variable(gm.model)

        _IM.relaxation_product(gm.model, pi, r, rpi)
        _IM.relaxation_product(gm.model, pj, r, rpj)
        _add_constraint!(gm, n, :compressor_ratio_value_ne1, k, JuMP.@constraint(gm.model, rpi <= pj + (2 - y - z) * i_pmax^2 * max_ratio^2))
        _add_constraint!(gm, n, :compressor_ratio_value_ne2, k, JuMP.@constraint(gm.model, rpi >= pj - (2 - y - z) * j_pmax^2))
        _add_constraint!(gm, n, :compressor_ratio_value_ne3, k, JuMP.@constraint(gm.model, rpj <= pi + (1 + y - z) * j_pmax^2 * max_ratio^2))
        _add_constraint!(gm, n, :compressor_ratio_value_ne4, k, JuMP.@constraint(gm.model, rpj >= pi - (1 + y - z) * i_pmax^2))
    else
        rpi = JuMP.@variable(gm.model)
        _IM.relaxation_product(gm.model, pi, r, rpi)
        _add_constraint!(gm, n, :compressor_ratio_value_ne1, k, JuMP.@constraint(gm.model, rpi <= pj + (1 - z) * i_pmax^2 * max_ratio^2))
        _add_constraint!(gm, n, :compressor_ratio_value_ne2, k, JuMP.@constraint(gm.model, rpi >= pj - (1 - z) * j_pmax^2))
    end
end


"Constraint: constrains the energy of the compressor"
function constraint_compressor_energy(gm::AbstractLRWPModel, n::Int, k, power_max, m, work)
    #TODO Linear convex hull equations in wp.jl
    #required variables
    r = var(gm, n, :rsqr, k)
    f = var(gm, n, :f_compressor, k)

    min_ratio = ref(gm, n, :compressor, k)["c_ratio_min"]
    max_ratio = ref(gm, n, :compressor, k)["c_ratio_max"]
    f_min = ref(gm, n, :compressor, k)["flow_min"]
    f_max = ref(gm, n, :compressor, k)["flow_max"]

    # variable_r_exp_lifted()
    if get(var(gm, n), :r_exp_comp_lifted) == false
        var(gm, n)[:r_exp_comp_lifted] = Dict()
    end
    var(gm, n, :r_exp_comp_lifted)[k] = JuMP.@variable(gm.model, base_name="r_exp_comp_lifted_$(k)")
    r_exp_lifted = var(gm, n, :r_exp_comp_lifted, k)


    # polyhedral relaxation for r_exp = r^(m/2)-1
        r_exp = x->x^(m/2)-1
        partition = [min_ratio,max_ratio]
        relaxation_data = add_lp_relaxation(r_exp, partition,gm.model,name="compressor_r_exp_lp_relaxation")
        vars_in_relax = relaxation_data[1];
        r_index = relaxation_data[2];
        r_exp_index = relaxation_data[3];
        _add_constraint!(gm, n, :comp_energy_1, k, JuMP.@constraint(gm.model, r == vars_in_relax[r_index]))
        _add_constraint!(gm, n, :comp_energy_2, k, JuMP.@constraint(gm.model, r_exp_lifted == vars_in_relax[r_exp_index]))

        r_exp_l = min_ratio^(m/2)-1;
        r_exp_u = max_ratio^(m/2)-1;

        _add_constraint!(gm, n, :comp_energy_3, k, JuMP.@constraint(gm.model, fr_lifted <= power_max / work))
        _add_constraint!(gm, n, :comp_energy_3_mc1, k, JuMP.@constraint(gm.model, fr_lifted - f*r_exp_l - f_min*r_exp_lifted + f_min*r_exp_l >= 0))
        _add_constraint!(gm, n, :comp_energy_3_mc2, k, JuMP.@constraint(gm.model, fr_lifted - f*r_exp_l - f_max*r_exp_lifted + f_max*r_exp_l <= 0))
        _add_constraint!(gm, n, :comp_energy_3_mc3, k, JuMP.@constraint(gm.model, fr_lifted - f*r_exp_u - f_max*r_exp_lifted + f_max*r_exp_u >= 0))
        _add_constraint!(gm, n, :comp_energy_3_mc4, k, JuMP.@constraint(gm.model, fr_lifted - f*r_exp_u - f_min*r_exp_lifted + f_min*r_exp_u <= 0))
end


"Constraint: constrains the energy of the compressor"
function constraint_compressor_energy_ne(gm::AbstractLRWPModel, n::Int, k, power_max, m, work)
    #TODO Linear convex hull equations in wp.jl
end
