# Define LRDWP implementations of Gas Models


"Variables needed for modeling flow in LRDWP models"
function variable_flow(gm::AbstractLRDWPModel, n::Int = nw_id_default; bounded::Bool = true, report::Bool = true)
    variable_mass_flow(gm, n; bounded = bounded, report = report)
    variable_connection_direction(gm, n; report = report)
end


"Variables needed for modeling flow in LRDWP models"
function variable_flow_ne(gm::AbstractLRDWPModel, n::Int = nw_id_default; bounded::Bool = true, report::Bool = true)
    variable_mass_flow_ne(gm, n; bounded = bounded, report = report)
    variable_connection_direction_ne(gm, n; report = report)
end

######################################################################################################
## Constraints
######################################################################################################

"Constraint: Weymouth equation--not applicable for LRDWP models"
function constraint_pipe_weymouth(gm::AbstractLRDWPModel, n::Int, k, i, j, f_min, f_max, w, pd_min, pd_max)
    pi = var(gm, n, :psqr, i)
    pj = var(gm, n, :psqr, j)
    f = var(gm, n, :f_pipe, k)

    var(gm, n)[:fmodf_pipe_lifted] = Dict()
    var(gm, n, :fmodf_pipe_lifted)[k] = JuMP.@variable(gm.model, base_name = "fmodf_pipe_lifted_$(k)")
    fmodf_lifted = var(gm, n, :fmodf_pipe_lifted, k)

    _add_constraint!(gm, n, :weymouth1, k, JuMP.@constraint(gm.model, w*(pi - pj) == fmodf_lifted))

    # relaxation for fmodf
        fmodf = x->x*abs(x)
        if(f_min<=0 && f_max>=0)
            partition = [f_min,0,f_max]
        else
            partition = [f_min, f_max]
        end
        relaxation_data = add_milp_relaxation(fmodf, partition,gm.model,name="pipe_fmodf_milp_relaxation")
        vars_in_relax = relaxation_data[1];
        f_index = relaxation_data[2];
        fmodf_index = relaxation_data[3];
        _add_constraint!(gm, n, :weymouth2, k, JuMP.@constraint(gm.model, f == vars_in_relax[f_index]))
        _add_constraint!(gm, n, :weymouth3, k, JuMP.@constraint(gm.model, fmodf_lifted == vars_in_relax[fmodf_index]))
end

"Auxillary function: relaxation for pressure square at junctions incdident on resistors"
function constraint_resistor_junction_pressure_squared_relaxation(gm::AbstractLRDWPModel, n::Int, i::Int, p_min, p_max)
    p = var(gm, n, :p, i)
    psqr = var(gm, n, :psqr, i)

    var(gm, n)[:pp_lifted] = Dict()
    var(gm, n, :pp_lifted)[i] = JuMP.@variable(gm.model, base_name="pp_lifted_$(i)")
    pp_lifted = var(gm, n, :pp_lifted, i)

    _add_constraint!(gm, n, :sqrd_pressure_relaxation_1, i, JuMP.@constraint(gm.model, psqr == pp_lifted))

    # relaxation for p*p
    pp = x->x*x
    partition = [p_min, p_max]

    relaxation_data = add_milp_relaxation(pp, partition,gm.model,name="pp_milp_relaxation")
    vars_in_relax = relaxation_data[1];
    p_index = relaxation_data[2];
    pp_index = relaxation_data[3];
    _add_constraint!(gm, n, :sqrd_pressure_relaxation_2, k, JuMP.@constraint(gm.model, p == vars_in_relax[p_index]))
    _add_constraint!(gm, n, :sqrd_pressure_relaxation_3, k, JuMP.@constraint(gm.model, pp_lifted == vars_in_relax[pp_index]))

end

"Constraint: Darcy-Weisbach equation--not applicable for LRDWP models"
function constraint_resistor_darcy_weisbach(gm::AbstractLRDWPModel, n::Int, k, i, j, f_min, f_max, w, pd_min, pd_max)
    p_i = var(gm, n, :p, i)
    p_j = var(gm, n, :p, j)
    f = var(gm, n, :f_resistor, k)

    var(gm, n)[:fmodf_resistor_lifted] = Dict()
    var(gm, n, :fmodf_resistor_lifted)[k] = JuMP.@variable(gm.model, base_name="fmodf__resistor_lifted_$(k)")
    fmodf_lifted = var(gm, n, :fmodf__resistor_lifted, k)

    _add_constraint!(gm, n, :darcy_weisbach_1, k, JuMP.@constraint(gm.model, (1.0/w)* (p_i - p_j) == fmodf_lifted))

    #relaxation for fmodf
        fmodf = x->x*(abs(x))
        if(f_min<=0 && f_max>=0)
            partition = [f_min,0,f_max]
        else
            partition = [f_min, f_max]
        end

        relaxation_data = add_milp_relaxation(fmodf, partition,gm.model,name="resistor_fmodf_milp_relaxation")
        vars_in_relax = relaxation_data[1];
        f_index = relaxation_data[2];
        fmodf_index = relaxation_data[3];
        _add_constraint!(gm, n, :darcy_weisbach_2, k, JuMP.@constraint(gm.model, f == vars_in_relax[f_index]))
        _add_constraint!(gm, n, :darcy_weisbach_3, k, JuMP.@constraint(gm.model, fmodf_lifted == vars_in_relax[fmodf_index]))

        constraint_resistor_junction_pressure_squared_relaxation(gm,n,i,ref(gm, n, :junction, i)["p_min"],ref(gm, n, :junction, i)["p_max"])
        constraint_resistor_junction_pressure_squared_relaxation(gm,n,j,ref(gm, n, :junction, j)["p_min"],ref(gm, n, :junction, j)["p_max"])

end


"Constraint: Define pressures across a resistor"
function constraint_resistor_pressure(gm::AbstractLRDWPModel, n::Int, k::Int, i::Int, j::Int, pd_min::Float64, pd_max::Float64)
end


"Constraint: Constraints which define pressure drop across a loss resistor"
function constraint_loss_resistor_pressure(gm::AbstractLRDWPModel, n::Int, k::Int, i::Int, j::Int, pd::Float64)
    p_i = var(gm, n, :p, i)
    p_j = var(gm, n, :p, j)
    f = var(gm, n, :f_loss_resistor, k)

    var(gm, n)[:y_loss_resistor] = Dict()
    var(gm, n, :y_loss_resistor)[k] = JuMP.@variable(gm.model, binary=true, base_name = "y_loss_resistor_$(k)")
    y = var(gm, n, :y_loss_resistor, k)

    constraint_resistor_junction_pressure_squared_relaxation(gm,n,i,ref(gm, n, :junction, i)["p_min"],ref(gm, n, :junction, i)["p_max"])
    constraint_resistor_junction_pressure_squared_relaxation(gm,n,j,ref(gm, n, :junction, j)["p_min"],ref(gm, n, :junction, j)["p_max"])

    _add_constraint!(gm, n, :loss_resistor_1, k, JuMP.@constraint(gm.model, p_i -p_j == (2*y -1)*pd))
    f_min = ref(gm, n, :loss_resistor, k)["f_min"]
    f_max = ref(gm, n, :loss_resistor, k)["f_max"]
    _add_constraint!(gm, n, :loss_resistor_2, k, JuMP.@constraint(gm.model, (1-y)*f_min <= f))
    _add_constraint!(gm, n, :loss_resistor_3, k, JuMP.@constraint(gm.model, f<= y*f_max))
end


"Constraint: Weymouth equation"
function constraint_pipe_weymouth_ne(gm::AbstractLRDWPModel, n::Int, k, i, j, w, f_min, f_max, pd_min, pd_max)
    #TODO Linear convex hull of the weymouth equations in crdwp.jl
end


"Constraint: constrains the ratio to be ``p_i \\cdot \\alpha = p_j``"
function constraint_compressor_ratio_value(gm::AbstractLRDWPModel, n::Int, k, i, j, type, i_pmax, j_pmax, max_ratio)
    pi = var(gm, n, :psqr, i)
    pj = var(gm, n, :psqr, j)
    r = var(gm, n, :rsqr, k)
    y = var(gm, n, :y_compressor, k)

    if type == 0
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
function constraint_compressor_ratio_value_ne(gm::AbstractLRDWPModel, n::Int, k, i, j, type, i_pmax, j_pmax, max_ratio)
    pi = var(gm, n, :psqr, i)
    pj = var(gm, n, :psqr, j)
    r = var(gm, n, :rsqr_ne, k)
    y = var(gm, n, :y_ne_compressor, k)
    z = var(gm, n, :zc, k)

    if type == 0
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
function constraint_compressor_energy(gm::AbstractLRDWPModel, n::Int, k, power_max, m, work)
    r = var(gm, n, :rsqr, k)
    f = var(gm, n, :f_compressor, k)


    var(gm, n)[:r_exp_comp_lifted] = Dict()
    var(gm, n, :r_exp_comp_lifted)[k] = JuMP.@variable(gm.model, base_name="r_exp_comp_lifted_$(k)")
    r_exp_lifted = var(gm, n, :r_exp_comp_lifted, k)

    var(gm, n)[:fr_comp_lifted] = Dict()
    var(gm, n, :fr_comp_lifted)[k] = JuMP.@variable(gm.model, base_name="fr_comp_lifted_$(k)")
    fr_lifted = var(gm, n, :fr_comp_lifted, k)

    min_ratio = ref(gm, n, :compressor, k)["c_ratio_min"]
    max_ratio = ref(gm, n, :compressor, k)["c_ratio_max"]
    f_min = ref(gm, n, :compressor, k)["flow_min"]
    f_max = ref(gm, n, :compressor, k)["flow_max"]

    # polyhedral relaxation for r_exp = r^(m/2)-1
        r_exp = x->x^(m/2)-1
        partition = [min_ratio,max_ratio]
        relaxation_data = add_milp_relaxation(r_exp, partition,gm.model,name="compressor_r_exp_milp_relaxation")
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
function constraint_compressor_energy_ne(gm::AbstractLRDWPModel, n::Int, k, power_max, m, work)
    #TODO - lrdwp relaxation
end
