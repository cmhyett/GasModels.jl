"variables associated with density (transient)"
function variable_density(
    gm::AbstractGasModel,
    nw::Int = gm.cnw;
    bounded::Bool = true,
    report::Bool = true,
)
    rho =
        var(gm, nw)[:density] = JuMP.@variable(
            gm.model,
            [i in ids(gm, nw, :junction)],
            base_name = "$(nw)_rho",
            start = comp_start_value(
                ref(gm, nw, :junction),
                i,
                "rho_start",
                ref(gm, nw, :junction, i)["p_min"],
            )
        )

    if bounded
        for (i, junction) in ref(gm, nw, :junction)
            JuMP.set_lower_bound(rho[i], junction["p_min"])
            JuMP.set_upper_bound(rho[i], junction["p_max"])
        end
    end

    report && _IM.sol_component_value(gm, nw, :junction, :density, ids(gm, nw, :junction), rho)
    report && _IM.sol_component_value(gm, nw, :junction, :pressure, ids(gm, nw, :junction), rho)
end

"variables associated with compressor mass flow (transient)"
function variable_compressor_flow(
    gm::AbstractGasModel,
    nw::Int = gm.cnw;
    bounded::Bool = true,
    report::Bool = true,
)
    max_mass_flow = ref(gm, nw, :max_mass_flow)
    f =
        var(gm, nw)[:compressor_flow] = JuMP.@variable(
            gm.model,
            [i in ids(gm, nw, :compressor)],
            base_name = "$(nw)_f_compressor",
            start = comp_start_value(ref(gm, nw, :compressor), i, "f_compressor_start", 0.0)
        )

    if bounded
        for (i, compressor) in ref(gm, nw, :compressor)
            JuMP.set_lower_bound(f[i], -max_mass_flow)
            JuMP.set_upper_bound(f[i], max_mass_flow)
        end
    end

    report && _IM.sol_component_value(gm, nw, :compressor, :flow, ids(gm, nw, :compressor), f)
end

"variables associated with pipe flux (transient)"
function variable_pipe_flux(
    gm::AbstractGasModel,
    nw::Int = gm.cnw;
    bounded::Bool = true,
    report::Bool = true,
)
    phi =
        var(gm, nw)[:pipe_flux] = JuMP.@variable(
            gm.model,
            [i in ids(gm, nw, :pipe)],
            base_name = "$(nw)_flux_pipe",
            start = comp_start_value(
                ref(gm, nw, :pipe),
                i,
                "phi_start",
                ref(gm, nw, :pipe, i)["flux_min"],
            )
        )

    if bounded
        for (i, pipe) in ref(gm, nw, :pipe)
            JuMP.set_lower_bound(phi[i], pipe["flux_min"])
            JuMP.set_upper_bound(phi[i], pipe["flux_max"])
        end
    end

    report && _IM.sol_component_value(gm, nw, :pipe, :flux, ids(gm, nw, :pipe), phi)
    sol_f = Dict(i => phi[i] * ref(gm, nw, :pipe, i)["area"] for i in ids(gm, nw, :pipe))
    report && _IM.sol_component_value(gm, nw, :pipe, :flow, ids(gm, nw, :pipe), sol_f)
end

"variables associated with compression ratio (transient)"
function variable_c_ratio(
    gm::AbstractGasModel,
    nw::Int = gm.cnw;
    bounded::Bool = true,
    report::Bool = true,
)
    c_ratio =
        var(gm, nw)[:compressor_ratio] = JuMP.@variable(
            gm.model,
            [i in ids(gm, nw, :compressor)],
            base_name = "$(nw)_c_ratio",
            start = comp_start_value(ref(gm, nw, :compressor), i, "c_ratio_start", 1.0)
        )

    if bounded
        for (i, compressor) in ref(gm, nw, :compressor)
            JuMP.set_lower_bound(c_ratio[i], compressor["c_ratio_min"])
            JuMP.set_upper_bound(c_ratio[i], compressor["c_ratio_max"])
        end
    end

    report &&
    _IM.sol_component_value(gm, nw, :compressor, :c_ratio, ids(gm, nw, :compressor), c_ratio)
end

