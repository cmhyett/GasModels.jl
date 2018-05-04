function check_pressure_status(sol, gm)
    for (idx,val) in sol["junction"]
        @test val["p"] <= gm.ref[:nw][gm.cnw][:junction][parse(Int64,idx)]["pmax"]
        @test val["p"] >= gm.ref[:nw][gm.cnw][:junction][parse(Int64,idx)]["pmin"]
    end
end

function check_ratio(sol, gm)
    for (idx,val) in sol["connection"]
        k = parse(Int64,idx)
        connection = gm.ref[:nw][gm.cnw][:connection][parse(Int64,idx)]
        if connection["type"] == "compressor" || connection["type"] == "control_valve"          
            @test val["ratio"] <= connection["c_ratio_max"] + 1e-6
            @test val["ratio"] >= connection["c_ratio_min"] - 1e-6
        end
    end
end


#Check the second order code model
@testset "test misocp gf" begin
    @testset "gaslib 40 case" begin
        result = run_gf("../test/data/gaslib-40.json", MISOCPGasModel, misocp_solver)
        @test result["status"] == :LocalOptimal || result["status"] == :Optimal
        @test isapprox(result["objective"], 0; atol = 1e-6)
        data = GasModels.parse_file("../test/data/gaslib-40.json")  
        gm = GasModels.build_generic_model(data, MINLPGasModel, GasModels.post_gf)        
        check_pressure_status(result["solution"], gm)
        check_ratio(result["solution"], gm)             
    end      
    @testset "gaslib 135 case" begin
        result = run_gf("../test/data/gaslib-135.json", MISOCPGasModel, misocp_solver)
        @test result["status"] == :LocalOptimal || result["status"] == :Optimal
        @test isapprox(result["objective"], 0; atol = 1e-6)
        data = GasModels.parse_file("../test/data/gaslib-135.json")  
        gm = GasModels.build_generic_model(data, MINLPGasModel, GasModels.post_gf)                  
        check_pressure_status(result["solution"], gm) 
        check_ratio(result["solution"], gm)          
    end
end

@testset "test minlp gf mathematical program" begin
    data = GasModels.parse_file("../test/data/gaslib-582.json")
    gm = GasModels.build_generic_model(data, MINLPGasModel, GasModels.post_gf)
    @test length(gm.var[:nw][gm.cnw][:p])  == 610
    @test length(gm.var[:nw][gm.cnw][:f])  == 637
    @test length(gm.var[:nw][gm.cnw][:yp]) == 637
    @test length(gm.var[:nw][gm.cnw][:yn]) == 637
    @test haskey(gm.var[:nw][gm.cnw],:l)   == false
    @test length(gm.var[:nw][gm.cnw][:v])  == 72  
       
    ref = gm.con[:nw][gm.cnw][:junction_flow_balance][100]  
    c = gm.model.linconstr[ref.idx]
    @test isapprox(c.ub, 0.0555264; atol = 1e-4)
    @test JuMP.sense(c) == :(==)
    @test length(c.terms.coeffs) == 1
    @test isapprox(c.terms.coeffs[1], -1.0; atol = 1e-4)
    @test c.terms.vars[1] == gm.var[:nw][gm.cnw][:f][128]
    
    # -f[426] - f[77] + f[78] == 0
    ref = gm.con[:nw][gm.cnw][:junction_flow_balance][306]  
    c = gm.model.linconstr[ref.idx]
    @test isapprox(c.ub, 0.0; atol = 1e-4)
    @test JuMP.sense(c) == :(==)
    @test length(c.terms.coeffs) == 3
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:f][426]
            @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:f][77]
            @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)       
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:f][78]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)
        else
            @test true == false
        end
    end

    # -f[360] - f[269] == -45.4464
    ref = gm.con[:nw][gm.cnw][:junction_flow_balance][26]  
    c = gm.model.linconstr[ref.idx]      
    @test isapprox(c.ub, -45.4464; atol = 1e-4)
    @test JuMP.sense(c) == :(==)
    @test length(c.terms.coeffs) == 2
                
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:f][360]
            @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:f][269]
            @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)       
         else
            @test true == false
         end
    end
    
    # "yp[360] + yp[269] >= 1"    
    ref = gm.con[:nw][gm.cnw][:source_flow][26]  
    c = gm.model.linconstr[ref.idx]      
    @test isapprox(c.lb, 1.0; atol = 1e-4)
    @test JuMP.sense(c) == :>=
    @test length(c.terms.coeffs) == 2
                
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:yp][360]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:yp][269]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)       
        else
            @test true == false
        end
    end
    
    #"yp[483] >= 1"
    ref = gm.con[:nw][gm.cnw][:sink_flow][112]  
    c = gm.model.linconstr[ref.idx]      
    @test isapprox(c.lb, 1.0; atol = 1e-4)
    @test JuMP.sense(c) == :>=
    @test length(c.terms.coeffs) == 1
    @test isapprox(c.terms.coeffs[1], 1.0; atol = 1e-4)
    @test c.terms.vars[1] == gm.var[:nw][gm.cnw][:yp][483]
  
    #  "yn[0] + yn[294] + yn[327] + yn[295] + yn[293] + yn[296] + yp[248] + yp[275] >= 1"
    ref = gm.con[:nw][gm.cnw][:sink_flow][32]  
    c = gm.model.linconstr[ref.idx]      
    @test isapprox(c.lb, 1.0; atol = 1e-4)
    @test JuMP.sense(c) == :>=
    @test length(c.terms.coeffs) == 8
                    
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] in [gm.var[:nw][gm.cnw][:yn][0], gm.var[:nw][gm.cnw][:yn][294], gm.var[:nw][gm.cnw][:yn][327], gm.var[:nw][gm.cnw][:yn][295], gm.var[:nw][gm.cnw][:yn][293], gm.var[:nw][gm.cnw][:yn][296], gm.var[:nw][gm.cnw][:yp][248], gm.var[:nw][gm.cnw][:yp][275] ]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)
        else
            @test true == false
        end
    end
    
    # yn[239] - yp[238] == 0
    ref = gm.con[:nw][gm.cnw][:conserve_flow1][523]  
    c = gm.model.linconstr[ref.idx]
              
    @test isapprox(c.ub, 0.0; atol = 1e-4)
    @test JuMP.sense(c) == :(==)
    @test length(c.terms.coeffs) == 2                  
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:yn][239]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:yp][238]
            @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)                
        else
            @test true == false
        end
    end
            
    # yn[221] - yn[514] == 0
    ref = gm.con[:nw][gm.cnw][:conserve_flow1][496] 
    c = gm.model.linconstr[ref.idx]              
    @test isapprox(c.ub, 0.0; atol = 1e-4) 
    @test JuMP.sense(c) == :(==) 
    @test length(c.terms.coeffs) == 2
                    
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:yn][221]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:yn][514]
            @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)                
        else
            @test true == false
        end
    end
          
    # yp[77] - yp[76] == 0
    ref = gm.con[:nw][gm.cnw][:conserve_flow1][305]   
    c = gm.model.linconstr[ref.idx]            
    @test isapprox(c.ub, 0.0; atol = 1e-4)
    @test JuMP.sense(c) == :(==)
    @test length(c.terms.coeffs) == 2                    
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:yp][77]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:yp][76]
            @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)                
        else
            @test true == false
        end
    end
        
    # yp[178] + yn[178] == 1
    ref = gm.con[:nw][gm.cnw][:flow_direction_choice][178]        
    c = gm.model.linconstr[ref.idx]      
    @test isapprox(c.lb, 1.0; atol = 1e-4)
    @test JuMP.sense(c) == :(==)
    @test length(c.terms.coeffs) == 2                
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:yp][178]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:yn][178]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)       
        else
            @test true == false
        end
    end
    
    # 5038.8285000000005 yp[161] - p[503] + p[415] <= 5038.8285000000005
    # p[503] - p[415] + 5038.8285000000005 yn[161] <= 5038.8285000000005
    ref = gm.con[:nw][gm.cnw][:on_off_pressure_drop2][161]  
    c = gm.model.linconstr[ref.idx]      
    @test isapprox(c.ub, 5038.8285000000005; atol = 1e-4)
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 3                
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:yn][161]
            @test isapprox(c.terms.coeffs[i], 5038.8285000000005; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][503]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)       
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][415]
            @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)          
        else
            @test true == false
        end
    end
                
    ref = gm.con[:nw][gm.cnw][:on_off_pressure_drop1][161]  
    c = gm.model.linconstr[ref.idx]      
    @test isapprox(c.ub, 5038.8285000000005; atol = 1e-4)
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 3                

    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:yp][161]
            @test isapprox(c.terms.coeffs[i], 5038.8285000000005; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][503]
            @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)       
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][415]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)          
        else
            @test true == false
       end
    end
   
    # 19.877719475620037 yp[186] - f[186] <= 19.877719475620037
    # f[186] + 19.877719475620037 yn[186] <= 19.877719475620037
    ref = gm.con[:nw][gm.cnw][:on_off_pipe_flow_direction1][186]
    c = gm.model.linconstr[ref.idx]      
    @test isapprox(c.ub, 19.877719475620037; atol = 1e-4)
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 2                
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:yp][186]
            @test isapprox(c.terms.coeffs[i], 19.877719475620037; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:f][186]
            @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)       
        else
            @test true == false
        end
    end
                
    ref = gm.con[:nw][gm.cnw][:on_off_pipe_flow_direction2][186]
    c = gm.model.linconstr[ref.idx]      
    @test isapprox(c.ub, 19.877719475620037; atol = 1e-4)
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 2                
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:yn][186]
            @test isapprox(c.terms.coeffs[i], 19.877719475620037; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:f][186]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)       
        else
            @test true == false
        end
    end
    
    # "6.124894594 * (p[498] - p[129]) - (f[222] ^ 2.0 - (1.0 - yp[222]) * 325.31057760000004 ^ 2.0) >= 0"
    # "6.124894594 * (p[498] - p[129]) - (f[222] ^ 2.0 + (1.0 - yp[222]) * 325.31057760000004 ^ 2.0) <= 0"                
    # "6.124894594 * (p[129] - p[498]) - (f[222] ^ 2.0 + (1.0 - yn[222]) * 325.31057760000004 ^ 2.0) <= 0"                
    # "6.124894594 * (p[129] - p[498]) - (f[222] ^ 2.0 - (1.0 - yn[222]) * 325.31057760000004 ^ 2.0) >= 0"                                
    ref = gm.con[:nw][gm.cnw][:weymouth1][222]
    c = gm.model.nlpdata.nlconstr[ref.idx]
    @test JuMP.sense(c) == :>=    
    @test isapprox(c.lb, 0.0; atol = 1e-4)
    @test length(c.terms.nd) == 17 
    #@test string(ref) == "6.124894594 * (p[498] - p[129]) - (f[222] ^ 2.0 - (1.0 - yp[222]) * 325.31057760000004 ^ 2.0) >= 0"

    ref = gm.con[:nw][gm.cnw][:weymouth2][222]
    c = gm.model.nlpdata.nlconstr[ref.idx]  
    @test JuMP.sense(c) == :<=    
    @test isapprox(c.ub, 0.0; atol = 1e-4)
    @test length(c.terms.nd) == 17 
    #@test string(ref) == "6.124894594 * (p[498] - p[129]) - (f[222] ^ 2.0 + (1.0 - yp[222]) * 325.31057760000004 ^ 2.0) <= 0"
    
    ref = gm.con[:nw][gm.cnw][:weymouth4][222]
    c = gm.model.nlpdata.nlconstr[ref.idx]
    @test JuMP.sense(c) == :<=    
    @test isapprox(c.ub, 0.0; atol = 1e-4)
    @test length(c.terms.nd) == 17
    #@test string(ref) == "6.124894594 * (p[129] - p[498]) - (f[222] ^ 2.0 + (1.0 - yn[222]) * 325.31057760000004 ^ 2.0) <= 0"
    
    ref = gm.con[:nw][gm.cnw][:weymouth3][222]
    c = gm.model.nlpdata.nlconstr[ref.idx]
    @test JuMP.sense(c) == :>=    
    @test isapprox(c.lb, 0.0; atol = 1e-4)
    @test length(c.terms.nd) == 17 
    #@test string(ref) == "6.124894594 * (p[129] - p[498]) - (f[222] ^ 2.0 - (1.0 - yn[222]) * 325.31057760000004 ^ 2.0) >= 0"
                    
    # p[302] - p[83] == 0
    ref = gm.con[:nw][gm.cnw][:short_pipe_pressure_drop][423]  
    c = gm.model.linconstr[ref.idx]      
    @test isapprox(c.ub, 0.0; atol = 1e-4)
    @test JuMP.sense(c) == :(==)
    @test length(c.terms.coeffs) == 2                

    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][302]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][83]
            @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)                            
        else
            @test true == false
        end
   end
    
    # 162.6552888 yp[321] - f[321] <= 325.31057760000004
    # f[321] + 162.6552888 yn[321] <= 325.31057760000004
   ref = gm.con[:nw][gm.cnw][:on_off_short_pipe_flow_direction1][321]
   c = gm.model.linconstr[ref.idx]
   @test isapprox(c.ub, 162.6552888; atol = 1e-4)
   @test JuMP.sense(c) == :<=
   @test length(c.terms.coeffs) == 2                
   for i = 1:length(c.terms.vars)
       if c.terms.vars[i] == gm.var[:nw][gm.cnw][:yp][321]
           @test isapprox(c.terms.coeffs[i], 162.6552888; atol = 1e-4)
       elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:f][321]
           @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)       
        else
           @test true == false
        end
    end
                
    ref = gm.con[:nw][gm.cnw][:on_off_short_pipe_flow_direction2][321]
    c = gm.model.linconstr[ref.idx]      
    @test isapprox(c.ub, 162.6552888; atol = 1e-4)
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 2                

    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:yn][321]
            @test isapprox(c.terms.coeffs[i], 162.6552888; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:f][321]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)       
        else
            @test true == false
        end
    end
   
    # 162.6552888 yp[549] - f[549] <= 162.6552888
    # f[549] + 162.6552888 yn[549] <= 162.6552888
    ref = gm.con[:nw][gm.cnw][:on_off_compressor_flow_direction1][549]
    c = gm.model.linconstr[ref.idx]      
    @test isapprox(c.ub, 162.6552888; atol = 1e-4)
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 2                

    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:yp][549]
            @test isapprox(c.terms.coeffs[i], 162.6552888; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:f][549]
            @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)       
        else
            @test true == false
        end
    end

    ref = gm.con[:nw][gm.cnw][:on_off_compressor_flow_direction2][549]
    c = gm.model.linconstr[ref.idx]      
    @test isapprox(c.ub, 162.6552888; atol = 1e-4)
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 2                

    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:yn][549]
            @test isapprox(c.terms.coeffs[i], 162.6552888; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:f][549]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)       
        else
            @test true == false
        end
    end
    
    # p[2200560] - 25 p[560] + 7075.038568175958 yp[551] <= 7075.038568175958   
    # p[560] - 25 p[2200560] + 5042.8816755625 yn[551] <= 5042.8816755625
    # p[560] - p[2200560] + 5042.8816755625 yp[551] <= 5042.8816755625
    # p[2200560] - p[560] + 7075.038568175958 yn[551] <= 7075.038568175958
    ref = gm.con[:nw][gm.cnw][:on_off_compressor_ratios4][551]
    c = gm.model.linconstr[ref.idx]      
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 3 
    
    @test isapprox(c.ub, 7075.038568175958; atol = 1e-4)    
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][2200560]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][560]
            @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)       
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:yn][551]
            @test isapprox(c.terms.coeffs[i], 7075.038568175958; atol = 1e-4)       
        else
            @test true == false
        end
    end

    ref = gm.con[:nw][gm.cnw][:on_off_compressor_ratios1][551]
    c = gm.model.linconstr[ref.idx]      
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 3               
    @test isapprox(c.ub, 7075.038568175958; atol = 1e-4)
    
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][2200560]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][560]
            @test isapprox(c.terms.coeffs[i], -25.0; atol = 1e-4)       
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:yp][551]
            @test isapprox(c.terms.coeffs[i], 7075.038568175958; atol = 1e-4)       
        else
            @test true == false
        end
    end
   
    ref = gm.con[:nw][gm.cnw][:on_off_compressor_ratios3][551]
    c = gm.model.linconstr[ref.idx]      
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 3               
    @test isapprox(c.ub, 5042.8816755625; atol = 1e-4)
    
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][560]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)
         elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][2200560]
            @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)       
         elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:yn][551]
            @test isapprox(c.terms.coeffs[i], 5042.8816755625; atol = 1e-4)       
         else
            @test true == false
         end
    end

    ref = gm.con[:nw][gm.cnw][:on_off_compressor_ratios2][551]
    c = gm.model.linconstr[ref.idx]      
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 3               
    @test isapprox(c.ub, 5042.8816755625; atol = 1e-4)
        
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][560]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][2200560]
            @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)       
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:yp][551]
            @test isapprox(c.terms.coeffs[i], 5042.8816755625; atol = 1e-4)       
        else
            @test true == false
       end
    end
    
    # f[558] - 162.6552888 v[558] <= 0
    # -162.6552888 v[558] - f[558] <= 0
    # 162.6552888 yp[558] - f[558] <= 162.6552888
    # f[558] + 162.6552888 yn[558] <= 162.6552888
    ref = gm.con[:nw][gm.cnw][:on_off_valve_flow_direction2][558]   
    c = gm.model.linconstr[ref.idx]      
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 2                
    @test isapprox(c.ub, 162.6552888; atol = 1e-4)
        
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:f][558]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:yn][558]
            @test isapprox(c.terms.coeffs[i], 162.6552888; atol = 1e-4)       
        else
            @test true == false
        end
    end

    ref = gm.con[:nw][gm.cnw][:on_off_valve_flow_direction1][558]   
    c = gm.model.linconstr[ref.idx]      
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 2                
    @test isapprox(c.ub, 162.6552888; atol = 1e-4)

    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:f][558]
            @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:yp][558]
            @test isapprox(c.terms.coeffs[i], 162.6552888; atol = 1e-4)       
        else
            @test true == false
       end
    end
        
    ref = gm.con[:nw][gm.cnw][:on_off_valve_flow_direction3][558]   
    c = gm.model.linconstr[ref.idx]      
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 2                
    @test isapprox(c.ub, 0.0; atol = 1e-4)    
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:f][558]
            @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:v][558]
            @test isapprox(c.terms.coeffs[i], -162.6552888; atol = 1e-4)       
        else
            @test true == false
        end
    end
                
    ref = gm.con[:nw][gm.cnw][:on_off_valve_flow_direction4][558]   
    c = gm.model.linconstr[ref.idx]      
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 2                
    @test isapprox(c.ub, 0.0; atol = 1e-4)
 
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:f][558]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:v][558]
            @test isapprox(c.terms.coeffs[i], -162.6552888; atol = 1e-4)       
        else
            @test true == false
        end
    end
            
    # p[164] + 7398.2791755625 v[571] - p[170] <= 7398.2791755625
    # p[170] - p[164] + 7398.2791755625 v[571] <= 7398.2791755625
    ref = gm.con[:nw][gm.cnw][:on_off_valve_pressure_drop2][571]
    c = gm.model.linconstr[ref.idx]      
    @test isapprox(c.ub, 7398.2791755625; atol = 1e-4)
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 3                      
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][164]
            @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:v][571]
            @test isapprox(c.terms.coeffs[i], 7398.2791755625; atol = 1e-4)       
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][170]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)          
        else
            @test true == false
        end
    end

    ref = gm.con[:nw][gm.cnw][:on_off_valve_pressure_drop1][571]
    c = gm.model.linconstr[ref.idx]      
    @test isapprox(c.ub, 7398.2791755625; atol = 1e-4)
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 3                
   
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][164]
           @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:v][571]
           @test isapprox(c.terms.coeffs[i], 7398.2791755625; atol = 1e-4)       
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][170]
           @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)          
        else
           @test true == false
        end
    end
    
    # -162.6552888 v[591] - f[591] <= 0
    # f[591] + 325.31057760000004 yn[591] <= 162.6552888
    # 162.6552888 yp[591] - f[591] <= 162.6552888
    # f[591] - 162.6552888 v[591] <= 0
    ref = gm.con[:nw][gm.cnw][:on_off_control_valve_flow_direction4][591]
    c = gm.model.linconstr[ref.idx]      
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 2                
    @test isapprox(c.ub, 0.0; atol = 1e-4)
    
    #    ok = true
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:v][591]
            @test isapprox(c.terms.coeffs[i], -162.6552888; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:f][591]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)       
        else
            @test true == false
        end
    end
   
    ref = gm.con[:nw][gm.cnw][:on_off_control_valve_flow_direction1][591]
    c = gm.model.linconstr[ref.idx]      
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 2                
    @test isapprox(c.ub, 162.6552888; atol = 1e-4)
   
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:yp][591]
            @test isapprox(c.terms.coeffs[i], 162.6552888; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:f][591]
            @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)       
         else
            @test true == false
         end
    end

    ref = gm.con[:nw][gm.cnw][:on_off_control_valve_flow_direction2][591]
    c = gm.model.linconstr[ref.idx]      
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 2                
    @test isapprox(c.ub, 162.6552888; atol = 1e-4)
    
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:yn][591]
            @test isapprox(c.terms.coeffs[i], 162.6552888; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:f][591]
            @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)       
        else
            @test true == false
        end
    end
        
    ref = gm.con[:nw][gm.cnw][:on_off_control_valve_flow_direction3][591]
    c = gm.model.linconstr[ref.idx]      
    @test JuMP.sense(c) == :<=
    @test length(c.terms.coeffs) == 2                
    @test isapprox(c.ub, 0.0; atol = 1e-4)
    for i = 1:length(c.terms.vars)
        if c.terms.vars[i] == gm.var[:nw][gm.cnw][:v][591]
            @test isapprox(c.terms.coeffs[i], -162.6552888; atol = 1e-4)
        elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:f][591]
            @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)       
        else
            @test true == false
        end
    end
    
    # 0_p[2600217] - 0_p[217] + 7398.2791755625 0_yp[585] + 7398.2791755625 0_v[585] <= 14796.558351125 
    # - p[2600217] + 7398.2791755625 yp[585]   + 7398.2791755625 v[585]   <= 14796.558351125
    # 0_p[2600217] - 0_p[217] + 7398.2791755625 0_yn[585] + 7398.2791755625 0_v[585] <= 14796.558351125
    # 0_p[217] - 0_p[2600217] + 7398.2791755625 0_yn[585] + 7398.2791755625 0_v[585] <= 14796.558351125
    ref = gm.con[:nw][gm.cnw][:on_off_control_valve_pressure_drop2][585]
    c = gm.model.linconstr[ref.idx]      
    @test JuMP.sense(c) == :<=
    @test isapprox(c.ub, 14796.558351125; atol = 1e-4) && length(c.terms.coeffs) == 3
    
    for i = 1:length(c.terms.vars)
       if c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][2600217]
           @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)       
       elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:yp][585]
           @test isapprox(c.terms.coeffs[i], 7398.2791755625; atol = 1e-4)       
       elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:v][585]
           @test isapprox(c.terms.coeffs[i], 7398.2791755625; atol = 1e-4)          
       else
           @test true == false
       end
   end

   ref = gm.con[:nw][gm.cnw][:on_off_control_valve_pressure_drop4][585]
   c = gm.model.linconstr[ref.idx]      
   @test JuMP.sense(c) == :<=
   @test isapprox(c.ub, 14796.558351125; atol = 1e-4) && length(c.terms.coeffs) == 4
    
   for i = 1:length(c.terms.vars)
       if c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][217]
           @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)
       elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][2600217]
           @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)       
       elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:yn][585]
           @test isapprox(c.terms.coeffs[i], 7398.2791755625; atol = 1e-4)       
       elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:v][585]
           @test isapprox(c.terms.coeffs[i], 7398.2791755625; atol = 1e-4)          
       else
           @test true == false
       end
   end

   ref = gm.con[:nw][gm.cnw][:on_off_control_valve_pressure_drop1][585]
   c = gm.model.linconstr[ref.idx]      
   @test JuMP.sense(c) == :<=
   @test isapprox(c.ub, 14796.558351125; atol = 1e-4) && length(c.terms.coeffs) == 4
    
   for i = 1:length(c.terms.vars)
       if c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][217]
           @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)
       elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][2600217]
           @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)       
       elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:yp][585]
           @test isapprox(c.terms.coeffs[i], 7398.2791755625; atol = 1e-4)       
       elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:v][585]
           @test isapprox(c.terms.coeffs[i], 7398.2791755625; atol = 1e-4)          
       else
           @test true == false
       end
   end
  
   ref = gm.con[:nw][gm.cnw][:on_off_control_valve_pressure_drop3][585]
   c = gm.model.linconstr[ref.idx]      
   @test JuMP.sense(c) == :<=    
   @test isapprox(c.ub, 14796.558351125; atol = 1e-4) && length(c.terms.coeffs) == 4
   
   for i = 1:length(c.terms.vars)
       if c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][217]
           @test isapprox(c.terms.coeffs[i], -1.0; atol = 1e-4)
       elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:p][2600217]
           @test isapprox(c.terms.coeffs[i], 1.0; atol = 1e-4)       
       elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:yn][585]
           @test isapprox(c.terms.coeffs[i], 7398.2791755625; atol = 1e-4)       
       elseif c.terms.vars[i] == gm.var[:nw][gm.cnw][:v][585]
           @test isapprox(c.terms.coeffs[i], 7398.2791755625; atol = 1e-4)          
       else
           @test true == false
       end
   end

end

