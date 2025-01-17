@testset "transient parsing" begin
    data = parse_transient("../test/data/transient/time-series-case-6b.csv")
    @test length(data) == 505
    @test all(length(item) == 5 for item in data)

    mn_data = parse_files("../test/data/matgas/case-6.m", "../test/data/transient/time-series-case-6b.csv", spatial_discretization = 1e4, additional_time = 0.0)
    @test length(mn_data["nw"]) == 25
    @test isapprox(mn_data["time_step"] * mn_data["base_time"], 3600.0; atol = 1e-3)

    mn_data = parse_files("../test/data/matgas/case-6.m", "../test/data/transient/time-series-case-6a.csv", spatial_discretization = 1e4, additional_time = 0.0)
    @test length(mn_data["nw"]) == 25
    @test isapprox(mn_data["time_step"] * mn_data["base_time"], 3600.0; atol = 1e-3)

    ss_data = parse_file("../test/data/matgas/case-6.m")
    pipe_ids = collect(keys(ss_data["pipe"]))

    @test length(first(mn_data["nw"]).second["original_pipe"]) == length(pipe_ids)
    @test all(k in pipe_ids for k in keys(first(mn_data["nw"]).second["original_pipe"]))
end

@testset "transient (steady state replicate) case" begin
    mn_data = parse_files("../test/data/matgas/case-6.m", "../test/data/transient/time-series-case-6b.csv", spatial_discretization = 1e4, additional_time = 0.0)
    result = run_transient_ogf(mn_data, WPGasModel, nlp_solver)
    @test result["termination_status"] == LOCALLY_SOLVED
    @test isapprox(result["objective"], -5160.91, atol = 1e-1)
    make_si_units!(result["solution"])
    @test isapprox(result["solution"]["nw"]["1"]["receipt"]["1"]["injection"], 94.80, atol = 1e-1)
    @test isapprox(result["solution"]["nw"]["2"]["receipt"]["1"]["injection"], 94.80, atol = 1e-1)
    @test isapprox(result["solution"]["nw"]["2"]["compressor"]["1"]["power_var"], 3e6, atol = 10)
end

@testset "transient time-periodic withdrawal case" begin
    mn_data = parse_files("../test/data/matgas/case-6.m", "../test/data/transient/time-series-case-6a.csv", spatial_discretization = 1e4, additional_time = 0.0)
    result = run_transient_ogf(mn_data, WPGasModel, nlp_solver)
    @test result["termination_status"] == LOCALLY_SOLVED
    @test isapprox(result["objective"], -5002, atol = 1)
    make_si_units!(result["solution"])
    @test isapprox(result["solution"]["nw"]["1"]["receipt"]["1"]["injection"], 87.75, atol = 1e-1)
    @test isapprox(result["solution"]["nw"]["2"]["receipt"]["1"]["injection"], 93.27, atol = 1e-1)
    @test isapprox(result["solution"]["nw"]["2"]["compressor"]["1"]["power_var"], 3e6, atol = 10)
end

@testset "transient (steady state replicate) case with storage" begin
    mn_data = parse_files("../test/data/matgas/case-6-storage.m", "../test/data/transient/time-series-case-6b.csv", spatial_discretization = 1e4, additional_time = 3600.0)
    result = run_transient_ogf(mn_data, WPGasModel, nlp_solver)
    @test result["termination_status"] == LOCALLY_SOLVED
    make_si_units!(result["solution"])
    @test isapprox(result["solution"]["nw"]["1"]["receipt"]["1"]["injection"], 44, atol = 1)
    @test isapprox(result["solution"]["nw"]["2"]["receipt"]["1"]["injection"], 44, atol = 1)
    @test isapprox(result["solution"]["nw"]["1"]["storage"]["1"]["withdrawal"], -86, atol = 1)
    @test isapprox(result["solution"]["nw"]["2"]["storage"]["1"]["withdrawal"], -86, atol = 1)
end

@testset "transient time-periodic withdrawal case with storage" begin
    mn_data = parse_files("../test/data/matgas/case-6-storage.m", "../test/data/transient/time-series-case-6a.csv", spatial_discretization = 1e4, additional_time = 3600.0)
    result = run_transient_ogf(mn_data, WPGasModel, nlp_solver)
    @test result["termination_status"] == LOCALLY_SOLVED
    @test isapprox(result["objective"], -9282, atol = 1)
    make_si_units!(result["solution"])
    @test isapprox(result["solution"]["nw"]["1"]["receipt"]["1"]["injection"], 0, atol = 1e-1)
    @test isapprox(result["solution"]["nw"]["2"]["receipt"]["1"]["injection"], 0.0, atol = 1e-1)
    @test isapprox(result["solution"]["nw"]["1"]["storage"]["1"]["withdrawal"], -83, atol = 1)
    @test isapprox(result["solution"]["nw"]["2"]["storage"]["1"]["withdrawal"], -83, atol = 1)
end
