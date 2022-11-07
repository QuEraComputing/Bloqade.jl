using Unitful:
    μm,
    m,
    rad,
    s,
    @u_str
using BloqadeSchema:
    BloqadeSchema,
    parse_unit,
    convert_units_recursive,
    get_device_capabilities,
    get_device_capabilities_SI
using Test

@testset "parse_units" begin
    tests = [
        ("m",:(m)),
        ("s",:(s)),
        ("m/s",:(m/s)),
        ("m*s^2",:(m*s^2))
    ]
    for (str,expr) in tests
        str_val =  Base.eval(@__MODULE__,parse_unit(@__MODULE__, str))
        expr_val = Base.eval(@__MODULE__,expr)
        @test str_val ===expr_val
    end



end

@testset "convert_units_recursive" begin

    @testset "base units" begin
    # main units that are targeted for conversion, 
    # m should become μm
    # s should become μs

    # create one dictionary for units
    unit_dict = Dict(
        "a" => "m",
        "b" => "s"
    )
    # another dictionary with the actual values
    value_dict = Dict(
        "a" => 2.0,
        "b" => 5.1
    )

    # expected results
    expected_dict = Dict(
        "a" => 2.0e6,
        "b" => 5.1e6
    )

    @test convert_units_recursive(value_dict, unit_dict) == expected_dict
    end
    
    @testset "compound units" begin
    # Test compound units, 
    # rad/s should become rad/μs
    # rad/s^2 should become rad/μs^2
    # rad*m^6/s should become rad*um^6 / μs
    unit_dict = Dict(
        "a" => "rad/s",
        "b" => "rad/s^2",
        "c" => "rad*m^6/s"
    )

    value_dict = Dict(
        "a" => 1.58E7,
        "b" => 2.5E14,
        "c" => 5.42E-24,
    )

    expected_dict = Dict(
        "a" => 15.8,
        "b" => 2.5e2,
        "c" => 5.42e6
    )

    @test convert_units_recursive(value_dict, unit_dict) == expected_dict
    end

    @testset "mixed units" begin
        # units that should be converted should be converted,
        # units that shouldn't be are left alone
        unit_dict = Dict(
            "a" => "rad",
            "b" => "s",
            "c" => "m"
        )

        value_dict = Dict(
            "a" => 5.0E-7,
            "b" => 1.0E-9,
            "c" => 4.0E-6
        )

        expected_dict = Dict(
            "a" => 5.0e-7,
            "b" => 1.0e-3,
            "c" => 4.0
        )

        @test convert_units_recursive(value_dict, unit_dict) == expected_dict
    end

    @testset "no units" begin
        # anything specified with NoUnits should be left alone
        unit_dict = Dict(
            "a" => "NoUnits",
            "b" => "NoUnits",
            "c" => "NoUnits"
        )

        value_dict = Dict(
            "a" => 100,
            "b" => 1000,
            "c" => 100
        )

        @test convert_units_recursive(value_dict, unit_dict) == value_dict
    end

    @testset "nested_dictionary" begin
        # should be able to handle dictionaries inside dictionaries
        unit_dict = Dict(
            "a" => Dict(
                "a1" => "NoUnits",
                "a2" => "NoUnits",
                "a3" => Dict(
                    "a3.1" => "rad/s",
                    "a3.2" => "rad/s"
                )
            ),
            "b" => "rad*m^6/s",
            "c" => Dict(
                "c1" => "rad/s^2",
                "c2" => Dict(
                    "c2.1" => "s",
                    "c2.2" => "m"
                )
            )
        )

        value_dict = Dict(
            "a" => Dict(
                "a1" => 1,
                "a2" => 100,
                "a3" => Dict(
                    "a3.1" => 0.0,
                    "a3.2" => 1.58E7
                )
            ),
            "b" => 6.21E-24,
            "c" => Dict(
                "c1" => 3.1E15,
                "c2" => Dict(
                    "c2.1" => 1.0E-9,
                    "c2.2" => 7.5E-5
                )
            )
        )

        expected_dict = Dict(
            "a" => Dict(
                "a1" => 1,
                "a2" => 100,
                "a3" => Dict(
                    "a3.1" => 0.0,
                    "a3.2" => 15.8 
                )
            ),
            "b" => 6.21e6,
            "c" => Dict(
                "c1" => 3.1e3,
                "c2" => Dict(
                    "c2.1" => 1.0e-3,
                    "c2.2" => 75
                )
            )
        )

        @test convert_units_recursive(value_dict, unit_dict) == expected_dict
    end 
end


@testset "get_device_capabilities" begin

    @test typeof(get_device_capabilities()) === BloqadeSchema.DeviceCapabilities 
    @test typeof(get_device_capabilities_SI()) === BloqadeSchema.DeviceCapabilities

end