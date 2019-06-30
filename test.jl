using Test

function test_mpi()
    @testset "O Exemplo Perdido" begin

    c = [-2.; -1; 0; 0; 0]
    A = [1 8/3 1 0 0;
         1 1 0 1 0;
         2 0 0 0 1]
    b = [4.; 2; 3]

    x₀ = [.5; 1; 5/6; .5; 2] # x1 e x2 obtidos visualmente e os demais substituídos
    y₀ = [-120.; -100; -150] # tentativa e erro
    z₀ = [518.; 419; 120; 100; 150]

    x, y, z =  mpi_basic!(c, A, b, x₀, y₀, z₀)

    @test all(x .>= 0)
    @test all(z .>= 0)
    @test dot(x, c) ≈ -3.5 atol = 1e-4

    x, y, z =  mpi_mehrotra!(c, A, b, x₀, y₀, z₀)

    @test all(x .>= 0)
    @test all(z .>= 0)
    @test dot(x, c) ≈ -3.5 atol = 1e-4

    end
end
