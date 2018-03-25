using Base.Test, BasicArnoldi

function test(n = 100_000)
    # Creates a sparse matrix A and a random vector v1 for testing that the properties of H and V hold
    A = spdiagm((fill(-1.0, n-1), fill(2.0, n), fill(-1.1, n-1)), (-1,0,1))
    v1 = rand(n)
    A, basic_arnoldi(A, v1, 30)
end

@testset "Testing" begin
    A, (V, H) = test()
    @test vecnorm(V' * V -  I) < 1e-10 #This should equal to 0
    @test vecnorm(A * V[:,1:30] - V*H) < 1e-10 #This should equal to 0
end