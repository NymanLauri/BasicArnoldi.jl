using BenchmarkTools

function bench(n = 100_000)
    A = spdiagm((fill(-1.0, n-1), fill(2.0, n), fill(-1.1, n-1)), (-1,0,1))
    v1 = rand(n)
    @benchmark basic_arnoldi($A, $v1, 30)
end

bench()