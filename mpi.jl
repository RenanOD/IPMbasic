using SparseArrays, LinearAlgebra

function mpi_basic!(c, A, b, x, y, z, σ = .2)
    m, n = size(A)
    iter = 0
    tired = false
    solved = false
    dx, dz, dy = zeros(n), zeros(n), zeros(m)

    while !(solved || tired)

        μ = dot(x,z)/n
        dy .= (A*(x.*z.^(-1).*A'))\(A*(x - σ*μ*z.^(-1)))
        dz .= -A'*dy
        dx .= -x + σ*μ*z.^(-1) - z.^(-1).*x.*dz

        α = 1.0
        for i in 1:n
          if (dx[i] < 0 && -x[i]/dx[i] < α) α = -x[i]/dx[i] end
          if (dz[i] < 0 && -z[i]/dz[i] < α) α = -z[i]/dz[i] end
        end
        α = 0.9*α

        x += α*dx
        y += α*dy
        z += α*dz
        iter += 1

        solved = dot(x,z) < 1e-5
        tired = iter > 1e3
    end
    return x, y, z
end

function mpi_mehrotra!(c, A, b, x, y, z)
    m, n = size(A)
    iter = 0
    tired = false
    solved = false
    dx, dy, dz = zeros(n), zeros(m), zeros(n)
    dxaff, dyaff, dzaff = zeros(n), zeros(m), zeros(n)

    while !(solved || tired)

    	μ = dot(x,z)/n

       # predictive phase

    	F = lu(A*(x.*z.^(-1).*A'))

        dyaff .= F\(A*x)
        dzaff .= -A'*dyaff
        dxaff .= -x - z.^(-1).*x.*dzaff

        αaff = 1.0
        for i in 1:n
          if (dxaff[i] < 0 && -x[i]/dxaff[i] < αaff) αaff = -x[i]/dxaff[i] end
          if (dzaff[i] < 0 && -z[i]/dzaff[i] < αaff) αaff = -z[i]/dzaff[i] end
        end
        αaff *= 0.9

        μaff = dot(x + αaff*dxaff , z + αaff*dzaff)/n

        σ = (μaff/μ)^3

        # corrective phase

        dy .= F\(A*(x - σ*μ*z.^(-1) + z.^(-1).*dxaff.*dzaff))
        dz .= -A'*dy
        dx .= -x + σ*μ*z.^(-1) - z.^(-1).*x.*dz - z.^(-1).*dxaff.*dzaff

        α = 1.0
        for i in 1:n
          if (dx[i] < 0 && -x[i]/dx[i] < α) α = -x[i]/dx[i] end
          if (dz[i] < 0 && -z[i]/dz[i] < α) α = -z[i]/dz[i] end
        end
        α *= 0.9

        x += α*dx
        y += α*dy
        z += α*dz
        iter += 1

        solved = dot(x,z) < 1e-4
        tired = iter > 1e3
    end
    return x, y, z
end
