# Tests regularized ipm

include("LDLFactorizations.jl/src/LDLFactorizations.jl")
include("RipQP.jl/src/ipm.jl")

using QPSReader, LinearAlgebra, NLPModelsIpopt, QuadraticModels, Test, Main.LDLFactorizations, QuadraticModels, LinearOperators

function getnlp(qps::QPSData)
  return nlp = QuadraticModel(qps.c, qps.Q, opHermitian(qps.Q), qps.A, qps.lcon, qps.ucon, qps.lvar, qps.uvar, c0=qps.c0)
end

function ipoptqp(qps::QPSData)
  return ipopt(getnlp(qps), print_level=0)
end

function ripqp(qps::QPSData)
  return ipm(getnlp(qps), itmax=400)
end

function test_ripqp()
  include("problems.jl")

  for problem in problems
    println("Problem $(problem):")
    qps = readqps("qpss/$(problem).qps")

    output = ripqp(qps)

    @test output[11] == true
  end
end