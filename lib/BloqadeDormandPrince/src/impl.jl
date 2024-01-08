
# implement DormandPrince API for BloqadeSolver
DormandPrince.get_current_state(solver::BloqadeDPSolver) = solver.reg
DormandPrince.integrate_core!(solver::BloqadeDPSolver, time::Real) = integrate_core!(solver.dp_solver, time)

# implement BloqadeExpr API for BloqadeSolver
function BloqadeExpr.emulate!(prob::DormandPrinceProblem) 
    integrate!(prob.solver, prob.tend)
    return prob
end
