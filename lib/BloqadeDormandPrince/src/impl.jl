
# implement DormandPrince API for BloqadeSolver
DormandPrince.get_current_state(solver::BloqadeDPSolver) = solver.reg
DormandPrince.integrate_core!(solver::BloqadeDPSolver, time::Real) = integrate_core!(solver.dp_solver, time)

function register(solver::BloqadeDPSolver)
    return get_current_state(solver)
end


# implement BloqadeExpr API for BloqadeSolver
function BloqadeExpr.emulate!(prob::BloqadeDPSolver) 
    integrate!(prob.dp_solver, prob.tend)
    return prob
end
