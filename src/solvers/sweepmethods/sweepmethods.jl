function solve_wavenumber(solver::Union{BoundaryIntegralMethod,DLP_kress,DLP_kress_global_corners,CFIE_kress,CFIE_kress_corners,CFIE_kress_global_corners,CFIE_kress_composite_solver},basis::Ba,billiard::Bi,k,dk;multithreaded::Bool=true,use_krylov::Bool=false,which::Symbol=:det_argmin) where {Ba<:AbsBasis,Bi<:BilliardGeometry.AbsBilliard}
    pts=evaluate_points(solver,billiard,k)
    f(kk)=solve(solver,basis,pts,kk;multithreaded=multithreaded,use_krylov=use_krylov,which=which)
    res=Optim.optimize(f,k-0.5*dk,k+0.5*dk,Optim.Brent();rel_tol=1e-14,abs_tol=1e-14)
    return Optim.minimizer(res),Optim.minimum(res)
end