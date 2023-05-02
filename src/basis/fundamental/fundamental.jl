include("fundamentalhankel.jl")


@inline function greens_fun(basis::AbsFundamentalBasis, k::T, pts::AbstractArray, source_pts::AbstractArray) where {T<:Real}
    let symmetries = basis.symmetries
        M =  length(pts)
        N = length(source_pts)
        G = zeros(basis.return_type,M,N)
        Threads.@threads for j in 1:N
            pts_loc = [pt - source_pts[j] for pt in pts]
            arg = [k*hypot(pt[1], pt[2]) for pt in pts_loc]
            G[:,j] .= basis_fun(basis, arg)
            if ~isnothing(symmetries)
                for sym in symmetries
                    pts_sym = sym.sym_map.(pts)
                    pts_loc = [pt - source_pts[j] for pt in pts_sym]
                    arg = [k*hypot(pt[1], pt[2]) for pt in pts_loc]
                    G[:,j] .+= sym.parity .* basis_fun(basis, arg)
                end
            end
        end
        if ~isnothing(symmetries)
            return G./(length(symmetries)+1.0)
        else
            return G
        end 
    end
end

@inline function greens_fun(basis::AbsFundamentalBasis, k::T, pts::AbstractArray) where {T<:Real}
    let symmetries = basis.symmetries
        #M =  length(pts)
        N = length(pts)
        G = zeros(basis.return_type,N,N)
        #fill upper triangle
        Threads.@threads for j in 1:N
            for i in 1:j-1
                pt = pts[i] - pts[j]
                arg = k*hypot(pt[1], pt[2])
                G[i,j] = basis_fun(basis, arg)
                if ~isnothing(symmetries)
                    for sym in symmetries
                        pt_sym = sym.sym_map.(pts[i]) - pts[j]
                        arg = k*hypot(pt_sym[1], pt_sym[2])
                        G[i,j] += sym.parity * basis_fun(basis, arg)
                    end
                end
            end
        end
        G = G + transpose(G) #lower triangle is equal to transpose of upper triangle
        #diagonal terms
        Threads.@threads for i in 1:N
            pt = pts[i] - pts[i]
            arg = k*hypot(pt[1], pt[2])
            G[i,i] = basis_fun(basis, arg)
            if ~isnothing(symmetries)
                for sym in symmetries
                    pt_sym = sym.sym_map.(pts[i]) - pts[i]
                    arg = k*hypot(pt_sym[1], pt_sym[2])
                    G[i,i] += sym.parity * basis_fun(basis, arg)
                end
            end
        end

        if ~isnothing(symmetries)
            return G./(length(symmetries)+1.0)
        else
            return G
        end 
    end
end


@inline function greens_gradient(basis::AbsFundamentalBasis, k::T, pts::AbstractArray, source_pts::AbstractArray) where {T<:Real}
    let symmetries = basis.symmetries
        M =  length(pts)
        N = length(source_pts)
        dG_dx = zeros(basis.return_type,M,N)
        dG_dy = zeros(basis.return_type,M,N)
        Threads.@threads for j in 1:N
            pts_loc = [pt - source_pts[j] for pt in pts]
            r = [hypot(pt[1], pt[2]) for pt in pts_loc]
            x = getindex.(pts_loc, 1)
            y = getindex.(pts_loc, 2) 
            arg = k.*r
            dg = derivative_fun(basis, arg)
            dG_dx[:,j] .= @. k*x/r*dg
            dG_dy[:,j] .= @. k*y/r*dg
            if ~isnothing(symmetries)
                for sym in symmetries
                    pts_sym = sym.sym_map.(pts)
                    pts_loc = [pt - source_pts[j] for pt in pts_sym]
                    arg = k.*r
                    dg = sym.parity * derivative_fun(basis, arg)
                    dG_dx[:,j] .= @. k*x/r*dg
                    dG_dy[:,j] .= @. k*y/r*dg
                end
            end
        end
        if ~isnothing(symmetries)
            n = (length(symmetries)+1.0)
            return dG_dx./n, dG_dy./n
        else
        #println(size(s))
            return dG_dx, dG_dy
        end 
    end
end

@inline function greens_gradient(basis::AbsFundamentalBasis, k::T, pts::AbstractArray; return_diagonal=false) where {T<:Real}
    let symmetries = basis.symmetries
        #M =  length(pts)
        N = length(pts)
        dG_dx = zeros(basis.return_type,N,N)
        dG_dy = zeros(basis.return_type,N,N)
        #fill upper triangle
        for j in 1:N #Threads.@threads
            for i in 1:j-1
                pt = pts[i] - pts[j]
                r = hypot(pt[1], pt[2])
                x = pt[1]
                y = pt[2]
                arg = k.*r
                dg = derivative_fun(basis, arg)
                dG_dx[i,j] = k*x/r*dg
                dG_dy[i,j] = k*y/r*dg
                if ~isnothing(symmetries)
                    for sym in symmetries
                        pt_sym = sym.sym_map(pts[i]) - pts[j]
                        r = hypot(pt_sym[1], pt_sym[2])
                        x = pt[1]
                        y = pt[2]
                        arg = k.*r
                        dg = derivative_fun(basis, arg)
                        dG_dx[i,j] = k*x/r*dg #does this have to be abs(x)?
                        dG_dy[i,j] = k*y/r*dg
                    end
                end
            end
        end
        dG_dx = dG_dx + transpose(dG_dx) #lower triangle is equal to transpose of upper triangle
        dG_dy = dG_dy + transpose(dG_dy)
        #diagonal terms
        if return_diagonal
            Threads.@threads for i in 1:N
                pt = pts[i] - pts[i]
                r = hypot(pt[1], pt[2])
                x = pt[1]
                y = pt[2]
                arg = k.*r
                dg = derivative_fun(basis, arg)
                dG_dx[i,i] = k*x/r*dg
                dG_dy[i,i] = k*y/r*dg
                if ~isnothing(symmetries)
                    for sym in symmetries
                        pt_sym = sym.sym_map.(pts[i]) - pts[j]
                        r = hypot(pt_sym[1], pt_sym[2])
                        x = pt[1]
                        y = pt[2]
                        arg = k.*r
                        dg = derivative_fun(basis, arg)
                        dG_dx[i,i] = k*x/r*dg #does this have to be abs(x)?
                        dG_dy[i,i] = k*y/r*dg
                    end
                end
            end
        end

        if ~isnothing(symmetries)
            n = (length(symmetries)+1.0)
            return dG_dx./n, dG_dy./n
        else
        #println(size(s))
            return dG_dx, dG_dy
        end 
    end
end





#=
using Bessels
-im * 0.25 * Bessels.hankelh1(0, 0.0) 
-im * 0.25 * Bessels.hankelh1(1, 0.0)
a = Bessels.hankelh1(0, 0.0)
-a * im * 0.25
using Random
G = rand(Complex{Float64},(5,5))

transpose(G)

G = zeros(Complex{Float64},5,5)
Threads.@threads for j in 1:size(G,1)
    for i in 1:j-1
        G[i,j] = rand(Complex{Float64})
    end
end
G

G = G + transpose(G)
=#