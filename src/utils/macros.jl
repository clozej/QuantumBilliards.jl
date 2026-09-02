
const MAX_BLAS_THREADS=Sys.CPU_THREADS

"""
    use_threads(args...)
    
The macro expects either:
(a) A keyword argument "multithreading" followed by a loop expression, or
b) A lone loop expression (in which case multithreading defaults to true).
NOTE: Already @inbounds 
"""
macro use_threads(args...)
    if length(args)>=2 && args[1] isa Expr && args[1].head==:(=) && args[1].args[1]==:multithreading
        if length(args)!=2
            error("Usage: @use_threads multithreading=[true|false] for ...")
        end
        # Extract the provided multithreading value and the loop expression.
        multithreading_val=args[1].args[2]
        loop_expr=args[2]
        # If the provided value is literally true or false, we branch accordingly at compile time.
        if multithreading_val===true
            return esc(:(@inbounds Threads.@threads $loop_expr))
        elseif multithreading_val===false
            return esc(:(@inbounds $loop_expr))
        else
            # If not a literal (i.e. it's some expression that will be evaluated at runtime),
            # we generate code that conditionally selects the threaded version.
            return esc(quote
                @inbounds begin
                    if $multithreading_val
                        Threads.@threads $loop_expr
                    else
                        $loop_expr
                    end
                end
            end)
        end
    elseif length(args)==1
        # No keyword argument provided. Default behavior is multithreading=true.
        loop_expr=args[1]
        return esc(:(@inbounds Threads.@threads $loop_expr))
    else
        error("Usage: @use_threads [multithreading=[true|false]] for ...")
    end
end

"""
    @benchit [timeit=true|false] [\"label\"] expr

The `@benchit` macro is a benchmarking tool that allows you to measure the execution time, memory allocation, and garbage collection time of a given expression.
"""
macro benchit(args...)
    isempty(args) && error("@benchit needs arguments")
    timeit_expr=:(false)
    label_expr="\"benchmark\""
    i=1
    if args[1] isa Expr && args[1].head== :(=) && args[1].args[1]==:timeit
        timeit_expr=args[1].args[2]
        i+=1
    end
    i>length(args) && error("@benchit is missing a body")
    if args[i] isa String || (args[i] isa Expr && args[i].head != :(=))
        label_expr=args[i]
        i+=1
    end
    i>length(args) && error("@benchit is missing a body")
    body=
        if i==length(args)
            args[i]
        else
            Expr(:block,args[i:end]...)
        end

    return esc(quote
        if $timeit_expr
            local _stats=Base.@timed $body
            local _t=_stats.time
            local _bytes=_stats.bytes
            local _gctime=_stats.gctime
            local _alloc_mb=_bytes / 1024^2
            println(
                $label_expr,": ",
                round(_t; digits=6)," s, ",
                round(_alloc_mb; digits=3)," MiB alloc, ",
                round(_gctime; digits=6)," s gc")
            _stats.value
        else
            $body
        end
    end)
end

"""
    @blas_threads n expr

Temporarily set BLAS threads to `n` for `expr`, then restore the previous value.
"""
macro blas_multi(n,expr)
    quote
        _old_blas_threads=LinearAlgebra.BLAS.get_num_threads()
        LinearAlgebra.BLAS.set_num_threads($(esc(n)))
        $(esc(expr))
        LinearAlgebra.BLAS.set_num_threads(_old_blas_threads)    
    end
end

"""
    @blas_1 expr

Temporarily set BLAS threads to `1` for `expr`.
"""
macro blas_1(expr)
    quote
        LinearAlgebra.BLAS.set_num_threads(1)
        $(esc(expr))
    end
end

"""
    @blas_multi_then_1 n expr

Set BLAS threads to `n` for `expr`, then set it to 1 afterward.
"""
macro blas_multi_then_1(n,expr)
    quote
        LinearAlgebra.BLAS.set_num_threads($(esc(n)))
        $(esc(expr))
        LinearAlgebra.BLAS.set_num_threads(1)
    end
end

"""
    svd_or_det_solve(A,use_krylov,which,blas_threads)

Solve for either the smallest singular value or the determinant of A, using either a Krylov method or a direct BLAS call, depending on `use_krylov`. The `which` argument specifies whether to compute the determinant (`:det`) or the smallest singular value (`:svd`). The `blas_threads` argument controls the number of threads used for BLAS operations when not using Krylov methods.
"""
macro svd_or_det_solve(A,use_krylov,which,blas_threads)
    quote
        if $(esc(use_krylov))
            if $(esc(which))===:det
                @blas_1 mu,_,_,_=svdsolve($(esc(A)),1,:SR)
                return mu[1]
            elseif $(esc(which))===:det_argmin
                @blas_1 mu,_,_,_=svdsolve($(esc(A)),1,:SR)
                return mu[1]
            elseif $(esc(which))===:svd
                @blas_1 mu,_,_,_=svdsolve($(esc(A)),1,:SR)
                return mu[1]
            else
                error("Invalid option for `which`. Use :det, :det_argmin, or :svd.")
            end
        else
            if $(esc(which))===:det
                @blas_multi_then_1 $(esc(blas_threads)) d=det($(esc(A)))
                return d
            elseif $(esc(which))===:det_argmin
                @blas_multi_then_1 $(esc(blas_threads)) d=det($(esc(A)))
                return abs(d)
            elseif $(esc(which))===:svd
                @blas_multi_then_1 $(esc(blas_threads)) _,S,_=LAPACK.gesvd!('N','N',$(esc(A)))
                return S[end]
            else
                error("Invalid option for `which`. Use :det, :det_argmin, or :svd.")
            end
        end
    end
end

# helper for determining the byte size of object "a"
memory_size(a)=Base.format_bytes(Base.summarysize(a)) 

"""
    try_MKL!()

Tries to use the MKL on x86_64 architecture if possible. Otherwise it defaults to the stock BLAS backend :lbt.
macOS currently not supported due to Accelerate.jl causing issues with non-even FFTW. MKL.jl does not have this issue (for now).
"""
function try_MKL!()
    if Sys.ARCH==:x86_64
        try
            @eval using MKL
            println(BLAS.get_config())
        catch e
            println(e)
            @warn "Install Math Kernel Library (MKL) via MKL.jl"
            @info "Defaulting to stock BLAS backend: $(BLAS.vendor())"
        end
    else
        @info "Not on x86_64 architecture ($(Sys.ARCH)), defaulting to stock BLAS backend: $(BLAS.vendor())"
    end
end