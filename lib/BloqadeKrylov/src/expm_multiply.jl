## loop up table for now, move to const stack later
struct θ_table
    ms::Vector{Int}
    θs::Vector{Float64}
    max_m::Int
    function θ_table()
        ms = Vector{Int}([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
            17,18,19,20,21,22,23,24,25,26,27,28,29,30,
            35,40,45,50,55])
        ts = Vector{Float64}([2.29e-16,
                            2.58e-8,
                            1.39e-5,
                            3.40e-4,
                            2.40e-3,
                            9.07e-3,
                            2.38e-2,
                            5.00e-2,
                            8.96e-2,
                            1.44e-1,
                            2.14e-1,
                            3.00e-1,
                            4.00e-1,
                            5.14e-1,
                            6.41e-1,
                            7.81e-1,
                            9.31e-1,
                            1.09,
                            1.26,
                            1.44,
                            1.62,
                            1.82,
                            2.01,
                            2.22,
                            2.43,
                            2.64,
                            2.86,
                            3.08,
                            3.31,
                            3.54,
                            4.7,
                            6.0,
                            7.2,
                            8.5,
                            9.9])


        return new(ms,ts,55)
    end
end
function get_theta(θ_tbl::θ_table, m::Int)
    if !(m in θ_tbl.ms)
        error("m is not in the table")
    end
    return θ_tbl.θs[indexin(m, θ_tbl.ms)[1]]
end

function get_optimal_sm(t::Number, A::AbstractMatrix, m_max::Int=55, ell::Int=2)
    """
        given an matrix A and t, find the optimal s and m_star for expm_multiply.
        [return] m_star, s, μ, A_1norm, As
        mu is the shifted constant μ
        A_1norm is the 1-norm of A
        As is the shifted A, As = A - mu*I
    """
     # 1) shift the A:
    traceA = LinearAlgebra.tr(A)
    n = size(A,1)
    μ = traceA/size(A,1) 
    
    As = A - μ*LinearAlgebra.I(n)
    A_1norm = opnorm(As,1)

    m_star::Int = 1
    s::Int = 1
    if t*A_1norm != 0
        m_star, s = _calc_optimal_sm(t*As, t*A_1norm, m_max, ell)
    end

    return m_star, s, μ, A_1norm, As
end

function expm_multiply(t::Number,
                       A,
                       v::AbstractVector{T},
                       tol = nothing) where {T}
    v_prom = similar(v, promote_type(eltype(A), T, typeof(t)))
    copyto!(v_prom, v)

    out = similar(v_prom)
    expm_multiply!(out, t, A, v_prom, tol)
    return out
end


function expm_multiply!(w::AbstractVector{T}, 
                       t::Number,
                       A,
                       v::AbstractVector{T},
                       tol = nothing
                       ) where {T}
    if size(v, 1) != size(A, 2)
        error("dimension mismatch")
    end


    if tol === nothing
        # infer tol base on the size of T
        size_T = sizeof(real(T))
        if size_T == 8
            tol = 2^(-53)
        elseif size_T == 4
            tol = 2^(-24)
        else
            error("Unsupported type T of vec{T} for auto set tol")
        end
    end
    
    # 1) get the optimal s and m_star, and As--shifted A
    m_star, s, μ, A_1norm, As = get_optimal_sm(t,A)

    println(m_star, " " , s, "endl")
    # 3) here lies the impl call:
    _expm_multiply_impl!(w, t, As, v, μ, s, m_star, tol)
    return w 
end 
                        
# n0 = 1
function _condition_3_13(θ_tbl, tA_1norm, m_max::Int, ell::Int)
    ## compute p_max:
    sqrt_m_max = sqrt(m_max)
    p_low = Int(floor(sqrt_m_max))
    p_high = Int(ceil(sqrt_m_max+1))
    p_max = maximum([p for p in p_low:1:p_high if p*(p-1) <= m_max + 1]) 

    a = 2*ell * p_max * (p_max+3)
    b = get_theta(θ_tbl,m_max) / float(m_max)

    return tA_1norm <= a*b, p_max

end


function calc_d(A,p::Int,ell::Int)
    #@warn "using explicit onenormest is slow!"
    exp = onenormest_explicit(A,p)
    println(exp)
    est = onenormest(A,p)
    println(est)
    return est^(1.0/p)
end 

function _calc_optimal_sm(As, tA_1norm, m_max::Int=55, ell::Int=2) 
    θ_tbl = θ_table()
    if m_max > θ_tbl.max_m
        error("m_max > ",θ_tbl.max_m, " is not supported")
    end

    m_star = nothing
    s_star = nothing

    flag, p_max = _condition_3_13(θ_tbl, tA_1norm, m_max, ell)
    if flag
        for (i, θ) in enumerate(θ_tbl.θs)
            m = θ_tbl.ms[i]
            s = Int(ceil(tA_1norm/θ))
            #println(tA_1norm, " ", θ, " ", m*s)
            if m_star === nothing
                m_star = m
                s_star = s
            elseif m * s < m_star*s_star
                m_star = m
                s_star = s
            end
        end
    else
        # eq(3.11)
        d = calc_d(As, 2, ell)
        for p in 2:p_max
            # compute d_p+1
            d1 = calc_d(As, p+1, ell)
            αp = max(d,d1)
            d = d1
            for m in p*(p-1):m_max
                if m in θ_tbl.ms
                    s = Int(ceil(αp / get_theta(θ_tbl,m))) #compute_cost_div_m
                    if m_star === nothing
                        m_star = m
                        s_star = s
                    elseif m * s < m_star*s_star
                        m_star = m
                        s_star = s
                    end
                end
            end
        end
        s_star = max(s_star,1)  
    end   
    return m_star, s_star

end


"""
    _expm_multiply_impl!(w::AbstractVector{T}
                         t::Number,
                         As, 
                         v::AbstractVector{T},
                         mu::Real, #shift parameter
                         s::Int,
                         m_star::Int,
                         tol::Real,
                        ) where {T}
    Calculate matrix exponential acting on some vector, ``w = e^{tA}v``,
    using the Krylov subspace approximation.

    [Note] As is a sugeried matrix, i.e. As= A - mu*I

"""
function _expm_multiply_impl!(w::AbstractVector{T},
                             t::Number,
                             As, 
                             v::AbstractVector{T},
                             mu::Number, #shift parameter
                             s::Int,
                             m_star::Int,
                             tol::Real,
                            ) where {T}
    F = deepcopy(v)
    vs = deepcopy(v)
    η = exp(t*mu/s)
    for i in 1:s
        c1 = norm(vs,Inf)
        for j in 1:m_star
            coef = t / (s*j)
            mul!(vs, As, coef*vs)
            c2 = norm(vs,Inf)
            F.+=vs
            if c1+c2 <= tol*norm(F,Inf)
                break
            end
            c1 = c2
        end
        F.*=η
        copyto!(vs,F)
    end 
    copy!(w,F)
    return w
end