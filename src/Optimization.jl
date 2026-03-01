module Optimization

export local_fminbnd

"""
    local_fminbnd(fun, ax, cx, tol=0.001)

Finds the local minimum of function `fun` in interval `[ax, cx]` using Brent's method 
(Golden Section Search + Parabolic Interpolation).
Ported from standard numerical recipes / verified MATLAB implementation.
"""
function local_fminbnd(fun::Function, ax::Real, cx::Real, tol::Real=0.01)
    # Phi is the golden ratio conjugate
    phi = (3 - sqrt(5)) / 2
    
    # Initialize
    a = min(ax, cx)
    b = max(ax, cx)
    v = a + phi * (b - a)
    w = v
    x = v
    
    fv = fun(v)
    fw = fv
    fx = fv
    
    # Golden section step size
    d = 0.0
    e = 0.0
    
    iter = 0
    max_iter = 100
    
    while iter < max_iter
        iter += 1
        
        xm = 0.5 * (a + b)
        tol1 = tol * abs(x) + 1e-10
        tol2 = 2.0 * tol1
        
        # Check for convergence
        if abs(x - xm) <= (tol2 - 0.5 * (b - a))
            break
        end
        
        if abs(e) > tol1
            # Attempt parabolic interpolation
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = 2.0 * (q - r)
            if q > 0
                p = -p
            end
            q = abs(q)
            etemp = e
            e = d
            
            # Check if parabolic step is acceptable
            if abs(p) >= abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)
                # Reject parabolic, use Golden Section
                if x >= xm
                    e = a - x
                else
                    e = b - x
                end
                d = phi * e
            else
                # Accept parabolic step
                d = p / q
                u = x + d
                if (u - a) < tol2 || (b - u) < tol2
                   d = sign(xm - x) * tol1
                end
            end
        else
            # Golden Section step
            if x >= xm
                e = a - x
            else
                e = b - x
            end
            d = phi * e
        end
        
        # Perform step
        if abs(d) >= tol1
            u = x + d
        else
            u = x + sign(d) * tol1
        end
        
        fu = fun(u)
        
        # Update book-keeping variables
        if fu <= fx
            if u >= x
                a = x
            else
                b = x
            end
            v = w; fv = fw
            w = x; fw = fx
            x = u; fx = fu
        else
            if u < x
                a = u
            else
                b = u
            end
            if fu <= fw || w == x
                v = w; fv = fw
                w = u; fw = fu
            elseif fu <= fv || v == x || v == w
                v = u; fv = fu
            end
        end
    end
    
    return x, fx
end

end # module
