module ShubertPiyavskii

using DataStructures: BinaryMaxHeap

include("ErrorMessages.jl")
using .ErrorMessages.ShubertPiyavskii

export GlobalMaximizationResult, maximize_shubert


const RealT = Float64
const DEFAULT_TOL = 1e-5 # TODO: Perhaps it's worth making this larger (1e-2 or 1e-3)?


struct GlobalMaximizationResult
    maximizer::RealT
    maximum::RealT
    tol::RealT
end


struct _SearchInterval
    x_lower::RealT
    x_upper::RealT
    y_lower::RealT
    y_upper::RealT
    upper_bound::RealT
    
    function _SearchInterval(
        x_lower::RealT, x_upper::RealT, y_lower::RealT, y_upper::RealT, lipschitz::RealT,
    )
        upper_bound = (lipschitz * (x_upper - x_lower) + y_lower + y_upper) / 2
        return new(x_lower, x_upper, y_lower, y_upper, upper_bound)
    end
end

Base.isless(A::_SearchInterval, B::_SearchInterval) = A.upper_bound < B.upper_bound


# TODO: Add parallelization over subintervals of the initial optimization range
function maximize_shubert(
    objective::Function, x_lower::RealT, x_upper::RealT, lipschitz::RealT;
    tol::RealT=DEFAULT_TOL,
)
    x_lower < x_upper || throw(DomainError([x_lower, x_upper], OPTIM_RANGE_ERR))
    lipschitz > 0 || throw(DomainError(lipschitz, LIPSCHITZ_ERR))
    tol > 0 || throw(DomainError(tol, TOL_ERR))
    
    y_lower = objective(x_lower)
    y_upper = objective(x_upper)
    maximizer, maximum = (y_lower >= y_upper) ? (x_lower, y_lower) : (x_upper, y_upper)
    
    search_intervals = BinaryMaxHeap{_SearchInterval}()
    initial_interval = _SearchInterval(x_lower, x_upper, y_lower, y_upper, lipschitz)
    push!(search_intervals, initial_interval)
    
    while !isempty(search_intervals)
        search_interval = pop!(search_intervals)
        
        if search_interval.upper_bound - maximum >= tol
            child1, child2 = _split_search_interval(search_interval, objective, lipschitz)
            
            for child in (child1, child2)
                if child.upper_bound - maximum >= tol
                    push!(search_intervals, child)
                end
            end
            
            x_sample = child1.x_upper
            y_sample = child1.y_upper
            
            if maximum <= y_sample
                maximizer, maximum = (x_sample, y_sample)
            end
        end
    end
    
    return GlobalMaximizationResult(maximizer, maximum, tol)
end


@inline # TODO: Should benchmark to see if this actually helps
function _split_search_interval(
    search_interval::_SearchInterval, objective::Function, lipschitz::RealT,
)
    x_lower = search_interval.x_lower
    x_upper = search_interval.x_upper
    y_lower = search_interval.y_lower
    y_upper = search_interval.y_upper
    
    x_sample = clamp(
        (x_lower + x_upper + (y_upper - y_lower) / lipschitz) / 2,
        x_lower,
        x_upper,
    )
    y_sample = objective(x_sample)
    
    child1 = _SearchInterval(x_lower, x_sample, y_lower, y_sample, lipschitz)
    child2 = _SearchInterval(x_sample, x_upper, y_sample, y_upper, lipschitz)
    return child1, child2
end

end
