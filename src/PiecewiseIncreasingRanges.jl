module PiecewiseIncreasingRanges
export PiecewiseIncreasingRange, NoNearestSampleError, findnearest

using Base.Order

function combine_ranges{R<:Range}(ranges::Vector{R}, firstrg::R, firstrgidx::Int)
    newranges = R[]
    offsets = Int[1]

    step(firstrg) < 0 && throw(ArgumentError("ranges must be strictly monotonically increasing"))
    currg = firstrg

    for i = firstrgidx:length(ranges)
        newrg = ranges[i]
        isempty(newrg) && continue

        if step(newrg) == step(currg) && first(newrg) == last(currg)+step(currg)
            # Can extend current range
            currg = first(currg):step(currg):last(newrg)
        else
            if first(newrg) <= last(currg) || step(newrg) < 0
                throw(ArgumentError("ranges must be strictly monotonically increasing"))
            end

            # Need to make a new range
            push!(newranges, currg)
            push!(offsets, offsets[end]+length(currg))
            currg = newrg
        end
    end

    push!(newranges, currg)
    (newranges, offsets)
end

function combine_ranges{R<:FloatRange}(ranges::Vector{R}, firstrg::R, firstrgidx::Int)
    newranges = R[]
    offsets = Int[1]

    if signbit(firstrg.step) != signbit(firstrg.divisor)
        throw(ArgumentError("ranges must be strictly monotonically increasing"))
    end
    currg = firstrg

    for i = firstrgidx:length(ranges)
        newrg = ranges[i]
        isempty(newrg) && continue

        if step(newrg) == step(currg) && newrg.start*currg.divisor == (currg.start+currg.step*currg.len)*newrg.divisor
            currg = FloatRange(currg.start, currg.step, currg.len + newrg.len, currg.divisor)
        else
            if first(newrg) <= last(currg) || signbit(newrg.step) != signbit(newrg.divisor)
                throw(ArgumentError("ranges must be strictly monotonically increasing"))
            end

            # Need to make a new range
            push!(newranges, currg)
            push!(offsets, offsets[end]+length(currg))
            currg = newrg
        end
    end

    push!(newranges, currg)
    (newranges, offsets)
end

immutable PiecewiseIncreasingRange{T,R<:Range,S} <: AbstractVector{T}
    ranges::Vector{R}
    offsets::Vector{Int}
    divisor::S

    function PiecewiseIncreasingRange(ranges::Vector{R}, divisor)
        isempty(ranges) && return new(ranges, Int[])

        # Find first non-empty range
        firstrg = ranges[1]
        j = 0
        for j = 1:length(ranges)
            firstrg = ranges[j]
            !isempty(firstrg) && break
        end
        isempty(firstrg) && return new(R[], Int[], divisor)

        newranges, offsets = combine_ranges(ranges, firstrg, j+1)
        new(newranges, offsets, divisor)
    end
end
PiecewiseIncreasingRange{R<:Range}(ranges::Vector{R}, divisor) = PiecewiseIncreasingRange{typeof(inv(one(eltype(R)))),R,typeof(divisor)}(ranges, divisor)
PiecewiseIncreasingRange{R<:Range}(ranges::Vector{R}) = PiecewiseIncreasingRange{eltype(R),R,Nothing}(ranges, nothing)

# Avoid applying the divisor if it is one, to get types right
divide_divisor{T,R}(r::PiecewiseIncreasingRange{T,R,Nothing}, x) = x
multiply_divisor{T,R}(r::PiecewiseIncreasingRange{T,R,Nothing}, x) = x
divide_divisor{T,R,S}(r::PiecewiseIncreasingRange{T,R,S}, x) = x/r.divisor
multiply_divisor{T,R,S}(r::PiecewiseIncreasingRange{T,R,S}, x) = x*r.divisor

function Base.size(r::PiecewiseIncreasingRange)
    isempty(r.ranges) && return 0
    return (r.offsets[end]+length(r.ranges[end])-1,)
end

function Base.getindex(r::PiecewiseIncreasingRange, i::Integer)
    rgidx = searchsortedlast(r.offsets, i, Forward)
    divide_divisor(r, r.ranges[rgidx][i-r.offsets[rgidx]+1])
end

function Base.getindex{T,R,S}(r::PiecewiseIncreasingRange{T,R,S}, x::Range{Int})
    isempty(x) && return PiecewiseIncreasingRange{T,R,S}(R[], Int[], r.divisor)
    (first(x) >= 1 && last(x) <= length(r)) || throw(BoundsError())

    # Not type-stable!
    step(x) < 0 && return invoke(getindex, (typeof(r), AbstractVector{Int}), r, x)

    firstrgidx = searchsortedlast(r.offsets, first(x), Forward)
    lastrgidx = searchsortedlast(r.offsets, last(x), Forward)
    newrgs = Array(R, lastrgidx-firstrgidx+1)
    for irange = firstrgidx:lastrgidx
        elmax = min(last(x)-r.offsets[irange]+1, length(r.ranges[irange]))
        newrg = newrgs[irange-firstrgidx+1] = r.ranges[irange][first(x)-r.offsets[irange]+1:elmax]
        x = x[length(newrg)+1:end]
    end
    PiecewiseIncreasingRange(newrgs, r.divisor)
end

# searchsortedfirst, searchsortedlast
immutable PiecewiseIncreasingRangeFirstOrdering <: Ordering end
Base.Order.lt(o::PiecewiseIncreasingRangeFirstOrdering, a, b) = isless(first(a), first(b))
immutable PiecewiseIncreasingRangeLastOrdering <: Ordering end
Base.Order.lt(o::PiecewiseIncreasingRangeLastOrdering, a, b) = isless(last(a), last(b))

function Base.searchsortedfirst(r::PiecewiseIncreasingRange, x)
    isempty(r.ranges) && return 1
    xd = multiply_divisor(r, x)

    rgidx = searchsortedfirst(r.ranges, xd, PiecewiseIncreasingRangeLastOrdering())
    rgidx > length(r.ranges) && return length(r) + 1
    searchsortedfirst(r.ranges[rgidx], xd, Forward) + r.offsets[rgidx] - 1
end

function Base.searchsortedlast(r::PiecewiseIncreasingRange, x)
    isempty(r.ranges) && return 1
    xd = multiply_divisor(r, x)

    rgidx = searchsortedlast(r.ranges, xd, PiecewiseIncreasingRangeFirstOrdering())
    rgidx == 0 && return 0
    searchsortedlast(r.ranges[rgidx], xd, Forward) + r.offsets[rgidx] - 1
end

immutable NoNearestSampleError <: Exception end

function findnearest(r::PiecewiseIncreasingRange, x, within_half_step::Bool=false)
    isempty(r.ranges) && throw(NoNearestSampleError())
    xd = multiply_divisor(r, x)

    rgidx = searchsortedfirst(r.ranges, xd, PiecewiseIncreasingRangeLastOrdering())
    if rgidx > length(r.ranges)
        rgend = r.ranges[end]
        within_half_step && xd > rgend[end]+step(rgend)/2 && throw(NoNearestSampleError())
        return length(r)
    end

    rg = r.ranges[rgidx]
    idxinrg = searchsortedfirst(rg, xd, Forward)
    idx = idxinrg + r.offsets[rgidx] - 1
    d = rg[idxinrg] - xd

    if idxinrg == 1
        if rgidx == 1
            # First element of all elements
            within_half_step && d > step(rg)/2 && throw(NoNearestSampleError())
            return 1
        end

        # Could be closer to last element of preceding range
        rgprev = r.ranges[rgidx-1]
        dprev = xd - rgprev[end]
        within_half_step && d > step(rg)/2 && dprev > step(rgprev)/2 && throw(NoNearestSampleError())
        return idx - (dprev <= d)
    end

    ifelse(d >= step(rg)/2, idx-1, idx)
end

end # module
