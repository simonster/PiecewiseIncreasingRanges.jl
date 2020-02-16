module PiecewiseIncreasingRanges
export PiecewiseIncreasingRange, NoNearestSampleError, findnearest, resample

using Base.Order

constructrange(::Type{T}, start, step, stop) where T<:AbstractUnitRange = T(start, stop)
constructrange(::Type{T}, start, step, stop) where T<:StepRange = T(start, step, stop)
constructrange(::Type{T}, start, step, stop) where T<:StepRangeLen = start:step:stop

function combine_ranges(ranges::Vector{R}, firstrg::R, firstrgidx::Int) where R<:AbstractRange
    newranges = R[]
    offsets = Int[1]

    step(firstrg) < 0 && throw(ArgumentError("ranges must be strictly monotonically increasing"))
    curstart = first(firstrg)
    curstop = last(firstrg)
    curstep = step(firstrg)

    for i = firstrgidx:length(ranges)
        newrg = ranges[i]
        isempty(newrg) && continue

        if step(newrg) == curstep && first(newrg) == curstop+curstep
            # Can extend current range
            curstop = last(newrg)
        else
            if first(newrg) <= curstop || step(newrg) < 0
                throw(ArgumentError("ranges must be strictly monotonically increasing"))
            end

            # Need to make a new range
            currg = constructrange(R, curstart, curstep, curstop)
            push!(newranges, currg)
            push!(offsets, offsets[end]+length(currg))

            curstart = first(newrg)
            curstop = last(newrg)
            curstep = step(newrg)
        end
    end

    push!(newranges, constructrange(R, curstart, curstep, curstop))
    (newranges, offsets)
end

struct PiecewiseIncreasingRange{T,R<:AbstractRange,S} <: AbstractVector{T}
    ranges::Vector{R}
    offsets::Vector{Int}
    divisor::S

    function PiecewiseIncreasingRange(ranges::Vector{R}, divisor) where R<:AbstractRange
        ranges = convert(Vector{R}, ranges)
        if isempty(ranges)
            offsets = Int[]
        else
            # Find first non-empty range
            firstrg = ranges[1]
            j = 1
            while j <= length(ranges)
                firstrg = ranges[j]
                j += 1
                !isempty(firstrg) && break
            end
            if isempty(firstrg)
                ranges = R[]
                offsets = Int[]
            else
                ranges, offsets = combine_ranges(ranges, firstrg, j)
            end
        end
        new{eltype(R),R,typeof(divisor)}(ranges, offsets, divisor)
    end
end
PiecewiseIncreasingRange(ranges::Vector{R}) where R<:AbstractRange = PiecewiseIncreasingRange(ranges, nothing)

Base.convert(::Type{PiecewiseIncreasingRange{T,R,S}}, x::PiecewiseIncreasingRange{T,R,S}) where {T,R<:AbstractRange,S} = x
Base.convert(::Type{PiecewiseIncreasingRange{T,R,S}}, x::PiecewiseIncreasingRange) where {T,R,S} =
    PiecewiseIncreasingRange{T,R,S}(x.ranges, x.divisor)

# Avoid applying the divisor if it is one, to get types right
divide_divisor(r::PiecewiseIncreasingRange{T,R,Nothing}, x) where {T,R} = x
multiply_divisor(r::PiecewiseIncreasingRange{T,R,Nothing}, x) where {T,R} = x
divide_divisor(r::PiecewiseIncreasingRange{T,R,S}, x) where {T,R,S} = x/r.divisor
multiply_divisor(r::PiecewiseIncreasingRange{T,R,S}, x) where {T,R,S} = x*r.divisor

function Base.size(r::PiecewiseIncreasingRange)
    isempty(r.ranges) && return (0,)
    return (r.offsets[end]+length(r.ranges[end])-1,)
end

function Base.getindex(r::PiecewiseIncreasingRange, i::Integer)
    rgidx = searchsortedlast(r.offsets, i, Forward)
    divide_divisor(r, r.ranges[rgidx][i-r.offsets[rgidx]+1])
end

function Base.getindex(r::PiecewiseIncreasingRange{T,R,S}, x::AbstractRange{Int}) where {T,R,S}
    isempty(x) && return PiecewiseIncreasingRange{T,R,S}(R[], Int[], r.divisor)
    (first(x) >= 1 && last(x) <= length(r)) || throw(BoundsError())

    # Not type-stable!
    step(x) < 0 && return invoke(getindex, (typeof(r), AbstractVector{Int}), r, x)

    firstrgidx = searchsortedlast(r.offsets, first(x), Forward)
    lastrgidx = searchsortedlast(r.offsets, last(x), Forward)
    newrgs = Array{R}(undef, lastrgidx-firstrgidx+1)
    for irange = firstrgidx:lastrgidx
        elmax = min(last(x)-r.offsets[irange]+1, length(r.ranges[irange]))
        newrg = newrgs[irange-firstrgidx+1] = r.ranges[irange][first(x)-r.offsets[irange]+1:elmax]
        x = x[length(newrg)+1:end]
    end
    PiecewiseIncreasingRange(newrgs, r.divisor)
end

function resample(r::PiecewiseIncreasingRange{T,R,S}, ratio::Rational{Int}) where {T,R,S}
    excess = zero(eltype(R))
    curpt = 0
    newrgs = R[]
    for i = 1:length(r.ranges)
        endpt = last(r.ranges[i])*numerator(ratio)
        newrg = first(r.ranges[i])*numerator(ratio)+excess*step(r.ranges[i]):step(r.ranges[i])*denominator(ratio):endpt
        curpt += length(newrg)
        push!(newrgs, newrg)
        if i != length(r.ranges)
            # Need to linearly interpolate between endpoint and next range
            # nextind is the next index in r (= 1 + curpt/ratio) times numerator(ratio)
            nextind = numerator(ratio) + curpt*denominator(ratio)
            if nextind < r.offsets[i+1]*numerator(ratio)
                # Need to interpolate between this point and the next
                weight = mod(nextind, numerator(ratio))
                interpt = last(r.ranges[i])*(numerator(ratio)-weight) + first(r.ranges[i+1])*weight
                interpstep = (first(r.ranges[i+1]) - last(r.ranges[i]))*denominator(ratio)
                newrg = interpt:interpstep:first(r.ranges[i+1])*numerator(ratio)
                push!(newrgs, newrg)
                curpt += length(newrg)
                nextind += length(newrg)*denominator(ratio)
            end
            excess = nextind - r.offsets[i+1]*numerator(ratio)
        end
    end
    PiecewiseIncreasingRange(newrgs, multiply_divisor(r, numerator(ratio)))
end

# searchsortedfirst, searchsortedlast
struct PiecewiseIncreasingRangeFirstOrdering <: Ordering end
Base.Order.lt(o::PiecewiseIncreasingRangeFirstOrdering, a, b) = isless(first(a), first(b))
struct PiecewiseIncreasingRangeLastOrdering <: Ordering end
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

struct NoNearestSampleError <: Exception end

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
