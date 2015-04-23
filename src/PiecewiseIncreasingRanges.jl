module PiecewiseIncreasingRanges
export PiecewiseIncreasingRange, NoNearestSampleError, findnearest

using Base.Order

immutable PiecewiseIncreasingRange{T} <: AbstractVector{T}
    ranges::Vector{StepRange{T,T}}
    offsets::Vector{Int}

    function PiecewiseIncreasingRange(ranges::Vector{StepRange{T,T}})
        isempty(ranges) && return new(ranges, Int[])
        newranges = Range{T}[]
        offsets = Int[1]

        # Find first non-empty range
        currg = ranges[1]
        for i = 1:length(ranges)
            currg = ranges[i]
            !isempty(currg) && break
        end
        isempty(currg) && return PiecewiseIncreasingRange(newranges, offsets)
        step(currg) < 0 && throw(ArgumentError("ranges must be strictly monotonically increasing"))

        for i = 2:length(ranges)
            newrg = ranges[i]
            isempty(newrg) && continue

            if first(newrg) <= last(currg) || step(newrg) < 0
                throw(ArgumentError("ranges must be strictly monotonically increasing"))
            end

            if step(newrg) == step(currg) && first(newrg) == last(currg)+step(currg)
                # Can extend current range
                currg = first(currg):step(currg):last(newrg)
            else
                # Need to make a new range
                push!(newranges, currg)
                push!(offsets, offsets[end]+length(currg))
                currg = newrg
            end
        end
        push!(newranges, currg)

        new(newranges, offsets)
    end
end
PiecewiseIncreasingRange{T}(ranges::Vector{StepRange{T,T}}) = PiecewiseIncreasingRange{T}(ranges)

function Base.size(r::PiecewiseIncreasingRange)
    isempty(r.ranges) && return 0
    return (r.offsets[end]+length(r.ranges[end])-1,)
end

function Base.getindex(r::PiecewiseIncreasingRange, i::Integer)
    rgidx = searchsortedlast(r.offsets, i, Forward)
    r.ranges[rgidx][i-r.offsets[rgidx]+1]
end

# searchsortedfirst, searchsortedlast
immutable PiecewiseIncreasingRangeFirstOrdering <: Ordering end
Base.Order.lt(o::PiecewiseIncreasingRangeFirstOrdering, a, b) = isless(first(a), first(b))
immutable PiecewiseIncreasingRangeLastOrdering <: Ordering end
Base.Order.lt(o::PiecewiseIncreasingRangeLastOrdering, a, b) = isless(last(a), last(b))

function Base.searchsortedfirst(r::PiecewiseIncreasingRange, x)
    isempty(r.ranges) && return 1
    rgidx = searchsortedfirst(r.ranges, x, PiecewiseIncreasingRangeLastOrdering())
    rgidx > length(r.ranges) && return length(r) + 1
    searchsortedfirst(r.ranges[rgidx], x, Forward) + r.offsets[rgidx] - 1
end

function Base.searchsortedlast(r::PiecewiseIncreasingRange, x)
    isempty(r.ranges) && return 1
    rgidx = searchsortedlast(r.ranges, x, PiecewiseIncreasingRangeFirstOrdering())
    rgidx == 0 && return 0
    searchsortedlast(r.ranges[rgidx], x, Forward) + r.offsets[rgidx] - 1
end

immutable NoNearestSampleError <: Exception end

function findnearest(r::PiecewiseIncreasingRange, x, within_half_step::Bool=false)
    isempty(r.ranges) && throw(NoNearestSampleError())

    rgidx = searchsortedfirst(r.ranges, x, PiecewiseIncreasingRangeLastOrdering())
    if rgidx > length(r.ranges)
        rgend = r.ranges[end]
        within_half_step && x > rgend[end]+step(rgend)/2 && throw(NoNearestSampleError())
        return length(r)
    end

    rg = r.ranges[rgidx]
    idxinrg = searchsortedfirst(rg, x, Forward)
    idx = idxinrg + r.offsets[rgidx] - 1
    d = rg[idxinrg] - x

    if idxinrg == 1
        if rgidx == 1
            # First element of all elements
            within_half_step && d > step(rg)/2 && throw(NoNearestSampleError())
            return 1
        end

        # Could be closer to last element of preceding range
        rgprev = r.ranges[rgidx-1]
        dprev = x - rgprev[end]
        within_half_step && d > step(rg)/2 && dprev > step(rgprev)/2 && throw(NoNearestSampleError())
        return idx - (dprev <= d)
    end

    ifelse(d >= step(rg)/2, idx-1, idx)
end

end # module
