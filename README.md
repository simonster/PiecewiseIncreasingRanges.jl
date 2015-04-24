# PiecewiseIncreasingRanges

[![Build Status](https://travis-ci.org/simonster/PiecewiseIncreasingRanges.jl.svg?branch=master)](https://travis-ci.org/simonster/PiecewiseIncreasingRanges.jl)
[![Coverage Status](https://coveralls.io/repos/simonster/PiecewiseIncreasingRanges.jl/badge.svg?branch=master)](https://coveralls.io/r/simonster/PiecewiseIncreasingRanges.jl?branch=master)

PiecewiseIncreasingRanges.jl provides a PiecewiseIncreasingRange type that corresponds to a set of concatenated, monotonically increasing ranges. It implements indexing as well as `searchsortedfirst`, `searchsortedlast`, and a `findnearest` function. This type is useful for holding the potentially discontinuous timestamps corresponding to large quantities of sampled data, where creating and indexing a vector would be undesirably slow.

## Construction

You can construct a PiecewiseIncreasingRange as below:

```julia
julia> PiecewiseIncreasingRange(StepRange{Int,Int}[0:5:15, 18:2:20])
6-element PiecewiseIncreasingRanges.PiecewiseIncreasingRange{Int64,StepRange{Int64,Int64},Void}:
  0
  5
 10
 15
 18
 20
```

PiecewiseIncreasingRanges also accept an optional divisor, which is divided out of all elements:

```julia
julia> PiecewiseIncreasingRange(UnitRange{Int}[0:3, 15:16], 4)
6-element PiecewiseIncreasingRanges.PiecewiseIncreasingRange{Float64,UnitRange{Int64},Int64}:
 0.0 
 0.25
 0.5 
 0.75
 3.75
 4.0
```

## findnearest

`findnearest(rg, x, within_half_sample)` finds the index of the element of `rg` closest to `x`. If `x` is equidistant between two samples, it chooses the earlier one. If `within_half_sample` is true, `findnearest` throws a `NoNearestSampleError` if there is no sample within half of a step of `x`.