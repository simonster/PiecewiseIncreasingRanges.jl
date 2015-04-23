# PiecewiseIncreasingRanges

[![Build Status](https://travis-ci.org/simonster/PiecewiseIncreasingRanges.jl.svg?branch=master)](https://travis-ci.org/simonster/PiecewiseIncreasingRanges.jl)
[![Coverage Status](https://coveralls.io/repos/simonster/PiecewiseIncreasingRanges.jl/badge.svg?branch=master)](https://coveralls.io/r/simonster/PiecewiseIncreasingRanges.jl?branch=master)

PiecewiseIncreasingRanges.jl provides a PiecewiseIncreasingRange type that corresponds to a set of concatenated, monotonically increasing ranges. It implements indexing as well as `searchsortedfirst`, `searchsortedlast`, and a `findnearest` function. This type is useful for holding the potentially discontinuous timestamps corresponding to large quantities of sampled data, where creating and indexing a vector would be undesirably slow.

## Construction

You can construct a PiecewiseIncreasingRange as below:

```julia
PiecewiseIncreasingRange(StepRange{Rational{Int}, Rational{Int}}[0:1//20000:1, 5:1//20000:6])
```

At present, only StepRanges are supported.

## findnearest

`findnearest(rg, x, within_half_sample)` finds the index of the element of `rg` closest to `x`. If `x` is equidistant between two samples, it chooses the earlier one. If `within_half_sample` is true, `findnearest` throws a `NoNearestSampleError` if there is no sample within half of a step of `x`.