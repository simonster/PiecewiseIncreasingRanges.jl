using PiecewiseIncreasingRanges
using Base.Test

rgs = StepRange{Rational{Int},Rational{Int}}[0:1//2:5, 5+1//8:1//8:10, 10+1//4:1//4:14, 14+1//4:1//4:15, 17:1//2:18]
vcrgs = vcat(rgs...)
rg = PiecewiseIncreasingRange(rgs)
@test length(rg.ranges) == 4
@test vcrgs == rg

for i = 1:length(rg)
	@test searchsortedfirst(rg, rg[i]) == i
	@test searchsortedfirst(rg, rg[i]-1//17) == i
	@test searchsortedfirst(rg, rg[i]+1//17) == i+1

	@test searchsortedlast(rg, rg[i]) == i
	@test searchsortedlast(rg, rg[i]-1//17) == i-1
	@test searchsortedlast(rg, rg[i]+1//17) == i

	for within_half_sample in (false, true)
		@test findnearest(rg, rg[i], within_half_sample) == i
		@test findnearest(rg, rg[i]+1//16, within_half_sample) == i
		@test findnearest(rg, rg[i]-1//17, within_half_sample) == i
		@test findnearest(rg, rg[i]+1//17, within_half_sample) == i
	end
end

@test findnearest(rg, -1) == 1
@test_throws NoNearestSampleError findnearest(rg, -1, true)
@test findnearest(rg, 15.9) == searchsortedfirst(rg, 15)
@test findnearest(rg, 16.1) == searchsortedfirst(rg, 17)
@test_throws NoNearestSampleError findnearest(rg, 15.9, true)
@test_throws NoNearestSampleError findnearest(rg, 16.1, true)
@test findnearest(rg, 19) == length(rg)
@test_throws NoNearestSampleError findnearest(rg, 19, true)

# Empty test
rg = PiecewiseIncreasingRange(StepRange{Rational{Int},Rational{Int}}[])
@test searchsortedfirst(rg, 1) == 1
@test searchsortedlast(rg, 1) == 1
@test_throws NoNearestSampleError findnearest(rg, 1)

@test_throws ArgumentError PiecewiseIncreasingRange(StepRange{Rational{Int},Rational{Int}}[0:1//2:1, 3//4:1//2:4])
