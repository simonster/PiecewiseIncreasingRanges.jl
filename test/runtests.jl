using PiecewiseIncreasingRanges, Test, Interpolations

function test(rgs, divisor...)
    vcrgs = vcat(rgs...)
    !isempty(divisor) && (vcrgs /= divisor[1])
    rg = PiecewiseIncreasingRange(rgs, divisor...)
    @test length(rg.ranges) == 4
    @test vcrgs == rg

    # yi = InterpGrid(convert(Vector{Float64}, vcrgs), BCnil, InterpLinear)
    itp = interpolate(convert(Vector{Float64}, vcrgs), BSpline(Linear()))
    @test resample(rg, 3//7) ≈ itp(1:7//3:length(vcrgs))
    @test resample(rg, 7//3) ≈ itp(1:3//7:length(vcrgs))
    @test resample(rg, 5//9) ≈ itp(1:9//5:length(vcrgs))
    @test resample(rg, 61//3) ≈ itp(1:3//61:length(vcrgs))

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

        for j = i:length(rg)
            @test vcrgs[i:j] == rg[i:j]
            for k = 1:i-j+1
                @test vcrgs[i:k:j] == rg[i:k:j]
            end
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
end

test([0:4:40, 41:1:80, 82:2:112, 114:2:120, 136:4:144], 8)
test([0:1//2:5, 5+1//8:1//8:10, 10+1//4:1//4:14, 14+1//4:1//4:15, 17:1//2:18])
test([0:1//2:5., 5+1//8:1//8:10., 10+1//4:1//4:14., 14+1//4:1//4:15., 17:1//2:18.])

# Empty test
rg = PiecewiseIncreasingRange(StepRange{Rational{Int},Rational{Int}}[])
@test searchsortedfirst(rg, 1) == 1
@test searchsortedlast(rg, 1) == 1
@test_throws NoNearestSampleError findnearest(rg, 1)
rg = PiecewiseIncreasingRange(UnitRange{Int}[0:-1])
@test isempty(rg)

# Test with non-monotonic range
@test_throws ArgumentError PiecewiseIncreasingRange([0:1//2:1, 3//4:1//2:4])
@test_throws ArgumentError PiecewiseIncreasingRange([0:-1//2:-1, 10:1//2:8])
@test_throws ArgumentError PiecewiseIncreasingRange([0:1//2:1, 10:-1//2:8])
@test_throws ArgumentError PiecewiseIncreasingRange([0:0.5:1, 0.75:0.5:4])
@test_throws ArgumentError PiecewiseIncreasingRange([0:-0.5:-1, 10:0.5:8])
@test_throws ArgumentError PiecewiseIncreasingRange([0:0.5:1, 10:-0.5:8])

# Test with empty ranges interspersed with non-empty ranges
rg = PiecewiseIncreasingRange(UnitRange{Int}[0:-1, 1:0, 0:-1, 1:3, 0:-1, 5:10])
@test rg.ranges == UnitRange{Int}[1:3, 5:10]
