function sampleVariance(x, n, df, scale)
    return (dot(x,x) + df*scale)/rand(Chisq(n+df))
end