### allelefreq - calculate allele frequencies

# allele = 0, 2 for missing or heterozygous genotype frequencies
#          1, 3 for allele A, B including AB/2, excluding missing
# by     = dimension (1 = per individual, 2 = per SNP)
function allelefreq(x::SnpData, allele::Integer = 1, by::Integer = 2)
    if by < 1 || by > 2
        error("by out of range")
    end
    n = size(x, by)
    m = size(x, 3 - by)
    f = Vector{Real}(n)
    for (i in 1:n)
        if (by == 2)
            s = x[:,i]
        else
            s = x[i,:]
        end
        if allele in [0 2]
            # missing and heterozygous frequencies
            f[i] = length(findin(s, allele)) / m # count or mapreduce take more time
        else
            # diploid allele
            nmis = length(findin(s, 0))
            nall = length(findin(s, allele))
            nhet = length(findin(s, 2))/2
            f[i] = (nall + nhet) / (m - nmis)
        end
    end
    return(f)
end

function allelefreq(x::SnpMatrix, allele::Integer = 1, by::Integer = 2)
    allelefreq(x.data, allele, by)
end


### minor - calculate the minor allele frequency

# x = {Vector with each} value in [0;1], for diploidy
function minor{ T<:Union{Real, Vector{Real}} }(x::T)
    maf = min(x, 1 - x) # min works element wise on arrays
    return(maf)
end
