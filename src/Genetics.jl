module Genetics

    import Base.size

    include("types.jl")
    include("io.jl")
    include("attributes.jl")
    include("stats.jl")

    export

        # types
        SnpMatrix,
        LocusInfo,

        # io functions
        read_vcf,

        # object attribute functions
        size,

        # stats functions
        allelefreq,
        minor

end # module
