module IO

using ..Parameters: OUTDIR
using Printf

export save_table_tsv, save_all_tsv

function format_tsv_value(value)
    if value === missing || (value isa AbstractFloat && !isfinite(value))
        return "NA"
    elseif value isa AbstractFloat
        return @sprintf("%.10g", value)
    end
    return replace(string(value), '\t' => ' ', '\n' => ' ', '\r' => ' ')
end

function save_table_tsv(path::String, rows)
    isempty(rows) && error("Cannot save an empty table to $path")
    headers = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(headers), '\t'))
        for row in rows
            println(io, join((format_tsv_value(getproperty(row, header)) for header in headers), '\t'))
        end
    end
    return path
end

function save_all_tsv(results)
    degree_dir = joinpath(OUTDIR, "degreeResults")
    community_dir = joinpath(OUTDIR, "communityMetrics")
    mkpath(degree_dir)
    mkpath(community_dir)

    paths = (
        consumer_results=save_table_tsv(
            joinpath(degree_dir, "consumer_results.tsv"),
            results.consumer_results
        ),
        degree_replicate=save_table_tsv(
            joinpath(degree_dir, "degree_replicate_summary.tsv"),
            results.degree_replicate
        ),
        degree_summary=save_table_tsv(
            joinpath(degree_dir, "degree_summary.tsv"),
            results.degree_summary
        ),
        community_results=save_table_tsv(
            joinpath(community_dir, "community_results.tsv"),
            results.community_results
        )
    )
    return paths
end

end
