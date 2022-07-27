

function write_phenotype_batch(phenotypes, batch, batch_start, batch_end)
    batch_file = batch_filename(batch)
    open(batch_file, "w") do io
        for i in batch_start:batch_end
            println(io, phenotypes[i])
        end
    end
end

batch_filename(batch) = string("phenotypes_batch_", batch, ".csv")

function tmle_phenotypes_batches(parsed_args)
    phenotypes = filter(
        !=("SAMPLE_ID"), 
        split(readline(open(parsed_args["phenotypes-file"])), ",")
    )
    n_phenotypes = size(phenotypes, 1)
    batch_size = parsed_args["batch-size"]
    batch_size = batch_size == "max" ? size(phenotypes, 1) : parse(Int, batch_size)
    
    batch = 1
    batch_start = 1
    batch_end = batch_start + batch_size - 1
    while batch_end <= n_phenotypes
        write_phenotype_batch(phenotypes, batch, batch_start, batch_end)
        batch += 1
        batch_start = batch_end + 1
        batch_end = batch_start + batch_size - 1
    end

    if batch_start <= n_phenotypes
        write_phenotype_batch(phenotypes, batch, batch_start, n_phenotypes)
    end

end