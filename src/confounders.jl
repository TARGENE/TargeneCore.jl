# ADAPTING, MERGING AND FILTERING PRIOR TO PCA
#####################################################################


merge_beds(input_prefix, output) = merge_plink(input_prefix; des = output)

issnp(x::String) = length(x) == 1

all_batches_ok(row, batchcols) = all(x == 1 for x in row[batchcols])

bounds(pos, lb, ub) =  (min(lb, pos - 10^4), max(ub, pos + 10^4))


function notin_ldblocks(row, ldblocks)
    if (chr=row.chromosome,) ∉ keys(ldblocks)
        return true
    else
        group = ldblocks[(chr=row.chromosome,)]
        inblock = any(b.lower_bound < row.position < b.upper_bound for b in eachrow(group))
        return ~inblock
    end
end

ld_blocks_filter(data, ld_blocks_file::Nothing) = data

function ld_blocks_filter(data, ld_blocks_file)
    ld_blocks = CSV.File(
        ld_blocks_file;
        header=["rsid", "chr", "pos", "LDblock_lower", "LDblock_upper", "LDblock_length", "lower_bound", "upper_bound"]) |> DataFrame
    ld_blocks.chr = string.(ld_blocks.chr)
    transform!(ld_blocks,
        :,
        [:pos, :lower_bound, :upper_bound] => ByRow(bounds) => [:lower_bound, :upper_bound],
        )
    ld_blocks = groupby(ld_blocks, :chr)
    return filter(x -> notin_ldblocks(x, ld_blocks), data)
end

ukb_qc_filter(data, qcfile::Nothing) = data

function ukb_qc_filter(data, qcfile)
    qc_df = CSV.File(qcfile) |> DataFrame
    fully_genotyped_snps = innerjoin(
        data, 
        qc_df, 
        on = :snpid => :rs_id,
        makeunique = true
    )
    # Assayed in both genotyping arrays
    return filter(:array => ==(2), fully_genotyped_snps)
end


"""
    filter_chromosome(input, output, traits_file; 
    maf_threshold=0.01, 
    ld_block_file=nothing, 
    qc_file=nothing)

The purpose of this method is to filter SNPS and Individuals before applying a 
dimensionality reduction technique such as PCA that is later used
to extract population stratification compomemts.
We filter SNPs using quality control metrics from the following resource:
    - https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=1955
"""
function filter_chromosome(input, output, traits_file; 
    maf_threshold=0.01, 
    ld_block_file=nothing, 
    qc_file=nothing
    )
    snp_data = SnpData(input)

    # Remove SNP's with MAF < maf-threshold
    maf_threshold = maf_threshold
    snp_data.snp_info[!, "MAF"] = SnpArrays.maf(snp_data.snparray)
    mafpassed = filter(:MAF => >=(maf_threshold), snp_data.snp_info)
    
    # Remove LD regions specified by ld_blocks
    ld_pruned = ld_blocks_filter(mafpassed, ld_block_file)

    # The QC file contains information on fully genotyped SNPS
    # We only keep those
    qced = ukb_qc_filter(ld_pruned, qc_file)

    # If an RSID appears multiple times, it is because it has 
    # more than 2 possible alleles: we remove them 
    # (why? maybe because the PCA then cannot tackle them)
    duplicate_rsids = Set(qced.snpid[nonunique(qced, ["snpid"])])
    biallelic = filter(:snpid=>∉(duplicate_rsids), qced)

    # Keep only actual SNPs and not other kinds of variants
    actual_snps = subset(biallelic, :allele1 => ByRow(issnp), :allele2 => ByRow(issnp))

    # All batches pass QC 
    batch_cols = [x for x in names(actual_snps) if occursin("Batch", x)]
    final = filter(row -> all_batches_ok(row, batch_cols), actual_snps)
    
    rsids = Set(final.snpid)
    sample_ids = Set(CSV.read(traits_file, DataFrame, select=["SAMPLE_ID"], types=String)[!, 1])
    SnpArrays.filter(
        input; 
        des=output, 
        f_person = x -> x[:iid] ∈ sample_ids, 
        f_snp = x -> x[:snpid] ∈ rsids
    )
end


function adapt_flashpca(input, output)
    # I think FID and IID are just duplicates
    pcs = CSV.File(input, drop=["FID"]) |> DataFrame
    rename!(pcs, :IID => :SAMPLE_ID)
    CSV.write(output, pcs)
end
