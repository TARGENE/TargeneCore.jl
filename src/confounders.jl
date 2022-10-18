# ADAPTING, MERGING AND FILTERING PRIOR TO PCA
#####################################################################


function merge_beds(parsed_args)
    merge_plink(parsed_args["input"]; des = parsed_args["output"])
end

issnp(x::String) = length(x) == 1

all_batches_ok(row, batchcols) = all(x == 1 for x in row[batchcols])

bounds(pos, lb, ub) =  (min(lb, pos - 10^4), max(ub, pos + 10^4))

function in_ld_blocks(chr, position, ld_blocks)
    group = ld_blocks[(chr=chr,)]
    return any(b.lower_bound < position < b.upper_bound for b in eachrow(group))
end

function notin_ldblocks(row, ldblocks)
    if (chr=row.chromosome,) ∉ keys(ldblocks)
        return true
    else
        return ~in_ld_blocks(row.chromosome, row.position, ldblocks)
    end
end

function ld_blocks_by_chr(path)
    # Load and redefine LD bounds
    ld_blocks = CSV.File(
        path;
        header=["rsid","chr","pos","LDblock_lower","LDblock_upper","LDblock_length","lower_bound","upper_bound"]) |> DataFrame
    transform!(ld_blocks,
        :,
        [:pos, :lower_bound, :upper_bound] => ByRow(bounds) => [:lower_bound, :upper_bound],
        )
    return groupby(ld_blocks, :chr)
end

"""
    filter_chromosome(parsed_args)

The purpose of this method is to filter SNPS and Individuals before applying a 
dimensionality reduction technique such as PCA that is later used
to extract population stratification compomemts.
We filter SNPs using quality control metrics from the following resource:
    - https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=1955
"""
function filter_bed_file(parsed_args)

    qc_df = CSV.File(parsed_args["qcfile"]) |> DataFrame
    
    snp_data = SnpData(parsed_args["input"])
    snp_data.snp_info.chromosome = parse.(Int, snp_data.snp_info.chromosome)
    # Load and redefine LD bounds
    ld_blocks = ld_blocks_by_chr(parsed_args["ld-blocks"])

    # Remove SNP's with MAF < maf-threshold
    maf_threshold = parsed_args["maf-threshold"]
    snp_data.snp_info[!, "MAF"] = SnpArrays.maf(snp_data.snparray)
    mafpassed = filter(:MAF => >=(maf_threshold), snp_data.snp_info)
    
    # Remove LD regions specified by ld_blocks
    ld_pruned = filter(x -> notin_ldblocks(x, ld_blocks), mafpassed)

    # The QC file contains information on fully genotyped SNPS
    # We only keep those
    fully_genotyped_snps = innerjoin(
        ld_pruned, 
        qc_df, 
        on = :snpid => :rs_id,
        makeunique = true
    )

    # If an RSID appears multiple times, it is because it has 
    # more than 2 possible alleles: we remove them 
    # (why? maybe because the PCA then cannot tackle them)
    duplicate_rsids = Set(fully_genotyped_snps.snpid[nonunique(fully_genotyped_snps, ["snpid"])])
    biallelic = filter(:snpid=>∉(duplicate_rsids), fully_genotyped_snps)

    # Keep only actual SNPs and not other kinds of variants
    actual_snps = subset(biallelic, :allele1 => ByRow(issnp), :allele2 => ByRow(issnp))

    # All batches pass QC 
    batch_cols = [x for x in names(actual_snps) if occursin("Batch", x)]
    batches_ok = filter(row -> all_batches_ok(row, batch_cols), actual_snps)

    # Assayed in both genotyping arrays
    final = filter(:array => ==(2), batches_ok)
    
    rsids = Set(final.snpid)
    sample_ids = Set(CSV.read(parsed_args["traits"], DataFrame, select=["SAMPLE_ID"], types=String)[!, 1])
    SnpArrays.filter(
        parsed_args["input"]; 
        des=parsed_args["output"], 
        f_person = x -> x[:iid] ∈ sample_ids, 
        f_snp = x -> x[:snpid] ∈ rsids
    )
end


function adapt_flashpca(parsed_args)
    # I think FID and IID are just duplicates
    pcs = CSV.File(parsed_args["input"], drop=["IID"]) |> DataFrame
    rename!(pcs, :FID => :SAMPLE_ID)
    CSV.write(parsed_args["output"], pcs)
end


function filter_bgen_file(args)
    path = string(ukb_dir, "imputed/ukb_53116_chr", chr, ".bgen")
    sample_path = string(ukb_dir, "imputed/ukb_53116_chr", chr, ".sample")
    idx_path = string(ukb_dir, "imputed/ukb_53116_chr", chr, ".bgen.bgi")
    maf_threshold = 0.01

    bgen = Bgen(path; sample_path=sample_path, idx_path=idx_path)
    ld_blocks_by_chr(args["ld-blocks"])
    qc_df = CSV.File("/home/s2042526/UK-BioBank-53116/imputed/ukb_snp_qc.txt") |> DataFrame
    batch_cols = [x for x in names(qc_df) if occursin("Batch", x)]

    variant_mask = zeros(Bool, n_variants(bgen))

    genotyped_snps = Set(qc_df.rs_id)
    for (index, v) in enumerate(iterator(bgen))
        # Only genotyped SNPs
        v ∈ genotyped_snps || continue

        # SNP is not in LD blocks
        in_ld_blocks(parse(Int, v.chrom), pos(v), ld_blocks_by_chr) && continue

        ## Check SNP passes QC
        snp_qc = filter(x -> x.rs_id == v.rsid, qc_df)[1, :]
        # SNP assayed in both arrays
        snp_qc.array == 2 || continue
        # All batches pass QC 
        all(==(1), values(snp_qc[batch_cols])) || continue
        # Minor Allele frequency
        maf = mean(minor_allele_dosage!(bgen, v))
        maf >= maf_threshold || continue

        # All checks have passed
        variant_mask[index] = true
    end

    sample_ids = Set(CSV.read(parsed_args["traits"], DataFrame, select=["SAMPLE_ID"], types=String)[!, 1])
    sample_mask = [s ∈ sample_ids for s in samples(bgen)]

    BGEN.filter(args["out"], bgen, variant_mask, sample_mask)
end

function filter_ukb_genetic_file(args)
    if endswith(args["input"], "bgen")
        filter_bgen_file(args)
    else
        filter_bed_file(args)
    end
end