# ADAPTING, MERGING AND FILTERING PRIOR TO PCA
#####################################################################


function merge_beds(parsed_args)
    merge_plink(parsed_args["input"]; des = parsed_args["output"])
end

issnp(x) = length(x) == 1

all_batches_ok(row, batchcols) = all(x == 1 for x in row[batchcols])

"""
At least 2*10^4 long LD blocks.
"""
bounds(pos, lb, ub) =  (min(lb, pos - 10^4), max(ub, pos + 10^4))

function in_ld_blocks(chr, position, ld_blocks)
    (chr=chr, ) ∉ keys(ld_blocks) && return false
    group = ld_blocks[(chr=chr,)]
    return any(b.lower_bound < position < b.upper_bound for b in eachrow(group))
end

in_ld_blocks(row::DataFrameRow, ld_blocks) =
    in_ld_blocks(row.chromosome, row.position, ld_blocks)


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

    input_prefix = parsed_args["input"][1:end-4]
    output_prefix = parsed_args["output"][1:end-4]
    qc_df = CSV.File(parsed_args["qcfile"]) |> DataFrame
    
    snp_data = SnpData(input_prefix)
    snp_data.snp_info.chromosome = parse.(Int, snp_data.snp_info.chromosome)
    # Load and redefine LD bounds
    ld_blocks = ld_blocks_by_chr(parsed_args["ld-blocks"])

    # Remove SNP's with MAF < maf-threshold
    maf_threshold = parsed_args["maf-threshold"]
    snp_data.snp_info[!, "MAF"] = SnpArrays.maf(snp_data.snparray)
    mafpassed = filter(:MAF => >=(maf_threshold), snp_data.snp_info)
    
    # Remove LD regions specified by ld_blocks
    ld_pruned = filter(x -> ~in_ld_blocks(x, ld_blocks), mafpassed)

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
    actual_snps = subset(biallelic, :allele1_ref => ByRow(issnp), :allele2_alt => ByRow(issnp))

    # All batches pass QC 
    batch_cols = [x for x in names(actual_snps) if occursin("Batch", x)]
    batches_ok = filter(row -> all_batches_ok(row, batch_cols), actual_snps)

    # Assayed in both genotyping arrays
    final = filter(:array => ==(2), batches_ok)
    
    rsids = Set(final.snpid)
    sample_ids = Set(CSV.read(parsed_args["traits"], DataFrame, select=["SAMPLE_ID"], types=String)[!, 1])
    SnpArrays.filter(
        input_prefix; 
        des=output_prefix, 
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


function filter_bgen_file(parsed_args)
    maf_threshold = parsed_args["maf-threshold"]
    bgenpath = parsed_args["input"]
    bgi_path = bgenpath * ".bgi"
    sample_path = bgenpath[1:end-4] * "sample"

    # Load data
    bgen = Bgen(bgenpath; sample_path=sample_path, idx_path=bgi_path)
    ld_blocks = TargeneCore.ld_blocks_by_chr(parsed_args["ld-blocks"])
    qc_df = CSV.File(parsed_args["qcfile"]) |> DataFrame

    # Filter QC file so that only unique SNPs will be kept (uniqueness seems worthless but doesn't cost much)
    qc_df = unique(
        filter([:allele1_ref, :allele2_alt] => (ref, alt) -> issnp(ref) && issnp(alt), qc_df),
        :rs_id
        )
    # Retrieve Batch QC columns
    batch_cols = [x for x in names(qc_df) if occursin("Batch", x)]

    # Initialize mask
    variant_mask = zeros(Bool, n_variants(bgen))

    genotyped_snps = Set(qc_df.rs_id)
    for (index, v) in enumerate(iterator(bgen))
        # Only genotyped SNPs
        v.rsid ∈ genotyped_snps || continue

        # SNP is not in LD blocks
        in_ld_blocks(parse(Int, v.chrom), pos(v), ld_blocks) && continue

        ## Check SNP passes QC
        snp_qc = filter(x -> x.rs_id == v.rsid, qc_df)[1, :]
        # SNP assayed in both arrays
        snp_qc.array == 2 || continue
        # All batches pass QC 
        all(==(1), values(snp_qc[batch_cols])) || continue
        # Minor Allele frequency
        maf = mean(minor_allele_dosage!(bgen, v)) / 2
        maf >= maf_threshold || continue

        # All checks have passed
        variant_mask[index] = true
    end

    sample_ids = Set(CSV.read(parsed_args["traits"], DataFrame, select=["SAMPLE_ID"], types=String)[!, 1])
    sample_mask = [s ∈ sample_ids for s in samples(bgen)]

    outpath = parsed_args["output"]
    # Filter the BGEN file
    BGEN.filter(outpath, bgen, BitVector(variant_mask), BitVector(sample_mask))

    # Index the files
    run(`bgenix -index -g $outpath`)
end

"""
Filters a UK-Biobank genetic (.BGEN or .BED) file given by `input`.

The `output` is based on the following rules, are kept only:
    1. genotyped SNPs (in `qcfile`)
    2. SNPs not contained in the `ld-blocks` provided file
    3. SNPs that pass QC (info in `qcfile`):
        1. have been assayed in both arrays
        2. SNPs for which all batches pass QC
    5. SNPs with minor allele frequency >= `maf-threshold`
    6. bi-allelic SNPs
    6. Individuals present in the `traits` dataset

"""
function filter_ukb_genetic_file(args)
    if endswith(args["input"], "bgen")
        filter_bgen_file(args)
    else
        filter_bed_file(args)
    end
end