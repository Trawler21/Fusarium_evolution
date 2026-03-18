nextflow.enable.dsl=2

params.orthogroups_tsv = null
params.single_copy_txt = null
params.proteomes_dir = null
params.species_ref = null

params.expect_species_n = 193
params.min_coverage = 0.5
params.total_species = 193

params.outdir = 'results'
params.do_iqtree = false
params.iqtree_extra = ''
params.concat_script = "${projectDir}/concat_trimmed_alignments.py"

if (!params.orthogroups_tsv || !params.single_copy_txt || !params.proteomes_dir || !params.species_ref) {
    error "Missing required params: --orthogroups_tsv --single_copy_txt --proteomes_dir --species_ref"
}

process MAKE_SC_TABLE {
    tag 'make_sc_table'
    publishDir "${params.outdir}/tables", mode: 'copy'

    input:
    path orthogroups_tsv
    path single_copy_txt

    output:
    path 'sc_og_species_gene_wide.tsv', emit: wide
    path 'sc_og_species_gene_long.tsv', emit: long

    script:
    """
    python3 ${projectDir}/make_sc_og_species_gene_table.py \
      --og_tsv ${orthogroups_tsv} \
      --sc_file ${single_copy_txt} \
      --out_wide sc_og_species_gene_wide.tsv \
      --out_long sc_og_species_gene_long.tsv \
      --expect_species_n ${params.expect_species_n}
    """
}

process EXTRACT_OG_FASTA {
    tag 'extract_sc_ogs'
    publishDir "${params.outdir}/og_fasta", mode: 'copy'

    input:
    path long_table
    path proteomes_dir

    output:
    path 'singlecopy_og_fasta/*.faa', emit: og_fastas

    script:
    """
    python3 ${projectDir}/extract_sc_ogs_by_table.py \
      --table_long ${long_table} \
      --proteomes_dir ${proteomes_dir} \
      --out_dir singlecopy_og_fasta \
      --min_coverage ${params.min_coverage} \
      --total_species ${params.total_species}
    """
}

process CHECK_OG_COVERAGE {
    tag 'check_og_coverage'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path species_ref
    path og_fastas

    output:
    path 'og_coverage_dup_report.tsv'

    script:
    """
    python3 ${projectDir}/check_og_coverage_and_dup.py \
      --ref ${species_ref} \
      --dir . \
      --pattern "*.faa" \
      --report og_coverage_dup_report.tsv
    """
}

process MAFFT_ALIGN {
    tag { og_faa.baseName }
    publishDir "${params.outdir}/aln", mode: 'copy'

    input:
    path og_faa

    output:
    path "${og_faa.baseName}.aln"

    script:
    """
    mafft --auto ${og_faa} > ${og_faa.baseName}.aln
    """
}

process TRIMAL_ALIGN {
    tag { aln.baseName }
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    path aln

    output:
    path "${aln.baseName}.trimmed.aln"

    script:
    """
    trimal -in ${aln} -out ${aln.baseName}.trimmed.aln -automated1
    """
}

process CHECK_TRIMMED_LENGTHS {
    tag 'check_trimmed'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path trimmed_files

    output:
    path 'trimmed_length_report.tsv'

    script:
    """
    python3 ${projectDir}/check_trimmed_aln_lengths.py \
      --dir . \
      --pattern "*.trimmed.aln" \
      --report trimmed_length_report.tsv
    """
}

process CONCAT_ALIGNMENTS {
    tag 'concat'
    publishDir "${params.outdir}/concat", mode: 'copy'

    input:
    path trimmed_files
    path species_ref

    output:
    path 'concatenated_alignment.fasta', emit: concat_fa
    path 'partition.tsv', emit: partition

    script:
    """
    python3 ${params.concat_script} \
      --aln_dir . \
      --pattern "*.trimmed.aln" \
      --ref ${species_ref} \
      --out_fasta concatenated_alignment.fasta \
      --out_partitions partition.tsv
    """
}

process IQTREE {
    tag 'iqtree'
    publishDir "${params.outdir}/iqtree", mode: 'copy'

    input:
    path concat_fasta

    output:
    path 'iqtree_out/*'

    script:
    """
    mkdir -p iqtree_out
    iqtree2 -s ${concat_fasta} -m MFP -bb 1000 -alrt 1000 ${params.iqtree_extra} --prefix iqtree_out/concat
    """
}

workflow {
    ch_og_tsv = Channel.fromPath(params.orthogroups_tsv, checkIfExists: true)
    ch_sc_txt = Channel.fromPath(params.single_copy_txt, checkIfExists: true)
    ch_proteomes = Channel.fromPath(params.proteomes_dir, checkIfExists: true, type: 'dir')
    ch_species_ref = Channel.fromPath(params.species_ref, checkIfExists: true)

    sc_tables = MAKE_SC_TABLE(ch_og_tsv, ch_sc_txt)
    extracted = EXTRACT_OG_FASTA(sc_tables.long, ch_proteomes)

    // 用于 coverage 检查：保留为整体文件列表
    extracted.og_fastas
        .collect()
        .set { all_og_fastas }

    CHECK_OG_COVERAGE(ch_species_ref, all_og_fastas)

    // 用于 mafft：打散成单个 faa 文件逐个比对
    aligned = MAFFT_ALIGN(extracted.og_fastas.flatten())
    trimmed = TRIMAL_ALIGN(aligned)

    trimmed_list = trimmed.collect()

    CHECK_TRIMMED_LENGTHS(trimmed_list)
    concat_out = CONCAT_ALIGNMENTS(trimmed_list, ch_species_ref)

    if (params.do_iqtree) {
        IQTREE(concat_out.concat_fa)
    }
}
