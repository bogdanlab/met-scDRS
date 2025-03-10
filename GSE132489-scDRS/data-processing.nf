#!/usr/bin/env nextflow
/*
 * pipeline for methylation data processing:
 */

// define process:
process data_extraction {

        input:
        path executable
        val data_dir
        val output_dir

        output:
        val "${output_dir}full_100_kb_mch_count.csv", emit: mch_100kb_count
        val "${output_dir}full_100_kb_mch_coverage.csv", emit: mch_100kb_coverage
        val "${output_dir}full_100_kb_mcg_count.csv", emit: mcg_100kb_count
        val "${output_dir}full_100_kb_mcg_coverage.csv", emit: mcg_100kb_coverage
        val "${output_dir}full_gene_mch_count.csv", emit: mch_gene_count
        val "${output_dir}full_gene_mch_coverage.csv", emit: mch_gene_coverage
        val "${output_dir}full_gene_mcg_count.csv", emit: mcg_gene_count
        val "${output_dir}full_gene_mcg_coverage.csv", emit: mcg_gene_coverage

        script:
        """
        python $params.executable --input_dir=$data_dir --output_dir=$output_dir
        """
}

process mch_qc_processing {

        input:
        path executable
        val mch_100kb_count
        val mch_100kb_coverage
        val mcg_100kb_count
        val mcg_100kb_coverage
        val mch_gene_count
        val mch_gene_coverage
        val mcg_gene_count
        val mcg_gene_coverage
        path meta_data
        val output_path

        output:
        val "${output_path}", emit: mch_fraction

        script:
        """
        Rscript ${executable} \
                --count ${mch_gene_count} \
                --coverage ${mch_gene_coverage} \
                --meta ${meta_data} \
                --output ${output_path}

        """
}

process mcg_qc_processing {

        input:
        path executable
        val mch_100kb_count
        val mch_100kb_coverage
        val mcg_100kb_count
        val mcg_100kb_coverage
        val mch_gene_count
        val mch_gene_coverage
        val mcg_gene_count
        val mcg_gene_coverage
        path meta_data
        val output_path

        output:
        val "${output_path}", emit: mcg_fraction

        script:
        """
        Rscript ${executable} \
                --count ${mcg_gene_count} \
                --coverage ${mcg_gene_coverage} \
                --meta ${meta_data} \
                --output ${output_path}
        """
}

process concatenate_mch {

        input:
        val pass
        path file_paths
        val concat_file

        output:
        val concat_file

        script:
        """
        cat ${file_paths} > ${concat_file}
        """
}

process concatenate_mcg {

        input:
        val pass
        path file_paths
        val concat_file

        output:
        val concat_file

        script:
        """
        cat ${file_paths} > ${concat_file}
        """
}

process remove_duplication_mch {

        input:
        path concatenated_file
        path header_file, stageAs: 'header.csv'
        val unique_file

        output:
        val unique_file

        script:
        """
        sort -t',' -k1,1 -n ${concatenated_file} | uniq | tail -n +2 > "tmp-unique-file.txt"
        head -n1 header.csv | sed '1s/^"",/cell,/' > "header-line.txt"
        cat header-line.txt tmp-unique-file.txt > ${unique_file}
        rm header-line.txt
        rm tmp-unique-file.txt

        """
}

process remove_duplication_mcg {

        input:
        path concatenated_file
        path header_file, stageAs: 'header.csv'
        val unique_file

        output:
        val unique_file

        script:
        """
        sort -t',' -k1,1 -n ${concatenated_file} | uniq | tail -n +2 > "tmp-unique-file.txt"
        head -n1 header.csv | sed '1s/^"",/cell,/' > "header-line.txt"
        cat header-line.txt tmp-unique-file.txt > ${unique_file}
        rm header-line.txt
        rm tmp-unique-file.txt

        """
}

process spot_check_mch {
        input:
        path executable
        path count_path
        path coverage_path
        path fraction_path
        path meta_path
        path randomize_path

        output:
        path fraction_path

        script:
        """
        # grab out a randomized file:
        tail -n +2 ${fraction_path} > without_header.csv
        shuf without_header.csv > random_fraction.csv
        head -n1 ${fraction_path} > header.txt
        cat header.txt random_fraction.csv > ${randomize_path}

        # using this randomized file, spot check fractions:
        Rscript ${executable} \
                --count ${count_path} \
                --coverage ${coverage_path} \
                --fraction ${randomize_path} \
                --meta ${meta_path} \
                --lines 1000

        # file clean up:
        rm without_header.csv
        rm random_fraction.csv
        rm header.txt
        """

}

process spot_check_mcg {
        input:
        path executable
        path count_path
        path coverage_path
        path fraction_path
        path meta_path
        path randomize_path

        output:
        path fraction_path

        script:
        """
        # grab out a randomized file:
        tail -n +2 ${fraction_path} > without_header.csv
        shuf without_header.csv > random_fraction.csv
        head -n1 ${fraction_path} > header.txt
        cat header.txt random_fraction.csv > ${randomize_path}

        # using this randomized file, spot check fractions:
        Rscript ${executable} \
                --count ${count_path} \
                --coverage ${coverage_path} \
                --fraction ${randomize_path} \
                --meta ${meta_path} \
                --lines 1000

        # file clean up:
        rm without_header.csv
        rm random_fraction.csv
        rm header.txt
        """

}

// define parameters:
// data extraction:
params.executable = "/u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE132489-scDRS/data-extraction.py"
params.data_dir = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/"
params.output_dir = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/extracted/"

// R processing for getting the fraction
params.processing_script = "/u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE132489-scDRS/data-QC.R"
params.meta_data = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/41586_2020_3182_MOESM9_ESM.xlsx"
params.mch_gene_output_path = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/extracted/chunk/mch_gene_qced_fraction"
params.mcg_gene_output_path = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/extracted/chunk/mcg_gene_qced_fraction"

// concatenate all the fraction that we have generated:
params.mch_concat_file = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/mch_gene_fraction.csv"
params.mcg_concat_file = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/mcg_gene_fraction.csv"
params.mch_header = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/extracted/chunk/mch_gene_qced_fraction-chunk0.csv"
params.mcg_header = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/extracted/chunk/mcg_gene_qced_fraction-chunk0.csv"

// spot check to make sure fraction are ok:
params.spot_check_script = "/u/home/l/lixinzhe/project-github/scDRS-applications/code/GSE132489-scDRS/data-validata.R"
params.mch_count = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/extracted/full_gene_mch_count.csv"
params.mch_coverage = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/extracted/full_gene_mch_coverage.csv"
params.mcg_count = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/extracted/full_gene_mcg_count.csv"
params.mcg_coverage = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/extracted/full_gene_mcg_coverage.csv"
params.mch_randomized = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/randomized_mch_gene_fraction.csv"
params.mcg_randomized = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/randomized_mcg_gene_fraction.csv"

// remove redundant line:
params.unique_mch_destination = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/unique_mch_gene_fraction.csv"
params.unique_mcg_destination = "/u/project/geschwind/lixinzhe/data/Liu_et_al_2021_methylation_gse132489/unique_mcg_gene_fraction.csv"

// define workflow
workflow {
        // extract file:    
        extracted_files = data_extraction(params.executable, params.data_dir, params.output_dir)
        mch_file = mch_qc_processing(params.processing_script, extracted_files, params.meta_data, params.mch_gene_output_path)
        mcg_file = mcg_qc_processing(params.processing_script, extracted_files, params.meta_data, params.mcg_gene_output_path)

        // start concatenation:
        mch_file_channels = Channel.fromPath("${params.mch_gene_output_path}*").collect()
        mcg_file_channels = Channel.fromPath("${params.mcg_gene_output_path}*").collect()
        concat_mch = concatenate_mch(mch_file, mch_file_channels, params.mch_concat_file)
        concat_mcg = concatenate_mcg(mcg_file, mcg_file_channels, params.mcg_concat_file)

        // remove that one duplicated line:
        unique_mch = remove_duplication_mch(concat_mch, params.mch_header, params.unique_mch_destination)
        unique_mcg = remove_duplication_mcg(concat_mcg, params.mcg_header, params.unique_mcg_destination)

        // fraction spot check:
        mch_check = spot_check_mch(params.spot_check_script, params.mch_count, params.mch_coverage, unique_mch, params.meta_data, params.mch_randomized)
        mcg_check = spot_check_mcg(params.spot_check_script, params.mcg_count, params.mcg_coverage, unique_mcg, params.meta_data, params.mcg_randomized)
}
