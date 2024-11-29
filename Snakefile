# Snakefile

# Define configuration file that lists samples and metadata
# That is the file you should edit prior to analysis
configfile: "config.yaml"

# Get the software container from config file
singularity: config["singularity"]

# Get list of sample names from config file
SAMPLES = config["samples"]

# Get list of environment factors from config file
ENVFACTOR_NAMES = config["envfactor_names"]

# Get parameters from config file
WZA_WINDOW_SIZE = config["parameters"]["WZA_window_size"]
WZA_FDR = config["parameters"]["WZA_fdr"]

# Get resources from config file
RESOURCES = config["resources"]

rule all:
    input:
        expand("res/{sample}/tmp_{envfactor}.csv", sample=SAMPLES, envfactor=ENVFACTOR_NAMES),
        expand("dat/WZA/{sample}_WZA_input.csv", sample=SAMPLES),
        expand("res/WZA_res/{sample}_{envfactor}_WZA_output.csv", sample=SAMPLES, envfactor=ENVFACTOR_NAMES),
        expand("res/WZA_res/{sample}_WZA_manhattan_plots.png", sample=SAMPLES)
    resources:
        runtime=RESOURCES["all"]["runtime"],
        mem_mb=RESOURCES["all"]["mem_mb"],
        slurm_partition=RESOURCES["all"]["slurm_partition"]

# This rule prepares intermediate files required for the WZA (Windowed Z-score Analysis).
# It processes LFMM result files for each environmental factor, extracts the relevant columns,
# removes headers, adds new headers with the environmental factor names,
# and merges the temporary CSV files into a single file for each sample.
rule gather_pvalues:
    resources:
        runtime=RESOURCES["gather_pvalues"]["runtime"],
        mem_mb=RESOURCES["gather_pvalues"]["mem_mb"],
        slurm_partition=RESOURCES["gather_pvalues"]["slurm_partition"]
    output:
        WZA_pre_input = temp("res/{sample}/tmp_{envfactor}.csv"),
    shell:
        """
        for file in res/{wildcards.sample}/Full_analysis_K_1/environment_{wildcards.envfactor}/LFMM_AllResults_env_{wildcards.envfactor}_K1.csv; do
            cut -d',' -f3 $file \
                    | sed '1d' \
                    | sed '1i {wildcards.envfactor}' > {output.WZA_pre_input}
        done
        """

rule merge_tmp_files:
    input:
        expand("res/{sample}/tmp_{envfactor}.csv", sample=SAMPLES, envfactor=ENVFACTOR_NAMES),
    output:
        merged_tmp = temp("res/{sample}/merged_tmp.csv"),
    resources:
        runtime=RESOURCES["merge_tmp_files"]["runtime"],
        mem_mb=RESOURCES["merge_tmp_files"]["mem_mb"],
        slurm_partition=RESOURCES["merge_tmp_files"]["slurm_partition"]
    shell:
        """
        paste -d',' res/{wildcards.sample}/tmp_*.csv > res/{wildcards.sample}/merged_tmp.csv
        """

# This rule prepares the input file for the WZA (Windowed Z-score Analysis).
# It extracts the chromosome and position columns from the Allele Frequency Table,
# calculates the mean MAF (Minor Allele Frequency) for each SNP,
# and generates the input file for the WZA analysis by combining the chromosome, position, MAF, and merged LFMM results.
rule prepare_WZA_input:
    input:
        frequency_table = "dat/{sample}_AlleleFrequencyTable.txt",
        merged_tmp = "res/{sample}/merged_tmp.csv",
    params:
        window_size = WZA_WINDOW_SIZE,
    output:
        WZA_input = "dat/WZA/{sample}_WZA_input.csv",
    resources:
        runtime=RESOURCES["prepare_WZA_input"]["runtime"],
        mem_mb=RESOURCES["prepare_WZA_input"]["mem_mb"],
        slurm_partition=RESOURCES["prepare_WZA_input"]["slurm_partition"]
    shell:
        """
        cut -d'\t' -f1 {input.frequency_table} \
            | awk -F'_' 'BEGIN{{OFS=","}} {{print $1, $2, $3}}' \
            | cut -d',' -f2,3 \
            | sed '1d' | sed '1i CHR,POS' > "res/{wildcards.sample}/CHR_POS.csv"
        
        awk -f scripts/Calculate_Mean_MAF.sh {input.frequency_table} | cut -d',' -f3 > "res/{wildcards.sample}/tmpMAF.csv"

        paste -d',' res/{wildcards.sample}/CHR_POS.csv res/{wildcards.sample}/merged_tmp.csv res/{wildcards.sample}/tmpMAF.csv \
            | awk -f scripts/Make_Genomic_Windows.sh -v window_size={params.window_size} > {output.WZA_input}

        rm res/{wildcards.sample}/tmpMAF.csv
        rm res/{wildcards.sample}/CHR_POS.csv
        """

# Perform the Windowed Z-score Analysis (WZA) using the input file generated in the previous step.
rule WZA:
    input:
        WZA_input = "dat/WZA/{sample}_WZA_input.csv",
    resources:
        runtime=RESOURCES["WZA"]["runtime"],
        mem_mb=RESOURCES["WZA"]["mem_mb"],
        slurm_partition=RESOURCES["WZA"]["slurm_partition"]
    output:
        WZA_output = protected("res/WZA_res/{sample}_{envfactor}_WZA_output.csv"),
    shell:
        """
        python3 scripts/general_WZA_script.py \
            --correlations {input.WZA_input} \
            --summary_stat {wildcards.envfactor} \
            --window window_id \
            --MAF MAF \
            --sep "," \
            --retain POS \
            --output {output.WZA_output}
        """

# Generates diagnostic plots (P-value distributions, QQplots) for the WZA results. 
# Also generates manhattan plots for each environmental factor.
rule WZA_diagnostics:
    input:
        envfactor_names = "dat/{sample}_efile_envfactor_names",
        WZA_output = expand("res/WZA_res/{{sample}}_{envfactor}_WZA_output.csv", envfactor=ENVFACTOR_NAMES),
    params:
        prefix = "res/WZA_res/{sample}_",
        FDR_level = WZA_FDR,
    resources:
        runtime=RESOURCES["WZA_diagnostics"]["runtime"],
        mem_mb=RESOURCES["WZA_diagnostics"]["mem_mb"],
        slurm_partition=RESOURCES["WZA_diagnostics"]["slurm_partition"]
    output:
        WZA_manhattan_plots = "res/WZA_res/{sample}_WZA_manhattan_plots.png",
    script: "scripts/WZA_diagnostics.R"