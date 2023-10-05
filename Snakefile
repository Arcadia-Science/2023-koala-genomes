import pandas as pd
        
samples_df = pd.read_table('inputs/samples.tsv').set_index("AWSFileName",drop=False)
SAMPLE = list(samples_df['AWSFileName'])
GENES= ['RUNX2', 'ZIC2', 'FOXL2', 'ARX']

gene_data = {
    'RUNX2': {"coordinates": "MSTS01000055.1:12711350-12714427", "pept": "VAAQ\+E\+A\+"},
    'ZIC2': {"coordinates": "MSTS01000150.1:3308891-3309559", "pept": "LSPA\+"},
    'FOXL2': {"coordinates": "MSTS01000024.1:18722964-18774112", "pept": "CQMA\+"},
    'ARX': {"coordinates": "MSTS01000037.1:11243615-11249999", "pept": "LQGA\+"}
}

rule all:
    input:
        expand("outputs/counts/{sample}_{gene}_count.txt", sample=SAMPLE, gene = GENES),
        expand("outputs/counts/{sample}_RUNX2_count_q.txt", sample=SAMPLE),
        "outputs/counts/finaltable.txt"

rule download_extract:
    """
    This rule downloads bams for 430 koala genomes (only keep one at a time),  extracts the reads from the original bam that mapped to the genes of interest, and removes the bam file after.
    """
    output:
        "outputs/fastq/{sample}_mapped_{gene}.fastq"
    conda:"envs/samtools.yml"
    params:
        folder = lambda wildcards: samples_df.loc[wildcards.sample]["AWSFolderName"],
        coords = lambda wildcards: gene_data[wildcards.gene]["coordinates"],
        bam = "outputs/downloadbams/{wildcards.sample}.bam"
    shell: """
        export AWS_REGION=ap-southeast-2
        aws s3 cp "s3://koalagenomes/{params.folder}/bam/{wildcards.sample}.bam" {params.bam}
        samtools view -b {params.bam} "{params.coords}" > outputs/bams/{wildcards.sample}_mapped_{wildcards.gene}.bam
        bedtools bamtofastq -i outputs/bams/{wildcards.sample}_mapped_{wildcards.gene}.bam -fq {output}
        rm {params.bam}
        """

rule assemble:
    """
    This rule assembles the kmer-baited reads and mapped reads for the region of interest
    """
    input:
        mapped="outputs/fastq/{sample}_mapped_{gene}.fastq"
    output:
        "outputs/assembled/{sample}_{gene}/{sample}.fasta"
    conda:"envs/spades.yml"
    shell: """
        spades.py --isolate -s "{input.mapped}" -o "outputs/assembled/{wildcards.sample}_{wildcards.gene}"
        mv outputs/assembled/{wildcards.sample}_{wildcards.gene}/contigs.fasta {output}
        """

rule orf_call:
    """
    This rule translates ORFs and then pulls out the peptide of interest
    """
    input:
        "outputs/assembled/{sample}_{gene}/{sample}.fasta"
    output:
        "outputs/peptides/{sample}/{sample}_peptide_filt_{gene}.fa"
    conda:"envs/orfipy.yml"
    params:
        pept = lambda wildcards: gene_data[wildcards.gene]["pept"]
    shell: """
        orfipy {input} --pep orf_peptides_{wildcards.gene}.fa --outdir outputs/peptides/{wildcards.sample}/ --partial-3
        awk '/^>/ {{if (seq) print seq; print; seq=""; next;}} {{seq=seq $0;}} END {{if (seq) print seq;}}' outputs/peptides/{wildcards.sample}/orf_peptides_{wildcards.gene}.fa | grep -B 1 '{params.pept}' > outputs/peptides/{wildcards.sample}/{wildcards.sample}_peptide_filt_{wildcards.gene}.fa
        """

rule count_expansion:
    """
    This rule finds the expansion and determines its length
    """
    input:
        "outputs/peptides/{sample}/{sample}_peptide_filt_{gene}.fa"
    output:
        "outputs/counts/{sample}_{gene}_count.txt"
    shell:
        """
        longest_stretch=$(grep -v ">" {input} | tr -d '\n' | grep -o 'A\\+' | awk '{{ print length($0) }}' | sort -n | tail -n 1)
        echo -e "{wildcards.sample}\t${{longest_stretch}}\t{wildcards.gene}_a" > "{output}"
        """

rule count_q:
    """
    This rule finds the runx2 Q expansion and determines its length
    """
    input:
        "outputs/peptides/{sample}/{sample}_peptide_filt_RUNX2.fa"
    output:   
        "outputs/counts/{sample}_RUNX2_count_q.txt"
    shell:
        """
        longest_stretch=$(grep -v ">" {input} | tr -d '\n' | grep -o 'Q\\+' | awk '{{ print length($0) }}' | sort -n | tail -n 1)
        echo -e "{wildcards.sample}\t${{longest_stretch}}\tRUNX2_q" > "outputs/counts/{wildcards.sample}_RUNX2_count_q.txt"
        """

rule combine_runx2:
    """
    # Rule for combining RUNX2 files
    """
    input: 
        files=expand("outputs/counts/{sample}_RUNX2_count_q.txt", sample=SAMPLE)
    output: 
        temp("outputs/counts/runx2_combined.txt")
    run:
        dfs = [pd.read_csv(f, sep='\t', header=None) for f in input.files]
        merged = pd.concat(dfs)
        merged.to_csv(output[0], sep='\t', header=False, index=False)

rule combine_gene:
    """
    # Rule for combining gene files
    """
    input: 
        files=expand("outputs/counts/{sample}_{gene}_count.txt",sample=SAMPLE, gene = GENES)
    output: 
        temp("outputs/counts/gene_combined.txt")
    run:
        dfs = [pd.read_csv(f, sep='\t', header=None) for f in input.files]
        merged = pd.concat(dfs)
        merged.to_csv(output[0], sep='\t', header=False, index=False)

rule combine:
    """
    Combine the consolidated files
    """
    input:
        runx2="outputs/counts/runx2_combined.txt",
        gene="outputs/counts/gene_combined.txt"
    output:
        "outputs/counts/finaltable.txt"
    run:
        dfs = [pd.read_csv(f, sep='\t', header=None) for f in [input.runx2, input.gene]]
        merged = pd.concat(dfs)
        merged.to_csv(output[0], sep='\t', header=False, index=False)
