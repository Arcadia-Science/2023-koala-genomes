import pandas as pd
        
samples_df = pd.read_table('inputs/samples.tsv').set_index("AWSFileName",drop=False)
SAMPLE = list(samples_df['AWSFileName'])
GENES= ['RUNX2', 'ZIC2', 'FOXL2', 'ARX']

gene_data = {
    'RUNX2': {"coordinates": "MSTS01000055.1:12711350-12714427", "pept": "VAAQ\+E\+A\+"},
    'ZIC2': {"coordinates": "MSTS01000150.1:3308891-3309559", "pept": "LSPA\+"},
    'FOXL2': {"coordinates": "MSTS01000024.1:18722964-18774112", "pept": "CQMA\+"},
    'ARX': {"coordinates": "MSTS01000037.1:11244865-11245057", "pept": "LQGA\+"}
}

rule all:
    input:
        expand("outputs/counts/{sample}_{gene}_count.txt", sample=SAMPLE, gene = GENES),
        expand("outputs/counts/{sample}_RUNX2_count_q.txt", sample=SAMPLE),
        "outputs/counts/finaltable.txt"

rule download_bams:
      """
      This rule downloads bams for 430 koala genomes (only keep one at a time)
      """
      output:
        temp("outputs/downloadbams/{sample}.bam"),
        temp("outputs/downloadbams/{sample}.bam.bai")
      shell: """
        export AWS_REGION=ap-southeast-2
        s5cmd cp "s3://koalagenomes/*/bam/{wildcards.sample}.bam" "{output[0]}"
        s5cmd cp "s3://koalagenomes/*/bam/{wildcards.sample}.bam.bai" "{output[1]}"
        """

rule get_unmapped:
    """
    This rule gets unmapped reads from bam file
    """
    input: "outputs/downloadbams/{sample}.bam"
    output: 
        out1="outputs/bams/unmap1_{sample}.bam",
        out2="outputs/bams/unmap2_{sample}.bam",
        out3="outputs/bams/unmap3_{sample}.bam"
    conda:"envs/samtools.yml"
    threads: 10
    shell: """
        samtools view -u -f 4 -F 264 {input} > {output.out1}
        samtools view -u -f 8 -F 260 {input} > {output.out2}
        samtools view -u -f 12 -F 256 {input} > {output.out3}
        """

rule merge_unmapped:
    """
    This rule gets unmapped reads from bam file                              
    """
    input:
        out1="outputs/bams/unmap1_{sample}.bam",
        out2="outputs/bams/unmap2_{sample}.bam",
        out3="outputs/bams/unmap3_{sample}.bam"
    conda:"envs/samtools.yml"
    output:
        "outputs/bams/unmapped_{sample}.bam"
    threads: 10
    shell: """
        samtools merge -u {output} {input}
        """

rule unmapped_to_reads:
    """
    This rule sorts the unmapped reads bam and creates fastq files
    """
    input:
        "outputs/bams/unmapped_{sample}.bam"
    output:
        "outputs/fastq/unmapped_r1_{sample}.fastq",
        "outputs/fastq/unmapped_r2_{sample}.fastq"
    conda:"envs/samtools.yml"
    threads: 10
    shell: """
        samtools sort -n {input} -o outputs/bams/unmapped_sort_{wildcards.sample}.bam
        bedtools bamtofastq -i outputs/bams/unmapped_sort_{wildcards.sample}.bam -fq {output[0]} -fq2 {output[1]}
        """

rule create_kmer_file:
    """
    This rule generates kmers from genes of interest
    """
    input: "inputs/{gene}_koala.fasta"
    output:
        "outputs/{gene}_kmers.fasta"
    conda:"envs/bbmap.yml"
    shell: """
        kmercountexact.sh in={input} out={output} k=21
        """

rule kmer_bait:
    """
    This rule filters reads based on kmers of interest
    """
    input:
        r1="outputs/fastq/unmapped_r1_{sample}.fastq",
        r2="outputs/fastq/unmapped_r2_{sample}.fastq",
        kmers="outputs/{gene}_kmers.fasta"
    output:
        r1="outputs/fastq/{sample}_unmapped_baited1_{gene}.fastq",
        r2="outputs/fastq/{sample}_unmapped_baited2_{gene}.fastq",
    conda:"envs/bbmap.yml"
    shell: """
        bbmap.sh ref={input.kmers} in={input.r1} in2={input.r2} outm1={output.r1} outm2={output.r2} nodisk 
        """

rule mapped_reads:
    """
    This rule extracts the reads from the original bam that mapped to the gene of interest.
    """
    input:
        bam="outputs/downloadbams/{sample}.bam"
    output:
        "outputs/fastq/{sample}_mapped_{gene}.fastq",
    conda:"envs/samtools.yml"
    params:
        coords = lambda wildcards: gene_data[wildcards.gene]["coordinates"]
    shell: """
        samtools view -b {input.bam} "{params.coords}" > outputs/bams/{wildcards.sample}_mapped_{wildcards.gene}.bam
        bedtools bamtofastq -i outputs/bams/{wildcards.sample}_mapped_{wildcards.gene}.bam -fq {output}
        """

rule assemble:
    """
    This rule assembles the kmer-baited reads and mapped reads for the region of interest
    """
    input:
        mapped="outputs/fastq/{sample}_mapped_{gene}.fastq",
        baited1="outputs/fastq/{sample}_unmapped_baited1_{gene}.fastq",
        baited2="outputs/fastq/{sample}_unmapped_baited2_{gene}.fastq",
    output:
        "outputs/assembled/{sample}_{gene}/{sample}.fasta",
    conda:"envs/spades.yml"
    shell: """
        run_spades() {{
            local fastq1=$1
            local fastq2=$2
            local fastq_s=$3
            local output_dir=$4

            if [[ -s "$fastq1" ]]; then
                input1="-1 $fastq1"
            else
                input1=""
            fi

            if [[ -s "$fastq2" ]]; then
                input2="-2 $fastq2"
            else
                input2=""
            fi

            if [[ "$input1" || "$input2" ]]; then
                spades.py --isolate $input1 $input2 -s "$fastq_s" -o "$output_dir" -m 1024
            else
                spades.py --isolate -s "$fastq_s" -o "$output_dir" -m 1024
            fi
        }}

        run_spades "{input.baited1}" "{input.baited2}" "{input.mapped}" "outputs/assembled/{wildcards.sample}_{wildcards.gene}"
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

checkpoint combine_runx2:
    """
    # Checkpoint for combining RUNX2 files
    """
    input: 
        files=expand("outputs/counts/{sample}_RUNX2_count_q.txt", sample=SAMPLE)
    output: 
        temp("outputs/counts/runx2_combined.txt")
    run:
        dfs = [pd.read_csv(f, sep='\t', header=None) for f in input.files]
        merged = pd.concat(dfs)
        merged.to_csv(output[0], sep='\t', header=False, index=False)

checkpoint combine_gene:
    """
    # Checkpoint for combining gene files
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
