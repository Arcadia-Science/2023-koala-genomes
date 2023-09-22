import pandas as pd
        
samples_df = pd.read_table('inputs/samples.tsv').set_index("AWSFileName",drop=False)
SAMPLE = list(samples_df['AWSFileName'])
GENES= ['RUNX2', 'ZIC2', 'FOXL2', 'ARX']
rule all:
    input:
        "outputs/counts/stretch_lengths.txt",
        expand("outputs/counts/{sample}_{gene}_count.txt", sample=SAMPLE, gene = GENES)
     
rule get_unmapped:
    """
    This rule gets unmapped reads from bam file
    """
    input: "outputs/bams/{sample}.bam"
    output: 
        out1="outputs/bams/{sample}_unmap1.bam",
        out2="outputs/bams/{sample}_unmap2.bam",
        out3="outputs/bams/{sample}_unmap3.bam"
    conda:"envs/samtools.yml"
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
        "outputs/bams/{sample}_unmap1.bam", 
        "outputs/bams/{sample}_unmap2.bam",
        "outputs/bams/{sample}_unmap3.bam"
    conda:"envs/samtools.yml"
    output:
        "outputs/bams/{sample}_unmapped.bam"
    shell: """
        samtools merge -u {output} {input}
        """

rule unmapped_to_reads:
    """
    This rule sorts the unmapped reads bam and creates fastq files
    """
    input:
        "outputs/bams/{sample}_unmapped.bam"
    output:
        "outputs/fastq/{sample}_unmapped_r1.fastq",
        "outputs/fastq/{sample}_unmapped_r2.fastq"
    conda:"envs/samtools.yml"
    shell: """
        samtools sort -n {input} -o outputs/bams/{wildcards.sample}_unmapped_sort.bam
        bedtools bamtofastq -i outputs/bams/{wildcards.sample}_unmapped_sort.bam -fq {output[0]} -fq2 {output[1]}
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
        r1="outputs/fastq/{sample}_unmapped_r1.fastq",
        r2="outputs/fastq/{sample}_unmapped_r2.fastq",
        kmers="outputs/{gene}_kmers.fasta",
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
        bam="outputs/bams/{sample}.bam"
    output:
        "outputs/fastq/{sample}_mapped_{gene}.fastq",
    conda:"envs/samtools.yml"
    shell: """
        samtools view -b {input.bam} "MSTS01000055.1:12711350-12714427" > outputs/bams/{wildcards.sample}_mapped_runx2.bam
        bedtools bamtofastq -i outputs/bams/{wildcards.sample}_mapped_{gene}.bam -fq {output}
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
        "outputs/assembled/{sample}_runx2/{sample}.fasta",
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

        run_spades "outputs/fastq/{wildcards.sample}_unmapped_baited1_runx2.fastq" "outputs/fastq/{wildcards.sample}_unmapped_baited2_runx2.fastq" "outputs/fastq/{wildcards.sample}_mapped_runx2.fastq" "outputs/assembled/{wildcards.sample}_runx2"
        mv outputs/assembled/{wildcards.sample}_runx2/contigs.fasta outputs/assembled/{wildcards.sample}_runx2/{wildcards.sample}.fasta

        run_spades "outputs/fastq/{wildcards.sample}_unmapped_baited1_zic2.fastq" "outputs/fastq/{wildcards.sample}_unmapped_baited2_zic2.fastq" "outputs/fastq/{wildcards.sample}_mapped_zic2.fastq" "outputs/assembled/{wildcards.sample}_zic2"
        mv outputs/assembled/{wildcards.sample}_zic2/contigs.fasta outputs/assembled/{wildcards.sample}_zic2/{wildcards.sample}.fasta

        run_spades "outputs/fastq/{wildcards.sample}_unmapped_baited1_foxl2.fastq" "outputs/fastq/{wildcards.sample}_unmapped_baited2_foxl2.fastq" "outputs/fastq/{wildcards.sample}_mapped_foxl2.fastq" "outputs/assembled/{wildcards.sample}_foxl2"
        mv outputs/assembled/{wildcards.sample}_foxl2/contigs.fasta outputs/assembled/{wildcards.sample}_foxl2/{wildcards.sample}.fasta

        run_spades "outputs/fastq/{wildcards.sample}_unmapped_baited1_arx.fastq" "outputs/fastq/{wildcards.sample}_unmapped_baited2_arx.fastq" "outputs/fastq/{wildcards.sample}_mapped_arx.fastq" "outputs/assembled/{wildcards.sample}_arx"
        mv outputs/assembled/{wildcards.sample}_arx/contigs.fasta outputs/assembled/{wildcards.sample}_arx/{wildcards.sample}.fasta
        """

rule orf_call:
    """
    This rule translates ORFs and then pulls out the peptide of interest
    """
    input:
        "outputs/assembled/{sample}_runx2/{sample}.fasta",
        "outputs/assembled/{sample}_zic2/{sample}.fasta", 
        "outputs/assembled/{sample}_foxl2/{sample}.fasta",
        "outputs/assembled/{sample}_arx/{sample}.fasta"
    output:
        "outputs/peptides/{sample}/{sample}_peptide_filt_runx2.fa",
        "outputs/peptides/{sample}/{sample}_peptide_filt_zic2.fa",
        "outputs/peptides/{sample}/{sample}_peptide_filt_foxl2.fa",
        "outputs/peptides/{sample}/{sample}_peptide_filt_arx.fa"
    conda:"envs/orfipy.yml"
    params:
        sequence_start_runx2="MRIPVDP",
        sequence_start_zic2="VHESSPQ",
        sequence_start_foxl2="MMASYPE",
        sequence_start_arx="LVAVHGT"
    shell: """
        orfipy outputs/assembled/{wildcards.sample}_runx2/{wildcards.sample}.fasta --pep orf_peptides_runx2.fa --outdir outputs/peptides/{wildcards.sample}/ --partial-3
        awk -v seq_start={params.sequence_start_runx2} 'BEGIN{{RS=">"}} $0 ~ seq_start {{print ">" $0}}' outputs/peptides/{wildcards.sample}/orf_peptides_runx2.fa > outputs/peptides/{wildcards.sample}/{wildcards.sample}_peptide_filt_runx2.fa
        orfipy outputs/assembled/{wildcards.sample}_zic2/{wildcards.sample}.fasta --pep orf_peptides_zic2.fa --outdir outputs/peptides/{wildcards.sample}/ --partial-3
        awk -v seq_start={params.sequence_start_zic2} 'BEGIN{{RS=">"}} $0 ~ seq_start {{print ">" $0}}' outputs/peptides/{wildcards.sample}/orf_peptides_zic2.fa > outputs/peptides/{wildcards.sample}/{wildcards.sample}_peptide_filt_zic2.fa
        orfipy outputs/assembled/{wildcards.sample}_foxl2/{wildcards.sample}.fasta --pep orf_peptides_foxl2.fa --outdir outputs/peptides/{wildcards.sample}/ --partial-3
        awk -v seq_start={params.sequence_start_foxl2} 'BEGIN{{RS=">"}} $0 ~ seq_start {{print ">" $0}}' outputs/peptides/{wildcards.sample}/orf_peptides_foxl2.fa > outputs/peptides/{wildcards.sample}/{wildcards.sample}_peptide_filt_foxl2.fa
        orfipy outputs/assembled/{wildcards.sample}_arx/{wildcards.sample}.fasta --pep orf_peptides_arx.fa --outdir outputs/peptides/{wildcards.sample}/ --partial-3
        awk -v seq_start={params.sequence_start_arx} 'BEGIN{{RS=">"}} $0 ~ seq_start {{print ">" $0}}' outputs/peptides/{wildcards.sample}/orf_peptides_arx.fa > outputs/peptides/{wildcards.sample}/{wildcards.sample}_peptide_filt_arx.fa        
        """

rule init_count_file:
    output:
        "outputs/counts/stretch_lengths.txt"
    shell:
        """
        echo -e "SampleName\tLongestStretch\tGene" > {output}
        """

rule count_expansion:
    """
    This rule finds the expansion and determines its length
    """
    input:
        runx2="outputs/peptides/{sample}/{sample}_peptide_filt_runx2.fa",
        zic2="outputs/peptides/{sample}/{sample}_peptide_filt_zic2.fa",
        foxl2="outputs/peptides/{sample}/{sample}_peptide_filt_foxl2.fa",
        arx="outputs/peptides/{sample}/{sample}_peptide_filt_arx.fa"
    output:
        "outputs/counts/{sample}_{gene}_count.txt"
    params:
        output_table="outputs/counts/stretch_lengths.txt"
    shell:
        """
        longest_stretch=$(grep -v ">" {input} | tr -d '\n' | grep -o 'A\\+' | awk '{{ print length($0) }}' | sort -n | tail -n 1)
        echo -e "{wildcards.sample}\t${{longest_stretch}}\t{wildcards.gene}_a" >> {params.output_table}
        echo -e "{wildcards.sample}\t${{longest_stretch}}\t{wildcards.gene}_a" > "{output}"
        """
