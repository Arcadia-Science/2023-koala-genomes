import pandas as pd
        
samples_df = pd.read_table('inputs/samples.tsv').set_index("AWSFileName",drop=False)
SAMPLE = list(samples_df['AWSFileName'])

rule all:
    input:
        "outputs/counts/stretch_lengths.txt",
        expand("outputs/counts/{sample}count.txt", sample=SAMPLE)
     
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
        samtools merge -u outputs/bams/{wildcards.sample}_unmapped.bam outputs/bams/{wildcards.sample}_unmap1.bam outputs/bams/{wildcards.sample}_unmap2.bam outputs/bams/{wildcards.sample}_unmap3.bam        
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
        samtools sort -n outputs/bams/{wildcards.sample}_unmapped.bam -o outputs/bams/{wildcards.sample}_unmapped_sort.bam
        bedtools bamtofastq -i outputs/bams/{wildcards.sample}_unmapped_sort.bam -fq outputs/fastq/{wildcards.sample}_unmapped_r1.fastq -fq2 outputs/fastq/{wildcards.sample}_unmapped_r2.fastq
        """

rule create_kmer_file:
    """
    This rule generates kmers from genes of interest
    """
    output:
        "outputs/kmers.fasta"
    conda:"envs/bbmap.yml"
    shell: """
        kmercountexact.sh in=inputs/RUNX2_koala.fasta  out=outputs/kmers.txt k=21
        mv outputs/kmers.txt outputs/kmers.fasta
        """

rule kmer_bait:
    """
    This rule filters reads based on kmers of interest
    """
    input:
        "outputs/fastq/{sample}_unmapped_r1.fastq",
        "outputs/fastq/{sample}_unmapped_r2.fastq",
        "outputs/kmers.fasta"
    output:
        "outputs/fastq/{sample}_unmapped_baited1.fastq",
        "outputs/fastq/{sample}_unmapped_baited2.fastq"
    conda:"envs/bbmap.yml"
    shell: """
        bbmap.sh ref=outputs/kmers.fasta in=outputs/fastq/{wildcards.sample}_unmapped_r1.fastq in2=outputs/fastq/{wildcards.sample}_unmapped_r2.fastq outm1=outputs/fastq/{wildcards.sample}_unmapped_baited1.fastq outm2=outputs/fastq/{wildcards.sample}_unmapped_baited2.fastq nodisk 
        """

rule mapped_reads:
    """
    This rule extracts the reads from the original bam that mapped to the gene of interest.
    """
    input:
        bam="outputs/bams/{sample}.bam"
    output:
        "outputs/fastq/{sample}_mapped.fastq"
    conda:"envs/samtools.yml"
    shell: """
        samtools view -b {input.bam} "MSTS01000055.1:12711350-12714427" > outputs/bams/{wildcards.sample}_mapped.bam
        bedtools bamtofastq -i outputs/bams/{wildcards.sample}_mapped.bam -fq outputs/fastq/{wildcards.sample}_mapped.fastq
        """

rule assemble:
    """
    This rule assembles the kmer-baited reads and mapped reads for the region of interest
    """
    input:
        "outputs/fastq/{sample}_mapped.fastq",
        "outputs/fastq/{sample}_unmapped_baited1.fastq",
        "outputs/fastq/{sample}_unmapped_baited2.fastq"
    output:
        "outputs/assembled/{sample}/{sample}.fasta"
    conda:"envs/spades.yml"
    shell: """
        spades.py --isolate -1 outputs/fastq/{wildcards.sample}_unmapped_baited1.fastq -2 outputs/fastq/{wildcards.sample}_unmapped_baited2.fastq -s outputs/fastq/{wildcards.sample}_mapped.fastq -o outputs/assembled/{wildcards.sample} -m 1024
        mv outputs/assembled/{wildcards.sample}/contigs.fasta outputs/assembled/{wildcards.sample}/{wildcards.sample}.fasta
        """

rule orf_call:
    """
    This rule translates ORFs and then pulls out the peptide of interest
    """
    input:
        "outputs/assembled/{sample}/{sample}.fasta"
    output:
        "outputs/peptides/{sample}/{sample}_peptide_filt.fa"
    conda:"envs/orfipy.yml"
    params:
        sequence_start="MRIPVDP"
    shell: """
        orfipy outputs/assembled/{wildcards.sample}/{wildcards.sample}.fasta --pep orf_peptides.fa --outdir outputs/peptides/{wildcards.sample}/ --partial-3
        awk -v seq_start={params.sequence_start} 'BEGIN{{RS=">"}} $0 ~ seq_start {{print ">" $0}}' outputs/peptides/{wildcards.sample}/orf_peptides.fa > outputs/peptides/{wildcards.sample}/{wildcards.sample}_peptide_filt.fa        
        """

rule init_count_file:
    output:
        "outputs/counts/stretch_lengths.txt"
    shell:
        """
        echo -e "SampleName\tLongestQStretch" > {output}
        """

rule count_expansion:
    """
    This rule finds the expansion and determines its length
    """
    input:
        peptide="outputs/peptides/{sample}/{sample}_peptide_filt.fa"
    output:
        "outputs/counts/{sample}count.txt"
    params:
        output_table="outputs/counts/stretch_lengths.txt"
    shell:
        """
        longest_stretch=$(grep -v ">" {input.peptide} | tr -d '\n' | grep -o 'Q\\+' | awk '{{ print length($0) }}' | sort -n | tail -n 1)
        echo -e "{wildcards.sample}\t${{longest_stretch}}" >> {params.output_table}
        echo -e "{wildcards.sample}\t${{longest_stretch}}" > "outputs/counts/{wildcards.sample}count.txt"
        """
