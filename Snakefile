import pandas as pd
        
samples_df = pd.read_table('inputs/samples.tsv').set_index("AWSFileName",drop=False)
SAMPLE = list(samples_df['AWSFileName'])

rule all:
    input:
        expand("outputs/bams/{sample}.bam",sample=SAMPLE)

rule download_bams:
    """
    This rule downloads bams for 430 koala genomes
    """
    output: 
        "outputs/bams/{sample}.bam.bai",
        "outputs/bams/{sample}.bam"
    shell: """
        export AWS_REGION=ap-southeast-2
        s5cmd run s5cmds.txt
        """
