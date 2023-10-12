import subprocess
import sys
import json

def run_command(cmd):
    subprocess.run(cmd, shell=True, executable="/bin/bash")

def main(sample, folder, gene_data_json):
    # Load gene_data from the JSON string
    gene_data = json.loads(gene_data_json)

    bam = f"outputs/downloadbams/{sample}.bam"
    bai = f"outputs/downloadbams/{sample}.bam.bai"

    # Download files
    run_command(f"s5cmd cp --show-progress s3://koalagenomes/{folder}/bam/{sample}.bam {bam}")
    run_command(f"s5cmd cp --show-progress s3://koalagenomes/{folder}/bam/{sample}.bam.bai {bai}")

    # Process each gene
    for gene, info in gene_data.items():
        command1 = f"samtools view -b -@ 16 {bam} {info['coordinates']} > outputs/bams/{sample}_mapped_{gene}.bam"
        command2 = f"bedtools bamtofastq -i outputs/bams/{sample}_mapped_{gene}.bam -fq outputs/fastq/{sample}_mapped_{gene}.fastq"
        run_command(command1)
        run_command(command2)

    # Cleanup
    run_command(f"rm {bam}")
    run_command(f"rm {bai}")

if __name__ == "__main__":
    sample = sys.argv[1]
    folder = sys.argv[2]
    gene_data_json = sys.argv[3]
    main(sample, folder, gene_data_json)
