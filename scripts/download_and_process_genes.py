import subprocess
import sys
import json

def run_command(cmd):
    subprocess.run(cmd, shell=True, executable="/bin/bash")

def main(sample, folder, gene_data_file):
    # Load gene_data from the JSON file
    with open(gene_data_file, "r") as f:
        gene_data = json.load(f)

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
    sample, folder, gene_data_file = sys.argv[1], sys.argv[2], sys.argv[3]
    main(sample, folder, gene_data_file)
