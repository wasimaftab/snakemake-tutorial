# ### Step 1: Mapping reads
# rule bwa_map:
#     input:
#         "data/genome.fa",
#         "data/samples/A.fastq"
#     output:
#         "mapped_reads/A.bam"
#     shell:
#         "bwa mem {input} | samtools view -Sb - > {output}"


SAMPLES = ["A", "B"]

# ### Step 2: Generalizing the read mapping rule
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"


### Step 3: Sorting read alignments
rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} -O bam {input} > {output}"


### Step 4: Indexing read alignments and visualizing the DAG of jobs
rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"


### Step 5: Calling genomic variants
rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"