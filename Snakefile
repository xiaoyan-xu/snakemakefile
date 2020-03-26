
configfile: "config.yaml"
	

rule all:
  input:
    "count/merge.count.matrix"

rule index:
  input:
    "reference/genome.fna"
  output:
    "reference/genome.1.ht2"
  shell:
    "hisat2-build {input} reference/genome"

rule quality_check:
  input:
    lambda wildcards: config["samples"][wildcards.sample]
  output:
    "qc/"
  shell:
    "fastqc -o {output} {input}"

rule hisat2:
  input:
    fastq = lambda wildcards: config["samples"][wildcards.sample],
    index = "reference/genome.1.ht2"
  output:
    "bwa/{sample}.sam"
  params:
    IND= "reference/genome"
  shell:
    "hisat2 -x {params.IND} -U {input.fastq} -S {output}"

rule bam:
  input:
    index = "reference/genome.1.ht2",
    sam = "bwa/{sample}.sam"
  output:
    "bam/{sample}.bam"
  shell:
    "samtools view -Sb {input.sam} > {output}"

rule sort:
  input:
    "bam/{sample}.bam"
  output:
    "sort/{sample}_sort.bam"
  shell:
    "samtools sort {input} -o {output}"


rule count:
  input:
    sample="sort/{sample}_sort.bam",
    annotation=config["annotation"]
  output:
    "count/{sample}.count"
  shell:
    "htseq-count -t CDS -f bam -i name -s no -m intersection-nonempty {input.sample}  {input.annotation} > {output}"

rule matrix:
  input:
    expand("count/{sample}.count",sample=config["samples"])
  output:
    "count/merge.count.matrix"
  run:
    """
    awk 'BEGIN {{OFS="\t"}} NR==FNR{{if(FNR==1){{f="Gene\t"FILENAME}};a[FNR]=$1;gene[$1]=$2;next}} {{if(FNR==1){{f=f"\t"FILENAME}};gene[$1]=gene[$1]"\t"$2}} END{{print f;for (i in a) print a[i]"\t"gene[a[i]]}}' count/*.count > {output}
    """
