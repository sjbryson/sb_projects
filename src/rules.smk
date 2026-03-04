
from snakemake.utils import multiext

## -- Minimap2 Rules -- ##

rule PreIndexMinimap2SR:
    input:
        fasta = "{sample}.fasta"
    output:
        mmi = "{sample}.mmi"
    threads: 4
    shell:
        "minimap2 -k 21 -w 11 -t {threads} -d {output.mmi} {input.fasta}"


rule PreIndexMinimap2LR:
    input:
        fasta = "{sample}.fasta"
    output:
        mmi = "{sample}.mmi"
    threads: 4
    shell:
        "minimap2 -H -t {threads} -d {output.mmi} {input.fasta}"


rule Minimap2MapSRtoBam:
    input:
        mmi = "ref.mmi",
        r1 = "{sample}_R1.fastq",
        r2 = "{sample}_R2.fastq"
    output:
        bam = "{sample}.bam"
    threads: 8
    shell:
        "minimap2 -ax sr -t {threads} {input.mmi} {input.r1} {input.r2} | "
        "samtools view -bS - > {output.bam}"

rule Minimap2MapSRtoSortedBam:
    input:
        mmi = "ref.mmi",
        r1 = "{sample}_R1.fastq",
        r2 = "{sample}_R2.fastq"
    output:
        bam = "{sample}.bam"
    threads: 8
    shell:
        "minimap2 -ax sr -t {threads} {input.mmi} {input.r1} {input.r2} | "
        "samtools sort -bS - > {output.bam}"

rule Minimap2SRHumanDepletion:
    input:
        mmi = "human_ref.mmi",
        r1 = "{sample}_R1.fastq",
        r2 = "{sample}_R2.fastq"
    output:
        bam = "{sample}_human_mapped.bam"
    threads: 8
    shell:
        "minimap2 -ax sr --secondary=no -t {threads} {input.mmi} {input.r1} {input.r2} | "
        "samtools view -bS - > {output.bam}"


rule Minimap2MapLR5:
    input:
        mmi = "ref.mmi",
        reads = "{sample}.fastq"
    output:
        bam = "{sample}.bam"
    threads: 8
    shell:
        "minimap2 -ax asm5 -t {threads} {input.mmi} {input.reads} | samtools view -bS - > {output.bam}"


## -- Samtools Rules -- ##


rule SamtoolsSort:
    input:
        bam = "{sample}.bam"
    output:
        sorted_bam = "{sample}.sorted.bam"
    threads: 4
    shell:
        "samtools sort -@ {threads} -o {output.sorted_bam} {input.bam}"


rule SamtoolsIndex:
    input:
        sorted_bam = "{sample}.sorted.bam"
    output:
        bai = "{sample}.sorted.bam.bai"
    shell:
        "samtools index {input.sorted_bam}"


rule SamtoolsMpileup:
    input:
        fasta = "ref.fasta",
        bam = "{sample}.sorted.bam"
    output:
        pileup = "{sample}.pileup"
    shell:
        "samtools mpileup -f {input.fasta} {input.bam} > {output.pileup}"


rule SamtoolsStats:
    input:
        bam = "{sample}.sorted.bam"
    output:
        txt = "{sample}.stats.txt"
    shell:
        "samtools stats {input.bam} > {output.txt}"


rule GetProperMappedIDs:
    input:
        bam = "{sample}.sorted.bam"
    output:
        ids = "{sample}.proper_ids.txt"
    threads: 2
    shell:
        "samtools view -f 1 -F 12 -@ {threads} {input.bam} | cut -f1 | sort -u > {output.ids}"


## -- Blast Rules -- ##


rule MakeBlastnDB:
    input:
        fasta = "{db_name}.fasta"
    output:
        db = multiext("{db_name}", ".nhr", ".nin", ".nsq")
    params:
        db_path = "{db_name}"
    shell:
        "makeblastdb -in {input.fasta} -input_type fasta -dbtype nucl -out {params.db_path} -parse_seqids -title {wildcards.db_name}"


rule Blastn:
    input:
        fasta = "{sample}.fasta",
        db = multiext("refs/my_db", ".nhr", ".nin", ".nsq")
    output:
        tsv = "{sample}.blast.tsv"
    threads: 4
    params:
        db_prefix = "refs/my_db"
    shell:
        "blastn -db {params.db_prefix} -query {input.fasta} -out {output.tsv} -outfmt 6 -num_threads {threads}"


## -- SB_Tools & Utils -- ##


#rule SBFilterFastq:
#    input:
#        r1 = "{sample}_R1.fastq",
#        r2 = "{sample}_R2.fastq",
#        filter_list = "{sample}.ids.txt"
#    output:
#        r1_filt = "{sample}_filtered_R1.fastq",
#        r2_filt = "{sample}_filtered_R2.fastq"
#    params:
#        prefix = "{sample}_filtered",
#        gz = "--gz", # Or empty string if False
#        invert = ""   # Or "--invert" if True
#    shell:
#        "sb_filter_fastq {params.gz} {params.invert} --r1 {input.r1} --r2 {input.r2} "
#        "--filter {input.filter_list} --out-prefix {params.prefix}"