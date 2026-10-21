if RNA_INPUT_MODE == "fastq":
    rule fastp_pe:
        input:
            preflight=INPUT_PREFLIGHT_COMPLETE,
            read1=lambda wildcards: str(_sample_record(wildcards).read1),
            read2=lambda wildcards: str(_sample_record(wildcards).read2),
        output:
            read1=f"{READS_DIR}/PE/{{sample}}/{{sample}}.clean.R1.fastq.gz",
            read2=f"{READS_DIR}/PE/{{sample}}/{{sample}}.clean.R2.fastq.gz",
            json=f"{READS_DIR}/PE/{{sample}}/{{sample}}.fastp.json",
            html=f"{READS_DIR}/PE/{{sample}}/{{sample}}.fastp.html",
        log:
            f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/fastp/PE/{{sample}}.log",
        threads: workflow.cores
        conda:
            "../../envs/fastp.yaml"
        shell:
            r"""
            fastp \
              --in1 {input.read1:q} --in2 {input.read2:q} \
              --out1 {output.read1:q} --out2 {output.read2:q} \
              --json {output.json:q} --html {output.html:q} \
              --report_title {wildcards.sample:q} \
              --detect_adapter_for_pe --thread {threads} \
              > {log:q} 2>&1
            """
    rule fastp_se:
        input:
            preflight=INPUT_PREFLIGHT_COMPLETE,
            read1=lambda wildcards: str(_sample_record(wildcards).read1),
        output:
            read1=f"{READS_DIR}/SE/{{sample}}/{{sample}}.clean.fastq.gz",
            json=f"{READS_DIR}/SE/{{sample}}/{{sample}}.fastp.json",
            html=f"{READS_DIR}/SE/{{sample}}/{{sample}}.fastp.html",
        log:
            f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/fastp/SE/{{sample}}.log",
        threads: workflow.cores
        conda:
            "../../envs/fastp.yaml"
        shell:
            r"""
            fastp \
              --in1 {input.read1:q} --out1 {output.read1:q} \
              --json {output.json:q} --html {output.html:q} \
              --report_title {wildcards.sample:q} --thread {threads} \
              > {log:q} 2>&1
            """
    rule hisat2_build:
        input:
            preflight=INPUT_PREFLIGHT_COMPLETE,
            genome=lambda wildcards: config["genome"],
        output:
            HISAT2_INDEX,
        params:
            prefix=lambda wildcards, output: str(output[0])[:-len(".1.ht2")],
        log:
            f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/hisat2/build.log",
        threads: workflow.cores
        conda:
            "../../envs/alignment.yaml"
        shell:
            r"""
            hisat2-build -p {threads} {input.genome:q} {params.prefix:q} \
              > {log:q} 2>&1
            """
    rule hisat2_align_pe:
        input:
            index=HISAT2_INDEX,
            read1=f"{READS_DIR}/PE/{{sample}}/{{sample}}.clean.R1.fastq.gz",
            read2=f"{READS_DIR}/PE/{{sample}}/{{sample}}.clean.R2.fastq.gz",
        output:
            bam=f"{ALIGNMENT_DIR}/PE/{{sample}}/{ASSEMBLY_ID}.{{sample}}.sorted.bam",
            csi=f"{ALIGNMENT_DIR}/PE/{{sample}}/{ASSEMBLY_ID}.{{sample}}.sorted.bam.csi",
            summary=f"{ALIGNMENT_DIR}/PE/{{sample}}/{ASSEMBLY_ID}.{{sample}}.hisat2.summary.txt",
            flagstat=f"{ALIGNMENT_DIR}/PE/{{sample}}/{ASSEMBLY_ID}.{{sample}}.flagstat.txt",
            idxstats=f"{ALIGNMENT_DIR}/PE/{{sample}}/{ASSEMBLY_ID}.{{sample}}.idxstats.tsv",
        params:
            index_prefix=lambda wildcards, input: str(input.index[0])[:-len(".1.ht2")],
            strand=_hisat2_strand_option,
            align_threads=lambda wildcards, threads: max(1, threads - 1),
            sort_threads=lambda wildcards, threads: max(0, threads - 1),
        log:
            f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/hisat2/PE/{{sample}}.stderr.log",
        threads: workflow.cores
        conda:
            "../../envs/alignment.yaml"
        shell:
            r"""
            set -o pipefail
            if [ {threads} -eq 1 ]; then
              # A one-CPU allocation cannot run an aligner and sorter concurrently.
              sam={output.bam:q}.unsorted.sam
              trap 'rm -f -- "$sam"' EXIT
              hisat2 --dta {params.strand} -p 1 \
                -x {params.index_prefix:q} -1 {input.read1:q} -2 {input.read2:q} \
                --summary-file {output.summary:q} -S "$sam" 2> {log:q}
              samtools sort -@ 0 -m 1G -T {output.bam:q}.sort -o {output.bam:q} "$sam"
            else
              hisat2 --dta {params.strand} -p {params.align_threads} \
                -x {params.index_prefix:q} -1 {input.read1:q} -2 {input.read2:q} \
                --summary-file {output.summary:q} 2> {log:q} \
              | samtools sort -@ 0 -m 1G \
                -T {output.bam:q}.sort -o {output.bam:q} -
            fi
            samtools quickcheck -v {output.bam:q}
            samtools index -@ {params.sort_threads} -c {output.bam:q} {output.csi:q}
            samtools flagstat -@ {params.sort_threads} {output.bam:q} > {output.flagstat:q}
            samtools idxstats {output.bam:q} > {output.idxstats:q}
            """
    rule hisat2_align_se:
        input:
            index=HISAT2_INDEX,
            read1=f"{READS_DIR}/SE/{{sample}}/{{sample}}.clean.fastq.gz",
        output:
            bam=f"{ALIGNMENT_DIR}/SE/{{sample}}/{ASSEMBLY_ID}.{{sample}}.sorted.bam",
            csi=f"{ALIGNMENT_DIR}/SE/{{sample}}/{ASSEMBLY_ID}.{{sample}}.sorted.bam.csi",
            summary=f"{ALIGNMENT_DIR}/SE/{{sample}}/{ASSEMBLY_ID}.{{sample}}.hisat2.summary.txt",
            flagstat=f"{ALIGNMENT_DIR}/SE/{{sample}}/{ASSEMBLY_ID}.{{sample}}.flagstat.txt",
            idxstats=f"{ALIGNMENT_DIR}/SE/{{sample}}/{ASSEMBLY_ID}.{{sample}}.idxstats.tsv",
        params:
            index_prefix=lambda wildcards, input: str(input.index[0])[:-len(".1.ht2")],
            strand=_hisat2_strand_option,
            align_threads=lambda wildcards, threads: max(1, threads - 1),
            sort_threads=lambda wildcards, threads: max(0, threads - 1),
        log:
            f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/hisat2/SE/{{sample}}.stderr.log",
        threads: workflow.cores
        conda:
            "../../envs/alignment.yaml"
        shell:
            r"""
            set -o pipefail
            if [ {threads} -eq 1 ]; then
              # A one-CPU allocation cannot run an aligner and sorter concurrently.
              sam={output.bam:q}.unsorted.sam
              trap 'rm -f -- "$sam"' EXIT
              hisat2 --dta {params.strand} -p 1 \
                -x {params.index_prefix:q} -U {input.read1:q} \
                --summary-file {output.summary:q} -S "$sam" 2> {log:q}
              samtools sort -@ 0 -m 1G -T {output.bam:q}.sort -o {output.bam:q} "$sam"
            else
              hisat2 --dta {params.strand} -p {params.align_threads} \
                -x {params.index_prefix:q} -U {input.read1:q} \
                --summary-file {output.summary:q} 2> {log:q} \
              | samtools sort -@ 0 -m 1G \
                -T {output.bam:q}.sort -o {output.bam:q} -
            fi
            samtools quickcheck -v {output.bam:q}
            samtools index -@ {params.sort_threads} -c {output.bam:q} {output.csi:q}
            samtools flagstat -@ {params.sort_threads} {output.bam:q} > {output.flagstat:q}
            samtools idxstats {output.bam:q} > {output.idxstats:q}
            """
else:
    rule import_rna_bam:
        input:
            code=source_inputs(workflow.basedir, "rna_import"),
            preflight=INPUT_PREFLIGHT_COMPLETE,
            manifest=f"{INPUT_PREFLIGHT_DIR}/input_manifest.tsv",
            samples=f"{INPUT_PREFLIGHT_DIR}/normalized_samples.tsv",
            bam=lambda wildcards: str(_sample_record(wildcards).bam),
            csi=lambda wildcards: str(_sample_record(wildcards).csi),
        output:
            bam=f"{ALIGNMENT_DIR}/{{layout}}/{{sample}}/{ASSEMBLY_ID}.{{sample}}.sorted.bam",
            csi=f"{ALIGNMENT_DIR}/{{layout}}/{{sample}}/{ASSEMBLY_ID}.{{sample}}.sorted.bam.csi",
            audit=f"{ALIGNMENT_DIR}/{{layout}}/{{sample}}/{ASSEMBLY_ID}.{{sample}}.sorted.bam.import.json",
        params:
            repo_root=str(workflow.basedir),
        log:
            f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/rna_import/{{layout}}/{{sample}}.log",
        conda:
            "../../envs/preflight.yaml"
        script:
            "../scripts/import_rna_bam.py"

    rule rna_import:
        input:
            bams=SAMPLE_BAMS,
            csis=SAMPLE_CSIS,
            audits=RNA_IMPORT_AUDITS,

rule stringtie_sample:
    input:
        bam=_sample_bam,
        csi=_sample_csi,
    output:
        gtf=f"{TRANSCRIPTS_DIR}/samples/{{sample}}/{ASSEMBLY_ID}.{{sample}}.stringtie.gtf",
    params:
        strand=_stringtie_strand_option,
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/stringtie/{{sample}}.log",
    threads: workflow.cores
    conda:
        "../../envs/stringtie.yaml"
    shell:
        r"""
        stringtie {params.strand} -p {threads} -l {wildcards.sample:q} \
          -o {output.gtf:q} {input.bam:q} > {log:q} 2>&1
        test -s {output.gtf:q}
        awk -F '\t' '$0 !~ /^#/ && $3 == "transcript" {{found=1; exit}} END {{exit !found}}' \
          {output.gtf:q}
        """
rule stringtie_merge:
    input:
        gtfs=expand(
            f"{TRANSCRIPTS_DIR}/samples/{{sample}}/{ASSEMBLY_ID}.{{sample}}.stringtie.gtf",
            sample=SAMPLES,
        ),
    output:
        gtf=f"{TRANSCRIPTS_DIR}/merged/{ASSEMBLY_ID}.stringtie.merged.gtf",
        merge_list=f"{TRANSCRIPTS_DIR}/merged/{ASSEMBLY_ID}.stringtie.merge_list.txt",
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/stringtie/merge.log",
    threads: workflow.cores
    conda:
        "../../envs/stringtie.yaml"
    shell:
        r"""
        printf '%s\n' {input.gtfs:q} > {output.merge_list:q}
        stringtie --merge -p {threads} -o {output.gtf:q} {output.merge_list:q} \
          > {log:q} 2>&1
        test -s {output.gtf:q}
        awk -F '\t' '$0 !~ /^#/ && $3 == "transcript" {{found=1; exit}} END {{exit !found}}' \
          {output.gtf:q}
        """
if RNA_INPUT_MODE == "fastq":
    rule fastp_version:
        output:
            f"{TRANSCRIPTS_DIR}/software/fastp.version.txt",
        log:
            f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/software/fastp_version.log",
        conda:
            "../../envs/fastp.yaml"
        shell:
            r"""
            fastp --version > {output:q} 2>&1
            test -s {output:q}
            echo 'fastp version capture PASS' > {log:q}
            """


    rule alignment_versions:
        output:
            f"{TRANSCRIPTS_DIR}/software/alignment_versions.txt",
        log:
            f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/software/alignment_versions.log",
        conda:
            "../../envs/alignment.yaml"
        shell:
            r"""
            {{
              hisat2 --version
              samtools --version
            }} > {output:q} 2>&1
            test -s {output:q}
            echo 'alignment version capture PASS' > {log:q}
            """


rule stringtie_version:
    output:
        f"{TRANSCRIPTS_DIR}/software/stringtie.version.txt",
    log:
        f"{EXECUTION_ROOT}/{EXECUTION_ID}/logs/software/stringtie_version.log",
    conda:
        "../../envs/stringtie.yaml"
    shell:
        r"""
        stringtie --version > {output:q} 2>&1
        test -s {output:q}
        echo 'StringTie version capture PASS' > {log:q}
        """
rule rna_evidence:
    input:
        merged=f"{TRANSCRIPTS_DIR}/merged/{ASSEMBLY_ID}.stringtie.merged.gtf",
        mapping_versions=(
            [f"{TRANSCRIPTS_DIR}/software/fastp.version.txt",
             f"{TRANSCRIPTS_DIR}/software/alignment_versions.txt"]
            if RNA_INPUT_MODE == "fastq" else []
        ),
        import_audits=RNA_IMPORT_AUDITS,
        stringtie_version=f"{TRANSCRIPTS_DIR}/software/stringtie.version.txt",
        bams=expand(
            f"{ALIGNMENT_DIR}/PE/{{sample}}/{ASSEMBLY_ID}.{{sample}}.sorted.bam",
            sample=PE_SAMPLES,
        ) + expand(
            f"{ALIGNMENT_DIR}/SE/{{sample}}/{ASSEMBLY_ID}.{{sample}}.sorted.bam",
            sample=SE_SAMPLES,
        ),
