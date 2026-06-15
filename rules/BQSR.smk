rule BaseRecalibrator:
    input:
        bam = config['bam_dir']+"/{sample}.bam",
        R = working_dir + "/genome_prepare/"+ref_basename,
        fai = working_dir+"/genome_prepare/"+ref_basename+".fai",
        dict = working_dir+"/genome_prepare/"+ref_name+".dict",
        vcf = config['dbSNP']
    
    output:
        working_dir+"/BQSR/{sample}.table",
        
    log:
        log_dir+"/BQSR/BaseRecalibrator/{sample}.log",
    
    params:
        outdir = working_dir+"/BQSR"
        
    threads: 1
    
    conda:
        "../envs/gatk.yml"

    shell:
        "mkdir -p {params.outdir} && "
        "gatk --java-options '-Xmx20G' "
        "BaseRecalibrator "
        "-O {output} "
        "-R {input.R} "
        "-I {input.bam} "
        "--known-sites {input.vcf} "        
        ">{log} 2>&1"
       
        
rule ApplyBQSR:
    input:
        bam = config['bam_dir']+"/{sample}.bam",
        R = working_dir + "/genome_prepare/"+ref_basename,
        table = working_dir+"/BQSR/{sample}.table",
    
    output:
        working_dir+"/BQSR/{sample}_recal.bam"
        
    log:
        log_dir+"/BQSR/ApplyBQSR/{sample}.log"
        
        
    params:
        outdir = working_dir+"/BQSR"
        
    threads: 1
    
    conda:
        "../envs/gatk.yml"
    
    shell:        
        "gatk --java-options '-Xmx20G' "
        "ApplyBQSR "
        "-O {output} "
        "-R {input.R} "
        "-I {input.bam} "
        "--bqsr-recal-file {input.table} "
        ">{log} 2>&1"
