rule GenomicsDBImport:
    input:
        expand(working_dir+"/HaplotypeCaller/{sample}/{idx}.g.vcf.gz", sample=samples.keys(), idx=IDX),

    
    output:
        directory("/dev/shm/pallares_lab/database_{idx}"),
    
    log:
        log_dir + "/JointCallSNPs/GenomicsDBImport_{idx}.log",   
    
    params:
        outdir = "/dev/shm/pallares_lab",
        input_params = lambda wildcards, input: " ".join(["-V "+ i for i in input]),
        interval = lambda wildcards: " ".join([f"-L {chrom}:{start+1}-{end}" for chrom, start, end in INTERVALS[int(wildcards.idx)]]),
        
    threads: 4
    
    conda:
        "../envs/gatk.yml"
    
    shell:
        "mkdir -p {params.outdir} && "
        "gatk --java-options '-Xmx100G' "
        "GenomicsDBImport "
        "--genomicsdb-workspace-path {output} "
        "{params.interval} "
        "{params.input_params} "
        "--batch-size 50 "
        "--reader-threads {threads} "
        ">{log} 2>&1 "

rule GenotypeGVCFs:
    input:
        R = working_dir + "/genome_prepare/"+ref_basename,
        db = "/dev/shm/pallares_lab/database_{idx}",
    
    output:
        working_dir + "/JointCallSNPs/{idx}.vcf.gz",
    
    log:
        log_dir + "/JointCallSNPs/GenotypeGVCFs_{idx}.log",     
        
    threads :1
    
    conda:
        "../envs/gatk.yml"
    
    shell:     
        "gatk --java-options '-Xmx100G' "
        "GenotypeGVCFs "
        "-R {input.R} "
        "-V gendb://{input.db} "
        "-O {output} && "
        "rm -rf {input.db} "
        ">{log} 3>&1 "


rule GatherVcfs:
    input:
       expand(working_dir + "/JointCallSNPs/{idx}.vcf.gz", idx=IDX)
    output:
        working_dir + "/JointCallSNPs/all_samples.vcf.gz"
    log:
        log_dir + "/JointCallSNPs/GatherVcfs.log",     
    
    params:
        input_params = lambda wildcards,input:" ".join(["I= "+ f for f in input])
        
    conda:
        "../envs/gatk.yml"

    shell: 
        "rm -rf /dev/shm/pallares_lab/ && "    
        "gatk GatherVcfs "
        "{params.input_params} O={output} && "
        "rm {input} && "
        "gatk IndexFeatureFile -I {output} "
        ">{log} 3>&1 " 
     