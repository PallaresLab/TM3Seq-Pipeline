rule VariantRecalibrator:
    input:
        V = working_dir + "/JointCallSNPs/all_samples.vcf.gz",
        R = working_dir + "/genome_prepare/"+ref_basename,
        L = config['bed_file'] ,
        vcf = config['dbSNP']
    
    output:
        recal = working_dir + "/VQSR/SNP.recal",
        tranches = working_dir + "/VQSR/output.tranches",
        r = working_dir + "/VQSR/recalibration_plots.R",
    
    log:
        log_dir + "/VQSR/VariantRecalibrator.log",
   
    params:
        outdir = working_dir+"/VQSR"
        
    conda:
        "../envs/vqsr.yml"

    shell:
        "mkdir -p {params.outdir} && "
        "gatk --java-options '-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=8' "
        "VariantRecalibrator "
        "-V {input.V} "
        "-R {input.R} "
        "-L {input.L} "
        "-resource:dbSNP_Nex_Sep28.19,known=false,training=true,truth=true,prior=10.0 {input.vcf} "
        "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 "
        "-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP "
        "-mode SNP "
        "--tranches-file {output.tranches} "
        "--rscript-file {output.r} "
        "--trust-all-polymorphic "
        "-O {output.recal} "
        ">{log} 2>&1 "
    
        
                
rule ApplyVQSR:
    input:
        V = working_dir + "/JointCallSNPs/all_samples.vcf.gz",
        R = working_dir + "/genome_prepare/"+ref_basename,
        recal = working_dir + "/VQSR/SNP.recal",
        tranches = working_dir + "/VQSR/output.tranches",
    
    output:
        working_dir + "/VQSR/SNP_recal.vcf.gz",

    
    log:
        log_dir + "/VQSR/ApplyVQSR.log",
   
    
    params:
        outdir = working_dir+"/VQSR"
        
    conda:
        "../envs/vqsr.yml"

    shell:
        "gatk --java-options '-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=8' "
        "ApplyVQSR "
        "-V {input.V} "
        "-R {input.R} "
        "--tranches-file {input.tranches} "
        "-mode SNP "
        "--recal-file {input.recal} "
        "--truth-sensitivity-filter-level 99.5 "
        "-O {output} "
        ">{log} 2>&1"