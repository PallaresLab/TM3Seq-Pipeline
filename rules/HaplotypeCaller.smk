rule HaplotypeCaller:
  input:
    R = working_dir+"/genome_prepare/"+ref_basename,
    table = working_dir+"/BQSR/{sample}.table",
    recal_bam = working_dir+"/BQSR/{sample}_recal.bam",
    
  
  output:
    working_dir+"/HaplotypeCaller/{sample}/{idx}.g.vcf.gz",
    
  log:
    log_dir+"/HaplotypeCaller/{sample}/{idx}.log",
    
  params:
    outdir=working_dir+"/HaplotypeCaller/{sample}",
    interval = lambda wildcards: " ".join([f"-L {chrom}:{start+1}-{end}" for chrom, start, end in INTERVALS[int(wildcards.idx)]])
    
  threads :1
  
  conda:
    "../envs/gatk.yml"
  
  shell:
    "mkdir -p {params.outdir} && "
    "gatk --java-options '-Xmx20G' "
    "HaplotypeCaller "
    "-O {output} "
    "-R {input.R} "
    "-I {input.recal_bam} "
    "{params.interval} "
    "-ERC GVCF "
    ">{log} 2>&1 "
    
    


