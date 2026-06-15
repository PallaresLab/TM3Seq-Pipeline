shell.executable("bash")
import os

from snakemake.utils import min_version
import multiprocessing
min_version("8.10.7")

configfile: "config.yaml"
log_dir = config['output_dir'] + "/logs"
working_dir = config['output_dir'] + "/working"
ref_basename=os.path.basename(config['reference_genome'])
ref_name=os.path.splitext(os.path.basename(config['reference_genome']))[0]
sample_files = snakemake.utils.listfiles(config["bam_dir"]+"/{sample}.bam")
samples = dict((y[0], x) for x, y in sample_files)

#singularity: "docker://broadinstitute/gatk"

def process_intervals(file_path, max_size=500000):
    result = []
    small_group = []
    small_total = 0

    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 3:
                continue  # skip lines

            chrom, start, end = parts
            start, end = int(start), int(end)
            size = end - start

            if size >= max_size:
                # Split into chunks of size max_size
                while start < end:
                    chunk_end = min(start + max_size, end)
                    result.append([[chrom, start, chunk_end]])
                    start = chunk_end
            else:
                # Add to small group for merging
                small_group.append([chrom, start, end])
                small_total += size

                if small_total >= max_size:
                    # Flush the current group
                    result.append(small_group)
                    small_group = []
                    small_total = 0

    # Add any remaining small intervals
    if small_group:
        result.append(small_group)

    return result


INTERVALS= process_intervals(config['bed_file'])
IDX=list(range(len(INTERVALS)))
print(IDX)
"""
INTERVALS = []

for n in result:
    l = []
    for i in n:
        chrom, start, end = i
        l.append(f" -L {chrom}:{start+1}-{end}")  # Format string
    INTERVALS.append(l)




INTERVALS = [f"-L {name}:{start+1}-{end}" for name, length in gd.items() for start, end in length]
big_stuff={"2L","2R","3L","3R","X","Y","4"}
bigints=[interv for interv in INTERVALS if interv.split(":")[0].split(" ")[1] in big_stuff]
smallints=[interv for interv in INTERVALS if interv.split(":")[0].split(" ")[1] not in big_stuff]
bigints.append(" ".join(smallints))

print(len(INTERVALS),len(bigints),len(smallints))
INTERVALS=bigints
realints=[x for x in INTERVALS]

"""
rule all:
    input:
        working_dir + "/JointCallSNPs/all_samples.vcf.gz"


include: "rules/genome_prepare.smk"
include: "rules/BQSR.smk"
include: "rules/HaplotypeCaller.smk"
include: "rules/JointCallSNPs.smk"
#include: "rules/VQSR.smk"




