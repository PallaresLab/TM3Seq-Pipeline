import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main(file1, file2, if_SE, output):
    df1 = pd.read_table(file1, index_col=0)
    df2 = pd.read_table(file2, index_col=0)
    df = pd.merge(df1, df2, on="Sample")
    assert df.shape[0] == df2.shape[0] and df.shape[1] == df1.shape[1] + df2.shape[1]
    df_table = pd.DataFrame(index=df.index)
    if if_SE=="True":
        df_table['input_reads'] = df['input_reads'].astype(int)

    elif if_SE=="False":
        df_table['input_read_pairs'] = df['input_read_pairs'].astype(int)
        df_table['forward_only_surviving'] = df['forward_only_surviving'].astype(int)
        df_table['reverse_only_surviving'] = df['reverse_only_surviving'].astype(int)

    df_table['after_trimmed'] = df['surviving'].astype(int)
    
    if "STAR_mqc-generalstats-star-uniquely_mapped" in df.columns:
        df_table['uniquely_mapped'] = df['STAR_mqc-generalstats-star-uniquely_mapped'].astype(int)
    elif "uniquely_mapped" in df.columns:
        df_table['uniquely_mapped'] = df['uniquely_mapped'].astype(int)
    else:
        raise KeyError(f"Missing required column: uniquely_mapped")
    
    if "featureCounts_mqc-generalstats-featurecounts-Assigned" in df.columns:
        df_table['assigned'] = df['featureCounts_mqc-generalstats-featurecounts-Assigned'].astype(int)
    elif "Assigned" in df.columns:
        df_table['assigned'] = df['Assigned'].astype(int)
    else:
        raise KeyError(f"Missing required column: assigned")
    if if_SE=="True":
        df_table=df_table[~df_table.index.str.contains('R2')]
    df_table.to_csv(f'{output}_count.csv', index=True)
 
    plt.figure(figsize=(8, 6))
    df_percent=pd.DataFrame(index=df_table.index)
    df_percent["after_trimmed/input_reads"]=df_table["after_trimmed"]/df_table['input_reads']
    df_percent["uniquely_mapped/after_trimmed"]=df_table["uniquely_mapped"]/df_table['after_trimmed']
    df_percent["assigned/after_trimmed"]=df_table["assigned"]/df_table['after_trimmed']
    
    sns.boxplot(data=df_percent)
    plt.xticks(rotation=45)
    plt.title("Boxplot of Percentage")
    plt.ylabel("Percentage")
    plt.tight_layout()
    plt.savefig(f'{output}_percentage.pdf')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--folder', type=str, help='Path to the multiqc folder')
    parser.add_argument('-s', '--if_SE', type=str, help='if single end or paired end reads')
    parser.add_argument('-o', '--output', type=str, help='Path to the output file')

    args = parser.parse_args()
    main(f'{args.folder}/multiqc_general_stats.txt', f'{args.folder}/multiqc_trimmomatic.txt', args.if_SE, args.output)
