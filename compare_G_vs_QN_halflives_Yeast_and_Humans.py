
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
from matplotlib.patches import Patch
    
def main():

    output = open('TableS1_Yeast_and_Human_LCD_Protein_HalfLives.tsv', 'w')
    output.write('\t'.join( ['Organism', 'ORF Name', 'Common Name', 'LCD Class', 'Minimum LCD Length', 'Protein Half-life (mins for yeast, hrs for humans)'] ) + '\n')
    
    stats_file = open('Yeast_and_Human_ProteinHalfLife_Statistics.tsv', 'w')
    stats_file.write('\t'.join(['Organism', 'Minimum LCD Length', 'U-statistic', 'P-value']) + '\n')
        
    organisms = ['Yeast', 'Human']
    for organism in organisms:
        halflife_dict = pickle.load(open(organism + '_Protein_HalfLives.dat', 'rb'))
        df = {'LCD Class':[], 'Minimum LCD Length':[], 'Half-life':[]}

        grich_halflives = []
        qnrich_halflives = []
        for min_len in [20, 40]:
            for dataset in ['G', 'QN']:
                track_prots = set()
                if dataset == 'G':
                    h = open(organism + '_G_35percent_LCD-Composer_RESULTS.tsv')
                else:
                    h = open(organism + '_QN_50percent_LCD-Composer_RESULTS.tsv')
                for i in range(7):
                    h.readline()
                header = h.readline()
                
                common_names = []
                for line in h:
                    items = line.rstrip().split('\t')
                    if organism == 'Yeast':
                        accession, common_name, *junk = items[0].split(' ')
                    else:
                        junk, accession, common_name = items[0].split('|')
                        common_name, junk = common_name.split(' OS=')
                        common_name, *junk = common_name.split('_')
                        
                    seq = items[1]
                    
                    if len(seq) < min_len or accession not in halflife_dict:
                        continue
                        
                    if accession not in track_prots:
                        output.write('\t'.join( [organism, accession, common_name, dataset, str(min_len), str(halflife_dict[accession]) ] ) + '\n')
                            
                        df['LCD Class'].append(dataset)
                        if min_len == 20:
                            df['Minimum LCD Length'].append(0)
                        else:
                            df['Minimum LCD Length'].append(1)
                        df['Half-life'].append(float(halflife_dict[accession]))
                        track_prots.add(accession)
                h.close()

        plotting(grich_halflives, qnrich_halflives, 'Protein Half-life', df, organism)
        stats_file = run_stats(grich_halflives, qnrich_halflives, 'Protein Half-life', df, stats_file, organism)
        
    output.close()
    stats_file.close()

            
def plotting(grich_vals, qnrich_vals, data_type, df, organism):

    #PLOTTING
    if organism == 'Yeast':
        pal = ['#1f77b4', '#ff7f0e']
    else:
        pal=['#009292','#920000']

    ax = sns.boxplot(x='Minimum LCD Length', y='Half-life', hue='LCD Class', data=df, dodge=True, showfliers=False, palette=pal)
    ax = sns.stripplot(x='Minimum LCD Length', y='Half-life', hue='LCD Class', data=df, color='black', dodge=True, s=4, alpha=0.5)

    legend_elements = [Patch(facecolor=pal[0], edgecolor='0.2', label='G'),
                        Patch(facecolor=pal[1], edgecolor='0.2', label='Q/N')]
                        
    leg = plt.legend(handles=legend_elements, prop={'family':'Arial', 'size':14}, loc=2, bbox_to_anchor=(1,1), handletextpad=0.2)
    leg.set_title('LCD Class', prop={'family':'Arial', 'size':14})
    
    if organism == 'Yeast':
        plt.ylim(0, 3900)
        plt.ylabel(data_type + ' (mins)', fontname='Arial', fontsize=16)
    else:
        plt.ylim(0, 224)
        plt.ylabel(data_type + ' (hrs)', fontname='Arial', fontsize=16)
    
    plt.xlabel('Minimum LCD Length', fontname='Arial', fontsize=16)
    ax.set_xticklabels(['20aa', '40aa'])
    for label in ax.get_xticklabels():
        label.set_fontname('Arial')
        label.set_fontsize(14)
    for tick in ax.get_yticklabels():
        tick.set_fontname('Arial')
        tick.set_fontsize(14)
        
    plt.tight_layout()
    
    if organism == 'Yeast':
        plt.savefig('Fig1A - Yeast G-rich and QN-rich LCD Protein HalfLives.tiff', bbox_inches='tight', dpi=600)
    else:
        plt.savefig('Fig1B - Human G-rich and QN-rich LCD Protein HalfLives.tiff', bbox_inches='tight', dpi=600)

    plt.close()
    
            
def run_stats(grich_vals, qnrich_vals, data_type, df, stats_file, organism):

    #run Mann-Whitney U test
    df = pd.DataFrame.from_dict(df)
    for min_len_pos in [0, 1]:
        temp_df = df[df['Minimum LCD Length'] == min_len_pos]
        grich_vals = [temp_df.iloc[i]['Half-life'] for i in range(len(temp_df['Half-life'])) if temp_df.iloc[i]['LCD Class'] == 'G']
        qnrich_vals = [temp_df.iloc[i]['Half-life'] for i in range(len(temp_df['Half-life'])) if temp_df.iloc[i]['LCD Class'] in 'QN']
        results = stats.mannwhitneyu(grich_vals, qnrich_vals, alternative='two-sided')
        u_stat, pval = str(results).split(',')
        junk, u_stat = u_stat.split('=')
        junk, pval = pval[:-1].split('=')
        if min_len_pos == 0:
            stats_file.write('\t'.join( [organism, '20', u_stat, pval] ) + '\n')
        else:
            stats_file.write('\t'.join( [organism, '40', u_stat, pval] ) + '\n')
    
    return stats_file


if __name__ == '__main__':
    main()