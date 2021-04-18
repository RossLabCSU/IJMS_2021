
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import pandas as pd
from scipy import stats

deg_propensities = { 'A':0.2665627892708544, 'R':-0.7412501624697525, 'N':-1.1869752279573649, 'D':-0.728512694019163,
                    'C':0.4206842379683707, 'E':-0.7157286901525158, 'Q':-0.7233002187306328, 'G':-0.4440764148760815,
                    'H':-0.028047412555156075, 'I':0.8799390174848035, 'L':1.046555755328647, 'K':-0.8574121141055006,
                    'M':1.038969401500499, 'F':0.22392755762641842, 'P':-0.10366937834956014, 'S':-0.16794980890760758,
                    'T':-0.03867551503137706, 'W':0.12078754467668205, 'Y':0.25418799347663784, 'V':0.934778908124434 }
                    
organisms = ['Yeast', 'Human']
                
def main():

    stats_file = open('Yeast_Human_G-rich_vs_QN-rich_DegradationScoreStatistics.tsv', 'w')
    stats_file.write('Organism\tU Statistic\tP-value\n')
    for organism in organisms:
        df = {'Score':[], 'Category':[], 'Position':[]}
        leg_labels = []
        for res_set in [('G', '35'), ('QN', '50')]:
            res = res_set[0]
            comp = res_set[1]
            leg_labels.append( res )
            h = open(organism + '_' + res + '_' + comp + 'percent_LCD-Composer_RESULTS.tsv')
            for i in range(7):
                h.readline()
            header = h.readline()
            for line in h:
                items = line.rstrip().split('\t')
                id, *junk = items[0].split(' ')
                seq = items[1]
                if len(seq) < 40:
                    continue
                score = sum( [deg_propensities[aa] for aa in seq if aa not in res] )
                df['Score'].append(score)
                df['Category'].append(res)
                df['Position'].append(0)

            h.close()

        boxplot(df, leg_labels, organism)
        stats_file = run_stats(df, organism, stats_file)
        
    stats_file.close()
        
        
def boxplot(df, leg_labels, organism):

    if organism == 'Yeast':
        pal = ['#1f77b4', '#ff7f0e']
    else:
        pal=['#009292','#920000']
    pandas_df = pd.DataFrame.from_dict(df)
    ax = sns.boxplot(x='Position', y='Score', data=df, hue='Category', showfliers=False, palette=pal)
    
    sns.stripplot(x='Position', y='Score', data=df, hue='Category', dodge=True, color='0.2', alpha=0.5)
    plt.xticks([0], labels=[organism], fontname='Arial', fontsize=16)
    plt.ylabel('Degradation Score', fontname='Arial', fontsize=16)
    
    if organism == 'Yeast':
        plt.yticks([x for x in range(-15, 25, 5)], fontname='Arial', fontsize=16)
        ax.set_ylim(top=27)
    else:
        plt.yticks([x for x in range(-125, 75, 25)], fontname='Arial', fontsize=16)
        ax.set_ylim(top=100)
    
    legend_elements = [Patch(facecolor=pal[i], edgecolor='0.2', label=leg_labels[i].replace('QN', 'Q/N')) for i in range(len(leg_labels))]
    leg = plt.legend(handles=legend_elements, loc=9, prop={'family':'Arial', 'size':14}, labelspacing=0.1)
    leg.set_title('LCD Class', prop={'family':'Arial', 'size':16})
    
    fig = plt.gcf()
    fig.set_size_inches(3, 6)
    if organism == 'Yeast':
        plt.savefig('Fig 2A - Degradation Scores Yeast G-rich vs QN-rich LCDs.tiff', bbox_inches='tight', dpi=600)
    else:
        plt.savefig('Fig 2B - Degradation Scores Human G-rich vs QN-rich LCDs.tiff', bbox_inches='tight', dpi=600)
    plt.close()

    
def run_stats(df, organism, stats_file):

    temp_df = pd.DataFrame.from_dict(df)
    
    grich_vals = [temp_df.iloc[i]['Score'] for i in range(len(temp_df['Score'])) if temp_df.iloc[i]['Category'] == 'G']
    qnrich_vals = [temp_df.iloc[i]['Score'] for i in range(len(temp_df['Score'])) if temp_df.iloc[i]['Category'] in 'QN']
    
    results = stats.mannwhitneyu(grich_vals, qnrich_vals, alternative='two-sided')
    u_stat, pval = str(results).split(',')
    junk, u_stat = u_stat.split('=')
    junk, pval = pval[:-1].split('=')

    stats_file.write(organism + '\t' + u_stat + '\t' + pval + '\n')

    return stats_file
    
    
if __name__ == '__main__':
    main()
    