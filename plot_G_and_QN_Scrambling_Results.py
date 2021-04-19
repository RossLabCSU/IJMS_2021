
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import pandas as pd
from scipy import stats
import statistics
import pickle

combined_deg_propensities = { 'A':0.2665627892708544, 'R':-0.7412501624697525, 'N':-1.1869752279573649, 'D':-0.728512694019163,
                                'C':0.4206842379683707, 'E':-0.7157286901525158, 'Q':-0.7233002187306328, 'G':-0.4440764148760815,
                                'H':-0.028047412555156075, 'I':0.8799390174848035, 'L':1.046555755328647, 'K':-0.8574121141055006,
                                'M':1.038969401500499, 'F':0.22392755762641842, 'P':-0.10366937834956014, 'S':-0.16794980890760758,
                                'T':-0.03867551503137706, 'W':0.12078754467668205, 'Y':0.25418799347663784, 'V':0.934778908124434 }
                
def main():

    output = open('PairedTtests_Scrambled_G_vs_QN_LCDs_Yeast_and_Humans.tsv', 'w')
    output.write('Organism\tLCD Class\tComparison Group\tCluster Size\tt Statistic (paired t-test)\tP-value\n')
        
    for proteome in ['Yeast', 'Human']:
        
        for cluster_size in [4, 8]:
            scramble_score_df = pickle.load(open(proteome + '_Grich_and_QNrich_Scrambled_' + str(cluster_size) + 'aaCluster_Scores_10kScrambleIterations.dat', 'rb'))
                                
            df = {'Score':[], 'LCD Class':[], 'Position':[]}
            df = cluster_analysis(cluster_size, df, proteome, scramble_score_df)
            
            #MANUALLY ADD PRE-COMPUTED SCRAMBLING RESULTS TO MAIN df
            for res in ['G', 'QN']:
                for i in range(len(scramble_score_df[res])):
                    df['Score'].append( scramble_score_df[res][i] )
                    df['LCD Class'].append('Scrambled ' + res)
                    df['Position'].append(0)
            
            boxplot(df, proteome, cluster_size)
            output = run_stats(df, proteome, output, cluster_size)

    output.close()
        
    
def cluster_analysis(cluster_size, df, proteome, scramble_score_df):

    leg_labels = []
    for res_set in [('G', '35'), ('QN', '50')]:
        res = res_set[0]
        comp = res_set[1]
        leg_labels.append( res )
        
        h = open(proteome + '_' + res + '_' + comp + 'percent_LCD-Composer_RESULTS.tsv')
        for i in range(7):
            h.readline()
        header = h.readline()
        for line in h:
            items = line.rstrip().split('\t')
            seq = items[1]
            if len(seq) < 40:
                continue
            max_score = calc_max_cluster(seq, cluster_size, combined_deg_propensities, res)
            df['Score'].append(max_score)
            df['LCD Class'].append(res)
            df['Position'].append(0)

    return df
    
    
def calc_max_cluster(seq, cluster_size, deg_propensities, res):
    
    max_score = -100
    for i in range(0, len(seq)-cluster_size+1):
        window = seq[i:i+cluster_size]
        score = sum( [deg_propensities[aa] for aa in window if aa not in res] )
        if score > max_score:
            max_score = score
            
    return max_score
    

def boxplot(df, proteome, cluster_size):

    if proteome == 'Yeast':
        pal = ['#1f77b4', '#bdf6fe', '#ff7f0e', '#fed8b1']
    else:
        pal = ['#009292', '#90e4c1', '#920000', '#d9544d']

    ax = sns.boxplot(x='Position', y='Score', hue='LCD Class', data=df, palette=pal, dodge=True, showfliers=False, hue_order=['G', 'Scrambled G', 'QN', 'Scrambled QN'])
    ax = sns.stripplot(x='Position', y='Score', hue='LCD Class', data=df, dodge=True, color='black', alpha=0.5, hue_order=['G', 'Scrambled G', 'QN', 'Scrambled QN'])

    if proteome == 'Yeast':
        legend_elements = [Patch(facecolor='#1f77b4', edgecolor='0.2', label='G'),
                            Patch(facecolor='#bdf6fe', edgecolor='0.2', label='G Scr'), #pale sky blue
                            Patch(facecolor='#ff7f0e', edgecolor='0.2', label='Q/N'),
                            Patch(facecolor='#fed8b1', edgecolor='0.2', label='Q/N Scr')] #light orange
    else:
        legend_elements = [Patch(facecolor='#009292', edgecolor='0.2', label='G'),
                            Patch(facecolor='#90e4c1', edgecolor='0.2', label='G Scr'), #light teal
                            Patch(facecolor='#920000', edgecolor='0.2', label='Q/N'),
                            Patch(facecolor='#d9544d', edgecolor='0.2', label='Q/N Scr')] #pale red
                        
    leg = plt.legend(handles=legend_elements, prop={'family':'Arial', 'size':14}, loc=2, handletextpad=0.2, labelspacing=0.1)
    leg.set_title('LCD Class', prop={'family':'Arial', 'size':16})

    plt.ylabel('Maximum Cluster Score', fontname='Arial', fontsize=16)
    title='Clustering of Degradation Promoting Residues within G-rich and Q/N-rich Domains'
    if proteome == 'Yeast':
        ax.set_xticklabels([proteome])
    else:
        ax.set_xticklabels([proteome + 's'])

    for label in ax.get_xticklabels():
        label.set_fontname('Arial')
        label.set_fontsize(16)
    for tick in ax.get_yticklabels():
        tick.set_fontname('Arial')
        tick.set_fontsize(16)
        
    plt.tight_layout()
    fig = plt.gcf()
    fig.set_size_inches(4, 6)
    
    # CUSTOMIZED Y-LIMS
    if cluster_size == 4 and proteome == 'Yeast':
        ax.set_ylim(top=4.9)
        plt.savefig('Fig S2A - Max ' + str(cluster_size) + 'aa Score G vs QN LCDs_Yeast.tiff', bbox_inches='tight', dpi=600)
    elif cluster_size == 8 and proteome == 'Yeast':
        ax.set_ylim(top=7.0)
        plt.yticks([x for x in range(-1, 6, 1)], fontname='Arial', fontsize=16)
        plt.savefig('Fig 2C - Max ' + str(cluster_size) + 'aa Score G vs QN LCDs_Yeast.tiff', bbox_inches='tight', dpi=600)
    elif cluster_size == 4 and proteome == 'Human':
        ax.set_ylim(top=5.7)
        plt.yticks([x for x in range(0, 6, 1)], fontname='Arial', fontsize=16)
        plt.savefig('Fig S2B - Max ' + str(cluster_size) + 'aa Score G vs QN LCDs_Human.tiff', bbox_inches='tight', dpi=600)
    elif cluster_size == 8 and proteome == 'Human':
        ax.set_ylim(top=9.2)
        plt.yticks([x for x in range(-1, 7, 1)], fontname='Arial', fontsize=16)
        plt.savefig('Fig 2D - Max ' + str(cluster_size) + 'aa Score G vs QN LCDs_Yeast.tiff', bbox_inches='tight', dpi=600)
    
    plt.close()
    
        
def run_stats(df, proteome, output, cluster_size):

    temp_df = pd.DataFrame.from_dict(df)
    
    for res in ['G', 'QN']:
        obs_vals = [temp_df.iloc[i]['Score'] for i in range(len(temp_df['Score'])) if temp_df.iloc[i]['LCD Class'] == res]
        scrambled_vals = [temp_df.iloc[i]['Score'] for i in range(len(temp_df['Score'])) if temp_df.iloc[i]['LCD Class'] == 'Scrambled ' + res]
        
        test_plotting_df = {'Position':[], 'Value':[], 'Category':[]}
        for i in range(len(obs_vals)):
            test_plotting_df['Position'].append( 0 )
            test_plotting_df['Value'].append( obs_vals[i] )
            test_plotting_df['Category'].append( 'Observed' )
            
        for i in range(len(scrambled_vals)):
            test_plotting_df['Position'].append( 0 )
            test_plotting_df['Value'].append( scrambled_vals[i] )
            test_plotting_df['Category'].append( 'Scrambled' )
        
        results = stats.ttest_rel(obs_vals, scrambled_vals) #Two-sided paired t-test with DEPENDENT samples (since each scrambled scores is derived from each original sequence)
        u_stat, pval = str(results).split(',')
        junk, u_stat = u_stat.split('=')
        junk, pval = pval[:-1].split('=')
        
        output.write(proteome + '\t' + res + '\tScrambled '+res + '\t' + str(cluster_size) + '\t' + u_stat + '\t' + pval + '\n')

    return output
    

if __name__ == '__main__':
    main()