
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

combined_deg_propensities = { 'A':0.2665627892708544, 'R':-0.7412501624697525, 'N':-1.1869752279573649, 'D':-0.728512694019163,
                                'C':0.4206842379683707, 'E':-0.7157286901525158, 'Q':-0.7233002187306328, 'G':-0.4440764148760815,
                                'H':-0.028047412555156075, 'I':0.8799390174848035, 'L':1.046555755328647, 'K':-0.8574121141055006,
                                'M':1.038969401500499, 'F':0.22392755762641842, 'P':-0.10366937834956014, 'S':-0.16794980890760758,
                                'T':-0.03867551503137706, 'W':0.12078754467668205, 'Y':0.25418799347663784, 'V':0.934778908124434 }
                
def main():

    df = {'Score':[], 'Category':[], 'Position':[]}
    leg_labels = []
    comp = '40'
    
    for res in ['G', 'N', 'Q', 'S', 'T']:
        leg_labels.append( res )

        h = open('Yeast_' + res + '_' + comp + 'percent_LCD-Composer_RESULTS.tsv')
        for i in range(7):
            h.readline()
        header = h.readline()
        for line in h:
            items = line.rstrip().split('\t')
            id, *junk = items[0].split(' ')
            seq = items[1]
            if len(seq) < 40:
                continue
            score = sum( [combined_deg_propensities[aa] for aa in seq if aa not in res] )
            df['Score'].append(score)
            df['Category'].append(res)
            df['Position'].append(0)
            
        h.close()

    boxplot(df, leg_labels)
    
        
def boxplot(df, leg_labels):

    sns.boxplot(x='Position', y='Score', data=df, hue='Category', showfliers=False)
    
    sns.stripplot(x='Position', y='Score', data=df, hue='Category', dodge=True, color='0.2', alpha=0.5)
    plt.xticks([], labels=[])
    plt.yticks(fontname='Arial', fontsize=12)
    plt.ylabel('Degradation Scores\n(excluding AA of interest)', fontname='Arial', fontsize=12)
    
    pal = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    legend_elements = [Patch(facecolor=pal[i], edgecolor='0.2', label=leg_labels[i]) for i in range(len(leg_labels))]
    leg = plt.legend(handles=legend_elements, loc=2, bbox_to_anchor=(1.01, 1.027), prop={'family':'Arial', 'size':12})
    leg.set_title('LCD Class', prop={'family':'Arial', 'size':12})
    
    fig = plt.gcf()
    fig.set_size_inches(4.5, 4)
    plt.savefig('Fig 6A - PolarLCD DegradationScores.tiff', bbox_inches='tight', dpi=600)
    plt.close()

    
if __name__ == '__main__':
    main()
    