
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def main():

    df = {'Protein':[],
            'Score':[],
            'Position':[]}
    
    h = open('Grich_and_QNrich_LCD-Sup35_Fusion_Sequences_PAPA_RESULTS.tsv')
    header = h.readline()
    for line in h:
        prot, high_score, position, papa_scores, fold_index = line.rstrip().split('\t')
        papa_scores = papa_scores[1:-1].split(' ')
        papa_scores = [float(x) for x in papa_scores]
        i = 0
        for score in papa_scores:
            df['Protein'].append(prot)
            df['Score'].append(score)
            df['Position'].append( i )
            i += 1
            
    plotting(df)
        
        
def plotting(df):

    prots = ['Npl3', 'Ynl208W_Grich', 'Ynl208W_QNrich', 'Rnq1', 'Ure2']
    colors = ['0.85', '0.65', '0.45', '0.0']
    
    df = pd.DataFrame.from_dict(df)
    
    for prot in prots:
        color_index = 0
        min_val = 10    # PLACEHOLDER INITIALIZATION VALUE
        for variant in ['_WT', '_+2hyd', '_+4hyd', '_+6hyd']:
            temp_df = df[df['Protein'] == prot+variant]
            sns.lineplot(x='Position', y='Score', data=temp_df, color=colors[color_index])
            color_index += 1
            if min(temp_df['Score']) < min_val:
                min_val = min(temp_df['Score'])
        
        max_val = max(temp_df['Score'])
            
        plt.xlim(0, 56)
        plt.ylim(min_val-0.02, max_val+0.02)
        plt.xlabel('Position', fontname='Arial', fontsize=12)
        plt.ylabel('Prion Propensity', fontname='Arial', fontsize=12)
        plt.xticks(fontname='Arial', fontsize=10)
        plt.yticks(fontname='Arial', fontsize=10)
        title, *remainder = prot.split('_')
        plt.title(title, fontname='Arial', fontsize=14)
        
        fig = plt.gcf()
        fig.set_size_inches(6, 4)
        
        plt.savefig(prot + '_mPAPA_scores.tiff', bbox_inches='tight', dpi=600)
        plt.close()
    

if __name__ == '__main__':
    main()