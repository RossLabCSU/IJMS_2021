
"""
Description:
    CompositionPlotter is a composition-based method for identifying and visualizing low-complexity domains in protein sequences.
    Refer to Cascarina et. al. (2021) (INCLUDE LINK WHEN PUBLISHED) for a complete 
    description of algorithm testing and application.
    Documentation for ComopsitionPlotter is available at XXXX (INCLUDE LINK WHEN FINALIZED).
============================================================================================
License_info:
    CompositionPlotter is subject to the terms of the GPLv3.0 license. For a complete description of 
    license terms, please see the license at XXXX (INCLUDE LINK WHEN FINALIZED).
"""

__author__ = 'Sean M Cascarina'
__copyright__ = 'Copyright 2021'
__credits__ = ['Sean M Cascarina']
__license__ = 'GPLv3.0'
__version__ = '1.0'
__maintainer__ = 'Sean M Cascarina'
__email__ = 'Sean.Cascarina@colostate.edu'


all_aas = 'ACDEFGHIKLMNPQRSTVWY'

def main(args):

    win_size, amino_acids, comp_threshold, resolution, output_type, figure_length, figure_height, show_threshold, color_palette = get_params(args)
    print('Start time:', str(datetime.datetime.now()), 'Amino Acids:', amino_acids)
    
    file_info = args.fasta_file.split('.')
    filename = file_info[-2].replace('\\', '')
    h = open(args.fasta_file)

    for id, seq in fasta_parser(h):

        # WORKS WELL FOR UNIPROT FORMAT FASTA HEADERS, WHICH HAVE '|' THAT WILL THROW AN ERROR WHEN USED IN THE OUTPUT FILE NAME
        id, *remainder = id.replace('|', '-').split(' ')
        
        # PREP DICTIONARY WITH NECESSARY AMINO ACIDS
        df = {}
        if amino_acids == 'Auto-detect':
            for aa in all_aas:
                df[aa] = []
        else:
            for group in amino_acids:
                df[group] = []
        
        seq = check_stopcodons(id, seq)
        df = calc_composition(seq, df, amino_acids, win_size)
        plotting_aas = set_plotting_AAs(df, amino_acids, comp_threshold)
        
        plot_comp_disp(seq, id, df, plotting_aas, win_size, resolution, output_type, figure_length, figure_height, amino_acids, comp_threshold, color_palette)
        
        
def calc_composition(seq, df, amino_acids, win_size):

    for i in range(len(seq) - win_size+1):
        window = seq[i:i+win_size]
        
        #SKIP SEQUENCES SHORTER THAN THE WINDOW SIZE
        if len(window) < win_size:
            continue

        if amino_acids == 'Auto-detect':
            for aa in all_aas:
                df[aa].append( window.count(aa) / len(window) * 100 )
        else:
            for group in amino_acids:
                df[group].append( sum([window.count(aa) / len(window) * 100 for aa in group]) )
                
    return df
    
        
def check_stopcodons(id, seq):
    
    #REMOVE STOP CODON FROM C-TERMINUS AND WARN USERS IF A SEQUENCE CONTAINS MULTIPLE STOP CODONS
    if seq[-1] == '*':
        seq = seq[:-1]
    if '*' in seq:
        print('Runtime warning: protID ' + id + ' contains multiple stop codons, which will slightly affect the amino acid composition and separation calculations. Stop codons are automatically removed by the program from the C-terminus of each sequence but internal stop codons are not removed. Consider removing extra stop codons before evaluating, removing these sequences from analyses entirely, or evaluating sequences as-is.')
            
    return seq
    
        
def set_plotting_AAs(df, amino_acids, threshold):
    
    if amino_acids == 'Auto-detect':
        new_aas = ''
        for aa in all_aas:
            if max(df[aa]) >= threshold:
                new_aas += aa
        amino_acids = new_aas
        
    return amino_acids


def plot_comp_disp(seq, id, df, plotting_aas, win_size, resolution, output_type, figure_length, figure_height, amino_acids, comp_threshold, color_palette):

    if color_palette == 'Seaborn':
        pal = sns.color_palette('colorblind')
        pal += sns.color_palette('pastel')
    else:
        pal = color_palette
        
    linestyles = ['-']*10 + ['--']*10

    # PLOT DOTTED LINE AT COMPOSITION THRESHOLD
    if amino_acids == 'Auto-detect' or args.show_threshold:
        plt.plot((-1000, len(df[plotting_aas[0]])+1000), (comp_threshold, comp_threshold), linestyle='--', color='0.8')
        
    index = 0
    for aa in plotting_aas:
        plt.plot([i+(win_size/2)+1 for i in range(len(df[aa]))], [x for x in df[aa]], color=pal[index], linestyle=linestyles[index])
        index += 1
    
    xmargin = len(df[aa]) * 0.015
    plt.xlim(-xmargin, len(df[aa])+win_size+xmargin)
    plt.yticks([x for x in range(0, 120, 20)], labels=[0,20, 40, 60, 80, 100], fontname='Arial', fontsize=16)
    plt.xticks(fontname='Arial', fontsize=16)
    plt.ylabel('AA Composition', fontname='Arial', fontsize=18)
    plt.xlabel('Protein Position', fontname='Arial', fontsize=18)
    
    legend_elements = [Line2D([0], [0], color=pal[i], lw=4, label=plotting_aas[i], linestyle=linestyles[i]) for i in range(len(plotting_aas))]
    
    leg = plt.legend(handles=legend_elements, prop={'size':12, 'family':'Arial'}, loc=2, bbox_to_anchor=(1,1))
    leg.set_title('Amino Acid(s)', prop = {'size':12, 'family':'Arial'})

    fig = plt.gcf()
    
    fig.set_size_inches(figure_length, figure_height)

    plt.savefig(id + '_' + '-'.join(plotting_aas) + '_CompositionPlot.' + output_type, bbox_inches='tight', dpi=resolution)

    plt.close()
    
    
def fasta_parser(file):
    """Parses each instance in a FASTA formatted file into a gene id and a sequence.
    
    Yields:
        id, seq (both strings)
    """
    #INITIATES GENE AND SEQ FOR FIRST ITERATION
    gene, seq = '', []
    
    #LOOPING THROUGH FILE USING GENERATOR
    for line in file:
        line = line.rstrip()
        if len(line) == 0:
            continue
            
        #YIELDS GENE AND SEQ IF THE NEXT ID IS REACHED
        if line.startswith('>'):
            if gene != '':
                yield (gene[1:], ''.join(seq))
            gene, seq = line, []
        else:
            seq.append(line)
            
    #YIELDS FINAL INSTANCE IN FASTA FILE
    if gene != '':
        yield (gene[1:], ''.join(seq))

    
def get_params( args ):
    """Gather and define user parameters. Also generates an error message if an invalid reduced alphabet is passed in by the user.
    
    Returns:
        1) Window size (int)
        2) Amino acids of interest (str)
        3) Composition threshold (float)
        4) Resolution (int)
        5) Output file type (str)
        6) Figure length (float)
        7) Figure Height (float)
        8) Show Threshold (True/False)
    """
    
    comp_threshold = float(args.composition)-0.0000001
    
    #RUN GATEKEEPER CHECKS=========================
    if comp_threshold < 0-0.000001 or comp_threshold > 100+0.000001:
        print('\n Invalid composition threshold. The composition threshold must be a number between 0-100 (inclusive)\n')
        exit()
        
    if args.resolution not in range(150, 1201):
        print('\n Invalid resolution input. Resolution value must be an integer between 150-1200 dots per inch (dpi). 600dpi is the default and is recommended for publication quality.\n')
        exit()
        
    if args.output_type.lower() not in ['tiff', 'tif', 'jpg', 'png']:
        print('\n Invalid output file type. File type must be tiff, tif, jpg, or png.\n')
        exit()
        
    #GET USER-SPECIFIED PARAMETERS=================
    win_size = args.window_size
    if args.amino_acids.lower() != 'auto-detect':
        amino_acids = args.amino_acids.split('_')
        amino_acids = [x.upper() for x in amino_acids]
    else:
        amino_acids = args.amino_acids
        
    if args.color_palette != 'Seaborn':
        color_palette = args.color_palette.split('_')
        color_palette = ['#'+hex for hex in color_palette]
        matches = [re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', hex) for hex in color_palette]
        if False in matches:
            print('\n Invalid color palette. The color palette must be a series of valid color hex codes, each separated by an underscore ("_"), and must be of the same length as the number of specified amino acid groups.')
            exit()
    else:
        color_palette = args.color_palette

    return win_size, amino_acids, comp_threshold, int(args.resolution), args.output_type.lower(), args.figure_length, args.figure_height, args.show_threshold, color_palette


def get_args(arguments):
    parser = argparse.ArgumentParser(description='Plot the amino acid composition of a protein.', prog='CompositionPlotter')
    
    parser.add_argument('fasta_file', help="""Your sequence file (in FASTA format).""")

    parser.add_argument('-w', '--window_size', type=int, default=20, 
                        help="""Scanning window size used in the initial search. Default=20aa.
                        Only integers (whole numbers) between 5-10000 are valid window sizes""")
                        
    parser.add_argument('-c', '--composition', type=str, default='40',
                        help="""Composition threshold for defining X-rich regions, where X represents a particular amino acid or group of amino acids. Default = 40 (corresponding to 40% composition).
                        
                        Composition is calculated as the fraction of residues of interest in each window.
                        
                        Value must be between 0-100.
                        """)

    parser.add_argument('-a', '--amino_acids', type=str, default='Auto-detect',
                        help="""Amino acid(s) of interest.
                        
                        Simple amino acid criteria should be an unbroken string consisting of the single-letter abbreviation for a single amino acid or a group of amino acids.
                        e.g.
                        Q       (Q composition)
                        QN      (Combined composition of Q+N)
                        QNST    (Combined composition of Q+N+S+T)
                        
                        Complex amino acid criteria allow for specification of different individual amino acids or groups of amino acids (combined composition) to plot.
                        Separate amino acids or groups of amino acids should be separated by an underscore '_'.
                        Complex criteria use "AND" logic, meaning each amino acid or group separated by underscores will be considered separately.
                        e.g.
                        Q_P     (calculate and plot Q composition and P composition independently)
                        QN_ST   (calculated combined Q+N composition and combined S+T composition independently, and plot the two groups)
                        G_R_Y   (calculate and plot G composition, R composition, and Y composition independently)
                        DEKR_FWY_P    (calculate combined D+E+K+R composition, combined F+W+Y composition, and individual P composition and plot the three groups separately)
                        """)

    parser.add_argument('-r', '--resolution', type=int, default=600, 
                        help="""Resolution in dots per inch (dpi) for the resulting composition plot image.
                        Default is 600dpi, which generally works well for generating publication-quality images if the .tiff file option is used.
                        Allowable range is 150-1200dpi.
                        NOTE: lower dpi and .png file types will result in smaller file sizes but lower quality images (useful for large-scale exploratory analyses)""")

    parser.add_argument('-o', '--output_type', type=str, default='tiff', 
                        help="""Output file type.
                        Default is tiff file.
                        Additional options are:
                           1) jpg
                           2) png
                           3) tif""")
                           
    parser.add_argument('-l', '--figure_length', type=float, default=6, 
                        help="""Desired length of the plot (in inches).
                        Default is 6in.
                        DISCLAIMER: very large values with high resolution settings may result in failure to run or extremely large file sizes.""")
                        
    parser.add_argument('-e', '--figure_height', type=float, default=4, 
                        help="""Desired height of the plot (in inches).
                        Default is 4in.
                        DISCLAIMER: very large values with high resolution settings may result in failure to run or extremely large file sizes.""")

    parser.add_argument('-x', '--show_threshold', action='store_true', 
                        help="""Flag used without arguments to show a composition line at the specified composition threshold.""")

    parser.add_argument('-p', '--color_palette', default='Seaborn', 
                        help="""Color palette used for line colors.
                        Default color palette is the standard palette used in the Seaborn plotting package.
                        
                        You can specify your own color palette using color hex codes WITHOUT the hashtag (#), with each color separated by an underscore.
                        If you specify your own palette, the number of colors must match the number of amino acid groups that you've specified.
                        If you are using the "Auto-detect" mode (i.e. did not specify amino acids to plot) and you wish to specify a custom color palette, it is recommended that you choose a large palette so that all detected amino acids will have a unique color. Otherwise you may receive an error.
                        """)

                        
    args = parser.parse_args(arguments)
    
    return args

if __name__ == '__main__':
    import sys, argparse, datetime, re
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import seaborn as sns
    args = get_args(sys.argv[1:])
    main(args)