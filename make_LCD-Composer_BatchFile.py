
def main():

    output = open('RUN_LCD-Composer_Initial_G-QN_Search_Batch.bat', 'w')
    organisms = ['Yeast', 'Human']
    i = 0
    for file in ['orf_trans.fasta', 'uniprot-proteome_UP000005640_Human_ALL_ISOFORMS.fasta']:
        percents = ['35', '50']
        j = 0
        for aa_group in ['G', 'QN']:
            output.write('python LCD-Composer.py ' + file + ' ' + organisms[i] + '_' + aa_group + '_' + percents[j] + 'percent_LCD-Composer_RESULTS -a ' + aa_group + ' -c ' + percents[j] + '\n')
            j += 1
        i += 1

    output.close()

if __name__ == '__main__':
    main()