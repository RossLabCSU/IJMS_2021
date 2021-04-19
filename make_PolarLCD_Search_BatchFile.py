
def main():

    output = open('RUN_LCD-Composer_PolarLCD_Search_Batch.bat', 'w')

    for aa in ['G', 'N', 'Q', 'S', 'T']:
        output.write('python LCD-Composer.py orf_trans.fasta Yeast_' + aa + '_40percent_LCD-Composer_RESULTS -a ' + aa + '\n')

    output.close()

if __name__ == '__main__':
    main()