# Reproducing Computational Analyses in Cascarina *et al.* (2021) *IJMS*

### Instructions
1. Download LCD-Composer.py from https://github.com/RossLabCSU/LCD-Composer, as well as all files in the IJMS_2021 directory and place in the same folder.
2. Extract files from compressed folders in the same location as LCD-Composer.py
3. Download mPAPA.py from https://github.com/RossLabCSU/mPAPA and place in the same folder as the files from Steps 1 & 2.
4. Navigate to appropriate folder via command line.
5. Run the following commands in-sequence (NOTE: each run must be completed before issuing the next command):

```
python make_Initial_G-QN_Search_BatchFile.py
```

```
.\RUN_LCD-Composer_Initial_G-QN_Search_Batch.bat
```

```
python compare_G_vs_QN_halflives_Yeast_and_Humans.py
```

```
python plot_Grich_vs_QNrich_domain_DegradationScores.py
```

```
python plot_G_and_QN_Scrambling_Results.py
```

```
.\RUN_plot_ModelProtein_Compositions.bat
```

```
python make_PolarLCD_Search_BatchFile.py
```

```
.\RUN_LCD-Composer_PolarLCD_Search_Batch.bat
```

```
python plot_PolarLCD_domain_DegradationScores.py
```

```
python mPAPA.py Grich_and_QNrich_LCD-Sup35_Fusion_Sequences.txt -o Grich_and_QNrich_LCD-Sup35_Fusion_Sequences_PAPA_RESULTS.tsv --ignore_fold_index --verbose
```

```
python plot_PAPA_results.py
```

This series of commands generates Fig 1, Fig 2, Fig 3A, Fig 6A, Fig S2, Fig S4, Table S1, and all data appearing in Table S2.
