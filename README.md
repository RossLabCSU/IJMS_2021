# Reproducing Analyses Pertaining to Development of the Linear Dispersion Parameter

### Instructions
1. Download all files in the JBC_2021 directory.
2. Navigate to appropriate folder via command line.
3. Run the following commands in-sequence to generate Fig S1 and to validate determination of the minimum and maximum linear dispersion for benchmark sequences from 5aa to 30aa in length (NOTE: each run must be completed before issuing the next command):

```
python make_LCD-Composer_BatchFile.py
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

This series of commands generates Fig 1, Fig 2A-B, Table S1, and all data appearing in Table S2.
