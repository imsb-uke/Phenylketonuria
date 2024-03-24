Under development

# Phenylketonuria landscape analysis

Pipeline:
1. Apply Gaussian model and extract features
2. Perform QC (rmse - var)
3. Detect non-responders
4. Subpopulation detection

## Run pipeline

install packages
```
numpy
pandas
scipy
matplotlib
```
To run the pipeline use the script ```gm.py```, which requres experimant files as ```*.xlsm```.
Example 1) run the pipeline for one experimant
```
python gm.py --input_dir Data/experiments_all/20240305_EXP13.xlsm
```
Example 2) run the pipeline for all experimants in a directory (*e.g.* ```Data/experiments_all```)
```
python gm.py --input_dir Data/experiments_all
```
