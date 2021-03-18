#!/bin/bash


# This script will:
#   1. Take single .csv feature importance files that have the result for one seed
#   2. Combine them together to have the results for 100 seeds in one .csv file
#   3. We only keep the header of the first file.

# In the end, the combined_best file must be 101 lines. 1st line is the header and the 100 lines have the data of 100 files.
#             the combined_all file must have 100*(hyper-parameter number)+1 lines.
########################################################################################

SEARCH_DIR=results/runs
FINAL_DIR=results/

# Keep the first line of File1 and remove the first line of all the others and combine

for model in "glmnet" "rf" "svmRadial" "rpart2" "xgbTree"
do
  	head -1 $SEARCH_DIR/"$model"_100_feature-importance.csv  > $SEARCH_DIR/combined_feature-importance_"$model".csv; tail -n +2 -q $SEARCH_DIR/"$model"_*_feature-importance.csv >> $SEARCH_DIR/combined_feature-importance_"$model".csv
    mv $SEARCH_DIR/combined_feature-importance_"$model".csv $FINAL_DIR/combined_feature-importance_"$model".csv
done
