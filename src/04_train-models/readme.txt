#partition.r:
  #perform unique partition on genes across 12 traits before model fitting.
  #This avoids using the same gene multiple times for different traits thus remove overfitting bias.

#part_func.r:
#function called by  partition.r.

#Gen_partition_res.RData:
#results of partition.r

#trainModel.r:
#train xgboost model using 11 traits and predict to 1 trait left.