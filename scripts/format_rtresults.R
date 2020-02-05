# Script to format regression tree results to make gridplots of variable selection
# and variable importance in the tree construction.

# load regression tree models
tree <- readRDS("outputs/regressiontrees.RDS") 
#NOTE: not available due to stored bird abundance data.

# pull out variables selected for tree construction
vars = lapply(tree, function(x) levels(droplevels(x$frame$var, "<leaf>")))
# make dataframe with binary code for 1 = selected, 0 = excluded
df = data.frame(matrix(0, nrow = 55, ncol = 17), row.names = colnames(bird))
colnames(df) = c(colnames(land), "xcoord", "ycoord")
for(i in 1:length(vars)){df[i, vars[[i]]] = 1}
# write file
write.csv(df, "outputs/rt_variables.csv")

# pull out variable importance values
imp = lapply(tree, function(x) x$variable.importance)
# make dataframe with variable importance values
df2 = data.frame(matrix(NA, nrow = 55, ncol = 17),
                 row.names = colnames(bird))
colnames(df2) = colnames(df)
for(i in 1:length(imp)){
  vars = names(imp[[i]])
  for(n in 1:length(imp[[i]])){ df2[i, vars[n]] = unname(imp[[i]][n]) }
}
# write file
write.csv(df2, "outputs/rt_importance.csv")
