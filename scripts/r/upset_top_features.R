library(data.table)
library(svglite)
library(UpSetR)
library(ggplot2)

# Open top 5 feature df
feature_df <- fread(snakemake@input[['df']], header=TRUE, sep='\t')

# Select rows according to each model and then select only the feature column
log_reg <- subset(feature_df, model == 'log_reg')
feat_log_reg <- log_reg[['feature']]

svc <- subset(feature_df, model == 'svc')
feat_svc <- svc[['feature']]

gbm <- subset(feature_df, model == 'gbm')
feat_gbm <- gbm[['feature']]

knn <- subset(feature_df, model == 'knn')
feat_knn <- knn[['feature']]

# Create list input for UpSetR
listInput <- list(log_reg = feat_log_reg, svc = feat_svc, gbm = feat_gbm,
                  knn = feat_knn)

# Create the upset plot and order it by degree of intersection (intersection to all models, than 3 models, than 2 than 1)
fig <- upset(fromList(listInput), text.scale=c(2, 2, 2, 2, 2, 0), order.by='degree',
      sets.x.label="Number of top \npredictive features", point.size=6, line.size=1,
      matrix.color='black', main.bar.color='black', sets.bar.color='black')
svg(snakemake@output[['upset']], width=8, height=8)
print(fig)
dev.off()
