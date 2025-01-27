# ???????????????====================================
# ??????????????????
rm(list=ls())
# ??????R???
library(pheatmap)
# ??????????????????
setwd("C:/Users/linqingquan/Desktop/GSE_121787")

# ???????????????====================================
# ????????????
dataset <- read.table('diffgene_express.txt',header = TRUE, row.names = 1)
# ???????????????????????????????????????????????????
exp_ds = dataset[c(1:30),c(1:10)]
# ????????????????????????
cell_list=c(rep('cell_1',5),
            rep('cell_2',5))
annotation_c <- data.frame(cell_list)
rownames(annotation_c) <- colnames(exp_ds)

# ????????????=====================================
pheatmap(exp_ds, #????????????
         cluster_rows = T,#?????????
         cluster_cols = T,#?????????
         annotation_col =annotation_c, #??????????????????
         annotation_legend=TRUE, # ??????????????????
         show_rownames = T,# ????????????
         show_colnames = T,# ????????????
         scale = "row", #???????????????
         color =colorRampPalette(c("#8854d0", "#ffffff","#fa8231"))(100) # ??????????????????
)
