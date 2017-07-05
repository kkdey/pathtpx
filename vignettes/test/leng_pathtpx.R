

#################  Leng et al pathtpx  ##########################

devtools::install_github("kkdey/singleCellRNASeqHumanLengESC", force=TRUE)

library(singleCellRNASeqHumanLengESC)
data("HumanLengESC")
leng_gene_names <- Biobase::featureNames(HumanLengESC);

leng_data <- t(Biobase::exprs(HumanLengESC));
leng_metadata <- Biobase::pData(HumanLengESC)
leng_cell_state <- leng_metadata$cell_state;

table(leng_cell_state)

index_1 <- which(leng_cell_state=="G1");
index_2 <- which(leng_cell_state=="S");
index_3 <- which(leng_cell_state=="G2");

g1_cells <- leng_data[index_1, ]
s_cells <- leng_data[index_2,]
g2_cells <- leng_data[index_3,]

c1 <- colMeans(g1_cells)
plot(c1)

c2 <- colMeans(s_cells)
plot(c2)

c3 <- colMeans(g2_cells)
plot(c3)

index1 <- which(c1 > 1 & c2 < 0.2 & c3 < 0.2)
index2 <- which(c2 > 1 & c1 < 0.2 & c3 < 0.2)
index3 <- which(c3 > 1 & c1 < 0.2 & c2 < 0.2)

#index21 <- setdiff(index1, c(index2, index3))
#index22 <- setdiff(index2, c(index1, index3))
#index23 <- setdiff(index3, c(index1, index2))

pathways <- list()
pathways[[1]] <- leng_gene_names[index1]
pathways[[2]] <- leng_gene_names[index2]
pathways[[3]] <- leng_gene_names[index3]


topic_clus <- pathtopics(leng_data, pathways = pathways, 
                         tol = 0.01, ord=FALSE, tmax = 100)

theta <- topic_clus$theta
omega <- topic_clus$omega

omega[which(leng_cell_state=="G1"),]
omega[which(leng_cell_state=="S"),]
omega[which(leng_cell_state=="G2"),]

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(leng_cell_state,
                        levels = c("G1", "S", "G2", "H1") ) )


rownames(omega) <- annotation$sample_id;


CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Cell cycle phase",
                            order_sample = TRUE,
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))


c1 <- colSums(g1_cells)
plot(c1)

c2 <- colSums(s_cells)
plot(c2)

c3 <- colSums(g2_cells)
plot(c3)

index1 <- which(c1 > 100 & c2 < 10 & c3 < 10)
index2 <- which(c2 > 100 & c1 < 10 & c3 < 10)
index3 <- which(c3 > 100 & c1 < 10 & c2 < 10)

topic_clus <- pathtopics(leng_data, pathways = pathways, 
                         tol = 0.01, ord=FALSE, tmax = 100)

theta <- topic_clus$theta
omega <- topic_clus$omega

omega[which(leng_cell_state=="G1"),]
omega[which(leng_cell_state=="S"),]
omega[which(leng_cell_state=="G2"),]

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(leng_cell_state,
                        levels = c("G1", "S", "G2", "H1") ) )


rownames(omega) <- annotation$sample_id;


CountClust::StructureGGplot(omega = omega,
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Cell cycle phase",
                            order_sample = TRUE,
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))



counts <- leng_data
pathways <- pathways
shape=NULL
initopics=NULL
tol=0.1
bf=FALSE
kill=2
verb=1
admix=TRUE
nbundles=1
ord = FALSE
tmax=10000


wtol=10^{-4}
qn=100
nonzero=FALSE
dcut=-10
top_genes=150
burn_in=5

pathways = pathways_indices
alpha=shape
