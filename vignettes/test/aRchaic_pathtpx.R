

###############  aRchaic pathtpx  ##########################

moderns_data <- get(load("../../../ancient-damage/data/moderns_lite/moderns_lite.rda"))
sherpa_data <- get(load("../../../ancient-damage/data/Sherpa_data/Sherpa_data.rda"))

pooled_data <- rbind(moderns_data, sherpa_data)

features <- colnames(pooled_data)
muts <- substring(features, 2, 5)
pos <- sapply(features, function(x) return(as.numeric(tail(strsplit(x, "_")[[1]],1))))

damage_index <- names(which(muts=="C->T" & pos < 10))

pathway_names <- list()
pathway_names[[1]] <- damage_index
pathway_names[[2]] <- setdiff(colnames(pooled_data), pathway_names[[1]])

topic_clus <- pathtopics(pooled_data, pathways = pathway_names, tol = 0.01, ord=FALSE)
theta <- topic_clus$theta
omega <- topic_clus$omega

plot(theta[,1])
plot(theta[,2])
sort(theta[,2], decreasing = TRUE)[1:10]
sort(theta[,1], decreasing = TRUE)[1:10]

topic_clus_2 <- maptpx::topics(pooled_data, K=2, tol = 0.1)
omega2 <- topic_clus_2$omega
theta2 <- topic_clus_2$theta

omega2
plot(theta2[,1])
plot(theta2[,2])
sort(theta2[,2], decreasing = TRUE)[1:10]
sort(theta2[,1], decreasing = TRUE)[1:10]
