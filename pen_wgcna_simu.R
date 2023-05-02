#Simulation analysis
library(igraph)
library(effsize)
library(MASS)
library(WGCNA)
library(glmnet)
out_path <- "/Users/dokada/Desktop/work/lasso_simu_1228/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}

#Full code(NS)
#Create data X (sample by Gene)
n_genes <- 500
femData <- read.csv("/Users/dokada/Desktop/work/detarmin_data/FemaleLiver-Data/LiverFemale3600.csv")
datExpr0 = as.data.frame(t(femData[, -c(1:8)]));
names(datExpr0) = femData$substanceBXH;
rownames(datExpr0) = names(femData)[-c(1:8)];
nona_genes_idx <- complete.cases(t(datExpr0))
X_raw <- as.matrix(datExpr0[, nona_genes_idx])
vars <- apply(X_raw,2,var)
X_raw2 <-  X_raw[,order(vars, decreasing=T)[1:n_genes]]
covX <- cor(X_raw2)

#Replicate 1000 times and calculate the correlation coefficients and network concept
#noiseが大きすぎると係数が全てゼロになって、標準偏差がゼロですの警告
n_samples_can <- c(100, 500, 1000, 100, 500, 1000, 100, 500, 1000)
noise_can <- c(5, 5, 5, 10, 10, 10, 15, 15, 15)
n_type <- length(noise_can)
n_rep <- 100
res_mean_sd <- NULL
for(a in 1:n_type){

    #our_dir
    n_samples <- n_samples_can[a]
    noise <- noise_can[a]
    out_dir <- paste0("/Users/dokada/Desktop/work/lasso_simu_1228/n", n_samples, "_", noise, "/")
    if(file.exists(out_dir)){
        unlink(out_dir, recursive=TRUE)
        dir.create(out_dir)
    }else{
        dir.create(out_dir)
    }

    #analysis
    #options(warn=2) #for debug
    seed_can <- 1:n_rep
    labs <- c("lasso_et","lasso_ds", "lasso_cn", "lasso_cc", "lasso_ma", "lasso_gs", "lasso_wei", "lasso_tom",
    "elas_et", "elas_ds", "elas_cn", "elas_cc", "elas_ma", "elas_gs", "elas_wei", "elas_tom",
    "ridge_et","ridge_ds", "ridge_cn", "ridge_cc", "ridge_ma", "ridge_gs")
    res <- matrix(NA, n_rep, length(labs))
    colnames(res) <- labs
    for(s in 1:length(seed_can)){
        myseed <- seed_can[s]
        set.seed(myseed)

        #Generate X and WGCNA
        flg <- 0
        while(flg==0){
            mu <- rep(0, n_genes)
            X <- mvrnorm(n_samples, mu = mu, Sigma = covX)
            sft = pickSoftThreshold(X)
            #ft = suppressWarnings(pickSoftThreshold(X)) #for debug
            softPower <- sft$powerEstimate
            if(!is.na(softPower)) flg <- 1
        }
        adjacency = adjacency(as.data.frame(X), power=softPower)
        TOM = TOMsimilarity(adjacency)
        dissTOM = 1 - TOM
        tmp_res <- fundamentalNetworkConcepts(adjacency)
        Connectivity <- tmp_res[["Connectivity"]]
        ClusterCoef <- tmp_res[["ClusterCoef"]]
        MAR <- tmp_res[["MAR"]]


        #Generate phenotype
        n_tgs1 <- ncol(X)
        true_gs1 <- sample(colnames(X), n_tgs1, replace=FALSE)
        true_coef1 <- matrix(0, nrow=ncol(X), ncol=1)
        rownames(true_coef1) <- colnames(X)
        true_coef1[true_gs1,1] <- rnorm(n_tgs1)
        y_base1 <- X %*% true_coef1
        err <- matrix(rnorm(nrow(X), mean=0, sd=noise), nrow=nrow(X), ncol=1)
        y1 <- y_base1 + err

        #GS
        GS <- rep(NA, ncol(X))
        for(i in 1:ncol(X)){
            GS[i] <- abs(cor(X[,i], y1))
        }

        #Lasso
        alpha <- 1
        m <- cv.glmnet(x = X, y = y1, family = "gaussian", alpha = alpha)
        best.lambda <- m$lambda.min
        en.model <- glmnet(x = X, y = y1, family = "gaussian", lambda = best.lambda, alpha = alpha)
        est_beta <- en.model$beta
        est_y <- X %*% est_beta
        abs_eb <- abs(est_beta)

        #Example Visualization
        if(s==1){
            #Connectivity histograms
            png(paste0(out_dir, "ex_connecitivy_hist.png"), width=960, height=960)
            par(mar = c(9, 9, 6, 4)) ##bottom, left, top, right
            xlab = "Connecitivy"
            ylab = "Frequency"
            main = "Histogram of Connecitivy"
            hist(Connectivity, cex.axis=4, cex.lab=4, cex.main=4,xlab="", ylab="", main="")
            mtext(xlab, side=1, line=6, cex=4)
            mtext(ylab, side=2, line=6, cex=4)
            mtext(main, side=3, line=3, cex=4)
            dev.off()

            #Visualizayion in igraph
            geneTree = hclust(as.dist(dissTOM), method = "average")
            dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
            dynamicColors = labels2colors(dynamicMods)
            adj_viz <- adjacency
            cutoff <- 0.05
            adj_viz[adj_viz > cutoff] <- 1
            adj_viz[adj_viz <= cutoff] <- 0
            diag(adj_viz) <- 0
            adj_viz_sub2 <- adj_viz
            graph <- graph.adjacency(adj_viz_sub2, mode = "undirected", diag = FALSE)
            #V(graph)$color <- dynamicColors
            png(paste0(out_dir, "ex_network.png"), width=960, height=960)
            par(mar = c(9, 9, 6, 4)) ##bottom, left, top, right
            plot(graph,  vertex.size=3, vertex.label=NA, edge.width=6)
            mtext("Co-expression network", side=3, line=3, cex=4)
            dev.off()

            #Estimation accuracy
            png(paste0(out_dir, "ex_phenotype.png"), width=960, height=960)
            par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
            xlab = "Estimated phenotype"
            ylab = "Observed phenotype"
            main = ""
            plot(est_y, y1, xlab="", ylab="", cex.axis=4, cex=2, cex.lab=3, pch=16)
            mtext(xlab, side=1, line=6, cex=4)
            mtext(ylab, side=2, line=6, cex=4)
            mtext(main, side=3, line=3, cex=4, adj=0)
            dev.off()

            #Selected and Non-selected
            png(paste0(out_dir, "coef_gs_lasso.png"), width=960, height=960)
            par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
            xlab = "Abs(True Beta)"
            ylab = "GS"
            #corval <- cor(est_beta, true_coef1, method="spearman")
            main <- "LASSO"
            cols <- ifelse(est_beta!=0, "red", "black")
            plot(abs(true_coef1), GS, xlab="", ylab="", cex.axis=4, cex=2, cex.lab=3, col=cols, pch=16)
            mtext(xlab, side=1, line=6, cex=4)
            mtext(ylab, side=2, line=6, cex=4)
            mtext(main, side=3, line=3, cex=4, adj=0)
            dev.off()
        }

        #VS true beta
        corval <- cor(est_beta, true_coef1, method="spearman")
        res[s, "lasso_et"] <- cor(est_beta, true_coef1, method="spearman")
        res[s, "lasso_ds"] <- sum(sign(est_beta * true_coef1)==-1) / sum(est_beta != 0)

        #Node statistics
        mt <- "spearman"
        res[s, "lasso_cn"] <- cor(Connectivity, abs_eb, method=mt)
        res[s, "lasso_cc"] <- cor(ClusterCoef, abs_eb, method=mt)
        res[s, "lasso_ma"] <- cor(MAR, abs_eb, method=mt)
        res[s, "lasso_gs"] <- cor(GS, abs_eb, method=mt)

        #Edge statistics
        nonsel_idx <-  which(est_beta == 0)
        sel_idx <- which(est_beta != 0)

        #Weight (elected and non-selected subnetwork)            
        sel_adj <- adjacency[sel_idx, sel_idx]
        nonsel_adj <- adjacency[nonsel_idx, nonsel_idx]
        ssv <- sel_adj[upper.tri(sel_adj)]
        nsv <- nonsel_adj[upper.tri(nonsel_adj)]
        cohend <- cohen.d(ssv, nsv) #ssv - nsv
        res[s, "lasso_wei"] <- cohend$estimate

        #Violone plot
        sel_tom <- TOM[sel_idx, sel_idx]
        nonsel_tom <- TOM[nonsel_idx, nonsel_idx]
        ssv <- sel_tom[upper.tri(sel_tom)]
        nsv <- nonsel_tom[upper.tri(nonsel_tom)]
        cohend <- cohen.d(ssv, nsv)
        res[s, "lasso_tom"] <- cohend$estimate

        #Ridge Regression
        alpha <- 0
        m <- cv.glmnet(x = X, y = y1, family = "gaussian", alpha = alpha)
        best.lambda <- m$lambda.min
        en.model <- glmnet(x = X, y = y1, family = "gaussian", lambda = best.lambda, alpha = alpha)
        est_beta <- en.model$beta
        est_y <- X %*% est_beta
        abs_eb <- abs(est_beta)

        #VS true beta
        corval <- cor(est_beta, true_coef1, method="spearman")
        res[s, "ridge_et"] <- cor(est_beta, true_coef1, method="spearman")
        res[s, "ridge_ds"] <- sum(sign(est_beta * true_coef1) == -1) / sum(est_beta != 0)

        #Node statistics
        mt <- "spearman"
        res[s, "ridge_cn"] <- cor(Connectivity, abs_eb, method=mt)
        res[s, "ridge_cc"] <- cor(ClusterCoef, abs_eb, method=mt)
        res[s, "ridge_ma"] <- cor(MAR, abs_eb, method=mt)
        res[s, "ridge_gs"] <- cor(GS, abs_eb, method=mt)
        cat(s, "\n")

        #ElasticNet
        alpha <- 0.5
        m <- cv.glmnet(x = X, y = y1, family = "gaussian", alpha = alpha)
        best.lambda <- m$lambda.min
        en.model <- glmnet(x = X, y = y1, family = "gaussian", lambda = best.lambda, alpha = alpha)
        est_beta <- en.model$beta
        est_y <- X %*% est_beta
        abs_eb <- abs(est_beta)

        #VS true beta
        corval <- cor(est_beta, true_coef1, method="spearman")
        res[s, "elas_et"] <- cor(est_beta, true_coef1, method="spearman")
        res[s, "elas_ds"] <- sum(sign(est_beta * true_coef1)==-1) / sum(est_beta != 0)

        #Node statistics
        mt <- "spearman"
        res[s, "elas_cn"] <- cor(Connectivity, abs_eb, method=mt)
        res[s, "elas_cc"] <- cor(ClusterCoef, abs_eb, method=mt)
        res[s, "elas_ma"] <- cor(MAR, abs_eb, method=mt)
        res[s, "elas_gs"] <- cor(GS, abs_eb, method=mt)

        #Edge statistics
        nonsel_idx <-  which(est_beta == 0)
        sel_idx <- which(est_beta != 0)

        #Weight (elected and non-selected subnetwork)            
        sel_adj <- adjacency[sel_idx, sel_idx]
        nonsel_adj <- adjacency[nonsel_idx, nonsel_idx]
        ssv <- sel_adj[upper.tri(sel_adj)]
        nsv <- nonsel_adj[upper.tri(nonsel_adj)]
        cohend <- cohen.d(ssv, nsv)
        res[s, "elas_wei"] <- cohend$estimate

        #TOM
        sel_tom <- TOM[sel_idx, sel_idx]
        nonsel_tom <- TOM[nonsel_idx, nonsel_idx]
        ssv <- sel_tom[upper.tri(sel_tom)]
        nsv <- nonsel_tom[upper.tri(nonsel_tom)]
        cohend <- cohen.d(ssv, nsv)
        res[s, "elas_tom"] <- cohend$estimate

        #Selected vs Non-selected
        if(s==1){
            png(paste0(out_dir, "coef_gs_elastic.png"), width=960, height=960)
            par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
            xlab = "Abs(True Beta)"
            ylab = "GS"
            #corval <- cor(est_beta, true_coef1, method="spearman")
            main <- "ElasticNet"
            cols <- ifelse(est_beta!=0, "red", "black")
            plot(abs(true_coef1), GS, xlab="", ylab="", cex.axis=4, cex=2, cex.lab=3, col=cols, pch=16)
            mtext(xlab, side=1, line=6, cex=4)
            mtext(ylab, side=2, line=6, cex=4)
            mtext(main, side=3, line=3, cex=4, adj=0)
            dev.off()
        }

    }


    #Visualization
    l1 <- c("lasso_et","elas_et","ridge_et")
    l2 <- c("lasso_gs","elas_gs","ridge_gs")
    l3 <- c("lasso_cn","elas_cn","ridge_cn")
    l4 <- c("lasso_cc","elas_cc","ridge_cc")
    l5 <- c("lasso_ma","elas_ma","ridge_ma")
    l6 <- c("lasso_ds","elas_ds")
    l7 <- c("lasso_wei","elas_wei")
    l8 <- c("lasso_tom","elas_tom")
    labs_list <- list(l1, l2, l3, l4, l5, l6, l7, l8)
    labs_name <- c("Beta", "GS", "Connectivity", "ClusterCoef", "MAR", "Different Sign", "Weight", "TOV")
    ylabs <- c("Correlation", "Correlation", "Correlation", "Correlation", "Correlation", "Percentage", "Cohen's D", "Cohen's D")
    for(i in 1:length(labs_list)){
        png(paste0(out_dir, labs_name[i], ".png"), width=960, height=960)
        par(mar = c(9, 12, 9, 4)) ##bottom, left, top, right
        lab_sub <- labs_list[[i]]
        tab <- res[,lab_sub]
        if(ncol(tab)==3) colnames(tab) <- c("Lasso", "ElasticNet", "Ridge")
        if(ncol(tab)==2) colnames(tab) <- c("Lasso", "ElasticNet")
        boxplot(tab, cex.axis=4, cex.lab=4, main=labs_name[i], cex.main=4)
        mtext(ylabs[i], side=2, line=6, cex=4)
        dev.off()
    }

    #record
    write.csv(res, file=paste0(out_dir, "res.csv"))
    mus <- round(colMeans(res), 2)
    sds <- round(apply(res, 2, sd),2)
    ms_sd <- paste0(mus, "(", sds, ")")
    res_mean_sd <- cbind(res_mean_sd, ms_sd) 
}

#Record
res_mean_sd <- rbind(n_samples_can, noise_can, res_mean_sd)
rownames(res_mean_sd) <- c("n_samples", "noise", labs)
colnames(res_mean_sd) <- paste0("Condition", 1:ncol(res_mean_sd))
write.csv(res_mean_sd, file=paste0(out_path, "res_mean_sd.csv"))
