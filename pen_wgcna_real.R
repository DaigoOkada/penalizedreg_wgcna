##27k inplementation(Real data)
out_path <- "/Users/dokada/Desktop/work/lasso_real_1228/"
if(file.exists(out_path)){
    unlink(out_path, recursive=TRUE)
    dir.create(out_path)
}else{
    dir.create(out_path)
}

#analysis
library(data.table)
library(effsize)
library(MASS)
library(WGCNA)
library(glmnet)
dat <- read.csv("/Users/dokada/Desktop/work/detarmin_data/GSE41037_series_matrix_dat.csv", header=T, row.names=1)
cpgs <- rownames(dat)
annot <- read.csv("/Users/dokada/Desktop/work/detarmin_data/GSE41037_series_matrix_annot.csv", header=T) #all(colnames(x)==colnames(annot)), OK
ages <- sapply(annot[2,],function(x){as.numeric(strsplit(x,":")[[1]][2])})
sex <- sapply(annot[1,],function(x){gsub(" ","",strsplit(x,":")[[1]][2])})

#Additional File S3を読み込む
tab <- read.csv("/Users/dokada/Desktop/work/detarmin_data/13059_2013_3156_MOESM3_ESM.csv",header=T,stringsAsFactors=F)
horvath_cpgs <- tab[4:nrow(tab),1]

#Epigenetic Ageを計算
n_horvath <- length(horvath_cpgs)
int_mat <- NULL
for(i in 1:n_horvath){
    tmp_cpg_idx <- which(cpgs==horvath_cpgs[i])
    intensity <- as.numeric(dat[tmp_cpg_idx,])
    int_mat <- cbind(int_mat,intensity)
}
horvath_cpgs <- tab[4:nrow(tab),1]
coef_cpgs <- as.numeric(tab[4:nrow(tab),2])
intercept <- as.numeric(tab[3,2])
pred_age <- int_mat %*% matrix(coef_cpgs,ncol=1) + intercept
narm_idx <- which(!is.na(pred_age))
pred_age_narm <- pred_age[narm_idx] #Naあり
pred_ch_age <- rep(NA, length(pred_age_narm))
for(i in 1:length(pred_age_narm)){
    if(pred_age_narm[i] > 0){
        ch_age <- 21 * pred_age_narm[i] + 20
    }else{
        ch_age <- exp(pred_age_narm[i] + log(21)) - 1
    }
    pred_ch_age[i] <- ch_age
}

ages_narm <- ages[narm_idx]
png(paste0(out_path, "clock.png"), width=960, height=960)
par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
xlab = "Chlonological Age"
ylab = "Predicted Age"
main <- paste0("cor=", round(cor(ages_narm, pred_ch_age),2))
plot(ages_narm, pred_ch_age, xlab="", ylab="", cex.axis=4, cex=2, cex.lab=3)
abline(0, 1, col="red", lwd=3)
mtext(xlab, side=1, line=6, cex=4)
mtext(ylab, side=2, line=6, cex=4)
mtext(main, side=3, line=3, cex=4, adj=0)
dev.off()


#WGCNA analysis
n_genes <- 10000
X <- t(dat[order(apply(dat, 1, sd, na.rm=T), decreasing=T)[1:n_genes],])
sft = pickSoftThreshold(X)
softPower <- sft$powerEstimate
adjacency = adjacency(as.data.frame(X), power=softPower)
TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM
tmp_res <- fundamentalNetworkConcepts(adjacency)
Connectivity <- tmp_res[["Connectivity"]]
ClusterCoef <- tmp_res[["ClusterCoef"]]
MAR <- tmp_res[["MAR"]]


#Edge statistics
selflg <-  ifelse(colnames(X) %in% horvath_cpgs, "sel", "non")
sel_idx <- which(selflg == "sel")
nonsel_idx <- which(selflg == "non")
res <- matrix(NA, 3, 5)
colnames(res) <- c("Weight", "TOM", "Connectivity", "ClusterCoef", "MAR")
rownames(res) <- c("D", "P", "S-N")


#Weight and TOM
mat_list <- list(adjacency, TOM)
name_list <- list("Weight", "TOM")
for(i in 1:2){
    tmp_mat <- mat_list[[i]]
    sel_adj <- tmp_mat[sel_idx, sel_idx]
    nonsel_adj <- tmp_mat[nonsel_idx, nonsel_idx]
    ssv <- sel_adj[upper.tri(sel_adj)]
    nsv <- nonsel_adj[upper.tri(nonsel_adj)]
    png(paste0(out_path, name_list[[i]], ".png"), width=960, height=960)
    par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
    tmp_list <- list(ssv,nsv)
    names(tmp_list) <- c("Selected", "Non-selected")
    boxplot(tmp_list, cex.axis=4, cex.lab=4, main=name_list[[i]], cex.main=4)
    dev.off()

    #record
    res["D", name_list[[i]]] <- cohen.d(ssv, nsv)$estimate
    res["S-N", name_list[[i]]] <- mean(ssv) - mean(nsv)
    res["P", name_list[[i]]] <- t.test(ssv, nsv)$p.value
}      

#three betwork consept
feature_list <- list(Connectivity, ClusterCoef, MAR)
name_list <- list("Connectivity", "ClusterCoef", "MAR")
for(i in 1:3){
    feature <- feature_list[[i]]
    sel_val <- feature[sel_idx]
    non_val <- feature[nonsel_idx]
    png(paste0(out_path, name_list[[i]], ".png"), width=960, height=960)
    par(mar = c(9, 9, 9, 4)) ##bottom, left, top, right
    tmp_list <- list(sel_val, non_val)
    names(tmp_list) <- c("Selected", "Non-selected")
    boxplot(tmp_list, cex.axis=4, cex.lab=4, main=name_list[[i]], cex.main=4)
    dev.off()

    #record
    res["D", name_list[[i]]] <- cohen.d(sel_val, non_val)$estimate
    res["S-N", name_list[[i]]] <- mean(sel_val) - mean(non_val)
    res["P", name_list[[i]]] <- t.test(sel_val, non_val)$p.value
}

#out put
write.csv(res, file=paste0(out_path, "res.csv"))