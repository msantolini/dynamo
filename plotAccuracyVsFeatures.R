#!/usr/bin/env Rscript
require(RColorBrewer)
require(igraph)
library(car)


######################
## success features ##
######################
t <- read.csv('2017-02-24_BioModels-access-statistics.csv')

dnow <- as.Date('2017/02/24', format='%Y/%m/%d')
dall <- as.Date(t$Date.of.publication, format='%d/%m/%y')
days <- as.numeric(dnow-dall)

pdf('plotModelAccessVsAge.pdf')
plot.0(days,t$Access.hits,
       xlab='Model age (days)',
       ylab='Model access hits',
       log='y')
dev.off()

pdf('plotModelDownloadVsAge.pdf')
plot.0(days,t$Direct.downloads,
       xlab='Model age (days)',
       ylab='Model downloads',
       log='y')
dev.off()

pdf('plotModelDownloadVsAccess.pdf')
plot.0(t$Access.hits,t$Direct.downloads,
       xlab='Model access',
       ylab='Model downloads',
       log='xy')
dev.off()



#######################
## Network features  ##
#######################
dirs <- list.files('.','BIOM')

tconv <- read.delim('model_names.tsv')



dirs_keep <- NULL
model_size <- NULL
edge_density <- NULL
prop_reversible <- NULL
J_mean <- NULL
J_max <- NULL
betweenness_mean <- NULL
betweenness_max <- NULL
degree_max <- NULL
degree_out_max <- NULL
degree_in_max <- NULL
degree_out_mean <- NULL
degree_in_mean <- NULL
degree_mean <- NULL
degree_cv <- NULL
eigen_mean <- NULL
length_clusters_strong <- NULL
prop_clusters_strong <- NULL
nb_clusters_strong <- NULL
clustering <- NULL
structural_holes_mean <- NULL
structural_holes_max <- NULL
nb_infomap_communities <- NULL
access_hits <- NULL
direct_downloads <- NULL
model_age <- NULL
degree_cor_mean <- NULL
prop_infomap_communities <- NULL
for (dir in dirs){

    fjac <- paste0(dir,'/res_MICs.RData')
    if (!file.exists(fjac)) next
    fjac <- paste0(dir,'/jacobian.csv')

    if (!file.exists(fjac)) next
    J0 <- t(as.matrix(read.csv(fjac,row.names=1)))

    # remove satellite nodes
    mus <- apply(abs(J0),1,sum)
    inds_rm <- which(mus==0)
    J <- J0[-inds_rm,-inds_rm]

    if (is.null(dim(J))) next
    if (nrow(J) < 10) next

    model_size <- c(model_size, nrow(J))
    j <- c(J[upper.tri(J)], J[lower.tri(J)])
    edge_density <- c(edge_density, sum(j!=0) / (length(J)))

    jj <- j[j!=0]
    #coef_var_J <- c(coef_var_J, sd(jj) / mean(jj))
    J_mean <- c(J_mean, log10(mean(abs(jj))))
    J_max <- c(J_max, log10(max(abs(jj))))

    A <- abs(sign(J))
    cases <- A[upper.tri(A)] + A[lower.tri(A)]
    prop_reversible <- c(prop_reversible, sum(cases==2) / sum(cases>0))

    G = graph_from_adjacency_matrix(abs(t(sign(J))),weighted=T)
    betweenness_mean <- c(betweenness_mean, mean(edge_betweenness(G)))
    betweenness_max <- c(betweenness_max, max(edge_betweenness(G)))
    degree_max <- c(degree_max, max(degree(G)/(2 * length(V(G)))))
    degree_out_max <- c(degree_out_max, max(degree(G,mode='out')/(length(V(G)))))
    degree_in_max <- c(degree_in_max, max(degree(G,mode='in')/(length(V(G)))))
    degree_out_mean <- c(degree_out_mean, mean(degree(G,mode='out')/(length(V(G)))))
    degree_in_mean <- c(degree_in_mean, mean(degree(G,mode='in')/(length(V(G)))))
    degree_cv <- c(degree_cv, sd(degree(G))/ mean(degree(G)))
    degree_mean <- c(degree_mean, mean(degree(G)/(2 * length(V(G)))))
    eigen_mean <- c(eigen_mean, mean(eigen_centrality(G)$vector))
    length_clusters_strong <- c(length_clusters_strong, mean(clusters(G, mode='strong')$csize) / length(V(G)))
    prop_clusters_strong <- c(prop_clusters_strong, length(clusters(G, mode='strong')$csize) / length(V(G)))
    nb_clusters_strong <- c(nb_clusters_strong, length(clusters(G, mode='strong')$csize))
    clustering <- c(clustering, transitivity(G))
    structural_holes_max <- c(structural_holes_max, max(constraint(G)))
    structural_holes_mean <- c(structural_holes_mean, mean(constraint(G)))

    degree_cor_mean <- c(degree_cor_mean, assortativity_degree(G))
    
    ind <- match(dir,  t$BioModels.Identifier)
    access_hits <- c(access_hits, t$Access.hits[ind])
    direct_downloads <- c(direct_downloads, t$Direct.downloads[ind])
    model_age <- c(model_age, days[ind])
    

    nb_infomap_communities <- c(nb_infomap_communities, length(infomap.community(G)))
    prop_infomap_communities <- c(prop_infomap_communities, length(infomap.community(G))/length(V(G)))

    dirs_keep <- c(dirs_keep, dir)

}

dir.create('accuracy_vs_feature',showWarnings=F)

#for (method in c('kendall','spearman','MICs')){
for (method in c('spearman')){

    files <- list.files(dirs_keep,paste0('*res_cors_',method,'.RData'),full.names=T)

    Cors <- NULL
    names <- NULL
    for (file in files){
        names <- c(names, strsplit(file,'/')[[1]][1])
        ll <- load(file)
        if (method=='MICs'){
            Cors <- rbind(Cors, MICs)
        } else {
            Cors <- rbind(Cors, cors)
        }

    }


    #for (i in c(2,5)){
    for (i in c(2)){

        if (i==2){
            model <- 'Diffusion'
        } else if (i==5){
            model <- 'Distance'
        }


        #ylab = paste0(model,' model correlation (',method,')')
        ylab = 'Network model accuracy'

        pdf(paste0('accuracy_vs_feature/plotAccuracyVsLengthClustersStrong_',method,'_',model,'.pdf'))
        plot.cor(length_clusters_strong,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 fitline=TRUE,
                 ylab=ylab,
                 xlab='Mean length of strong connected components',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()

        pdf(paste0('accuracy_vs_feature/plotAccuracyVsPropClustersStrong_',method,'_',model,'.pdf'))
        plot.cor(prop_clusters_strong,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='x',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 fitline=TRUE,
                 fitlog=T,
                 ylab=ylab,
                 xlab='Proportion of strong connected components',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()


        pdf(paste0('accuracy_vs_feature/plotAccuracyVsNbClustersStrong_',method,'_',model,'.pdf'))
        plot.cor(nb_clusters_strong,
                 Cors[,i],
                 cut_by=0.05,
                 pch=19,
                 method='p',
                 log='x',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 fitline=F,
                 fitlog=F,
                 ylab=ylab,
                 xlab='Number of strongly connected components',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()

        pdf(paste0('accuracy_vs_feature/plotAccuracyVsNbClustersStrong_',method,'_',model,'_barplot.pdf'))
        barplot.split(nb_clusters_strong,
                      Cors[,i],
                      ylim=c(0,1),
                      ylab=ylab,
                      xlab='Number of strongly connected components',
                      main=paste(nrow(Cors),'models')
                      )
        dev.off()

        pdf(paste0('accuracy_vs_feature/plotAccuracyVsNbInfomapCommunities_',method,'_',model,'.pdf'))
        plot.cor(nb_infomap_communities,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 fitline=TRUE,
                 fitlog=T,
                 ylab=ylab,
                 xlab='Number of infomap communities',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()

        pdf(paste0('accuracy_vs_feature/plotAccuracyVsPropInfomapCommunities_',method,'_',model,'.pdf'))
        plot.cor(prop_infomap_communities,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 fitline=TRUE,
                 fitlog=T,
                 ylab=ylab,
                 xlab='Proportion of infomap communities',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()




        pdf(paste0('accuracy_vs_feature/plotAccuracyVsModelSize_',method,'_',model,'.pdf'))
        plot.cor(model_size,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 fitline=TRUE,
                 ylab=ylab,
                 xlab='Model size',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()


        pdf(paste0('accuracy_vs_feature/plotAccuracyVsEdgeDensity_',method,'_',model,'.pdf'))
        plot.cor(edge_density,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 #fitline=FALSE,
                 fitline=TRUE,
                 ylab=paste0(model, ' model correlation (',method,')'),
                 xlab='Edge density',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()


        pdf(paste0('accuracy_vs_feature/plotAccuracyVsPropReversible_',method,'_',model,'.pdf'))
        plot.cor(prop_reversible,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 #fitline=FALSE,
                 fitline=TRUE,
                 ylab=ylab,
                 xlab='Proportion of reversible reactions',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()

        pdf(paste0('accuracy_vs_feature/plotAccuracyVsJmean_',method,'_',model,'.pdf'))
        plot.cor(J_mean,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 #fitline=FALSE,
                 fitline=TRUE,
                 ylab=ylab,
                 xlab='Mean Jacobian (log10)',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()

        pdf(paste0('accuracy_vs_feature/plotAccuracyVsJmax_',method,'_',model,'.pdf'))
        plot.cor(J_max,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 #fitline=FALSE,
                 fitline=TRUE,
                 ylab=ylab,
                 xlab='Max Jacobian (log10)',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()


        pdf(paste0('accuracy_vs_feature/plotAccuracyVsBetweennessMean_',method,'_',model,'.pdf'))
        plot.cor(betweenness_mean,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 #fitline=FALSE,
                 fitline=TRUE,
                 ylab=ylab,
                 xlab='Mean edge betweenness',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()


        pdf(paste0('accuracy_vs_feature/plotAccuracyVsBetweennessMax_',method,'_',model,'.pdf'))
        plot.cor(betweenness_max,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 #fitline=FALSE,
                 fitline=TRUE,
                 ylab=ylab,
                 xlab='Max edge betweenness',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()

        pdf(paste0('accuracy_vs_feature/plotAccuracyVsDegreeMean_',method,'_',model,'.pdf'))
        plot.cor(degree_mean,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 #fitline=FALSE,
                 fitline=TRUE,
                 ylab=ylab,
                 xlab='Mean degree (normalized by number of nodes)',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()

        pdf(paste0('accuracy_vs_feature/plotAccuracyVsDegreeMax_',method,'_',model,'.pdf'))
        plot.cor(degree_max,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 #fitline=FALSE,
                 fitline=TRUE,
                 ylab=ylab,
                 xlab='Max degree (normalized by number of nodes)',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()

        pdf(paste0('accuracy_vs_feature/plotAccuracyVsEigenMean_',method,'_',model,'.pdf'))
        plot.cor(eigen_mean,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 #fitline=FALSE,
                 fitline=TRUE,
                 ylab=ylab,
                 xlab='Mean eigenvector centrality',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()

        pdf(paste0('accuracy_vs_feature/plotAccuracyVsDegreeCV_',method,'_',model,'.pdf'))
        plot.cor(degree_cv,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 #fitline=FALSE,
                 fitline=TRUE,
                 ylab=ylab,
                 xlab='Degree coefficient of variation',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()



        pdf(paste0('accuracy_vs_feature/plotAccuracyVsOutDegreeMax_',method,'_',model,'.pdf'))
        plot.cor(degree_out_max,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 #fitline=FALSE,
                 fitline=TRUE,
                 ylab=ylab,
                 xlab='Max out-degree (normalized by number of nodes)',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()

        pdf(paste0('accuracy_vs_feature/plotAccuracyVsInDegreeMax_',method,'_',model,'.pdf'))
        plot.cor(degree_in_max,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 #fitline=FALSE,
                 fitline=TRUE,
                 ylab=ylab,
                 xlab='Max in-degree (normalized by number of nodes)',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()

        pdf(paste0('accuracy_vs_feature/plotAccuracyVsOutDegreeMean_',method,'_',model,'.pdf'))
        plot.cor(degree_out_mean,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 #fitline=FALSE,
                 fitline=TRUE,
                 ylab=ylab,
                 xlab='Mean out-degree (normalized by number of nodes)',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()

        pdf(paste0('accuracy_vs_feature/plotAccuracyVsInDegreeMean_',method,'_',model,'.pdf'))
        plot.cor(degree_in_mean,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 #fitline=FALSE,
                 fitline=TRUE,
                 ylab=ylab,
                 xlab='Mean in-degree (normalized by number of nodes)',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()


        pdf(paste0('accuracy_vs_feature/plotAccuracyVsDegreeCorMean_',method,'_',model,'.pdf'))
        plot.cor(degree_cor_mean,
                 Cors[,i],
                 pch=19,
                 type='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 #fitline=FALSE,
                 fitline=TRUE,
                 ylab=paste0(model,' model correlation (',method,')'),
                 xlab='Degree assortativity',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()


        pdf(paste0('accuracy_vs_feature/plotAccuracyVsAccessHits_',method,'_',model,'.pdf'))
        plot.cor(log10(access_hits), 
                 Cors[,i],
                 pch=19,
                 type='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 fitline=TRUE,
                 #fitlog=TRUE,
                 ylab=paste0(model,' model correlation (',method,')'),
                 xlab='Model access hits on BioModels (log10)',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()


        pdf(paste0('accuracy_vs_feature/plotAccuracyVsDownloads_',method,'_',model,'.pdf'))
        plot.cor(log10(direct_downloads), 
                 Cors[,i],
                 pch=19,
                 type='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 fitline=TRUE,
                 #fitlog=TRUE,
                 ylab=paste0(model,' model correlation (',method,')'),
                 xlab='Model downloads on BioModels (log10)',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()


        pdf(paste0('accuracy_vs_feature/plotAccuracyVsModelAge_',method,'_',model,'.pdf'))
        plot.cor(model_age, 
                 Cors[,i],
                 pch=19,
                 type='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 fitline=TRUE,
                 #fitlog=TRUE,
                 ylab=paste0(model,' model correlation (',method,')'),
                 xlab='Model age on BioModels (days)',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()



        pdf(paste0('accuracy_vs_feature/plotAccuracyVsClustering_',method,'_',model,'.pdf'))
        plot.cor(clustering,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 #fitline=FALSE,
                 fitline=TRUE,
                 ylab=ylab,
                 xlab='Clustering',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()

        pdf(paste0('accuracy_vs_feature/plotAccuracyVsStructuralHolesMean_',method,'_',model,'.pdf'))
        plot.cor(structural_holes_mean,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 #fitline=FALSE,
                 fitline=TRUE,
                 ylab=ylab,
                 xlab='Structural holes (mean)',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()

        pdf(paste0('accuracy_vs_feature/plotAccuracyVsStructuralHolesMax_',method,'_',model,'.pdf'))
        plot.cor(structural_holes_max,
                 Cors[,i],
                 pch=19,
                 method='p',
                 log='',
                 mar=c(5,5,5,5),
                 #ylim=c(0,1),
                 #fitline=FALSE,
                 fitline=TRUE,
                 ylab=ylab,
                 xlab='Structural holes (max)',
                 main=paste(nrow(Cors),'models'),
                 #main=fullname,
                 cex=1)
        dev.off()

        df <- data.frame('id'=names,
                         'name'= tconv$Name[match(names,tconv$ID)],
                         'accuracy.diffusion.d.s' = Cors[,2],
                         'accuracy.diffusion.d' = Cors[,3],
                         'accuracy.diffusion.u' = Cors[,4],
                         'accuracy.distance.d' = Cors[,5],
                         'accuracy.distance.u' = Cors[,6],
                         'accuracy.first.neighbors.d.s' = Cors[,7],
                         'accuracy.first.neighbors.d' = Cors[,8],
                         'accuracy.first.neighbors.u' = Cors[,9],
                         'model.size'=model_size,
                         # 'length.clusters.strong'=length_clusters_strong,
                         'nb.clusters.strong'=nb_clusters_strong,
                         'prop.clusters.strong'=prop_clusters_strong,
                         'nb.infomap.communities'=nb_infomap_communities,
                         'prop.infomap.communities'=prop_infomap_communities,
                         'edge.density'=edge_density,
                         'prop.reversible.equations'=prop_reversible,
                         'J.mean.log10'=J_mean,
                         'J.max.log10'=J_max,
                         'betweenness.mean'=betweenness_mean,
                         'betweenness.max'=betweenness_max,
                         'degree.mean'=degree_mean,
                         'degree.max'=degree_max,
                         'eigenvector.mean'=eigen_mean,
                         'degree.cv'=degree_cv,
                         'degree.out.mean'=degree_out_mean,
                         'degree.in.mean'=degree_in_mean,
                         'degree.out.max'=degree_out_max,
                         'degree.in.max'=degree_in_max,
                         'structural.holes.mean'=-structural_holes_mean, # negative to get number of structural holes
                         'structural.holes.max'=-structural_holes_max
                         )
        write.csv(df[order(-df$accuracy.diffusion.d.s),],file=paste0('accuracy_vs_feature_summary.csv'),row.names=F)

        pdf('plotPairsAccuracies.pdf',9,9)
        #pairs(Cors[,-1],
              #cex=0.3,
              #xlim=c(0,1),
              #ylim=c(0,1),
              #main='Comparing accuracies across DYNAMOs')
        scatterplotMatrix(Cors[,-1],
                           pch=19,
                           smoother=F,
                           reg.line=F,
                           cex=0.3,
                           xlim=c(0,1),
                           ylim=c(0,1),
                           main='Comparing accuracies between DYNAMOs across 87 BioModels')
        dev.off()



        cc <- cor.0(df[,-(1:2)],method='s')
        cc_r <- cor.0(apply(df[,-(1:2)],2,sample),method='s')
        diag(cc_r) <- NA

        #m <- mean.0(cc_r)
        m <- 0
        s <- sd.0(cc_r)
        h <- m + 2*s

        pdf('plotHeatmapAccuracy.pdf')
        br <- c(seq(0,1,h)[-1],1)
        breaks <- c(-rev(br),br)
        heatmap.sym(cc,
                    cellnote=T,
                    notecex=0.4,
                    breaks=breaks
                    )
        dev.off()


        labels <- c('prop.clusters.strong'='Number of strong connected components',
                    'structural.holes.mean'='Number of structural holes',
                    'betweenness.mean'='Mean betweenness centrality',
                    'model.size'='Model size',
                    'prop.reversible.equations'='Number of reversible equations',
                    'eigenvector.mean'='Mean eigenvector centrality',
                    'J.mean.log10'='Mean Jacobian value',
                    'degree.mean'='Mean degree',
                    'edge.density'='Edge density'
                    )
        inds_keep <- match(names(labels), colnames(cc))

        cc1 <- cc[1,inds_keep]
        names(cc1) <- labels


        ################
        ## Error bars ##
        ################
        uiw <- tanh(atanh(cc1) + 1.96/sqrt(nrow(df)-3)) - cc1
        #uiw[cc1<0] <- NA
        liw <- cc1 - tanh(atanh(cc1) - 1.96/sqrt(nrow(df)-3))
        #liw[cc1>0] <- NA

        col <- rep('black',length(cc1))
        col[cc1 > h] <- 'darkred'
        col[cc1 < -h] <- 'darkblue'
        inds_ord <- order(cc1)

        mm <- max(abs(cc1+uiw),abs(cc1-liw))
        mm <- 1

        pdf('plotBarplotAccuracy.pdf',width=12,height=7)
        b <- barplot.0(cc1[inds_ord],
                       horiz=T,
                       xlab='Correlation with network model accuracy',
                       las=1,
                       cex=1.3,
                       #col=col[inds_ord],
                       col='white',
                       xlim=c(-mm,mm),
                       mar=c(5,22,5,5)
                       )
        polygon(c(-h,-h,h,h), c(min(b)-1,max(b)+1,max(b)+1,min(b)-1),
                col = "lightgrey", border = NA)
        #abline(v=h,lty=2,lwd=2,col='gray')
        #abline(v=-h,lty=2,lwd=2,col='gray')
        plotCI.0(cc1[inds_ord],
                 b,
                 uiw=uiw[inds_ord],
                 liw=liw[inds_ord],
                 col=col[inds_ord],
                 pch=19,
                 cex=1.5,
                 lwd=2,
                 barcol='black',
                 sfrac=0.0001,
                 err='x'
                 )
        dev.off()

    }
}



inds <- which(df$length.clusters.strong == 1)


li <- list(df$accuracy.diffusion[inds],
           df$accuracy.diffusion[-inds])
names(li) <- c(paste0('Fully strongly connected (',length(inds),')'),
               paste0('Other topologies (',nrow(df)-length(inds),')'))

pdf('plotStronglyConnected.pdf',height=8)
boxplot.0(li,
          las=3,
          mar=c(18,5,3,3),
          cex=1.2,
          ylim=c(0,1),
          ylab='Diffusion model accuracy')
dev.off()
