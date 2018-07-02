#!/usr/bin/env Rscript
require(MASS)
require(igraph)
require(RColorBrewer)
require(minerva)

dirs <- list.files('.','BIOM')

for (dir in dirs){

    print(dir)

    if (exists('Mbio')) rm('Mbio')

    fjac <- paste0(dir,'/jacobian.csv')
    if (!file.exists(fjac)) next

    J0 <- as.matrix(read.csv(fjac,row.names=1))

    # remove nodes with no out-degree
    mus <- apply(abs(J0),2,sum)
    inds_rm <- which(mus==0)
    if (length(inds_rm)>0){
        J <- J0[-inds_rm,-inds_rm]
    } else {
        J <- J0
    }

    if (is.null(dim(J))) next
    if (nrow(J) < 10) next

    G = graph_from_adjacency_matrix(t(sign(J)),weighted=T)

    speciesNames <- colnames(J)
    inds_conv <- match(speciesNames, V(G)$name)

    #################
    ## Biochemical ##
    #################

    Wp <- t(J)
    I <- diag(length(V(G)))
    try(Mbio <- solve(I-Wp),TRUE)

    if (!exists('Mbio')) next

    colnames(Mbio) = speciesNames
    rownames(Mbio) = speciesNames

    Mbio <- sapply(1:ncol(Mbio), function(i) Mbio[,i] / Mbio[i,i])

    if (sum(is.finite(Mbio)) == 0) next

    ###############
    ## Diffusion ##
    ###############
    alpha = .9

    #=================#
    # Directed signed #
    #=================#

    W <- t(sign(J))
    D1 <- matrix(0,length(V(G)),length(V(G)))
    D2 <- matrix(0,length(V(G)),length(V(G)))
    # out degree prob
    diag(D1) <- ifelse(sqrt(apply(abs(W),1,sum))==0,0,1/sqrt(apply(abs(W),1,sum)))
    # in degree prob
    diag(D2) <- ifelse(sqrt(apply(abs(W),2,sum))==0,0,1/sqrt(apply(abs(W),2,sum)))
    # Wp[i,j] is Wp divided by out degree of i and in degree of j
    Wp <- D1 %*% W %*% D2

    I <- diag(length(V(G)))
    Mprince_dir_sign <- solve(I-alpha * Wp) * (1-alpha)

    colnames(Mprince_dir_sign) = speciesNames
    rownames(Mprince_dir_sign) = speciesNames

    #==========#
    # Directed #
    #==========#

    W <- abs(t(sign(J)))
    D1 <- matrix(0,length(V(G)),length(V(G)))
    D2 <- matrix(0,length(V(G)),length(V(G)))
    # out degree prob
    diag(D1) <- ifelse(sqrt(apply(abs(W),1,sum))==0,0,1/sqrt(apply(abs(W),1,sum)))
    # in degree prob
    diag(D2) <- ifelse(sqrt(apply(abs(W),2,sum))==0,0,1/sqrt(apply(abs(W),2,sum)))
    # Wp[i,j] is Wp divided by out degree of i and in degree of j
    Wp <- D1 %*% W %*% D2

    I <- diag(length(V(G)))
    Mprince_dir <- solve(I-alpha * Wp) * (1-alpha)

    colnames(Mprince_dir) = speciesNames
    rownames(Mprince_dir) = speciesNames

    #============#
    # Undirected #
    #============#

    W <- t(abs(t(sign(J))) + abs(sign(J)))
    W[W>0] <- 1

    D1 <- matrix(0,length(V(G)),length(V(G)))
    D2 <- matrix(0,length(V(G)),length(V(G)))
    # out degree prob
    diag(D1) <- ifelse(sqrt(apply(abs(W),1,sum))==0,0,1/sqrt(apply(abs(W),1,sum)))
    # in degree prob
    diag(D2) <- ifelse(sqrt(apply(abs(W),2,sum))==0,0,1/sqrt(apply(abs(W),2,sum)))
    # Wp[i,j] is Wp divided by out degree of i and in degree of j
    Wp <- D1 %*% W %*% D2

    I <- diag(length(V(G)))
    Mprince_undir <- solve(I-alpha * Wp) * (1-alpha)

    colnames(Mprince_undir) = speciesNames
    rownames(Mprince_undir) = speciesNames

    ##############
    ## Distance ##
    ##############


    #==========#
    # Directed #
    #==========#

    distance_matrix <- matrix(NA,length(V(G)),length(V(G)))
    for (i in 1:length(V(G))){
        inds0 <- NULL
        for (j in 0:(length(V(G)-1))){
            inds <- which(speciesNames %in% names(ego(G,j,V(G)[i],'out')[[1]]))
            inds1 <- inds[!(inds %in% inds0)]
            distance_matrix[inds1,i] <- j
            inds0 <- c(inds0,inds)
        }
    }
    colnames(distance_matrix) = speciesNames
    rownames(distance_matrix) = speciesNames
    
    Mdist_dir <-  t(1 / (1 + distance_matrix))
    Mdist_dir[is.na(Mdist_dir)] <- 0
    
    #============#
    # Undirected #
    #============#

    distance_matrix <- matrix(NA,length(V(G)),length(V(G)))
    for (i in 1:length(V(G))){
        inds0 <- NULL
        for (j in 0:(length(V(G)-1))){
            inds <- which(speciesNames %in% names(ego(G,j,V(G)[i],'all')[[1]]))
            inds1 <- inds[!(inds %in% inds0)]
            distance_matrix[inds1,i] <- j
            inds0 <- c(inds0,inds)
        }
    }
    colnames(distance_matrix) = speciesNames
    rownames(distance_matrix) = speciesNames

    Mdist_undir <-  t(1 / (1 + distance_matrix))
    Mdist_undir[is.na(Mdist_undir)] <- 0

    #####################
    ## First neighbors ##
    #####################

    Madj_dir_sign <- t(sign(J))
    Madj_dir <- t(abs(sign(J)))
    Madj_undir <- t(abs(t(sign(J))) + abs(sign(J)))
    
    ###########
    ## Plot  ##
    ###########

    #=============#
    # correlation #
    #=============#

    method='spearman'

    
    # note that we take care of ties through jitter
    cors <- NULL
    cors <- c(cors, cor.0(as.numeric(Mbio), as.numeric(Mbio), method=method))
    cors <- c(cors, cor.0(abs(as.numeric(Mbio)), abs(as.numeric(Mprince_dir_sign)), method=method))
    cors <- c(cors, cor.0(abs(as.numeric(Mbio)), abs(as.numeric(Mprince_dir)), method=method))
    cors <- c(cors, cor.0(abs(as.numeric(Mbio)), abs(as.numeric(Mprince_undir)), method=method))
    cors <- c(cors, mean.0(sapply(1:10, function(i) cor.0(abs(as.numeric(Mbio)), jitter(abs(as.numeric(Mdist_dir)),amount=0.01), method=method))))
    cors <- c(cors, mean.0(sapply(1:10, function(i) cor.0(abs(as.numeric(Mbio)), jitter(abs(as.numeric(Mdist_undir)),amount=0.01), method=method))))
    cors <- c(cors, mean.0(sapply(1:10, function(i) cor.0(abs(as.numeric(Mbio)), jitter(abs(as.numeric(Madj_dir_sign)),amount=0.01), method=method))))
    cors <- c(cors, mean.0(sapply(1:10, function(i) cor.0(abs(as.numeric(Mbio)), jitter(abs(as.numeric(Madj_dir)),amount=0.01), method=method))))
    cors <- c(cors, mean.0(sapply(1:10, function(i) cor.0(abs(as.numeric(Mbio)), jitter(abs(as.numeric(Madj_undir)),amount=0.01), method=method))))



    cors_sign <- NULL
    cors_sign <- c(cors_sign, 1)
    s1 <- sign(as.numeric(Mbio)) 
    s2 <- sign(as.numeric(Mprince_dir_sign))
    inds_keep <-  1:length(s1)
    cors_sign <- c(cors_sign, sum(s1[inds_keep] == s2[inds_keep],na.rm=T) / length(inds_keep))
    s1 <- sign(as.numeric(Mbio)) 
    s2 <- sign(as.numeric(Madj_dir_sign))
    inds_keep <-  which(s1!=0 & s2!=0)
    inds_keep <-  1:length(s1)
    cors_sign <- c(cors_sign, sum(s1[inds_keep] == s2[inds_keep],na.rm=T) / length(inds_keep))

    cors_sign_r <- NULL
    for (ir in 1:10){
        cr <- NULL
        cr <- c(cr, 1)
        s1 <- sample(sign(as.numeric(Mbio)))
        s2 <- sign(as.numeric(Mprince_dir_sign))
        inds_keep <-  1:length(s1)
        cr <- c(cr, sum(s1[inds_keep] == s2[inds_keep],na.rm=T) / length(inds_keep))
        s2 <- sign(as.numeric(Madj_dir_sign))
        inds_keep <-  1:length(s1)
        cr <- c(cr, sum(s1[inds_keep] == s2[inds_keep],na.rm=T) / length(inds_keep))
        cors_sign_r <- rbind(cors_sign_r, cr)
    }
    cors_sign_r <- apply(cors_sign_r,2,mean.0)

    names(cors_sign) = c('Biochemical', 'Prince','Direct neighbors')

    cols <- NULL
    cols['Biochemical'] <- col2hex('black')
    cols['Propagation (d+s)'] <- brewer.pal(9, 'Blues')[8]
    cols['Propagation (d)'] <- rgb.0(brewer.pal(9, 'Blues')[8],0.7)
    cols['Propagation (u)'] <- rgb.0(brewer.pal(9, 'Blues')[8],0.4)
    cols['Distance (d)'] <- brewer.pal(9, 'Oranges')[7]
    cols['Distance (u)'] <- rgb.0(brewer.pal(9, 'Oranges')[7],0.7)
    cols['Adjacency (d+s)'] <- brewer.pal(9, 'Greys')[7]
    cols['Adjacency (d)'] <- rgb.0(brewer.pal(9, 'Greys')[7],0.7)
    cols['Adjacency (u)'] <- rgb.0(brewer.pal(9, 'Greys')[7],0.4)


    names(cors) <- names(cols)

    Sbio <- as.numeric(Mbio)
    cc_r <- sapply(1:100, function(i) cor.0(sample(Sbio),sample(Sbio), method=method))
    ss <- sd.0(cc_r)

    inds_ord <- 1:length(cors)

    cors[is.na(cors)] <- 0

    pdf(paste0(dir,'/plotBarplot_',method,'.pdf'))
    barplot.0(cors[inds_ord], 
              col=cols[inds_ord],
              mar=c(15,5,3,4),
              cex.names=1.2,
              ylab=paste0('Correlation (',method,')'),
              ylim=c(0,1),
              #main=fullname,
              cex=1)
    abline(h = 2 * ss, col='red', lty=2, lwd=3)
    dev.off()

    save(cors, cols,cors_sign, cors_sign_r,file=paste0(dir,'/res_cors_',method,'.RData'))
}

