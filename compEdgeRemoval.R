#!/usr/bin/env Rscript
require(igraph)
require(RColorBrewer)
require(gplots)

#############
# FUNCTION ##
#############
compCors <- function(Jfull, Nedges=2, Nr=20){

    set.seed(1)
    inds_edges <- which(abs(Jfull)>0, arr.ind=T)
    Cors <- NULL

    for (ir in 1:Nr){
        print(paste(ir,'/',Nr))
        rm('Mbio')

        inds_J <- inds_edges[sample(1:nrow(inds_edges),Nedges),]
        J <- Jfull
        J[inds_J] <- 0

        G = graph_from_adjacency_matrix(t(sign(J)),weighted=T)

        #################
        ## Biochemical ##
        #################

        Wp <- J
        try(Mbio <- solve(I-Wp),TRUE)
        if (!exists('Mbio')) next


        colnames(Mbio) = speciesNames
        rownames(Mbio) = speciesNames

        #diag(Mbio)=0
        Mbio <- sapply(1:ncol(Mbio), function(i) Mbio[,i] / Xss[i])
        Mbio <- sapply(1:ncol(Mbio), function(i) Mbio[,i] / Mbio[i,i])

        if (sum(is.finite(Mbio)) == 0) next

        diag(Mbio)=NA

        ###############
        ## Diffusion ##
        ###############

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
        Mprince_dir_sign <- t(solve(I-alpha * Wp) * (1-alpha))

        colnames(Mprince_dir_sign) = speciesNames
        rownames(Mprince_dir_sign) = speciesNames

        #diag(Mprince_dir_sign)=0
        Mprince_dir_sign <- sapply(1:ncol(Mprince_dir_sign), function(i) Mprince_dir_sign[,i] / Mprince_dir_sign[i,i])

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
        Mprince_dir <- t(solve(I-alpha * Wp) * (1-alpha))

        colnames(Mprince_dir) = speciesNames
        rownames(Mprince_dir) = speciesNames

        #diag(Mprince_dir)=0
        Mprince_dir <- sapply(1:ncol(Mprince_dir), function(i) Mprince_dir[,i] / Mprince_dir[i,i])


        #============#
        # Undirected #
        #============#


        W <- (abs(t(sign(J))) + abs(sign(J))) / 2
        D1 <- matrix(0,length(V(G)),length(V(G)))
        D2 <- matrix(0,length(V(G)),length(V(G)))
        # out degree prob
        diag(D1) <- ifelse(sqrt(apply(abs(W),1,sum))==0,0,1/sqrt(apply(abs(W),1,sum)))
        # in degree prob
        diag(D2) <- ifelse(sqrt(apply(abs(W),2,sum))==0,0,1/sqrt(apply(abs(W),2,sum)))
        # Wp[i,j] is Wp divided by out degree of i and in degree of j
        Wp <- D1 %*% W %*% D2

        I <- diag(length(V(G)))
        Mprince_undir <- t(solve(I-alpha * Wp) * (1-alpha))

        colnames(Mprince_undir) = speciesNames
        rownames(Mprince_undir) = speciesNames

        #diag(Mprince_undir)=0
        Mprince_undir <- sapply(1:ncol(Mprince_undir), function(i) Mprince_undir[,i] / Mprince_undir[i,i])


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
        Mdist_dir <- 0.9^distance_matrix
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
        Mdist_undir <- 0.9^distance_matrix
        Mdist_undir[is.na(Mdist_undir)] <- 0



        #####################
        ## First neighbors ##
        #####################

        Madj_dir_sign <- jitter(sign(J))

        Madj_dir <- jitter(abs(sign(J)))

        Madj_undir <- jitter((abs(t(sign(J))) + abs(sign(J))) / 2)

        method='s'

        cors <- NULL
        cors <- c(cors, cor.0(abs(as.numeric(Mbio_full)), abs(as.numeric(Mbio)), method=method))
        cors <- c(cors, cor.0(abs(as.numeric(Mbio_full)), abs(as.numeric(Mprince_dir_sign)), method=method))
        cors <- c(cors, cor.0(abs(as.numeric(Mbio_full)), abs(as.numeric(Mprince_dir)), method=method))
        cors <- c(cors, cor.0(abs(as.numeric(Mbio_full)), abs(as.numeric(Mprince_undir)), method=method))
        cors <- c(cors, cor.0(abs(as.numeric(Mbio_full)), jitter(abs(as.numeric(Mdist_dir))), method=method))
        cors <- c(cors, cor.0(abs(as.numeric(Mbio_full)), jitter(abs(as.numeric(Mdist_undir))), method=method))
        cors <- c(cors, cor.0(abs(as.numeric(Mbio_full)), jitter(abs(as.numeric(Madj_dir_sign))), method=method))
        cors <- c(cors, cor.0(abs(as.numeric(Mbio_full)), jitter(abs(as.numeric(Madj_dir))), method=method))
        cors <- c(cors, cor.0(abs(as.numeric(Mbio_full)), jitter(abs(as.numeric(Madj_undir))), method=method))
        Cors <- rbind(Cors, cors)

    }

    return(Cors)
}

##############
# Main loop ##
##############

dirs <- list.files('.','BIOM')

for (dir in dirs){

    print(dir)
    rm('Mbio')

    fjac <- paste0(dir,'/jacobian.csv')
    if (!file.exists(fjac)) next
        
    #f <- paste0(dir,'/res_perturbation.RData')
    #if (file.exists(f)) next

    J0 <- t(as.matrix(read.csv(fjac,row.names=1)))
    
    fss <- paste0(dir,'/steady_state.csv')
    if (!file.exists(fss)) next
    X_ss0 <- read.csv(paste0(dir,'/steady_state.csv'),header=F)$V1

    # remove satellite nodes
    mus <- apply(abs(J0),1,sum)
    inds_rm <- which(mus==0)
    if (length(inds_rm)>0){
        J <- J0[-inds_rm,-inds_rm]
        Xss <- X_ss0[-inds_rm]
    } else {
        J <- J0
        Xss <- X_ss0
    }


    if (is.null(dim(J))) next
    if (nrow(J) < 10) next

    Jfull <- J

    G = graph_from_adjacency_matrix(t(sign(J)),weighted=T)

    speciesNames <- colnames(J)
    inds_conv <- match(speciesNames, V(G)$name)

    alpha = 0.9

    #################
    ## Biochemical ##
    #################

    Wp <- Jfull

    I <- diag(length(V(G)))
    try(Mbio <- solve(I-Wp),TRUE)
    if (!exists('Mbio')) next

    colnames(Mbio) = speciesNames
    rownames(Mbio) = speciesNames

    Mbio <- sapply(1:ncol(Mbio), function(i) Mbio[,i] / Xss[i])
    Mbio <- sapply(1:ncol(Mbio), function(i) Mbio[,i] / Mbio[i,i])

    if (sum(is.finite(Mbio)) == 0) next

    diag(Mbio)=NA

    Mbio_full <- Mbio


    Nedges <- length(E(G))
    iseq <- c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99)

    # Let's remove 10 edges
    #
    cors_all <- list()
    for (i in iseq){
        print(i)
        Nrm <- floor(i * Nedges)
        #if (as.character(i) %in% names(cors_all)) next
        cors_all[[as.character(i)]] <- compCors(Jfull, Nrm, 5)
    }

    if (length(cors_all) != length(iseq)) next

    save(iseq,cors_all, file=paste0(dir,'/res_perturbation.RData'))

    cols <- NULL
    cols['Biochemical'] <- col2hex('black')
    cols['Diffusion (d+s)'] <- brewer.pal(9, 'Blues')[8]
    cols['Diffusion (d)'] <- rgb.0(brewer.pal(9, 'Blues')[8],0.7)
    cols['Diffusion (u)'] <- rgb.0(brewer.pal(9, 'Blues')[8],0.4)
    cols['Distance (d)'] <- brewer.pal(9, 'Oranges')[7]
    cols['Distance (u)'] <- rgb.0(brewer.pal(9, 'Oranges')[7],0.7)
    cols['Adjacency (d+s)'] <- brewer.pal(9, 'Greys')[7]
    cols['Adjacency (d)'] <- rgb.0(brewer.pal(9, 'Greys')[7],0.7)
    cols['Adjacency (u)'] <- rgb.0(brewer.pal(9, 'Greys')[7],0.4)

    mus_all <- NULL
    sds_all <- NULL
    for (i in as.numeric(names(cors_all))){
        print(i)
        Cors <- cors_all[[as.character(i)]]
        #plotCors(Cors,Nedge)
        mus <- apply(Cors,2,mean.0)
        sds <- apply(Cors,2,se.0)
        mus_all <- rbind(mus_all, mus)
        sds_all <- rbind(sds_all, sds)
    }


    Sbio <- as.numeric(Mbio_full)
    cc_r <- sapply(1:100, function(i) cor.0(sample(Sbio),sample(Sbio), method='s'))
    ss <- sd.0(cc_r)

    pdf(paste0(dir,'/plotAccuracyVsProportionEdgesRemoved.pdf'),width=10)
    par(mar=c(5,5,5,20))
    par(xpd=FALSE)
    inds_ord <- order(iseq)
    x <- iseq[inds_ord]
    x[1] <- x[2]/2 # for log only
    matplot(x, 
            mus_all[inds_ord,],
            type='o',
            log='x',
            lwd=3,
            lty=1,
            col=cols,
            pch=19,
            ylim=c(0,1),
            cex.lab=1.5,
            cex.axis=1.5,
            ylab='Spearman correlation (mean)',
            xlab='Proportion of edges removed')
    for (j in 1:ncol(mus_all)){
        plotCI.0(x,
                 t(mus_all[inds_ord,j]),
                 t(sds_all[inds_ord,j]),
                 pch=-1,
                 lwd=3,
                 sfrac=0.002,
                 col=cols[j],
                 add=T)
    }
    #random
    abline(h=2 * ss, lty=2, col='red', lwd=3)
    par(xpd=TRUE)
    #inds_sort <- order(-apply(mus_all[inds_ord[1:5],],2,mean.0))
    inds_sort <- order(-mus_all[1,])
    legend('topright',
           inset=c(-0.7,0),
           #leg=names(cols)[inds_sort],
           leg=names(cols),
           pt.cex=0.7,
           #col=cols[inds_sort],
           col=cols,
           lty=1,
           pch=19,
           lwd=3,
           bty='n',
           cex=1.3)
    dev.off()



}

