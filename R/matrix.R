library('som')

SENT.norm <- function(feats, feat.weights = NULL){
    s = as.matrix(feats) %*% matrix(1,nrow=dim(feats)[2],ncol=1);
    s = matrix(100/s,nrow=length(s),ncol=dim(feats)[2]);
    feats.norm = feats * s;
    rm(s);

    feats.norm[is.na(feats.norm)] = 0

    if (!is.null(feat.weights)){
      feats.norm =  feats.norm  * matrix(abs(feat.weights),ncol=length(feat.weights),nrow=dim(feats.norm)[1],byrow=T)
    }
    
    feats.norm;
}

SENT.prepare.matrix <- function(file.input, file.output, file.dict= NULL){
    feats = read.table(file.input, sep="\t", header=T, row.names=1,check.names=FALSE);
    
    if (!is.null(file.dict)){
        feats.weights = as.matrix(read.table(file=file.dict, sep="\t", row.names=1));
    }else {
        feats.weights = NULL;
    }

    good.words = apply(feats,2,sum) > 0
    feats = feats[,good.words]
    feats = SENT.norm(feats, feats.weights);

    write.table(file=file.output, feats, sep="\t", quote=FALSE)
}

SENT.join.results <- function(prefix){
    files.w <- Sys.glob(paste(prefix,'.matrix_w.*',sep=""))
    files.h <- Sys.glob(paste(prefix,'.matrix_h.*',sep=""))
    
    data.w <- NULL
    for (file in files.w){
        data <- read.table(file, sep="\t", header=T, row.names=1, check.names=FALSE)
        if (is.null(data.w)){
            data.w = data
        }else{
            data.w <- cbind(data.w,data)
        }
    }

    write.table(file=paste(prefix,'.features',sep=''),t(data.w), sep="\t", quote=FALSE, row.names = F)


    data.h <- NULL
    for (file in files.h){
        data <- read.table(file, sep="\t", header=T, row.names=1, check.names=FALSE)
        if (is.null(data.h)){
            data.h = data
        }else{
            data.h <- rbind(data.h,data)
        }
    }

    write.table(file=paste(prefix,'.profiles',sep=''),t(data.h), sep="\t", quote=FALSE,col.names=F)
}


SENT.analyze <- function(prefix, output, clusters = NULL, num.words = 15){
    profiles     <- read.table(paste(prefix, '.profiles',sep=""),sep="\t", row.names=1, check.names=F);
    features     <- read.table(paste(prefix, '.features',sep=""),sep="\t", header=T, check.names=F);

    # Assume 10 repetitions
    if (is.null(clusters)){
        clusters = dim(features)[1] / 10 ;
    }

    # Form a clustering
    fdist = dist(features)
    hfeatures <- hclust(fdist, method="ward");
    cfeatures <- cutree(hfeatures, k=clusters);

    coph <- cor(fdist,cophenetic(hfeatures));
    write(coph, file = paste(output, '.cophenetic',sep=""));



    # Average between clusters
    profiles.merged = vector();
    features.merged = vector();
    for (i in levels(factor(cfeatures))){
        profiles.merged = cbind(profiles.merged, apply(as.matrix(profiles[,cfeatures==i]),1,mean, trim=0.1));
        features.merged = rbind(features.merged, apply(as.matrix(features[cfeatures==i,]),2,mean, trim=0.1));
    }


    rownames(profiles.merged) <- rownames(profiles);
    colnames(features.merged) <- colnames(features);

    write.table(file=paste(output,'.merged.profiles',sep=''),profiles.merged, sep="\t", quote=FALSE,col.names=F)
    write.table(file=paste(output,'.merged.features',sep=''),t(features.merged), sep="\t", quote=FALSE,col.names=F)

    # Hard assign genes to features
    profiles.bin = profiles.merged
    for (i in 1:dim(profiles.bin)[1] ){
        m                   = sort(profiles.bin[i,],index.return = T,decreasing = T)$ix[1];
        profiles.bin[i,]       = 0;
        profiles.bin[i,m]      = 1;
    }

    profiles.sorted = c();
    profiles.bin.sorted = c();
    glabels=c();



    fgroups = cfeatures[unlist(dendrapply(as.dendrogram(hfeatures), function(e) attr(e, "label")))];
    flabels = sapply(seq(1,dim(features)[1]-1), function(i){ if(fgroups[i] != fgroups[i+1]){ '___'}else{''}});
    flabels = c(flabels,'');

    flabels[unlist(dendrapply(as.dendrogram(hfeatures), function(e) attr(e, "label")))] = flabels;

    order=unique(fgroups);
    for (i in order){
        profiles.sorted = rbind(profiles[profiles.bin[,i]==1,], profiles.sorted);
        if (sum(profiles.bin[,i]==1) == 0) next
        glabels = c(rep('',sum(profiles.bin[,i]==1)-1),glabels);
        glabels = c('___',glabels);
    }

    # Produce heatmap image
    bitmap(file=paste(output,'.jpg',sep=""),type='jpeg',res=75);
    heatmap(as.matrix(profiles),Rowv=NA,Colv=as.dendrogram(hfeatures),xlab="Factors from 10 factorizations", ylab="Genes", labRow=glabels, labCol=flabels, margins=c(4,4));


    # Produce heatmap image for hard assignment 
    bitmap(file=paste(output,'.hard.jpg',sep=""),type='jpeg',res=75);
    heatmap(as.matrix(profiles.sorted),Rowv=NA,Colv=as.dendrogram(hfeatures),xlab="Factors from 10 factorizations", ylab="Genes", labRow=glabels, labCol=flabels, margins=c(4,4));

    dev.off();

    features.merged.scores = apply(features.merged,2,function(x){ 
                                   sapply(x,function(v){ 
                                          (v - (sum(x) - v)/(length(x) - 1))
                                          })
                                   })

#    a = 0.1
#    features.merge.specificity = apply(features.merged, 2, function(x){ sapply(x, function(v){ v / mean(x)})})
#    features.merge.importance  = apply(features.merged, 1, function(x){ sapply(x, function(v){ v / mean(x)})})
#    features.scores = a * t(features.merge.importance) + (1-a) * features.merge.specificity





# Save Group Genes and Words
    g = 1
    for (i in order){
        genes = rownames(profiles)[profiles.bin[,i] == 1];
        cat(file=paste(output,g,'genes',sep="."),genes,sep="\n");
        words = names(sort(features.merged.scores[i,],decreasing=T))[1:num.words];
##features.t.test = apply(features, 2, function(x){ t.test(x[cfeatures == i],x[cfeatures != i])})
##words = names(sort(sapply(features.t.test,function(x) { x$p.value}),decreasing=F))[1:num.words]
#words = names(sort(features.scores[i,],decreasing=T))[1:num.words]
        cat(file=paste(output,g,'words',sep="."),words,sep="\n");
        g = g + 1
    }
}


