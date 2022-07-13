library(limma)
library(edgeR)
library(dplyr)

results.dir <- "../data/limma_results/"

Xall<- read.csv("../data/ALL_T_cells_pseudobulk.v2.csv")
Xall$type <- gsub("-","_",Xall$type)

Xall <- Xall[grepl("Infusion|CAR",Xall$type),]

for(subtype in c("CD4","CD8")) {
    X <- Xall[grepl(subtype,Xall$subtype),]
    counts <- t(X[,c(-1,-2,-3,-4,-5,-6,-7)])

    colnames(counts) <- paste(X$barcode,X$subtype,X$type,sep="_")
    
    dge <- DGEList(counts=counts)
    dge <- calcNormFactors(dge)

    condition <- paste(X$product,X$type,sep=".")

    design <- model.matrix(~0+condition  + X$gender)
    colnames(design) <- sub("condition","",colnames(design))

    colnames(design)[dim(design)[2]] <- "male"

    v <- voom(dge, design, plot=FALSE)
    corfit <- duplicateCorrelation(v,design,block=X$barcode)
    v <- voom(dge, design, plot=TRUE,block=X$barcode,correlation=corfit$consensus)
    
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    
    cm <- makeContrasts(
    IP_KvsY = (Kymriah.Infusion - Yescarta.Infusion),
    D7_KvsY = (Kymriah.D7_CAR_T - Yescarta.D7_CAR_T),
    Kymriah_IP_vs_D7 = (Kymriah.D7_CAR_T - Kymriah.Infusion ), 
    Yescarta_IP_vs_D7 = (Yescarta.D7_CAR_T - Yescarta.Infusion ),
    KvsY = (Kymriah.D7_CAR_T - Kymriah.Infusion ) - (Yescarta.D7_CAR_T - Yescarta.Infusion ),
    levels=design)

    fit2 <- contrasts.fit(fit, cm)
    fit2 <- eBayes(fit2)
    
    tests <- colnames(fit2$coefficients)
    
    res_list <- list()
    for(test in tests) {
        res <- topTable(fit2,test,n=Inf)
    
        # Save as standalone csv
        write.csv(res,paste0(results.dir,test,".",subtype,".csv"))
}
}
