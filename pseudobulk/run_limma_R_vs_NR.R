library(limma)
library(edgeR)
library(dplyr)

results.dir <- "../data/limma_results/limma_response_comparisons/"

Xall<- read.csv("../data/ALL_T_cells_pseudobulk.v2.csv")

Xall$type <- gsub("-","_",Xall$type)

# Remove patient with only PR
Xall <- Xall[Xall$response!='',]

for(subtype in c("CD4","CD8")) {
    for(product in c("Kymriah","Yescarta")) {
        for(type in c("Baseline","D7","Infusion","D7_CAR_T")) {
    X <- Xall[grepl(subtype,Xall$subtype)&(Xall$product==product)&(Xall$type==type),]
    counts <- t(X[,c(-1,-2,-3,-4,-5,-6,-7)])

    colnames(counts) <- paste(X$barcode,X$subtype,X$type,sep="_")
    
    dge <- DGEList(counts=counts)
    dge <- calcNormFactors(dge)

    condition <- X$response

    design <- model.matrix(~0+condition + X$gender)
    colnames(design) <- sub("condition","",colnames(design))

    colnames(design)[dim(design)[2]] <- "male"

    v <- voom(dge, design, plot=FALSE)
    
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    
    
    cm <- makeContrasts(
    resp = (R-NR),
    levels=design)

    fit2 <- contrasts.fit(fit, cm)
    fit2 <- eBayes(fit2)
     
    res <- topTable(fit2,"resp",n=Inf)

            
    # Save as standalone csv
    write.csv(res,paste0(results.dir,type,".",product,".",subtype,".csv"))
        }
    }
}


