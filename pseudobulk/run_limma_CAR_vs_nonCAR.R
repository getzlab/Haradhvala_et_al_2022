library(limma)
library(edgeR)
library(dplyr)

pseudobulk_data <- "../data/ALL_T_cells_pseudobulk_nosubtypesplit.v3.csv"

for(type in c("Infusion","D7")) {
        X <- read.csv(pseudobulk_data)
        X$type <- gsub("-","_",X$type)
        X <- X[grepl(type,X$type),]
    
        X$type <- sub(paste0(type,"_"),"",X$type)

        counts <- t(X[,c(-1,-2,-3,-4,-5,-6)])

        colnames(counts) <- paste(X$barcode,X$type,sep="_")
        
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

            Yescarta = (Yescarta.CAR_T - Yescarta.nonCAR ),
            Kymriah = (Kymriah.CAR_T - Kymriah.nonCAR ),
        levels=design)

        fit2 <- contrasts.fit(fit, cm)
        fit2 <- eBayes(fit2)
        
        for(product in c("Kymriah","Yescarta")) {
            res <- topTable(fit2,product,n=Inf)
            write.csv(res,paste0(results.dir,product,".",type,".csv"))
        }
        
}
