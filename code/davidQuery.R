library(RDAVIDWebService)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
source("C:/Users/15000/Desktop/taiji_/Rscripts/plotGO.R")
davidQuery <- function(fl, species="hs", is.trans=F, key="helloworld",
                       w=5, h=8, trailing=40, is.plot=T, david = david){
  geneL <- readLines(fl)
  if (species=="hs"){
    egids <- select(org.Hs.eg.db, geneL, "ENTREZID", "SYMBOL")[,2]
  }else if (species=="mm"){
    firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }
    geneL <- firstup(tolower(geneL))
    egids <- select(org.Mm.eg.db, geneL, "ENTREZID", "SYMBOL")[,2]
  }
  
  david <- DAVIDWebService$new(email='col003@ucsd.edu', url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  
  FG <- addList(david, egids, idType = "ENTREZ_GENE_ID", 
                listName = gsub(".txt","",fl), listType = "Gene")
  setAnnotationCategories(david, c("GOTERM_BP_DIRECT", "KEGG_PATHWAY"))
  FuncAnnotChart <- getFunctionalAnnotationChart(david)
  df <- as.data.frame(FuncAnnotChart@.Data)
  names(df) <- FuncAnnotChart@names
  
  # convert gene id to gene name
  conv <- data.frame("from"=egids, "to"=toupper(geneL))
  geneL <- strsplit(df$Genes, split = ", ")
  df$Genes <- sapply(geneL, function(x) paste0(conv[match(x, conv$from),"to"], collapse = ", "))
  
  write.table(df, file = paste0(gsub(".txt","",fl), "_kegg_go.txt"))
  if (is.plot){
    plotGO(df, fl = fl, is.trans = is.trans, key=key,w=w,h=h,trailing = trailing)
  }
  
}
