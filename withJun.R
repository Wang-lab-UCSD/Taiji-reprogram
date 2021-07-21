#### import packages-----------------------------------------------------------
library(gtools)
library(fmsb)
library(RColorBrewer)
library(scales)
library(kableExtra)
library(readxl)
library(tidyverse)
library(gridExtra)
library(webshot)
library(magick)
library(png)
setwd("C:/Users/15000/Desktop/Taiji/dataset/withJun/oe_laststep/recipes/")
source("C:/Users/15000/Desktop/taiji_/Rscripts/getRecipes.R")
source("C:/Users/15000/Desktop/taiji_/Rscripts/getRecipes_simplified.R")
source("C:/Users/15000/Desktop/taiji_/Rscripts/bubble.R")
source("C:/Users/15000/Desktop/taiji_/Rscripts/davidQuery.R")
options(max.print=1000000)
windowsFonts("Arial" = windowsFont("Arial"))


#### import Data----------------------------------------------------------------
# pagerank scores
Data <- read.table("../../pagerank_mine_withHMSC160_20200609.tsv", check.names = F)
Data <- Data[!rownames(Data) %in% c("DNMT1"),] # remove DNMT1

# group info
group_sorted <- read_excel("../../tissue_154_info.xlsx")
names(group_sorted)[4] <- "Group"
rownames(group_sorted) <- group_sorted$name
Data <- Data[,rownames(group_sorted)]

# rna
rna_seq <- read.table("../../expression_profile_with_6yue_20200609.tsv", check.names = F)
row.names(rna_seq) <- toupper(rownames(rna_seq))
rna_seq <- rna_seq[complete.cases(rna_seq),]
rna_seq <- rna_seq[,rownames(group_sorted)]

# update new TFnames
TFnames <- intersect(rownames(rna_seq),rownames(Data))
Data <- Data[TFnames,]
rna_seq <- rna_seq[TFnames,]
Data_normed <- as.data.frame(t(scale(t(as.matrix(Data)))))

#### generate repgramming recipes ####
s = "fibroblast"
t = c("stem","heart","myoblast","trophoblast","liver", "B.cell", "neural.cell", "kidney")
recipe_size = 3
rank_cut = 30
p = 0.001
method="prod"
m = "ratio.abs"

lapply(t, function(x) getRecipes(s, x, m=m,rank_cut = rank_cut))
getRecipes("B.cell", "monocyte", m=m,rank_cut = rank_cut)
getRecipes("liver", "neural.cell",m=m,rank_cut = rank_cut)
getRecipes("P13", "P6",m=m,rank_cut = rank_cut)

# get all the possible combinations
cellTypes <- unique(group_sorted$Group)
cbs <- as.data.frame(t(combn(cellTypes,2)))
mapply(function(x,y) getRecipes(soN = x, taN = y, rank_cut = rank_cut), 
       cbs$V1,cbs$V2)

tissues <- rownames(group_sorted)
cbs <- as.data.frame(t(combn(tissues,2)))
pts <- as.data.frame(gtools::permutations(n=160,r=2,v=tissues))
setwd("realThree_v2/")
mapply(function(x,y) getRecipes_simplified(name_s = x, name_t = y, rank_cut = rank_cut), 
       pts$V1,pts$V2)

getRecipes_simplified("R.emb.left.lung","R.emb.heart")

#### benckmark with experimental results ####
getRanks <- function(soN=s, taN, recipe_size, rank_cut=20, p=0.001, 
                     method="prod", m, reprog.type=c("oe"),
                     truth, No){
  sources <- rownames(group_sorted)[group_sorted$Group==soN]
  targets <-  rownames(group_sorted)[group_sorted$Group==taN]
  
  # choose top8 and correct reprog.type
  helper <- function(idx_s, idx_t){
    name_s <- sources[idx_s]
    name_t <- targets[idx_t]
    fn <- paste0("recipe_",name_s, "_to_",name_t,"_",method,"_",m,"_rank_",rank_cut,"_p_",p,"_size_",recipe_size,".txt")
    df <- read.delim(fn, sep = "")
    
    tmp <- df[which(df$type %in% reprog.type),]
    df2 <- tmp[1:min(8,nrow(tmp)),] 
    
    tp <- which(df2$TF %in% truth)
    tmp <- data.frame("count" = length(tp), "rank.sum"=sum(tp),
                      "genes" = paste(df2$TF[tp], collapse = ', '), "file"=fn)
    return(tmp)
  }
  
  output <- do.call(rbind,lapply(1:length(sources), function(x) do.call(rbind,lapply(1:length(targets), function(y) helper(x,y)))))
  output <- output[order(output$count, output$rank.sum, decreasing = c(T,F)),]
  write.table(output, paste0("benchmark_",soN, "_to_",taN,"_v",No,"_",
                             paste(reprog.type, collapse = ''),"_",m,"_",rank_cut,"_p_",p,"_size_",recipe_size,".txt"),
              row.names = F, quote = F)
}

## 1. fibroblast -> stem
getRanks(taN = t[1], truth = c("SOX2", "POU5F1", "KLF4","MYC"), recipe_size = 3, No = 1,m=m,rank_cut = rank_cut)
getRanks(taN = t[1], truth = c("SOX2", "POU5F1"), recipe_size = 3, No = 2,m=m,rank_cut = rank_cut)
getRanks(taN = t[1], truth = c("SOX2", "POU5F1", "NANOG"), recipe_size = 3, No = 3,m=m,rank_cut = rank_cut)

## 2. fibroblast -> heart
getRanks(taN = t[2], truth = c("GATA4", "MEF2C", "TBX5", "ESRRG"), recipe_size = 3, No = 1,m=m,rank_cut = rank_cut)
getRanks(taN = t[2], truth = c("GATA4", "MEF2C", "TBX5", "HAND2","NKX2.5"), recipe_size = 3, No = 2,m=m,rank_cut = rank_cut)
getRanks(taN = t[2], truth = c("GATA4", "MEF2C", "TBX5"), recipe_size = 3, No = 3,m=m,rank_cut = rank_cut)
getRanks(taN = t[2], truth = c("GATA4", "MEF2C", "TBX5", "ESRRG", "HAND2","NKX2.5"), recipe_size = 3, No = "summary",m=m,rank_cut = rank_cut)

## 3. fibroblast -> myoblast
getRanks(taN = t[3], truth = c("MYOD1"), recipe_size = 3, No = 1,m=m,rank_cut = rank_cut)

## 4. fibroblast -> trophoblast
getRanks(taN = t[4],truth = c("GATA3", "EOMES", "TFAP2C"), recipe_size = 3, No = 1,m=m,rank_cut = rank_cut)

## 5. fibroblast -> liver
getRanks(taN = t[5], truth = c("FOXA2", "HNF4A", "CEBPB"), recipe_size = 3, No = 1,m=m,rank_cut = rank_cut)
getRanks(taN = t[5], truth = c("HNF1A", "HNF4A", "ATF5", "PROX1", "CEBPA"), recipe_size = 3, No = 2,m=m,rank_cut = rank_cut)
getRanks(taN = t[5], truth = c("HNF1A", "HNF4A", "GATA4", "FOXA3"), recipe_size = 3, No = 3,m=m,rank_cut = rank_cut)
getRanks(taN = t[5], truth = c("HNF1A", "FOXA1", "FOXA2", "FOXA3"), recipe_size = 3, No = 4,m=m,rank_cut = rank_cut)
getRanks(taN = t[5], truth = c("HNF1A", "HNF4A", "ATF5", "PROX1", "CEBPA", "CEBPB", "FOXA1",  "FOXA2", "FOXA3", "GATA4"), recipe_size = 3, No = "summary",m=m,rank_cut = rank_cut)

## 6. fibroblast -> B.cell
getRanks(taN = t[6], truth = c("SPI1", "IRF8", "BATF3"), recipe_size = 3, No = 1,m=m,rank_cut = rank_cut)

## 7. fibroblast -> neural.cell
getRanks(taN = t[7], truth = c("FOXG1", "SOX2", "DLX5", "LHX6"), recipe_size = 3, No = 1,m=m,rank_cut = rank_cut)

## 8. fibroblast -> kidney
getRanks(taN = t[8], truth = c("EMX2", "HNF1B", "HNF4A", "PAX8"), recipe_size = 3, No = 1,m=m,rank_cut = rank_cut)

## 9. B.cell -> monocyte
getRanks(soN = "B.cell",taN = "monocyte", truth = c("CEBPA"), recipe_size = 3, No = 1,m=m,rank_cut = rank_cut)
getRanks(soN = "B.cell",taN = "monocyte", truth = c("CEBPA", "SPI1"), recipe_size = 3, No = 2,m=m,rank_cut = rank_cut)

#### visualization ####

#### 1. tables ####
setwd("C:/Users/15000/Desktop/Taiji/dataset/withJun/oe_laststep/tables/")
path <- "C:/Users/15000/Desktop/Taiji/dataset/withJun/benchmark_based_on_sheet32.xlsx"
recipes <- path %>%
  excel_sheets() %>%
  set_names() %>%
  map(read_excel,
      path = path)
str(recipes)

getKable <- function(no){
  dt <- recipes[[no]][3:6,1:9]
  dt$algorithms[4] <- "Taiji-reprogram"
  # remove the row if no prediction is made
  dt <- dt[rowSums(is.na(dt)) != ncol(dt)-1, ]
  
  truth <- unlist(recipes[[no]][7,2:9])
  truth <- toupper(truth[!is.na(truth)])
  dt$count <- unlist(lapply(1:nrow(dt), function(x) {length(which(dt[x,2:9]%in%truth))}))
  dt$rank.avg <- unlist(lapply(1:nrow(dt), function(x) {round(sum(which(dt[x,2:9]%in%truth))/length(which(dt[x,2:9]%in%truth)),3)}))
  # calculate scores as benchmark metric
  l <- length(truth)
  dt$scores <- round(dt$count/l * sum(1:l)/(l*dt$rank.avg),3)
  dt$scores <- round(((dt$count/dt$rank.avg)/(l^2/sum(1:l)))^0.5*100,1)
  dt$scores[which(!is.finite(dt$scores))] <- 0
  dt$rank.avg[which(!is.finite(dt$rank.avg))] <- 0
  write.table(dt,paste0(no,".txt"),row.names = F, quote = F, sep = "\t")
  
  dt[1] <- cell_spec(dt[[1]], bold = T,
                     color = ifelse(dt$scores==max(dt$scores), "white", "grey"),
                     background = ifelse(dt$scores==max(dt$scores), "#DC143C", "NULL"))
  dt[2:9] <- lapply(dt[2:9], function(x) {
    cell_spec(x,
              bold = ifelse(x%in%truth, "T", "F"),
              color = ifelse(x%in%truth, "white", "grey"),
              background = ifelse(x%in%truth, "	#32CD32", "NULL"))
  })
  dt %>%
    kbl(escape = F, format = "html",
        caption = paste0("<center><strong>",no, ": ", paste0(truth,collapse = ", "),"</strong></center>"),) %>%
    kable_styling("striped",full_width = F)  %>%
    save_kable(file=paste0(no,".png"), zoom = 5)
  
  
}

lapply(names(recipes),getKable)


#### 2. radar plot ####
read_plus <- function(flnm) {
  read_delim(flnm, delim = " ")[,c("algorithms","scores")] %>% 
    mutate(filename = sub("./", "", sub(".txt", "", flnm))) 
}

# add the name of subdirectory to the dataframe as a new column
dt <- list.files(path="./", pattern = "[0-9].txt",recursive=T, full.names = T) %>% 
  map_df(~read_plus(.)) %>% 
  filter(algorithms%in%c("Mogrify","Taiji-reprogram")) %>%
  spread(filename, scores) %>%
  column_to_rownames(var="algorithms") 

getRadarPlot <- function(df,fl,is.reverse=F,seg=4){
  
  rownames(df) <- c("Mogrify", "Taiji")
  step = (max(df)-min(df))/seg
  coul <- brewer.pal(5, "Set1")
  colors_border <- coul
  colors_in <- alpha(coul,0.3)
  # png(fl,width = 5, height = 5, res = 300, units = "in")
  pdf(fl,width = 5, height = 5)
  
  if(is.reverse){
    df <- rbind(rep(min(df), ncol(df)), rep(max(df),ncol(df)), df)
    radarchart( df, axistype=1, seg = seg,
                #custom polygon
                pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
                #custom the grid
                cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(max(df),min(df),-step), cglwd=0.8,
                #custom labels
                vlcex=0.8 
    )
  }else{
    df <- rbind(rep(max(df), ncol(df)), rep(min(df),ncol(df)), df)
    radarchart( df, axistype=1, seg = seg,
                #custom polygon
                pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
                #custom the grid
                cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,100,25), cglwd=0.8,
                #custom labels
                vlcex=0.8 
    )
  }
  
  legend(x=0.7, y=1.3, legend = rownames(df[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = colors_border, cex=1, pt.cex=1.5)
  
  dev.off()
}


getRadarPlot(dt,"ra_scores.pdf")


#### 4. arrange multiple images on the same page ####
rl = lapply(c(sprintf("liver_%i.png", 1:4),
              sprintf("stem_%i.png", 1:3),"myoblast_1.png",
              sprintf("monocyte_%i.png", 1:2),sprintf("heart_%i.png", 1:2),"bcell_1.png"), readPNG)
gl = lapply(rl, grid::rasterGrob)

png("summary.png",units="in", width=5, height=5, res=600)
grid.arrange(grobs = gl, ncol = 2, as.table = FALSE)
dev.off()

#### 5. bubble plot of predicted candidate TFs in source and target tissues ####
setwd("../figures/")
tmp <- recipes[["heart_1"]]
Tcell <- as.character(tmp[tmp$algorithms=="Taiji-v4",2:9])

tmp2 <- recipes[["liver_3"]]
Tcell2 <- as.character(tmp2[tmp2$algorithms=="Taiji_v4",2:9])
Tcell <- unique(c(Tcell2, Tcell))
Tcell <- Tcell[c(1:3, 5:7, 4,8, 9:14)]

s <- rownames(group_sorted)[group_sorted$Group=="fibroblast"]
t <- rownames(group_sorted)[group_sorted$Group=="heart"]
t <- rownames(group_sorted)[group_sorted$Group%in%c("liver","heart")]
log_exp <- log2(rna_seq+1)
g <- bubble(Tcell, log_exp[Tcell,c(t,s)], Data_normed[Tcell, c(t,s)], angle = 30)

# png("bb_heart_liver.png",units="in", width=9, height=4.5, res=300)
pdf("bb_heart_liver.pdf", width=9, height=4.5)
grid.draw(g)
dev.off()

#### 6. GO and KEGG pathway analysis ####
setwd("../network/")
L <- list.files(path="./", pattern = "heart.*_TFlist.txt")
davidQuery(fl=L[3], species = "mm",
           key = "heart|regeneration",trailing=50,is.plot = F)
plotGO(fl="KLF15_E.emb.heart_378_TFlist_kegg_go.txt",trailing = 50, w=5, h=8, key = "heart|regeneration")

davidQuery(fl="HNF4G_E.emb.liver_new.txt", species = "mm",
           key = "liver|regeneration",trailing=50,is.plot = F)
plotGO(fl="HNF4G_E.emb.liver_new_kegg_go.txt",trailing = 50, w=5, h=8, key = "liver|regeneration")


#### 290 TFs validatation ####
TF_290 <- read_excel("../../41587_2020_742_MOESM6_ESM.xlsx",skip = 1)
TF_target <- intersect(TF_290$`TF name`,rownames(Data))
TF_bg <- setdiff(rownames(Data),TF_target)
TF_candi <- c("ATOH1","NKX3-1","ETV2","SOX9","SOX6","ZNF536","NKX6.1","NKX6.2")

s = "stem"
t = setdiff(unique(group_sorted$Group),c("stem","P6","P8","P13"))
recipe_size = 1
rank_cut = nrow(Data)
p = 0.001
method="prod"
m = "ratio.abs"

sources = rownames(group_sorted)[group_sorted$Group==s][c(1:3,5)]
# targets = rownames(group_sorted)[group_sorted$Group=="neural.cell"]
targets = "E.emb.cardiac.muscle.cell"


# individual pair
checkCandidates <- function(idx_s,idx_t,candi,m="ratio.abs",check.type="oe"){
  
  name_s <- sources[idx_s]
  name_t <- targets[idx_t]
  message(name_s, " -> ", name_t)
  df = data.frame(diff=Data[,name_t]-Data[,name_s],
                  diff.abs=abs(Data[,name_t]-Data[,name_s]),
                  ratio=Data[,name_t]/Data[,name_s])
  if (sum(df$diff)==0){
    message("pagerank of source and target is exactly the same!")
    writeLines(paste0("pagerank of ",name_s," and ", name_t, "is the same!"),
               paste0("recipe_",name_s, "_to_",name_t,"_",method,"_",m,"_rank_",rank_cut,"_p_",p,"_size_",recipe_size,".txt"))
    return()
  }
  
  df$ratio.abs = unlist(lapply(df$ratio, function(x) if(x<1){1/x}else{x}))
  rownames(df) = rownames(Data)
  df_2 = head(df[order(df[[m]],decreasing = TRUE),], n = rank_cut)
  df_2$type <- ifelse(rna_seq[rownames(df_2),name_t]>rna_seq[rownames(df_2),name_s], "oe", "ko")
  
  # put check.type on top
  # df_2 <- 
  df_2$rank <- rep(1:nrow(df_2))
  message(candi,": ratio.abs=",df_2[candi,"ratio.abs"], "; rank: ", df_2[candi,"rank"])
  
  # get distribution
  dist1 <- df_2[TF_target,"rank"]
  dist2 <- df_2[TF_bg,"rank"]
  
  dt3 <- data.frame("group" = c(rep("target", length(dist1)), rep("background", length(dist2))),
                    "rank" = c(dist1, dist2))
  # return(df)
  # df$groups <- cut(df$Correlation, breaks=c(0,0.25,0.75,1))
  p <- ggplot(dt3, aes(x=group, y=rank, fill = group))+
    # geom_histogram(bins = 50, position = "identity", alpha=.75)+
    geom_boxplot()+
    ggtitle(paste0(name_s, " -> ", name_t))+
    scale_color_brewer(palette = "Set1")+
    scale_fill_brewer(palette = "Set1")+
    theme_light()
  # pdf(paste0("dist_",name_s,"_to_",name_t,".pdf"), width=5, height=5)
  # print(p)
  # dev.off()
  return(p)
  

}
checkCandidates(3,2,"ETV2")
mapply(function(x,y) checkCandidates(idx_s = x, idx_t = y, candi = "ETV2"),
       rep(1:length(sources),each=length(targets)),
       rep(1:length(targets),length(sources)))

# group to group
get_multiple_dist <- function(name_s="stem",name_t="heart",candi="SOX9",m="ratio.abs"){
  
  sources = rownames(group_sorted)[group_sorted$Group==name_s][c(1:3,5)]
  targets = rownames(group_sorted)[group_sorted$Group==name_t][4]
  height = length(targets)/2*10
  
  checkCandidates <- function(idx_s,idx_t,candi,m){
    
    name_s <- sources[idx_s]
    name_t <- targets[idx_t]
    message(name_s, " -> ", name_t)
    df = data.frame(diff=Data[,name_t]-Data[,name_s],
                    diff.abs=abs(Data[,name_t]-Data[,name_s]),
                    ratio=Data[,name_t]/Data[,name_s])
    if (sum(df$diff)==0){
      message("pagerank of source and target is exactly the same!")
      writeLines(paste0("pagerank of ",name_s," and ", name_t, "is the same!"),
                 paste0("recipe_",name_s, "_to_",name_t,"_",method,"_",m,"_rank_",rank_cut,"_p_",p,"_size_",recipe_size,".txt"))
      return()
    }
    
    df$ratio.abs = unlist(lapply(df$ratio, function(x) if(x<1){1/x}else{x}))
    rownames(df) = rownames(Data)
    df_2 = head(df[order(df[[m]],decreasing = TRUE),], n = rank_cut)
    df_2$rank <- seq(1:nrow(df_2))
    message(candi,": ",m,"=",df_2[candi,"ratio.abs"], "; rank: ", df_2[candi,"rank"])
    
    # get distribution
    dist1 <- df_2[TF_target,"rank"]
    dist2 <- df_2[TF_bg,"rank"]
    
    dt3 <- data.frame("group" = c(rep("target", length(dist1)), rep("background", length(dist2))),
                      "rank" = c(dist1, dist2))
    # return(df)
    # df$groups <- cut(df$Correlation, breaks=c(0,0.25,0.75,1))
    p <- ggplot(dt3, aes(x=group, y=rank, fill = group))+
      # geom_histogram(bins = 50, position = "identity", alpha=.75)+
      geom_boxplot()+
      ggtitle(paste0(name_s, " -> ", name_t))+
      scale_color_brewer(palette = "Set1")+
      scale_fill_brewer(palette = "Set1")+
      theme_light()
    # pdf(paste0("dist_",name_s,"_to_",name_t,".pdf"), width=5, height=5)
    # print(p)
    # dev.off()
    return(p)
    
    
  }
  
  ps <- mapply(function(x,y) checkCandidates(idx_s = x, idx_t = y, candi = candi, m=m),
               rep(1:length(sources),each=length(targets)),
               rep(1:length(targets),length(sources)),SIMPLIFY = F)
  p <- ggpubr::ggarrange(plotlist = ps,ncol = 2,nrow = ceiling(length(ps)/2),
                          common.legend = T)
  pdf(paste0("dist_",name_s,"_to_",name_t,"_",m,".pdf"),width = 10,height = height)
  print(p)
  dev.off()
  
}
get_multiple_dist()
get_multiple_dist(name_t = "heart",candi = "ETV2",m = "ratio")

# pick the best among all combinations
pick_best <- function(m){
  
  sources <- rownames(group_sorted)[group_sorted$Group=="stem"][c(1:3,5)]
  targets <- rownames(group_sorted)[group_sorted$Group%in%setdiff(unique(group_sorted$Group),c("stem","P6","P8","P13"))]
  
  best_rank_TF <- function(TF,m=m){
    
    rank_TF <- function(name_s,name_t,TF,m=m){
      
      # message(name_s, " -> ", name_t)
      df = data.frame(diff=Data[,name_t]-Data[,name_s],
                      diff.abs=abs(Data[,name_t]-Data[,name_s]),
                      ratio=Data[,name_t]/Data[,name_s])
      if (sum(df$diff)==0){
        message("pagerank of source and target is exactly the same!")
        return(1000)
      }
      
      df$ratio.abs = unlist(lapply(df$ratio, function(x) if(x<1){1/x}else{x}))
      rownames(df) = rownames(Data)
      df_2 = head(df[order(df[[m]],decreasing = TRUE),], n = rank_cut)
      df_2$rank <- seq(1:nrow(df_2))
      return(df_2[TF,"rank"])
      # message(TF,": ratio.abs=",df_2[TF,m], "; rank: ", df_2[TF,"rank"])
      
    }
    
    o <- mapply(function(x,y) rank_TF(name_s = x, name_t = y, TF = TF, m = m),
                rep(1:length(sources),each=length(targets)),
                rep(1:length(targets),length(sources)))
    return(min(o))
  }
  
  dist1 <- sapply(TF_target, function(x) best_rank_TF(TF = x, m = m))
  dist2 <- sapply(TF_bg, function(x) best_rank_TF(TF = x, m = m))
  dt3 <- data.frame("group" = c(rep("target", length(dist1)), rep("background", length(dist2))),
                    "rank" = c(dist1, dist2))
  
  stat <- t.test(dist1, dist2, "less")
  message("p.value=",stat$p.value)
  
  # box plot
  p1 <- ggplot(dt3, aes(x=group, y=rank, fill = group))+
    # geom_histogram(bins = 50, position = "identity", alpha=.75)+
    geom_boxplot()+
    # ggtitle(paste0(name_s, " -> ", name_t))+
    scale_color_brewer(palette = "Set1")+
    scale_fill_brewer(palette = "Set1")+
    theme_light()
  
  # dist plot
  p2 <- ggplot(dt3, aes(x=rank, fill = group))+
    geom_histogram(bins = 50, position = "identity", alpha=.75)+
    # geom_boxplot()+
    # ggtitle(paste0(name_s, " -> ", name_t))+
    scale_color_brewer(palette = "Set1")+
    scale_fill_brewer(palette = "Set1")+
    theme_light()
  
  # p <- ggpubr::ggarrange(plotlist = c(p1,p2),ncol = 2)
  p <- cowplot::plot_grid(p1,p2,ncol = 2)
  pdf(paste0("dist_pick_best_",m,".pdf"), width=10, height=5)
  print(p)
  dev.off()
}

pick_best("ratio")
pick_best("ratio.abs")

pick_top10 <- function(m){
  
  sources <- rownames(group_sorted)[group_sorted$Group=="stem"][c(1:3,5)]
  targets <- rownames(group_sorted)[group_sorted$Group%in%setdiff(unique(group_sorted$Group),c("stem","P6","P8","P13"))]
    
  top_TFs <- function(name_s,name_t,m=m,cut=10){
    
    # message(name_s, " -> ", name_t)
    df = data.frame(diff=Data[,name_t]-Data[,name_s],
                    diff.abs=abs(Data[,name_t]-Data[,name_s]),
                    ratio=Data[,name_t]/Data[,name_s])
    if (sum(df$diff)==0){
      message("pagerank of source and target is exactly the same!")
      return()
    }
    
    df$ratio.abs = unlist(lapply(df$ratio, function(x) if(x<1){1/x}else{x}))
    rownames(df) = rownames(Data)
    df_2 = head(df[order(df[[m]],decreasing = TRUE),], n = rank_cut)
    df_2$rank <- seq(1:nrow(df_2))
    return(rownames(df_2)[1:cut])
    
  }
    
  o <- mapply(function(x,y) top_TFs(name_s = x, name_t = y, m = m, cut=10),
              rep(1:length(sources),each=length(targets)),
              rep(1:length(targets),length(sources)))
  o <- unique(unlist(o))
  n1 <- length(intersect(o, TF_target))/length(TF_target)
  n2 <- length(intersect(o, TF_bg))/length(TF_bg)
  message("recover ", n1, " of target TFs")
  message("recover ", n2, " of background TFs")
  
  return(min(o))
  
}

#### response: find source-dependent samples ------------------------------------------
getRankFile <- function(source,target,id){
  
  fl <- paste0("recipe_",source, "_to_",target,"_prod_ratio.abs_rank_30_p_0.001_size_3.txt")
  recipe <- read.table(paste0("three/",fl),header = T)
  df <- Data[,c(source,target)]
  names(df) <- c("source","target")
  df$ratio <- df[,"target"]/df[,"source"]
  
  df <- df[order(df$source,decreasing = T),]
  df$source.rank <- c(1:nrow(df))
  df <- df[order(df$target,decreasing = T),]
  df$target.rank <- c(1:nrow(df))
  df <- df[order(df$ratio,decreasing = T),]
  df$ratio.rank <- c(1:nrow(df))
  
  write.table(df,paste0("../benchmark/rank_",id,"_",source,"_to_",target,".txt"),quote = F, sep = "\t")
  message("Success got file for recipe ", id)
  
}

# check heart_1 recipe
source <- "R.emb.skin.fibroblast"
target <- "E.adt.right.atrium.auricular.region"
id <-"heart_v1"
getRankFile(source = source, target = target, id = id)

# check liver_3 recipe
source <- "R.emb.fibroblast.of.skin.of.abdomen"
target <- "E.adt.right.lobe.of.liver"
id <- "liver_v3"
getRankFile(source = source, target = target, id = id)

# check stem_3 recipe
source <- "R.emb.fibroblast.of.skin.of.abdomen"
target <- "R.emb.H1.hESC"
id <- "stem_v3"
getRankFile(source = source, target = target, id = id)
