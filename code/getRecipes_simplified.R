getRecipes_simplified <- function(name_s, name_t, recipe_size=3,
                       rank_cut=30, p=0.001, 
                       method="prod",m="ratio.abs",
                       g=group_sorted,D=Data, R=rna_seq){
  
  getRecipeType <- function(o,s,t,c="TF"){
    o$exp_S <- rowMeans(R[match(o[[c]], rownames(R)),s,drop=F])
    o$exp_T <- rowMeans(R[match(o[[c]], rownames(R)),t,drop=F])
    with(o, ifelse(exp_T>exp_S, "oe", "ko"))
  }
  
  message(name_s, " -> ", name_t)
  df = data.frame(diff=D[,name_t]-D[,name_s],
                  diff.abs=abs(D[,name_t]-D[,name_s]),
                  ratio=D[,name_t]/D[,name_s])
  if (sum(df$diff)==0){
    message("pagerank of source and target is exactly the same!")
    writeLines(paste0("pagerank of ",name_s," and ", name_t, "is the same!"),
               paste0("recipe_",name_s, "_to_",name_t,"_",method,"_",m,"_rank_",rank_cut,"_p_",p,"_size_",recipe_size,".txt"))
    return()
  }
  
  df$ratio.abs = unlist(lapply(df$ratio, function(x) if(x<1){1/x}else{x}))
  rownames(df) = rownames(D)
  df_2 = head(df[order(df[[m]],decreasing = TRUE),], n = rank_cut)
  
  if (recipe_size == 1){
    df_2$rank <- rep(1:nrow(df_2))
    return(df_2)
  }else{
    tmp <- as.data.frame(t(combn(rownames(df_2),recipe_size)))
    tmp2 <- scale(apply(tmp,1,function(x) prod(df_2[x,"ratio.abs"])))
    tmp$agg <- tmp2
    tmp$count <- 1
    tmp <- tmp[tmp$agg>=qnorm(1-p),]
    if(nrow(tmp)>=1){
      
      output <- as.data.frame(table(unlist(tmp[,1:recipe_size])))
      output <- output[order(output$Freq, decreasing = T),]
      names(output) <- c("TF","freq")
      output$type <- getRecipeType(output, name_s, name_t)
      write.table(output, paste0("recipe_",name_s, "_to_",name_t,"_",method,"_",m,"_rank_",rank_cut,"_p_",p,"_size_",recipe_size,".txt"),
                  row.names = F, quote = F)
      # output <- tmp
      # names(output) <- c("TF1","TF2","TF3","z_score","freq")
      # output <- output[order(output[["z_score"]], decreasing = T),]
      # output$TF1_type <- ifelse(R[output$TF1,name_t]>R[output$TF1,name_s], "oe", "ko")
      # output$TF2_type <- ifelse(R[output$TF2,name_t]>R[output$TF2,name_s], "oe", "ko")
      # output$TF3_type <- ifelse(R[output$TF3,name_t]>R[output$TF3,name_s], "oe", "ko")
      # write.table(output, paste0("realRecipe_",name_s, "_to_",name_t,"_",method,"_",m,"_rank_",rank_cut,"_p_",p,"_size_",recipe_size,".txt"),
      #             row.names = F, quote = F)
      
    }else{
      message("No recipes pass the p-value cutoff")
      writeLines(paste0(name_s," -> ", name_t, ": No recipes pass the p-value cutoff of ",p),
                 paste0("recipe_",name_s, "_to_",name_t,"_",method,"_",m,"_rank_",rank_cut,"_p_",p,"_size_",recipe_size,".txt"))
    }
    # return(tmp)
    }
  
}
