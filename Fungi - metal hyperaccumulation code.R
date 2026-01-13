# Fungi - metal hyperaccumulation code

library(caper)
library(phytools)

setwd("C:/Users/jamie/OneDrive/Documents/Catherine_fungi")

tree = readNexus("1672taxa_290genes_bb_1 (1).treefile.nex.trees")
data = read.csv("hm_mel_data.csv")
data$HM = as.factor(data$HM)
data$Mel = as.factor(data$Mel)
data$combo = paste(data$HM, data$Mel)
table(data$combo)
head(data)
unique(data$combo)
colnames(data)

# D statistic
tree <- midpoint.root(tree)
comp <- comparative.data(phy = tree,
                         data = data,
                         names.col = "species",
                         vcv = TRUE,
                         warn.dropped = FALSE)

d_HM  <- phylo.d(comp, binvar = HM)
d_Mel <- phylo.d(comp, binvar = Mel)

# Hot-nodes
hot.nodes2 <- function(Sample, Tree, runs = 999, verbose = TRUE, return.null = FALSE) {
  if (!all(Tree$tip.label %in% names(Sample))) {
    stop("Some tip labels in the tree are missing in Sample names")
  }
  Sample <- Sample[Tree$tip.label]
  
  Table <- data.frame(Tree$tip.label, 1:length(Tree$tip.label))
  Tips_frac <- vector("list", length = length(Tree$tip.label) - 1)
  Start <- length(Tree$tip.label) + 1
  End <- Start + length(Tree$tip.label) - 2
  
  for (i in Start:End) {
    tips.i <- getDescendants(Tree, i)
    tips.i <- tips.i[tips.i <= length(Tree$tip.label)]
    Table.i <- Table[Table[, 2] %in% tips.i, ]
    Tips_frac[[i - (Start - 1)]] <- Table.i[, 1]
  }
  
  Trees_store <- vector("list", length = length(Tips_frac))
  for (j in seq_along(Tips_frac)) {
    if (length(Tips_frac[[j]]) > 1) {
      Trees_store[[j]] <- drop.tip(Tree, Tree$tip.label[!Tree$tip.label %in% Tips_frac[[j]]])
    } else {
      Trees_store[[j]] <- NULL
    }
  }
  
  Sample_1 <- Sample[Sample == 1]
  Tree_Sample_1 <- Tree$tip.label[Tree$tip.label %in% names(Sample_1)]
  Tree_Sample_0 <- Tree$tip.label[!Tree$tip.label %in% names(Sample_1)]
  Sample_All <- c(rep(1, length(Tree_Sample_1)), rep(0, length(Tree_Sample_0)))
  names(Sample_All) <- c(Tree_Sample_1, Tree_Sample_0)
  
  SES_Sample <- rep(NA, length(Trees_store))
  Pval <- rep(NA, length(Trees_store))
  Long_Lin <- rep(NA, length(Trees_store))
  Usadas <- rep(NA, length(Trees_store))
  Store_Store_NULL <- vector("list", length = length(SES_Sample))
  
  for (i in 2:length(Trees_store)) {
    Tree.i <- Trees_store[[i]]
    if (is.null(Tree.i)) next
    
    En_uso <- sum(Tree.i$tip.label %in% Tree_Sample_1)
    Usadas[i] <- En_uso
    Store_NULL <- numeric(runs)
    
    for (j in 1:runs) {
      Sample_NULL <- sample(Sample_All)
      names(Sample_NULL) <- names(Sample_All)
      Sample_NULL_1 <- Sample_NULL[Sample_NULL > 0]
      Store_NULL[j] <- sum(Tree.i$tip.label %in% names(Sample_NULL_1))
    }
    
    Store_Store_NULL[[i]] <- Store_NULL
    Long_Lin[i] <- length(Tree.i$tip.label)
    SES_Sample[i] <- (En_uso - mean(Store_NULL)) / sd(Store_NULL)
    Pval[i] <- sum(Store_NULL >= En_uso) / (runs + 1)
    
    if (verbose) {
      message(paste("node", i - 1, "out of", length(Trees_store) - 1, "completed"))
    }
  }
  
  Ntip <- length(Tree$tip.label)
  node_number <- 1:length(Trees_store)
  node_ape <- (Ntip + 1):(Ntip + length(Trees_store))
  
  df <- data.frame(node_number = node_number,
    node_ape = node_ape, Size = Long_Lin, Used = Usadas, SES = SES_Sample, p_value = Pval)
  
  df <- df[complete.cases(df), ]
  df <- df[, c("node_number", "node_ape", "Size", "Used", "SES", "p_value")]
  
  if (return.null) {
    return(list(Nodes = df, Null = Store_Store_NULL))
  } else {
    return(list(Nodes = df))
  }
}

hm_hotnodes = hot.nodes2(setNames(data$HM, data$species), tree, runs = 999, verbose = TRUE, return.null = FALSE)
mel_hotnodes = hot.nodes2(setNames(data$Mel, data$species), tree, runs = 999, verbose = TRUE, return.null = FALSE)
saveRDS(hm_hotnodes, "HM_hotnodes_custom_code.RDS")
saveRDS(mel_hotnodes, "Mel_hotnodes_custom_code.RDS")

hm_hot = hm_hotnodes$Nodes
mel_hot = mel_hotnodes$Nodes
hm_hot = hm_hot[hm_hot$Size <= 100,]
mel_hot = mel_hot[mel_hot$Size <= 100,]

sum(hm_hot$p_value <= 0.05)
sum(mel_hot$p_value <= 0.05)

hm_hot = hm_hot[hm_hot$p_value <= 0.05,]
mel_hot = mel_hot[mel_hot$p_value <= 0.05,]

write.csv(hm_hot, "hm_hot_nodes_significant_100.csv")
write.csv(mel_hot, "mel_hot_nodes_significant_100.csv")