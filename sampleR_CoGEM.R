# CoGEM analysis: communicability-guided random forest and module testing.
# Set `it` to a run index before execution; it is used in random seeds.

##################################################################################################
##################### 1. Parameters and input data
##################################################################################################
library(OmnipathR)
library(readr)
library(dplyr)

ntree <- 199              # Number of trees in each forest.
B <- 199                   # Number of shadow replicates for gene-level p-values.
Omega <- 199               # Number of phenotype permutations for the module p-value.
mtry <- 10                 # Number of candidate genes considered at each node.
kC <- 50                  # Maximum number of communicability-linked candidates per gene.
max_depth <- 8;            # Maximum tree depth.
min_leaf <- 5;             # Minimum samples allowed in each child node.
alpha <- 0.01              # Gene-selection cutoff applied to p_gene.

# Retrieve BioPlex protein-protein interactions via OmniPathR.
PPI0 <- OmnipathR::bioplex_all(unique=TRUE)

# Keep unique gene pairs.
PPI <- unique(data.frame(
  gene1=PPI0$SymbolA,
  gene2=PPI0$SymbolB
))

# Inputs: sample-by-gene expression matrix and matching binary phenotype labels.
X<-read.table("EXPRESSION.csv",sep=",")
y_true<-read.table("Phenotype.csv",sep=",")[,1]

Mgenes <- colnames(X)
p <- ncol(X); 
mtry <- min(mtry, p); 
kC <- min(kC, p - 1)

##################################################################################################
##################### 2. Network and communicability
##################################################################################################
# Restrict the interaction network to genes in the expression matrix.
PPI_M <- PPI[
  complete.cases(PPI) &
  PPI$gene1 %in% Mgenes &
  PPI$gene2 %in% Mgenes &
  PPI$gene1 != PPI$gene2,
]

PPI_M <- unique(PPI_M[,c("gene1","gene2")])

# Build an undirected adjacency matrix.
A <- matrix(0,p,p,dimnames=list(Mgenes,Mgenes))

if (nrow(PPI_M) > 0) {
  i <- match(PPI_M$gene1,Mgenes)
  j <- match(PPI_M$gene2,Mgenes)
  A[cbind(i,j)] <- 1
  A[cbind(j,i)] <- 1
}

diag(A) <- 0

# Degree matrix.
degree <- rowSums(A)
D <- diag(degree)

# Degree-normalized adjacency matrix.
D_inv_sqrt <- diag(ifelse(degree > 0, 1/sqrt(degree), 0))
A_tilde <- D_inv_sqrt %*% A %*% D_inv_sqrt

# Compute the matrix exponential.
eig <- eigen(A_tilde, symmetric=TRUE)
K <- sweep(eig$vectors, 2, exp(eig$values), "*") %*% t(eig$vectors)

# Remove self-communicability.
C <- K
diag(C) <- 0
C[C < 0 & abs(C) < 1e-12] <- 0
dimnames(C) <- list(Mgenes, Mgenes)


# Construct feature-sampling probabilities.
P <- matrix(0, p, p, dimnames=list(Mgenes, Mgenes))

for (i in seq_len(p)) {

  available <- setdiff(seq_len(p), i)

  candidate <- available[
    order(C[i,available], decreasing=TRUE)
  ]

  candidate <- head(candidate, kC)
  denominator <- sum(C[i,candidate])

  if (denominator > 0) {
    P[i,candidate] <- C[i,candidate]/denominator
  }
}


##################################################################################################
##################### 3. Communicability-guided random forest
##################################################################################################
y<-y_true
best_split <- function(j,id) {
  x <- X[id,j]; yy <- y[id]; o <- order(x)
  x <- x[o]; yy <- yy[o]; n <- length(yy)

  cut_id <- which(diff(x) > 0)
  cut_id <- cut_id[cut_id >= min_leaf & n-cut_id >= min_leaf]
  if (length(cut_id) == 0) return(NULL)

  left_n <- cut_id
  right_n <- n-cut_id
  left_1 <- cumsum(yy)[cut_id]
  right_1 <- sum(yy)-left_1

  parent_gini <- 2*mean(yy)*(1-mean(yy))
  left_gini <- 2*(left_1/left_n)*(1-left_1/left_n)
  right_gini <- 2*(right_1/right_n)*(1-right_1/right_n)

  gain <- parent_gini-
    (left_n/n)*left_gini-
    (right_n/n)*right_gini

  k <- cut_id[which.max(gain)]
  list(cut=(x[k]+x[k+1])/2,gain=max(gain))
}

grow_tree <- function(id,parent=NA,used=integer(),depth=0) {
  node_prob <- mean(y[id])
  available <- setdiff(seq_len(p),used)

  if (depth >= max_depth ||
      length(id) < 2*min_leaf ||
      node_prob %in% c(0,1) ||
      length(available) == 0) {
    return(list(leaf=TRUE,prob=node_prob))
  }

  if (is.na(parent)) {
    size <- min(mtry,length(available))
    candidate <- available[
      sample.int(length(available),size)
    ]
  } else {
    weight <- P[parent,available]
    eligible <- available[weight > 0]
    weight <- weight[weight > 0]

    if (length(eligible) == 0) {
      return(list(leaf=TRUE,prob=node_prob))
    }

    size <- min(mtry,length(eligible))
    candidate <- eligible[
      sample.int(
        length(eligible),
        size,
        replace=FALSE,
        prob=weight
      )
    ]
  }

  result <- lapply(
    candidate,
    function(j) best_split(j,id)
  )

  gain <- sapply(
    result,
    function(z) if (is.null(z)) -Inf else z$gain
  )

  if (all(!is.finite(gain)) || max(gain) <= 0) {
    return(list(leaf=TRUE,prob=node_prob))
  }

  q <- which.max(gain)
  split_gene <- candidate[q]
  split_cut <- result[[q]]$cut

  left_id <- id[X[id,split_gene] <= split_cut]
  right_id <- id[X[id,split_gene] > split_cut]

  list(
    leaf=FALSE,
    split=split_gene,
    gene=Mgenes[split_gene],
    cut=split_cut,
    left=grow_tree(
      left_id,
      parent=split_gene,
      used=c(used,split_gene),
      depth=depth+1
    ),
    right=grow_tree(
      right_id,
      parent=split_gene,
      used=c(used,split_gene),
      depth=depth+1
    )
  )
}

predict_tree <- function(tree,Xnew) {
  pred <- numeric(nrow(Xnew))

  walk <- function(node,id) {
    if (length(id) == 0) return(NULL)

    if (node$leaf) {
      pred[id] <<- node$prob
    } else {
      left <- Xnew[id,node$split] <= node$cut
      walk(node$left,id[left])
      walk(node$right,id[!left])
    }
  }

  walk(tree,seq_len(nrow(Xnew)))
  pred
}

# Fit the forest using bootstrap samples.
n <- nrow(X)
trees <- vector("list",ntree)
OOB <- vector("list",ntree)

for (t in seq_len(ntree)) {
  bootstrap_id <- sample(
    seq_len(n),
    size=n,
    replace=TRUE
  )

  OOB[[t]] <- setdiff(
    seq_len(n),
    unique(bootstrap_id)
  )

  trees[[t]] <- grow_tree(
    id=bootstrap_id,
    parent=NA,
    used=integer(),
    depth=0
  )

  if (t %% 10 == 0) {
    message("Tree ",t,"/",ntree)
  }
}

CoGEM_forest <- list(
  trees=trees,
  OOB=OOB,
  genes=Mgenes
)

# Compute full-forest predictions.
tree_prob <- sapply(
  CoGEM_forest$trees,
  function(tree) predict_tree(tree,X)
)

forest_prob <- rowMeans(tree_prob)
forest_class <- ifelse(forest_prob >= 0.5,1,0)

# Compute out-of-bag (OOB) predictions.
OOB_sum <- numeric(n)
OOB_count <- numeric(n)

for (t in seq_len(ntree)) {
  id <- CoGEM_forest$OOB[[t]]

  if (length(id) > 0) {
    OOB_sum[id] <- OOB_sum[id]+
      predict_tree(
        CoGEM_forest$trees[[t]],
        X[id,,drop=FALSE]
      )

    OOB_count[id] <- OOB_count[id]+1
  }
}

OOB_prob <- OOB_sum/OOB_count
OOB_prob[OOB_count == 0] <- NA
OOB_class <- ifelse(OOB_prob >= 0.5,1,0)

# Summarize OOB classification and prediction loss.
print(
  table(
    Observed=y,
    Predicted=OOB_class
  )
)

OOB_Brier <- mean(
  (y-OOB_prob)^2,
  na.rm=TRUE
)

cat("Number of trees:",length(CoGEM_forest$trees),"\n")
cat("OOB Brier loss:",OOB_Brier,"\n")


##################################################################################################
##################### 4. Observed gene importance
##################################################################################################
# Return genes used for splitting in a tree.
tree_genes <- function(node) {
  if (node$leaf) return(integer())

  unique(
    c(
      node$split,
      tree_genes(node$left),
      tree_genes(node$right)
    )
  )
}


I_obs <- numeric(p)
names(I_obs) <- Mgenes

for (t in seq_len(ntree)) {
  id <- CoGEM_forest$OOB[[t]]

  if (length(id) < 2) next

  tree <- CoGEM_forest$trees[[t]]
  X_OOB <- X[id,,drop=FALSE]
  y_OOB <- y[id]

  # OOB loss before permutation.
  original_prob <- predict_tree(tree,X_OOB)
  original_loss <- mean((y_OOB-original_prob)^2)

  # Evaluate genes used for splitting in this tree.
  used_in_tree <- tree_genes(tree)

  for (j in used_in_tree) {
    X_perm <- X_OOB

    # Permute gene j among OOB samples.
    X_perm[,j] <- sample(X_perm[,j])

    permuted_prob <- predict_tree(tree,X_perm)
    permuted_loss <- mean((y_OOB-permuted_prob)^2)

    I_obs[j] <- I_obs[j]+
      (permuted_loss-original_loss)
  }
}

# Average importance across trees.
I_obs <- I_obs/ntree

gene_importance <- data.frame(
  gene=Mgenes,
  importance=I_obs
)

gene_importance <- gene_importance[
  order(gene_importance$importance,decreasing=TRUE),
]

print(head(gene_importance,20))


##################################################################################################
##################### 5. Shadow-calibrated gene identification
##################################################################################################
X_original <- X

I_null <- matrix(
  0,
  nrow=p,
  ncol=B,
  dimnames=list(Mgenes,paste0("Shadow_",seq_len(B)))
)

for (b in seq_len(B)) {

  set.seed(1000+b+it)

  # Permute each gene independently to create a shadow matrix.
  X_shadow <- sapply(
    seq_len(p),
    function(j) {
      X_original[sample.int(nrow(X_original)),j]
    }
  )

  X_shadow <- as.matrix(X_shadow)
  colnames(X_shadow) <- Mgenes

  # The tree functions use the global expression matrix X.
  # Temporarily set X to the shadow matrix.
  X <- X_shadow

  shadow_trees <- vector("list",ntree)
  shadow_OOB <- vector("list",ntree)

  # Fit a forest to the shadow matrix.
  for (t in seq_len(ntree)) {

    bootstrap_id <- sample(
      seq_len(n),
      size=n,
      replace=TRUE
    )

    shadow_OOB[[t]] <- setdiff(
      seq_len(n),
      unique(bootstrap_id)
    )

    shadow_trees[[t]] <- grow_tree(
      id=bootstrap_id,
      parent=NA,
      used=integer(),
      depth=0
    )
  }

  # Compute gene importance for shadow replicate b.
  I_b <- numeric(p)

  for (t in seq_len(ntree)) {

    id <- shadow_OOB[[t]]

    if (length(id) < 2) next

    tree <- shadow_trees[[t]]
    X_OOB <- X_shadow[id,,drop=FALSE]
    y_OOB <- y[id]

    # OOB loss before permutation.
    original_prob <- predict_tree(tree,X_OOB)
    original_loss <- mean(
      (y_OOB-original_prob)^2
    )

    used_in_tree <- tree_genes(tree)

    for (j in used_in_tree) {

      X_perm <- X_OOB

      # Permute shadow gene j among OOB samples.
      X_perm[,j] <- sample(X_perm[,j])

      permuted_prob <- predict_tree(
        tree,
        X_perm
      )

      permuted_loss <- mean(
        (y_OOB-permuted_prob)^2
      )

      I_b[j] <- I_b[j]+
        (permuted_loss-original_loss)
    }
  }

  # Average importance across trees.
  I_null[,b] <- I_b/ntree

  if (b %% 10 == 0 || b == B) {
    message("Shadow replicate ",b,"/",B)
  }
}

# Restore the original expression matrix.
X <- X_original

# Compute gene-matched empirical p-values.
p_gene <- (
  1+rowSums(I_null >= I_obs)
)/(B+1)


##################################################################################################
##################### 6. Module-level phenotype association test
##################################################################################################
# Define the selected module.
candidate_genes <- Mgenes[p_gene < alpha]
module_id <- match(candidate_genes,Mgenes)
module_id <- module_id[!is.na(module_id)]

# Compute forest OOB predictions.
get_OOB_prob <- function(forest,Xnew) {
  n <- nrow(Xnew)
  pred_sum <- numeric(n)
  pred_n <- numeric(n)

  for (t in seq_along(forest$trees)) {
    id <- forest$OOB[[t]]
    if (length(id) == 0) next

    pred_sum[id] <- pred_sum[id]+
      predict_tree(
        forest$trees[[t]],
        Xnew[id,,drop=FALSE]
      )

    pred_n[id] <- pred_n[id]+1
  }

  prob <- pred_sum/pred_n
  prob[pred_n == 0] <- NA
  prob
}

# Compute the OOB Brier-loss increase after joint module permutation.
module_stat <- function(
  forest,X,y,module_id,
  H=20,seed=300+it
) {
  if (length(module_id) == 0) return(0)

  original_prob <- get_OOB_prob(forest,X)
  original_loss <- mean(
    (y-original_prob)^2,
    na.rm=TRUE
  )

  permuted_loss <- numeric(H)
  set.seed(seed)

  for (h in seq_len(H)) {
    X_perm <- X
    row_id <- sample.int(nrow(X))

    X_perm[,module_id] <-
      X[row_id,module_id,drop=FALSE]

    permuted_prob <- get_OOB_prob(
      forest,
      X_perm
    )

    permuted_loss[h] <- mean(
      (y-permuted_prob)^2,
      na.rm=TRUE
    )
  }

  mean(permuted_loss)-original_loss
}

# Fit a forest using the current global phenotype vector y.
make_forest <- function(seed) {
  set.seed(seed)

  trees <- vector("list",ntree)
  OOB <- vector("list",ntree)

  for (t in seq_len(ntree)) {
    bootstrap_id <- sample.int(
      n,n,replace=TRUE
    )

    OOB[[t]] <- setdiff(
      seq_len(n),
      unique(bootstrap_id)
    )

    trees[[t]] <- grow_tree(
      bootstrap_id,
      parent=NA,
      used=integer(),
      depth=0
    )
  }

  list(
    trees=trees,
    OOB=OOB,
    genes=Mgenes
  )
}

# 6.1 Observed module statistic

# Set the global phenotype vector to the observed phenotype.
y <- y_true

# Fit the evaluation forest on the observed phenotype.
evaluation_forest <- make_forest(
  seed=9000+it
)

# Compute the observed module statistic.
Q_obs <- module_stat(
  evaluation_forest,
  X,
  y_true,
  module_id,
  H=20,
  seed=300+it
)

# 6.2 Phenotype-permutation null distribution
# Keep the selected module fixed.

Q_null <- numeric(Omega)

for (r in seq_len(Omega)) {
  set.seed(5000+r+it)

  # Permute the observed phenotype.
  y <- sample(y_true)

  # Refit the forest using the permuted phenotype.
  forest_null <- make_forest(
    seed=10000+r+it
  )

  Q_null[r] <- module_stat(
    forest_null,
    X,
    y,
    module_id,
    H=20,
    seed=20000+r+it
  )

  if (r %% 10 == 0 || r == Omega) {
    message(
      "Phenotype permutation ",
      r,"/",Omega
    )
  }
}

# Restore the observed phenotype.
y <- y_true

# Compute the empirical module p-value.
p_module <- (
  1+sum(Q_null >= Q_obs)
)/(Omega+1)

cat("Module size:",length(candidate_genes),"\n")
cat("Observed module statistic:",Q_obs,"\n")
cat("Module p-value:",p_module,"\n")
