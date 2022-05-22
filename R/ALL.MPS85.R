#' @name MPS85
#' ALL Metabolic Typology Based on Transcriptomic Data
#' @description Based on the built-in training data metabolic types and the GSVA scores of 85 metabolic pathways, this package performs metabolic typing of the samples based on the GSVA results of the 85 metabolic pathways with the transcriptome data. Metabolic types include: MPS1, MPS2 and MPS3.
#' @param expr.data data.frame. Genes expression data.
#' @param gsva.methods character. The method of GSVA to estimate gene-set enrichment scores per sample. See \link{gsva}
#' @param ncores numeric. Number of threads used for GSVA analysis. See \link{gsva}
#' @param scale logical. Whether to scale the genes expression data.
#'
#' @import GSVA
#' @import pamr
#' @return MPS85.result
#' @export

MPS85 <- function(expr.data,
                  gsva.methods = "gsva",
                  ncores = 1,
                  scale = FALSE
){
  #GSVA分析
  message("Step1: Now we will perform GSVA analysis... \n Please wait a minite!")
  #data(genelist)
  gsvaresult<- gsva(expr = as.matrix(expr.data),
                     gset.idx.list = genelist,
                     method=gsva.methods,
                     verbose=T,
                     parallel.sz=ncores
  )

  #scale数据
  if(scale){
    message("Step2: Scale expression data for pamr...")
    names <- dimnames(expr.data)
    expr.data <- t(scale(t(expr.data)))
    dimnames(expr.data) <- names
  }else{
    message("Step2: Don't scale expression data...")
  }

  #pamr分型
  message("Step3: classify the ALL metabolism type...")
  set.seed(10086)
  #data(train.exp)
  #data(train.class)
  x <- list(train.Exp=train.exp,train.MPS=train.class)
  data.pamr <- list(x=x$train.Exp,
                    y=x$train.MPS$cluster,
                    genenames=rownames(x$train.Exp)
  )
  model.train <- pamr.train(data.pamr)
  cv <- pamr.cv(model.train, data=data.pamr)
  thr <- cv$threshold[which.min(cv$error)]
  new.data <- gsvaresult[rownames(train.exp),]
  Pred <- pamr.predict(model.train, newx=new.data,
                       threshold=thr, type="class")
  Pred <- as.data.frame(Pred)
  rownames(Pred) <- colnames(new.data)
  MPS85.result <- list(predict.class = Pred, gsvaresult = gsvaresult)

  message("MPS85 finished")
}

