FindVarBins <- function(object,assay.value='RNA',assay.save='integrated',ntop=50){
  GeneAves = AverageExpression(object,slot = 'data',assays = assay.value,group.by = 'project')
  exp = object@assays[[assay.value]]@data
  GeneAve = sort(GeneAves[[assay.value]][,1],decreasing = T)
  GeneAve = GeneAve[GeneAve > quantile(GeneAve,c(0.025,0.975))[1] &
                      GeneAve < quantile(GeneAve,c(0.025,0.975))[2]]
  bins = round(seq(1,length(GeneAve),length.out = 21))
  vargene = c()
  for (i in 1:(length(bins)-1)){
    genes = names(GeneAve[bins[i]:bins[i+1]])
    dispersion = sort(apply(exp[genes,],1,var)/apply(exp[genes,],1,mean),decreasing = T)
    print(dispersion[1:5])
    vargene = c(vargene,names(dispersion[1:ntop]))
  }
  object@assays[[assay.save]]@var.features = vargene
  return(object)
}