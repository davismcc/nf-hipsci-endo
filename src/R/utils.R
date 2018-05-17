load_data <- function(dataAll_0){
  genes = unlist(lapply(strsplit(dataAll_0$genes$ids,'[.]'), function(x)x[1]))
  cell_type = dataAll_0$variants$geno_id_map
  Iok = dataAll_0$qualityOK==TRUE & as.numeric(dataAll_0$qual$counts[,2])>100E3
  data_0 = t(dataAll_0$genes$counts)
  row.names(data_0) = genes
  colnames(data_0) = cell_type
  data_0 = data_0[,Iok]
  
  data_E = t(dataAll_0$ERCC$counts)
  row.names(data_E) = dataAll_0$ERCC$ids
  colnames(data_E) = cell_type
  data_E = data_E[,Iok]
  
  return(rbind(data_E, data_0))
}