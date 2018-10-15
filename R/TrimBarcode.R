TrimBarcode <- function(fq, length, seqtk) {
  nm <- sapply(strsplit(fq[1], '/'), function(x) rev(x)[1]);
  f1 <- paste('barcode_', nm, sep='');
  f2 <- paste('trimmed_', nm, sep='');
  
  if (nm != fq) {
    f1 <- sub(paste(nm, '$', sep=''), f1, fq);
    f2 <- sub(paste(nm, '$', sep=''), f2, fq);
  }
  
  c(trimmed=paste(seqtk, 'trimfq -L', length, fq, '| gzip >', f1), 
    barcode=paste(seqtk, 'trimfq -b', length, fq, '| gzip >', f2));
}