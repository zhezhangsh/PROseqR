RemoveBarcodeDup <- function(fin, fout, barcode) {
  # fin1        File with information about aligned reads, as an GRanges object
  # fout        Output file
  # barcode     Qname-barcode mapping as a named character vector or the file name with the vector
  
  require(PROseqR);
  require(GenomicRanges);
  require(GenomicAlignments);
  
  if (file.exists(barcode[1])) barcode <- readRDS(barcode); 
  
  gr <- readRDS(fin);
  od <- gr[rev(order(gr$mapq))];
  id <- paste(as.vector(seqnames(od)), start(od), end(od), barcode[od$qname], sep='_');
  du <- duplicated(id);
  rm <- table(table(id[du]));
  qn <- gr$qname[du];
  re <- gr[!(gr$qname %in% qn)];
  mt <- elementMetadata(re);
  mt <- mt[, !(colnames(mt) %in% 'qname'), drop=FALSE];
  
  elementMetadata(re) <- mt;
  saveRDS(re, fout);
  
  list(count=c(total=length(gr), remain=length(re)), duplicated=rm);
}