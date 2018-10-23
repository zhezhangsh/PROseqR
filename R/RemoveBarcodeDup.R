RemoveBarcodeDup <- function(fin, fout, barcode) {
  # fin1        File with information about aligned reads, as an GRanges object
  # fout        Output file
  # barcode     Qname-barcode mapping as a named character vector or the file name with the vector

  require(PROseqR);
  require(GenomicRanges);
  require(GenomicAlignments);

  if (file.exists(barcode[1])) barcode <- readRDS(barcode);

  dup <- readRDS(fin);
  loc <- paste(as.vector(seqnames(dup)), start(dup), end(dup), sep='_');
  rep <- loc[duplicated(loc)];

  du0 <- dup[loc %in% rep];
  du0 <- du0[rev(order(du0$mapq))];
  qn0 <- du0$qname;
  bc0 <- barcode[names(barcode) %in% qn0];
  bc0 <- as.vector(bc0[qn0]);
  lc0 <- paste(as.vector(seqnames(du0)), start(du0), end(du0), bc0, sep='_');
  ind <- which(duplicated(lc0));
  qn0 <- du0[ind]$qname;

  lc1 <- lc0[ind];
  ct1 <- table(table(lc1));
  du1 <- dup[!(dup$qname %in% qn0)];
  cl1 <- which(colnames(elementMetadata(du1))=='qname');
  elementMetadata(du1) <- elementMetadata(du1)[, -cl1, drop=FALSE];
  saveRDS(du1, fout);

  list(count=c(total=length(dup), remain=length(du1)), duplicated=ct1);
}
