LoadPairedEnd <- function(fp, region=GRanges(), primary.only=FALSE, min.mapq=1, simple.cigar=FALSE, sam.fields=c()) {
  require(GenomicAlignments);
  require(GenomicRanges);
  
  wht <- sam.fields[sam.fields %in% scanBamWhat()];
  flg <- scanBamFlag(isSecondaryAlignment=!primary.only);
  param <- ScanBamParam(flag=flg, simpleCigar = simple.cigar, mapqFilter=max(0, min.mapq), 
                        what=wht, which=region);
  
  ga  <- readGAlignmentPairs(fp, param = param);
  fst <- ga@first;
  lst <- ga@last;
  
  s1 <- as.vector(seqnames(fst));
  s2 <- as.vector(seqnames(lst));
  
  ind <- which(s1==s2);
  fst <- fst[ind];
  lst <- lst[ind];
  
  gr <- lapply(list(first=fst, last=lst), function(ga) {
    gr <- GRanges(as.vector(seqnames(ga)), IRanges(start(ga), end(ga)), strand=strand(ga));
    elementMetadata(gr) <- elementMetadata(ga)[, wht];
    if (length(gr) > 0) names(gr) <- 1:length(gr);
    gr; 
  }); 
}; 

# Load reads mapped to individual chromosomes with preset parameters
LoadPairedEndWrapper <- function(fin, fout, prefix, region, primary.only=TRUE, min.mapq=20, simple.cigar=TRUE, 
                                 sam.fields=c('mapq', 'qname'), max.width=2500, verbose=FALSE) {
  
  # fin     Input file of sorted bam file
  # fout    Directory of output file
  # prefix  Prefix of output file names
  
  if (!dir.exists(fout)) dir.create(fout, showWarnings = FALSE);
  if (identical(names(region), NULL)) names(region) <- 1:length(region);
  
  stat <- lapply(1:length(region), function(i) {
    if (verbose) cat('Loading', names(region)[i], '\n');
    
    loc <- region[i];
    gr  <- LoadPairedEnd(fin, region=loc, primary.only=primary.only, min.mapq=min.mapq, simple.cigar=TRUE, sam.fields=sam.fields);
    
    fst <- gr$first;
    lst <- gr$last;
    
    # stt1 <- start(fst);
    # stt2 <- start(lst);
    # end1 <- end(fst);
    # end2 <- end(lst);
    # strd <- as.vector(strand(fst));
    # 
    stt <- pmin(start(fst), start(lst));
    end <- pmax(end(fst), end(lst));
    gr0 <- GRanges(seqnames(fst), IRanges(stt, end), strand=strand(fst), mapq1=fst$mapq, mapq2=lst$mapq, qname=fst$qname);
    gr0 <- gr0[width(gr0)<=max.width];
    
    ttl <- c(number=length(gr0), mean_mapq1=mean(gr0$mapq1), mean_mapq1=mean(gr0$mapq2), mean_length=mean(width(gr0)));
    out <- list(overall=ttl, mapq=table(gr0$mapq1+gr0$mapq2), width=table(width(gr0)));
    
    fn1 <- paste(fout, '/', prefix, '_', names(region)[i], '.rds', sep='');
    fn2 <- paste(fout, '/', prefix, '_', names(region)[i], '_summary.rds', sep='');
    
    saveRDS(gr, fn1);
    saveRDS(out, fn2);
    
    out;
  }); 
  names(stat) <- names(region);
  
  stat;
}

LoadSingleEnd <- function(fs, region=GRanges(), primary.only=FALSE, min.mapq=1, simple.cigar=FALSE, sam.fields=c()) {
  require(GenomicAlignments);
  require(GenomicRanges);
  
  wht <- sam.fields[sam.fields %in% scanBamWhat()];
  flg <- scanBamFlag(isSecondaryAlignment=!primary.only);
  param <- ScanBamParam(flag=flg, simpleCigar = simple.cigar, mapqFilter=max(0, min.mapq), 
                        what=wht, which=region);
  
  ga <- readGAlignments(fs, param = param);
  gr <- GRanges(as.vector(seqnames(ga)), IRanges(start(ga), end(ga)), strand=strand(ga));
  elementMetadata(gr) <- elementMetadata(ga)[, wht];
  if (length(gr) > 0) names(gr) <- 1:length(gr);
  
  gr; 
};

# Load reads mapped to individual chromosomes with preset parameters
LoadSingleEndWrapper <- function(fin, fout, prefix, region, primary.only=TRUE, min.mapq=20, simple.cigar=TRUE, 
                                 sam.fields=c('mapq', 'qname'), verbose=FALSE) {
  
  # fin     Input file of sorted bam file
  # fout    Directory of output file
  # prefix  Prefix of output file names
  
  if (!dir.exists(fout)) dir.create(fout, showWarnings = FALSE);
  if (identical(names(region), NULL)) names(region) <- 1:length(region);
  
  stat <- lapply(1:length(region), function(i) {
    if (verbose) cat('Loading', names(region)[i], '\n');
    
    loc <- region[i];
    gr  <- LoadSingleEnd(fin, region=loc, primary.only=primary.only, min.mapq=min.mapq, simple.cigar=TRUE, sam.fields=sam.fields);
    ttl <- c(number=length(gr), mean_mapq=mean(gr$mapq), mean_length=mean(width(gr)));
    out <- list(overall=ttl, mapq=table(gr$mapq), width=table(width(gr)));
    
    fn1 <- paste(fout, '/', prefix, '_', names(region)[i], '.rds', sep='');
    fn2 <- paste(fout, '/', prefix, '_', names(region)[i], '_summary.rds', sep='');
    
    saveRDS(gr, fn1);
    saveRDS(out, fn2);
    
    out;
  }); 
  names(stat) <- names(region);
  
  stat;
}