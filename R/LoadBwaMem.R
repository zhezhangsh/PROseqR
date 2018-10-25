LoadPairedEnd <- function(fp, region=GRanges(), primary.only=FALSE, min.mapq=1, simple.cigar=FALSE, sam.fields=c()) {
  require(GenomicAlignments);
  require(GenomicRanges);

  flg <- scanBamFlag(isSecondaryAlignment=!primary.only);
  wht <- sam.fields[sam.fields %in% scanBamWhat()];
  if (length(wht)==0) wht <- character(0);
  param <- ScanBamParam(flag=flg, simpleCigar = simple.cigar, mapqFilter=max(0, min.mapq), what=wht, which=region);

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
LoadPairedEndWrapper <- function(fin, fout, prefix, region, min.mapq=20, simple.cigar=TRUE, max.width=2500, verbose=FALSE) {

  # fin     Input file of sorted bam file
  # fout    Directory of output file
  # prefix  Prefix of output file names

  if (!dir.exists(fout)) dir.create(fout, showWarnings = FALSE);
  if (identical(names(region), NULL)) names(region) <- 1:length(region);
  sam.fields <- c('flag', 'mapq', 'qname');
  primary.only=TRUE;

  stat <- lapply(1:length(region), function(i) {
    if (verbose) cat('Loading', names(region)[i], '\n');

    loc <- region[i];

    flg1   <- scanBamFlag(isSecondaryAlignment=!primary.only, isFirstMateRead=TRUE);
    param1 <- ScanBamParam(flag=flg1, simpleCigar = simple.cigar, mapqFilter=max(0, min.mapq), what=sam.fields, which=loc);

    flg2   <- scanBamFlag(isSecondaryAlignment=!primary.only, isSecondMateRead=TRUE);
    param2 <- ScanBamParam(flag=flg2, simpleCigar = simple.cigar, mapqFilter=max(0, min.mapq), what=sam.fields, which=loc);

    gr1 <- readGAlignments(fin, param = param1);
    gr2 <- readGAlignments(fin, param = param2);

    gr1 <- gr1[elementMetadata(gr1)$qname %in% elementMetadata(gr2)$qname];
    gr2 <- gr2[elementMetadata(gr2)$qname %in% elementMetadata(gr1)$qname];
    gr1 <- gr1[order(elementMetadata(gr1)$qname)];
    gr2 <- gr2[order(elementMetadata(gr2)$qname)];

    stt <- pmin(start(gr1), start(gr2));
    end <- pmax(end(gr1), end(gr2));
    gr0 <- GRanges(as.vector(seqnames(gr1)), IRanges(stt, end), strand=as.vector(strand(gr2)),
                   mapq1=elementMetadata(gr1)$mapq, mapq2=elementMetadata(gr2)$mapq, qname=elementMetadata(gr1)$qname);
    gr0 <- gr0[as.vector(seqnames(gr1))==as.vector(seqnames(gr2))];
    gr0 <- gr0[width(gr0)<=max.width];
    names(gr0) <- 1:length(gr0);

    mpq <- elementMetadata(gr0)[, c('mapq1', 'mapq2')];
    ttl <- c(number=length(gr0), mean_mapq=mean(mpq[,1]+mpq[,2]), mean_length=mean(width(gr0)));
    out <- list(overall=ttl, mapq=table(mpq[,1]+mpq[,2]), width=table(width(gr0)));

    fn1 <- paste(fout, '/', prefix, '_', names(region)[i], '.rds', sep='');
    fn2 <- paste(fout, '/', prefix, '_', names(region)[i], '_summary.rds', sep='');
    fn3 <- paste(fout, '/', prefix, '_', names(region)[i], '_dumped.rds', sep='');

    saveRDS(gr0, fn1);
    saveRDS(out, fn2);
    saveRDS(dmp, fn3);

    out;
  });
  names(stat) <- names(region);

  stat;
}

LoadSingleEnd <- function(fs, region=GRanges(), primary.only=FALSE, min.mapq=1, simple.cigar=FALSE, sam.fields=c()) {
  require(GenomicAlignments);
  require(GenomicRanges);

  flg <- scanBamFlag(isSecondaryAlignment=!primary.only);
  wht <- sam.fields[sam.fields %in% scanBamWhat()];
  if (length(wht)==0) wht <- character(0);
  param <- ScanBamParam(flag=flg, simpleCigar = simple.cigar, mapqFilter=max(0, min.mapq), what=wht, which=region);

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
    str <- c('-', '+')[as.integer(strand(gr))];
    strand(gr) <- str;

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
