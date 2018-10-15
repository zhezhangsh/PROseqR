SampleFastq <- function(fq1, fq2=NA, seqtk, n, seed=1234, output='', prefix='') {
  if (!identical(fq2[1], NA)) fq <- c(fq1[1], fq2[1]) else fq <- fq1[1];
  fnm <- sapply(strsplit(fq, '/'), function(x) rev(x)[1]);
  if (!identical(prefix, NA) & prefix[1]!='')  out <- paste(prefix, 'fnm', sep='_') else out <- fnm;
  if (!identical(output, NA) & output[1]!='') out <- paste(output, fnm, sep='/');
  
  paste(seqtk, 'sample -s', seed, fq, n, '>', out);
}