RunFlash2 <- function(fq1, fq2, flash2, name, output='', threads=1, 
                      max.overlap=NA, min.overlap=12, allow.outies=TRUE, min.overlap.outie=12, 
                      compress=TRUE, compress.prog='gzip', compress.suffix='gz') {
  
  cmmd <- c(flash2);
  
  if (identical(threads[1], NA)) threads <- 1;
  cmmd <- c(cmmd, paste('--threads=', threads, sep=''));
  
  if (identical(max.overlap[1], NA)) max.overlap <- nchar(readLines(fq1[1], n=2)[2]);
  cmmd <- c(cmmd, paste('--max-overlap=', max.overlap, sep=''));
  
  if (identical(min.overlap[1], NA)) min.overlap <- 10;
  cmmd <- c(cmmd, paste('--min-overlap=', min.overlap, sep=''));  
  
  if (allow.outies) {
    if (identical(min.overlap.outie[1], NA)) min.overlap.outie <- min.overlap;
    cmmd <- c(cmmd, '--allow-outies', paste('--min-overlap-outie=', min.overlap.outie, sep=''));
  };
  
  if (compress) {
    cmmd <- c(cmmd, '--compress', paste('--compress-prog=', compress.prog, sep=''), 
              paste('--output-suffix=', compress.suffix, sep=''));
  };
  
  if (output[1]=='' | identical(output[1], NA)) path <- name else path <- paste(output, name, sep='/');
  if (!dir.exists(path)) dir.create(path);
  cmmd <- c(cmmd, paste('--output-directory=', path, sep=''));
  
  cmmd <- c(cmmd, paste('--output-prefix=', name, sep='')); 
  
  cmmd <- c(cmmd, fq1, fq2);
  
  log  <- paste(path, '/', name, '.log', sep='');
  cmmd <- c(cmmd, paste('2>&1 | tee', log));
  
  paste(cmmd, collapse=' ');
}