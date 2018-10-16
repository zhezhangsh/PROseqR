RunBwaMemSingle <- function(fq, ref, bwa, samtools, output='', min.length=30, seed.length=19, threads=1, all=FALSE, sam=TRUE) {

  f1 <- sapply(strsplit(fq, '/'), function(x) rev(x)[1]);
  f1 <- sub('.gz$', '', f1);
  f1 <- sub('.fastq$', '', f1);
  f1 <- sub('.fq$', '', f1);
  f1 <- paste(f1, '.sam', sep='');
  if (output[1]!='' & !identical(output[1], NA)) f1 <- paste(output[1], f1, sep='/');
  f2 <- sub('.sam$', '.bam', f1);
  f3 <- sub('.bam$', '.sorted.bam', f2);

  cmd1 <- paste(bwa, 'mem', '-k', seed.length, '-T', min.length, '-t', threads);
  if (all) cmd1 <- paste(cmd1, '-a');
  if (sam) {
    cmd1 <- paste(cmd1, ref, fq, '>', f1);
    cmd2 <- paste(samtools, 'view -Sb', f1, '>', f2);
  } else {
    cmd1 <- paste(cmd1, ref, fq, '|', samtools, 'view -Sb >', f2);
    cmd2 <- '';
  }
  cmd3 <- paste(samtools, 'sort', f2, '-o', f3);
  cmd4 <- paste(samtools, 'index', f3);

  c(cmd1, cmd2, cmd3, cmd4);
}

#######################################################
RunBwaMemPair <- function(fq1, fq2, ref, bwa, samtools, output, min.length=30, seed.length=19, threads=1, all=FALSE) {

  f1 <- paste(output, '.sam', sep='');
  f2 <- sub('.sam$', '.bam', f1);
  f3 <- sub('.bam$', '.sorted.bam', f2);

  cmd1 <- paste(bwa, 'mem', '-k', seed.length, '-T', min.length, '-t', threads);
  if (all) cmd1 <- paste(cmd1, '-a');
  cmd1 <- paste(cmd1, ref, fq1, fq2, '>', f1);

  cmd2 <- paste(samtools, 'view -Sb', f1, '>', f2);
  cmd3 <- paste(samtools, 'sort', f2, '-o', f3);
  cmd4 <- paste(samtools, 'index', f3);

  c(cmd1, cmd2, cmd3, cmd4);
}
