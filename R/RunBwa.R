RunBwaMemSingle <- function(fq, ref, bwa, samtools, outpath='', prefix='', min.length=30, seed.length=19, threads=1, all=FALSE, sam=TRUE) {
  if (prefix[1]=='' | identical(prefix[1], NA)) {
    f1 <- paste(f1, '.sam', sep='');
    f1 <- sapply(strsplit(fq, '/'), function(x) rev(x)[1]);
    f1 <- sub('.gz$', '', f1);
    f1 <- sub('.fastq$', '', f1);
    f1 <- sub('.fq$', '', f1);
  } else f1 <- paste(prefix, '.sam', sep='');
  if (outpath[1]!='' & !identical(outpath[1], NA)) f1 <- paste(outpath[1], f1, sep='/');

  f2 <- sub('.sam$', '.bam', f1);
  f3 <- sub('.bam$', '.sorted.bam', f2);
  f4 <- sub('.sorted.bam$', '.sorted.primary.bam', f3);

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
  cmd5 <- paste(samtools, 'view -h -q 1 -F 4 -F 256', f3, '|', samtools, 'view -Sb |', samtools, 'sort >', f4);
  cmd6 <- paste(samtools, 'index', f4);

  c(cmd1, cmd2, cmd3, cmd4, cmd5, cmd6);
}

#######################################################
RunBwaMemPair <- function(fq1, fq2, ref, bwa, samtools, outpath, prefix, min.length=30, seed.length=19, threads=1, all=FALSE, sam=TRUE) {

  f1 <- paste(outpath, '/', prefix, '.sam', sep='');
  f2 <- sub('.sam$', '.bam', f1);
  f3 <- sub('.bam$', '.sorted.bam', f2);
  f4 <- sub('.sorted.bam$', '.sorted.primary.bam', f3);

  cmd1 <- paste(bwa, 'mem', '-k', seed.length, '-T', min.length, '-t', threads);
  if (all) cmd1 <- paste(cmd1, '-a');

  if (sam) {
    cmd1 <- paste(cmd1, ref, fq1, fq2, '>', f1);
    cmd2 <- paste(samtools, 'view -Sb', f1, '>', f2);
  } else {
    cmd1 <- paste(cmd1, ref, fq1, fq2, '|', samtools, 'view -Sb >', f2);
    cmd2 <- '';
  }
  cmd3 <- paste(samtools, 'sort', f2, '-o', f3);
  cmd4 <- paste(samtools, 'index', f3);
  cmd5 <- paste(samtools, 'view -h -q 1 -F 4 -F 256', f3, '|', samtools, 'view -Sb |', samtools, 'sort >', f4);
  cmd6 <- paste(samtools, 'index', f4);

  c(cmd1, cmd2, cmd3, cmd4, cmd5, cmd6);
}
