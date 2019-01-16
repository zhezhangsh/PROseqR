ParseDREG <- function(res, name, path=NA) {
  if (identical(path[1], NA)) out <- name else out <- paste(path, name, sep='/');
  dir.create(out, showWarnings = FALSE);
  
  lst <- readRDS(res);
  
  write.table(lst[[1]][, 1:3], paste0(out, '/', names(lst)[1], '.bed'), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE);
  write.table(lst[[1]], paste0(out, '/', names(lst)[1], '.txt'), sep='\t', row.names = TRUE, col.names = TRUE, quote = FALSE);
  
  write.table(lst[[2]][, 1:3], paste0(out, '/', names(lst)[2], '.bed'), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE);
  write.table(lst[[2]], paste0(out, '/', names(lst)[2], '.txt'), sep='\t', row.names = TRUE, col.names = TRUE, quote = FALSE);
 
  write.table(lst[[4]][, 1:3], paste0(out, '/', names(lst)[4], '.bed'), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE);
  write.table(lst[[4]], paste0(out, '/', names(lst)[4], '.txt'), sep='\t', row.names = TRUE, col.names = TRUE, quote = FALSE);
  
  write.table(lst[[5]][, 1:3], paste0(out, '/', names(lst)[5], '.bed'), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE);
  write.table(lst[[5]], paste0(out, '/', names(lst)[5], '.txt'), sep='\t', row.names = TRUE, col.names = TRUE, quote = FALSE);
  
  writeLines(as.character(lst[[3]]), paste0(out, '/', names(lst)[3], '.txt'));

  writeLines(paste(names(lst[[6]]), unlist(lst[[6]]), sep='='), paste0(out, '/', names(lst)[6], '.txt'));
  
  saveRDS(lst, paste0(out, '/dreg.rds'));
}