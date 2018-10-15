ParseFlash2 <- function(log=NA, hist=NA) {
  if (!file.exists(log[1])) mp1 <- NA else {
    l1 <- readLines(log[1]);
    l1 <- sub('^\\[FLASH\\] ', '', l1);
    l1 <- l1[(which(l1=='Read combination statistics:')+1):length(l1)];
    l1 <- gsub('\\s+', ' ', l1);
    l1 <- sub('^ ', '', l1);
    l1 <- l1[1:(which(l1=='')[1]-1)];
    l1 <- l1[grep('pairs:', l1)];
    m1 <- strsplit(l1, ' ');
    t1 <- sapply(m1, function(x) x[1]);
    n1 <- as.integer(sapply(m1, function(x) x[3]));
    names(n1) <- t1;
    mp1 <- n1;
  };
  
  if (!file.exists(hist[1])) mp1 <- NA else {
    l2 <- readLines(hist);
    m2 <- strsplit(l2, '\t');
    t2 <- as.integer(sapply(m2, function(x) x[1]));
    n2 <- as.integer(sapply(m2, function(x) x[2]));
    names(n2) <- t2;
    mp2 <- n2;
  };
  
  list(log=mp1, hist=mp2);
}