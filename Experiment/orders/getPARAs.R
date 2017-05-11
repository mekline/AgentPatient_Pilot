library(stringr)
order.files <- list.files(pattern='csv')

for(f in order.files) {
  tmp <- read.csv(paste(f, sep=''), header=T)
  tmp$TR.onset <- as.numeric(as.character(tmp$IntendedOnset))/2
  tmp$Condition <- str_trim(tmp$Condition)
  tmp$Cond.no <- 0
  tmp[as.character(tmp$Condition) == "Agent",]$Cond.no <- 1
  tmp[as.character(tmp$Condition) == "Patient",]$Cond.no <- 2
  tmp <- tmp[tmp$Cond.no != 0,]
  tmp <- tmp[, c('TR.onset', 'Cond.no')]
  fname <- paste(str_split(f,boundary('word'))[[1]][1], '.PARA', sep="")
  cat("#outputs\n\n", file=fname)
  write.table(tmp, fname, sep=" ", row.names = FALSE, col.names = FALSE, append=TRUE)
  cat("\n#names\nagt pat\n\n#durations\n2 2", file=fname, append=TRUE)
}
