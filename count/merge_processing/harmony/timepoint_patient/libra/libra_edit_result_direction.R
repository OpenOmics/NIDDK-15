library(stringr)
comparisons = c("2-3", "before-control", "after-control", "without-control", "after-before", "after-without", "before-without", "noflare-control", "flare-control", "flare-noflare")

for (filename in Sys.glob('*csv')) {
  #Get the comparison and alphabetize it (to match the default comparison order)
  comparison <- str_sort(strsplit(strsplit(grep('label', strsplit(filename, '_')[[1]], value=TRUE), '-')[[1]][[2]], 'vs')[[1]], numeric=TRUE)
  separated <- strsplit(filename, '_')[[1]]
  degs <- read.csv(filename, row.names=1)

  print(filename)
  if (paste(comparison, collapse='-') %in% comparisons) {
    separated[[2]] <- paste0('label-', paste(comparison, collapse='vs'))
    new_name <- paste(separated, collapse='_')  
    print("No Change")
    print(new_name)
    write.csv(degs, file.path('direction_adjusted', new_name))
  } else if (paste(rev(comparison), collapse='-') %in% comparisons) {
    separated[[2]] <- paste0('label-', paste(rev(comparison), collapse='vs'))
    new_name <- paste(separated, collapse='_')  
    print("Change")
    print(new_name)
    degs$avg_logFC <- -degs$avg_logFC
    write.csv(degs, file.path('direction_adjusted', new_name))
  } else{
    print("Comparison Not Found")
  }
  print("")
}
