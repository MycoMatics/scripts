# a function that takes a fastq file (with asci 33 encoding such as nanopore data)
# the function returns per read (1) the header information and (2) the average q-sccore for that read
# arguments are (1) the fastq file location (2) the number of seqs before progress echo (3) a percentage subset if needed

extract_avg_qscores <- function(fastq_file, print_every = 100000, subset_percent = 100) {
  # Estimate total number of lines in the file using 'wc -l'
  total_lines <- as.integer(system(sprintf("wc -l < '%s'", fastq_file), intern = TRUE))
  estimated_total_reads <- total_lines / 4
    # Adjust estimated total reads based on the subset percentage.
  if (subset_percent < 100) {
    estimated_total_reads <- floor(estimated_total_reads * subset_percent / 100)
  }
  message(sprintf("Estimated total reads to process: %d", estimated_total_reads))
  
  # Preallocate storage based on estimated total reads.
  fastq_ids <- character(estimated_total_reads)
  avg_qscores <- numeric(estimated_total_reads)
  
  counter <- 0 
  
  # Open file connection for reading.
  con <- file(fastq_file, "r")
  while (TRUE) {
    header <- readLines(con, n = 1)
    if (length(header) == 0) break  # End-of-file reached.
    
    sequence <- readLines(con, n = 1)
    separator <- readLines(con, n = 1)
    quality <- readLines(con, n = 1)
    
    counter <- counter + 1
    
    # Remove leading "@" from header to obtain the read ID.
    fastq_ids[counter] <- sub("^@", "", header)
    
    # Convert the quality string from ASCII to Phred scores.
    qs <- as.integer(charToRaw(quality)) - 33
    
    # Convert Q scores to error probabilities: p = 10^(-Q/10)
    p_vals <- 10^(-qs / 10)
    
    # Average the error probabilities over all bases in the read.
    avg_p <- mean(p_vals)
    
    # Convert the average error probability back into a Phred score.
    read_qscore <- -10 * log10(avg_p)
    avg_qscores[counter] <- read_qscore
    
    # Print progress message every 'print_every' reads.
    if (counter %% print_every == 0) {
      message(sprintf("Processed %d reads...", counter))
    }
    
    # Stop processing if the subset target is reached.
    if (counter >= estimated_total_reads) break
  }
  close(con)
  
  message(sprintf("Finished processing %d reads.", counter))
  
  # Trim the preallocated vectors to the actual number of processed reads.
  fastq_ids <- fastq_ids[1:counter]
  avg_qscores <- avg_qscores[1:counter]
  
  return(data.frame(fastq_id = fastq_ids, avg_qscore = avg_qscores))
}
