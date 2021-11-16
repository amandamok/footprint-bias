##################################################
### filter sam alignment file:
### report positionally-first read of alphabeetically-first transcript;
### assumes input sam alignment file is sorted by read name

library(optparse)

option_list <- list(make_option(c("-g", "--genome"), type="character", default="human_CoV",
                                help="bowtie index prefix", metavar="character"))
option_list <- list(make_option(c("-i", "--input"), type="character",
                                help="input sam alignment file, sorted by read name", metavar="character"),
                    make_option(c("-o", "--output"), type="character",
                                help="output filtered sam alignment file", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if(any(sapply(opt, is.null))) {
  cat("ERROR: missing input and/or output file")
  q(save="no")
}

process_alignment <- function(input_fname, output_fname) {
  # choose positionally-first alignment per read
  ## input_fname: character; file path to sam alignment file sorted by read name
  ## output_fname: character; file path to output sam alignment file
  input_con <- file(input_fname, "r")
  output_con <- file(output_fname, "w")
  first_read <- F
  while(TRUE) {
    line <- readLines(input_con, n=1)
    if(length(line) == 0) { # end of file
      # process last read group
      writeLines(choose_read(read_group), output_con)
      break
    } else {
      if(grepl("^@", line)) { # header line
        writeLines(line, output_con)
        next
      } else {
        if(!first_read) { # first alignment line
          last_read <- strsplit(line, split="\t")[[1]][1]
          read_group <- line
          first_read <- T
          next
        } else { # other alignment line
          this_read <- strsplit(line, split="\t")[[1]][1]
          if(this_read == last_read) { # same read, alternative alignment
            read_group <- c(read_group, line)
          } else { # new read
            # write positionally-first read of alphabetically-first transcript
            writeLines(choose_read(read_group), output_con)
            # start new read group
            last_read <- this_read
            read_group <- line
            next
          }
        }
      }
    }
  }
  close(input_con)
  close(output_con)
}

choose_read <- function(input_reads) {
  # choose positionally-first read of alphabetically-first transcript
  ## input_reads: character vector; from readLines() of sam alignment file
  input_reads_transcript <- sapply(input_reads, function(x) { strsplit(x, split="\t")[[1]][3] })
  first_transcript <- min(input_reads_transcript)
  if(sum(input_reads_transcript == first_transcript) == 1) { # one alignment to alphabetically first transcript
    which_read <- which(input_reads_transcript == first_transcript)
  } else { # read maps to multiple positions in ≥1 transcript, report positionally-first
    input_reads <- subset(input_reads,
                          input_reads_transcript == first_transcript)
    input_reads_position <- sapply(input_reads, function(x) { strsplit(x, split="\t")[[1]][4] })
    which_read <- which.min(input_reads_position)
  }
  return(input_reads[which_read])
}

process_alignment(opt$input, opt$output)