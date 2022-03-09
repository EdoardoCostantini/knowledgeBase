# Project:   knowledgeBase
# Objective: Quick guide to text mine (tm) R package
# Author:    Edoardo Costantini
# Created:   2022-03-04
# Modified:  2022-03-04
# Other ref: https://data.library.virginia.edu/reading-pdf-files-into-r-for-text-mining/
#            https://cran.r-project.org/web/packages/tm/vignettes/tm.pdf
#            https://www.r-bloggers.com/2011/10/reading-html-pages-in-r-for-text-processing/

# Set up -----------------------------------------------------------------------

  rm(list = ls())

# Load packages
  library(tm)

# Define the location of the pdfs
  pdfs_path <- "./pdfs/"

# Find what pdfs you have in the input folder
  files <- list.files(path = pdfs_path,
                      pattern = "pdf$")

# Create a corpus
  corp <- tm::Corpus(URISource(paste0(pdfs_path, files)),
                     readerControl = list(reader = readPDF))

# Simple example ---------------------------------------------------------------

# Individual documents can be accesed
  meta(corp[[2]], "id")
  inspect(corp[[2]])

# Create a term-document matrix (TDM)
  articles_tdm <- TermDocumentMatrix(corp,
                                     control =
                                       list(removePunctuation = TRUE,
                                            stopwords = TRUE,
                                            tolower = TRUE,
                                            stemming = TRUE,
                                            removeNumbers = TRUE,
                                            bounds = list(global = c(3, Inf))))

# Inspect the TDM matrix
  inspect(articles_tdm[1:10, ])

# Quickly find frequently occurring terms
  # find the most common terms
  highFreq_words <- findFreqTerms(articles_tdm, lowfreq = 100, highfreq = Inf)
  highFreq_words
  highFreq_table <- as.matrix(articles_tdm[highFreq_words, ])
  highFreq_table

  # summarize them
  sort(apply(highFreq_table, 1, sum), decreasing = TRUE)

# Transformations --------------------------------------------------------------

# Convert to lower case
  corp_lc <- tm_map(corp, content_transformer(tolower))
  inspect(corp[[2]])
  inspect(corp_lc[[2]])

# Remove stop words
  corp_nsp <- tm_map(corp, removeWords, stopwords("english"))
  inspect(corp[[2]])
  inspect(corp_nsp[[2]])

# Filters ----------------------------------------------------------------------
# This is possible only based on the tags store in the meta-data
  meta(corp[[1]])
  meta(corp[[2]])
  meta(corp[[3]])
  meta(corp[[4]])

  idx <- meta(corp, "origin") == "Arbortext Advanced Print Publisher 9.1.520/W"
  corp[idx]

# Creating Term-Document Matrices ----------------------------------------------

# Clean up the corp
  # Make lower case
  corp <- tm_map(corp, content_transformer(tolower))
  # Remove stop words
  corp <- tm_map(corp, removeWords, stopwords("english"))

# Create a TDM
  tdm <- TermDocumentMatrix(corp)
  inspect(tdm)

# Clean up the tdm
  inspect(removeSparseTerms(tdm, 0.4))

# Find frequent terms
  findFreqTerms(tdm, 100)

# Find words that are correlated
  findAssocs(tdm, "imputation", 0.8)

# Find a word of interest
  # list of words
  tdm$dimnames$Terms
  slam::row_sums(tdm[tdm$dimnames$Terms == "imputation", ])

# Find context for words used using
  q_corp <- quanteda::corpus(corp)
  toks <- quanteda::tokens(q_corp)
  quanteda::kwic(toks, pattern = "imputation*", valuetype = "glob", window = 5)
  quanteda::kwic(toks,
                 pattern = quanteda::phrase("multiple imputation*"),
                 window = 10)
