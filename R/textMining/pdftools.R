# Project:   knowledgeBase
# Objective: Read pdf files in with pdftools
# Author:    Edoardo Costantini
# Created:   2022-03-04
# Modified:  2022-03-04

# Set up -----------------------------------------------------------------------

  # Load packages
  library(pdftools)

# Read pdfs in -----------------------------------------------------------------

  # Define the location of the pdfs
  pdfs_path <- "./pdfs/"

  # Find what pdfs you have in the input folder
  files <- list.files(path = pdfs_path,
                      pattern = "pdf$")

  # Load the texts in your sessions by using pdftools::pdf_text()
  articles <- lapply(paste0(pdfs_path, files), pdf_text)

  # This creates a list of articles
  class(articles)
  length(articles)

  # Every article is a vector of the length of the pdf pages
  sapply(articles, length)

# You can then search for words in the texts with grep style functions

  grep("Gender", articles[[1]])
  # but this is inefficient! Try using the tm package instead