# Project:   knowledgeBase
# Objective: Quick guide to read html files
# Author:    Edoardo Costantini
# Created:   2022-03-04
# Modified:  2022-03-04

# Set up -----------------------------------------------------------------------

  # Load packages
  library(XML)

# Create a function to extract text from a local html file

  extractHTMLtext <- function (local_html){
    # Source: https://www.r-bloggers.com/2011/10/reading-html-pages-in-r-for-text-processing/
    # Read and parse HTML file
    doc.html <- htmlTreeParse(local_html, useInternal = TRUE)

    # Extract all the paragraphs (HTML tag is p, starting at the root of the document)
    doc.text <- xpathApply(doc.html, '//p', xmlValue)

    # Unlist flattens the list to create a character vector.
    doc.text <- unlist(doc.text)

    # Replace all \n by spaces
    doc.text <- gsub('\\n', ' ', doc.text)

    # Join all the elements of the character vector into a single
    # character string, separated by spaces
    doc.text <- paste(doc.text, collapse = ' ')

    # Return
    return(doc.text)
  }

  trial_text <- extractHTMLtext("./html/717103.html")

# Prepare for text mining analysis ---------------------------------------------

  # Define the location of the pdfs
  html_path <- "./html/"

  # Find what pdfs you have in the input folder
  files <- list.files(path = html_path,
                      pattern = "html$")

  # Load the texts in your sessions by using extractHTMLtext
  articles <- lapply(paste0(pdfs_path, files), extractHTMLtext)

  # Give names
  names(articles) <- paste0("article", 1:length(articles))

  # Create a corpus
  corp <- tm::Corpus(VectorSource(articles))

  # Use quanteda to extract contexts for words of interest
  q_corp <- quanteda::corpus(corp)
  toks <- quanteda::tokens(q_corp)
  quanteda::kwic(toks, pattern = "polit*", valuetype = "glob", window = 5)
  out <- quanteda::kwic(toks,
                        pattern = quanteda::phrase("multiple imputation*"),
                        window = 10)

  out <- quanteda::kwic(toks,
                          pattern = quanteda::phrase("missing values*"),
                          window = 10)
  out$docname