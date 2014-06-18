setwd('/Users/elizabethsweeney/Dropbox/Elizabeth_Sweeney_Documents/Current_Projects/Munich/code_for_nplarge')
install.packages('Matrix_1.1-4.tar.gz', repos = NULL, type="source")
install.packages('nplargedb_0.1-2.tar.gz', repos = NULL, type="source")
install.packages('nplargela_0.1-2.tar.gz', repos = NULL, type="source")

##############################################
## Explore S4 objects that are modelterms 
##############################################

showClass("modelterm")

