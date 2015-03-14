# About Phylofriend #

Phylofriend’s main purpose is to calculate genetic distances from Family Tree Y-DNA data. The results can be used as input for the [PHYLIP](http://evolution.genetics.washington.edu/phylip.html) program to create phylogenetic trees.

# Documentation #

  * [User Guide](https://phylofriend.googlecode.com/git/doc/phylofriend.pdf)
  * [Source Code](http://godoc.org/code.google.com/p/phylofriend)

# How to Create a Phylogenetic Tree #

  * Copy persons’ data from a Family Tree DNA project website into a spreadsheet. If the Y-STR values do not appear properly try inserting them into the spreadsheet as unformatted text.
  * Save the spreadsheet in CSV (comma separated values) format, for example _persons.csv_.
  * Start a terminal or command line interpreter and go to the directory where you stored _persons.csv_.
  * Create a matrix of genetic distances by typing `phylofriend -personsin persons.csv -phylipout infile`.
  * Use the PHYLIP program to create a tree in Newick format with `/usr/lib/phylip/bin/kitsch`. You will need to answer some questions. Usually the default values are good enough. The results will be two text files, one named _outtree_ which contains the tree in Newick format and another one named _outfile_ which contains a more human readable description.
  * Create an image of the tree by typing `/usr/lib/phylip/bin/drawgram`. Use _outtree_ as the input file name. The resulting tree will be stored in a file named _plotfile_.