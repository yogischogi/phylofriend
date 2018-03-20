# Phylofriend

Phylofriend’s main purpose is to calculate genetic
distances from Y-chromosome STR values. The results can be
used as input for the [PHYLIP](http://evolution.genetics.washington.edu/phylip.html)
program to create phylogenetic trees.

Phylofriend supports all 587 STR markers that are
reported by YFull (March 2018).


## Documentation

* [User Guide](https://github.com/yogischogi/phylofriend/blob/master/doc/phylofriend.pdf?raw=true)
* [Source Code](http://godoc.org/github.com/yogischogi/phylofriend)


## How to Create a Phylogenetic Tree

1. Copy persons’ data from a Family Tree DNA project website into a spreadsheet.
   If the Y-STR values do not appear properly try inserting them into the spreadsheet as unformatted text.
2. Save the spreadsheet in CSV (comma separated values) format, for example *persons.csv*.
3. Start a terminal or command line interpreter and go to the directory where you stored *persons.csv*.
4. Create a matrix of genetic distances by typing **phylofriend -personsin persons.csv -phylipout infile**.
5. Use the PHYLIP program to create a tree in Newick format with **/usr/lib/phylip/bin/kitsch**.
   You will need to answer some questions. Usually the default values are good enough.
   The results will be two text files, one named *outtree* which contains the tree in
   Newick format and another one named *outfile* which contains a more human readable description.
6. Create an image of the tree by typing **/usr/lib/phylip/bin/drawgram**.
   Use *outtree* as the input file name.
   The resulting tree will be stored in a file named *plotfile*. 


