# CutFree

CutFree is a tool created to design pools of DNA barcodes compatible with restriction enzymes.

Feel free to check out CutFree through the following link: https://jensenlab.shinyapps.io/cutfree/

# Instructions for Use

## In Julia

## In the Command Line
To use CutFree in the command line, follow the steps below.

### 1.) Download Julia
Download Julia from the following link: https://julialang.org/downloads/

When installing Julia, ensure you click "Add to Path" for proper functionality.

### 2.) Opening Julia in Command Prompt
Clone this repository.

```
git clone "https://github.com/wrickord/CutFree.git"
```

Open the repository in your command prompt window.

```
cd "ADD/PATH/HERE/CutFree"
```

### 3.) Running CutFree

To run CutFree, uuse the Main-CMD.jl file as seen below.

```
julia Main-CMD.jl --starting_oligo "NNNNNNNNNNNNNNNNNNNN" --restriction_sites "GGTCTC,GGCCGG" --min_blocks 1 --increase_diversity true
```

--starting_oligo (String): Input an oligo of any size to create a randomer of that size.
--restriction_sites (String): Input a vector of the sequences of the restriction enzymes you wish to block.
#### NOTE: The list of restiction sites must be entered as a single string with sites separated by commas (i.e., "GGTCTC,GGCCGG")
--min_block (Int): Input the minimum number of blocks you desire.
--increase_diversity (Bool): Input a boolean to tell the program whether or not it should find the most random output while maintaining its degeneracy.
