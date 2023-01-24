# CutFree

CutFree is a tool created to design pools of DNA barcodes compatible with restriction enzymes.

## CutFree GUI: https://jensenlab.shinyapps.io/cutfree/

## Instructions for Use

### Download Julia
Download Julia from the following link: https://julialang.org/downloads/

When installing Julia, ensure you click "Add to Path" for proper functionality.

### Clone Repository
Clone this repository.

```
git clone "https://github.com/wrickord/CutFree.git"
```

### CutFree in Julia
Straightforward.

### CutFree in the Command Line
To use CutFree in the command line, follow the steps below.

Open command prompt and navigate to the cloned repository.
```
cd "ADD/PATH/HERE/CutFree"
```

To run CutFree, use the Main-CMD.jl file as seen below.
```
julia Main-CMD.jl --starting_oligo "NNNNNNNNNNNNNNNNNNNN" --restriction_sites "GGTCTC,GGCCGG" --min_blocks 1 --increase_diversity true
```

#### Explanation of Arguments
--starting_oligo (String): Input an oligo of any size to create a randomer of that size.

--restriction_sites (String): Input a vector of the sequences of the restriction enzymes you wish to block.
##### NOTE: The list of restiction sites must be entered as a single string with sites separated by commas (i.e., "GGTCTC,GGCCGG")

--min_block (Int): Input the minimum number of blocks you desire.

--increase_diversity (Bool): Input a boolean to tell the program whether or not it should find the most random output while maintaining its degeneracy.
