# CutFree

CutFree is a tool created to design pools of DNA barcodes compatible with restriction enzymes.

## CutFree Paper: <br>https://pubmed.ncbi.nlm.nih.gov/28865135/

## CutFreeRL Paper: <br>https://pubmed.ncbi.nlm.nih.gov/35143615/

## CutFree GUI: <br>https://jensenlab.shinyapps.io/cutfree/

## Instructions for Use

### Download Julia
Download Julia from the following link: https://julialang.org/downloads/

When installing Julia, ensure you click "Add to Path" for proper functionality.

Add to path command for MacOS terminal (may need to change version number).
```
rm -f /usr/local/bin/julia
ln -s /Applications/Julia-1.8.app/Contents/Resources/julia/bin/julia /usr/local/bin/julia
```

### Download Gurobi
Download Gurobi from the following link: https://www.gurobi.com/downloads/

After installing Gurobi, make sure to activate your key in the command line using the instructions provided.

### Clone Repository
Clone this repository.
```
git clone "https://github.com/wrickord/CutFree.git"
```

### Changing the Model
In addition to using the Gurobi MILP optimization tool, CutFree's implementation makes it convenient to utilize open-source optimizers, such as GLPK (download included in cutfree-venv).

Simply install your open-source optimizer of choice, add "Using OPTIMIZER_NAME" (where OPTIMIZER_NAME=GLPK, for example) to the top of CutFree.jl, nagivate to **line 311 of CutFree.jl**, and change "Gurobi" to the name of your optimizer.

Example:
```
model = Model(GLPK.Optimizer)
```

### CutFree in Julia
Straightforward.

### CutFree in the Command Line
To use CutFree in the command line, follow the steps below.

Open command prompt/terminal and navigate to the cloned repository.
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
