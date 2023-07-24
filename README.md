# CutFree
CutFree constructs randomers free of user-specified restriction enzyme recognition sites to prevent improper DNA assembly in workflows utilizing restriction enzymes. The solution can be derived deterministically via mixed-integer linear programming (CutFree) or heuristically via reinforcement learning and stoichastic search (CutFreeRL). An algorithm classification model was developed by fine-tuning a pre-trained DistilBERT and multi-layer perceptron.

## CutFree GUI: <br>https://jensenlab.shinyapps.io/cutfree/

## CutFree Paper: <br>https://pubmed.ncbi.nlm.nih.gov/28865135/

## CutFreeRL Paper: <br>https://pubmed.ncbi.nlm.nih.gov/35143615/

# Instructions for Use
## Download Julia
Download Julia from the following link: https://julialang.org/downloads/

When installing Julia, ensure you click "Add to Path" for proper functionality.

Add to path command for MacOS terminal (may need to change version number).
```
rm -f /usr/local/bin/julia
ln -s /Applications/Julia-1.8.app/Contents/Resources/julia/bin/julia /usr/local/bin/julia
```

## Download Conda
Download Anaconda from the following link: https://www.anaconda.com/products/distribution#Downloads

To use Conda in Julia, it must be installed locally.

## Download Gurobi
Download Gurobi from the following link: https://www.gurobi.com/downloads/

After installing Gurobi, make sure to activate your key in the command line using the instructions provided.

## Clone Repository
Clone this repository.
```
git clone "https://github.com/wrickord/CutFree.git"
```

## Changing the Model
In addition to using the Gurobi MILP optimization tool, CutFree's implementation makes it convenient to utilize open-source optimizers, such as GLPK (download included in venv).

Simply install your open-source optimizer of choice, add "Using OPTIMIZER_NAME" (where OPTIMIZER_NAME=GLPK, for example) to the top of CutFree.jl, nagivate to **line 315 of CutFree.jl**, and change "Gurobi" to the name of your optimizer.

Example:
```
model = Model(GLPK.Optimizer)
```

## CutFree in the Command Line
To use CutFree in the command line, follow the steps below.

Open command prompt/terminal and navigate to the cloned repository.
```
cd "ADD/PATH/HERE/CutFree"
```

To run CutFree, use the Main.jl file as seen below.
```
julia Main.jl --starting_oligo "NNNNNNNNNNNNNNNNNNNN" --restriction_sites "GGTCTC,GGCCGG" --min_blocks 1 --increase_diversity true
```

### Explanation of Arguments
Starting Oligo:
- Starting DNA sequence that should be blocked from containing restriction sites. To generate a set of barcodes with the highest diversity, start with a string of N's the length of your oligo.

Restriction Enzyme Recognition Sites to Block:
- Sequences to block from the oligo pools. Separate multiple sequences by commas. Common names are accepted for many popular enzymes. For example, "BsaI" and "GGTCTC" are both valid inputs.

Minimum # of Blocks:
- Minimum number of blocks at each site, i.e. the minimum number of changes that need to be made for a cut site to appear anywhere in the oligo. Values greater than 1 are only valid for CutFree algorithm, not CutFreeRL.

Increase Diversity:
- Toggle to re-run the optimization to randomize the codes while maintaining the same number of oligos (i.e., same degeneracy). Only valid for CutFree algorithm, not CutFreeRL.

Algorithm Choice:
- Auto: Automatically chooses the best algorithm for the given input.
- CutFree: Uses the CutFree algorithm.
- CutFreeRL: Uses the CutFreeRL algorithm.