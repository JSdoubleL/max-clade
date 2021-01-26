# Max-Clade

Extracts maximal clades from a list of unrooted trees. In an unrooted context we define clades to be the tree on either half of any bipartition. Maximal clades are clades which do not contain duplications, however every clade that they are a subset of does.

## Dependencies

- Python 3
- [treeswift](https://github.com/niemasd/TreeSwift)

## Usage

```cmd
python max_clade -i <input_file> -o <ouput_file>
```

*Input*: File containing list of trees in newick format

*Output*: Maximal clades listed as trees in newick format

## Algorithm Description

For each tree in the input do the following:

- Arbitrarily root the tree (treeswift does this by default).
- *Postorder traversal:* label each edge with True/False indicating whether there's a duplication event under it (under defined with respect to the rooting). Store references to all clades without duplications.
- *Preorder traversal:* label each edge with True/False indicating whether there's a duplication event above it. Store references to all clades without duplications.
- Iterate through all duplication free clades. For each one, check all the edges which define a superset clade one layer up. If any of them are duplication free, then considered clade is not maximal, otherwise it is.
