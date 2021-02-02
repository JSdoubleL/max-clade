# Max-Clade

Extracts maximal clades from a list of unrooted trees. In an unrooted context we define clades to be the tree on either half of any bipartition. Maximal clades are clades which do not contain duplications, however every clade that they are a subset of does.

## Dependencies

- Python 3
- [treeswift](https://github.com/niemasd/TreeSwift)

## Usage

**Input**: File containing list of *multi-copy* trees in newick format

**Output**: Maximal clades listed as trees in newick format

```
python max_clade.py -i <input_file> -o <ouput_file> -t -d <delimiter>
```

### Arguments

- `-i`: Input tree list file
- `-o`: (optional) Output max clade list file
- `-t`: (optional) Include trivial clades in output (filtered by default).
- `-d`: (optional) Delimiter separating species name from rest of leaf label (default None).

### Example

```
python max_clade.py -i example/gtrees-mult.trees
```

## Algorithm Description

For each tree in the input do the following:

- If the tree is rooted, remove the root and connect the root's children with an edge. Then arbitrarily select internal vertex and call it the root. (The vertex chosen does not affect the result in any way; it is necessary, however, for defining destinct sides for the bipartitions).
- *Postorder traversal:* label each edge with True/False indicating whether there's a duplication event under it (*under* defined with respect to the rooting). Store references to all clades without duplications.
- *Preorder traversal:* label each edge with True/False indicating whether there's a duplication event above it. Store references to all clades without duplications.
- Iterate through all duplication free clades. For each one, check all the edges which define a superset clade one layer up. If any of them are duplication free, then considered clade is not maximal, otherwise it is.
