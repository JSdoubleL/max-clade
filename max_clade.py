import treeswift
import argparse
import copy
import types


def get_down_clades(tree):
    """
    Adapted from script "preproecess_multrees_v3.py" by Erin Molloy 
    https://github.com/ekmolloy/fastmulrfs 

    Gets dup free clades below each edge.

    Parameters
    ----------
    tree : treeswift tree object

    Returns list of vertices coresponding to clades without duplicates
    """
    clades = []

    for node in tree.traverse_postorder():
        node.dup_down = False
        if node.is_leaf():
            node.down = set([node.get_label()])
        else:
            node.down = set([])
            for child in node.child_nodes():
                expected_size = len(node.down) + len(child.down)
                node.down = node.down.union(child.down)
                if child.dup_down or len(node.down) < expected_size:
                    node.dup_down = True 
            # Check for duplications, if not add to list
            if not node.dup_down:
                clades.append(node)

    return clades


def get_up_clades(tree):
    """
    Adapted from script "preproecess_multrees_v3.py" by Erin Molloy 
    https://github.com/ekmolloy/fastmulrfs 

    Gets dup free clades above each edge. 

    NOTE: Must be called after 'get_down_clades'

    Parameters
    ----------
    tree : treeswift tree object

    Returns list of vertices coresponding to clades without duplicates
    """

    clades = []

    # Find root node
    for node in tree.traverse_preorder():
        if node.is_root():
            root = node
            root.skip = True
        #elif node.is_leaf():
        #    node.skip = True
        else:
            node.skip = False

    # Compute up profiles for children of root
    children_of_root = root.child_nodes()
    root.dup_down = False
    for node in children_of_root:
        # check if there's duplications below the root
        # in order to take care of edge case
        if node.dup_down:
            root.dup_down = True
        node.up = set([])
        node.skip = True
        node.dup_up = False
        for sibl in children_of_root:
            if node != sibl:
                expected_size = len(node.up) + len(sibl.down)
                node.up = node.up.union(sibl.down)
                if sibl.dup_down or len(node.up) < expected_size:
                    node.dup_up = True
        if not node.dup_up: 
            clades.append(node)
    children_of_root = set(children_of_root) # possible junk

    # Compute up profiles for remaining nodes
    for node in tree.traverse_preorder():
        label = node.get_label() # possible junk
        if not node.skip:            
            parent = node.get_parent()
            node.dup_up = parent.dup_up
            node.up = parent.up
            for sibl in parent.child_nodes():
                if node != sibl:
                    expected_size = len(node.up) + len(sibl.down)
                    node.up = node.up.union(sibl.down)
                    if sibl.dup_down or len(node.up) < expected_size:
                        node.dup_up = True
            if not node.dup_up:
                clades.append(node)

    return clades


def is_max_clade(node, up):
    """
    Calculates whether node corresponds to maximal clade by checking all possible
    superset clades.

    NOTE: Must be run after both 'get_up_clades' and 'get_down_clades'

    Parameters
    ----------
    node : vertex of a treeswift tree     
    up : boolean, when true we examine edge above vertex, else below

    Returns True if node corresponds to a maximal clade
    """
    if up:
        for child in node.child_nodes():
            if not child.dup_up and child != node:
                return False
    else:
        parent = node.get_parent()        
        if not parent.is_root() and not parent.get_parent().dup_down:
            return False
        for sibl in parent.child_nodes():
            if not sibl.dup_down and sibl != node:
                return False

    return True
    

def node_to_tree(tree, node, up):
    """
    Gets clade as tree given node

    Parameters
    ----------
    tree : original treeswift tree
    node : vertex of a treeswift tree     
    up : boolean, when true clade is above node, else below

    Returns treeswift subtree coresponding to clade
    """
    if up:
        parent = node.get_parent()
        c_tree = copy.copy(tree)
        c_tree.reroot(node)
        return c_tree.extract_subtree(parent)
    else:
        return tree.extract_subtree(node)


def find_max_clades(tree):
    """
    Find all the max clades in a tree.

    Parameters
    ----------
    tree : treeswift tree object

    Returns list containing all clades as treeswift tree objects
    """
    max_clades = []
    down_clades = get_down_clades(tree)
    up_clades = get_up_clades(tree)

    # if there's no duplication events below the root, then there's no duplications in tree
    if not tree.root.dup_down:
        max_clades.append(tree)
    else:
        for clade in down_clades:
            if is_max_clade(clade, False):
                max_clades.append(node_to_tree(tree, clade, False))
        for clade in up_clades:
            if is_max_clade(clade, True):
                max_clades.append(node_to_tree(tree, clade, True))

    return max_clades


def trivial(newick_str):
    """
    Determines if a newick string represents a trivial tree (tree containing now quartets).

    Parameters
    ----------
    newick_str: newick string

    Returns True if tree contains less than two '('
    """
    count = 0
    for c in newick_str:
        if c == '(':
            count += 1
        if count > 1:
            return False
    return True


def main(args):
    with open(args.input, 'r') as fi:
        #trees = treeswift.read_tree_newick(args.input)
        with open(args.output, 'w') as fo:
            for line in fi:
                tree = treeswift.read_tree_newick(line)
                max_clades = find_max_clades(tree)
                for c in max_clades:
                    c.suppress_unifurcations()
                    c.deroot() # return unrooted trees (also needed to make trivial() work)
                    c.is_rooted = False # for some reason deroot() does not do this
                    newk = c.newick()
                    if not trivial(newk):
                        fo.write(newk + '\n')


if __name__=="__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input gene tree file", required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Output max clade list file", required=True)

    main(parser.parse_args())