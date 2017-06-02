# -*- coding: utf-8 -*-
"""
Created on Fri May 27 14:02:14 2016

@author: steven
"""

def main():

    # check the command-line arguments:
    if len(sys.argv) != 2 or os.path.exists(sys.argv[1]) == False:
        print("Usage: %s newick_file" % sys.argv[0])
        sys.exit(1)
    newick_file = sys.argv[1]

    # the code below is based on examples on http://pythonhosted.org/ete2/tutorial/tutorial_trees.html
    # load a tree from a file. The root node is put into 't'
    t = Tree(newick_file, format=1) # format=1 means internal node names are loaded

    # iterate over all the nodes in the tree, and get the children of each node:
    for node in t.traverse("postorder"):
        # get the children of this node:
        children = node.get_children()
        for child in children:
            print "node %s has child %s" % (node.name, child.name)

    # iterate over all the nodes in the tree, and get the parent of each node:  
    for node in t.traverse("postorder"):
        # get the parent of this node:
        if not node.is_root(): # if it is not the root node
            print "node %s has parent %s" % (node.name, node.up.name)

    # iterate over all the nodes in the tree, and get the descendant species of each node:
    for node in t.traverse("postorder"):
        if not node.is_leaf(): # if it is not a leaf 
            # get the descendant species (leaves):
            leaves = node.get_leaves()
            leaf_names = [leaf.name for leaf in leaves]
            print "node %s has descendants %s" % (node.name, ','.join(leaf_names))

    # iterate over all the nodes in the tree, and get the ancestor nodes of each node, back 
    # to the root:
    for node in t.traverse("postorder"):
        # get the ancestors of this node: 
        if not node.is_root():
            ancestor = node.up
            ancestors = [ancestor]
            while not ancestor.is_root():
                ancestor = ancestor.up
                ancestors.append(ancestor)
            ancestor_names = [ancestor.name for ancestor in ancestors]
            print "node %s has ancestors %s" % (node.name, ','.join(ancestor_names))

    print("FINISHED")

#====================================================================#

if __name__=="__main__":
    main()

#====================================================================#