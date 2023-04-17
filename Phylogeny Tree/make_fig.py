from ete3 import Tree, TreeStyle, TextFace

t = Tree("phylogeny.nwk", format=9)
t.convert_to_ultrametric()
t.ladderize()
for node in t.traverse():
    node.img_style["size"] = 0
    node.img_style["vt_line_width"] = 3
    node.img_style["hz_line_width"] = 3
    face = TextFace(
        text=node.name.replace("_", " "), 
        ftype="Calibri",
        fstyle="italic",
        fsize=12,
    )
    face.margin_left = 6
    node.add_face(face, 0, position="aligned")
ts = TreeStyle()
ts.show_leaf_name = False
ts.show_scale = False
ts.scale = 40
ts.branch_vertical_margin = 10
t.render("phylogeny.pdf", w=200, tree_style=ts)
