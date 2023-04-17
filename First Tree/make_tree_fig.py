#! py -3.10

from ete3 import Tree, TreeStyle, TextFace
import random as rd

NEWICK_PATH = "./upto_600aa_clustalo_iter2.tree"
# Has to be a pdf since svg is bugged at the momment
PDF_PATH = "./tree.pdf"

ALIGN_LABELS = False
CIRCULAR_TREE = False
LABEL_FONT = "Calibri"
LABEL_FONT_SIZE = 9

SH_ALRT_THRESHOLD = 80
SUPPORT_THRESHOLD = 95


# Paints paths to sequences belonging to the same locus
# This helps to visualize loci that may be hybrid or poorly conserved
def paint_loci(t, loci):
    for loc, tups in loci.items():
        if len(tups) > 1:
            rd.seed(int(loc))
            # Creates random seeded color
            color = "#" + "".join(rd.choices('ABCDEF0123456789', k=6))
            prots = []
            for tup in tups:
                prots.append(tup[0])
                tup[1].fgcolor = color
                tup[2].fgcolor = color
                tup[3].fgcolor = color
            c_anc = t.get_common_ancestor(prots)
            for prot in prots:
                node = prot
                while node != c_anc:
                    node.img_style["hz_line_color"] = color
                    node = node.up


def make_layouts(t):
    # Creates backgrounds for these specific nodes  
    trim6 = t.get_common_ancestor("XP_037696217.1", "XP_044114634.1")
    trim34 = t.get_common_ancestor("XP_023598787.1", "NP_569074.2")
    trim22 = t.get_common_ancestor("XP_023501638.1", "XP_036077539.1")
    trim5 = t.get_common_ancestor("XP_037696204.1", "XP_026345387.1")
    trim6.img_style["bgcolor"] = "#D9D7F1"
    trim34.img_style["bgcolor"] = "#FFFDDE"
    trim22.img_style["bgcolor"] = "#E7FBBE"
    trim5.img_style["bgcolor"] = "#FFCBCB"

    accepted_clades = 0
    total_clades = 0
    loci = {}
    for node in t.traverse():
        if node.is_root():
            node.img_style["fgcolor"] = "Gray"
        elif node.is_leaf():
            node.img_style["size"] = 0
            sp_face = TextFace(text=node.species, fsize=LABEL_FONT_SIZE, ftype=LABEL_FONT, fstyle="italic", bold=True)
            sp_face.margin_left = 10
            trim_face = TextFace(text=node.trim, fsize=LABEL_FONT_SIZE, ftype=LABEL_FONT)
            trim_face.margin_left = 5
            trim_face.margin_right = 5
            id_face = TextFace(text=node.name, fsize=LABEL_FONT_SIZE, ftype=LABEL_FONT)
            id_face.margin_right = 10
            if ALIGN_LABELS:
                pos = "aligned"
            else:
                pos = "branch-right"
            node.add_face(sp_face, column=0, position=pos)
            node.add_face(trim_face, column=1, position=pos)
            node.add_face(id_face, column=2, position=pos)

            # Adds to list of proteins to be painted later
            loci.setdefault(node.gene_id, [])
            loci[node.gene_id].append((node, sp_face, trim_face, id_face))
        # Make circle according to whether the clade should be accepted as well supported (see: http://www.iqtree.org/doc/Frequently-Asked-Questions#how-do-i-interpret-ultrafast-bootstrap-ufboot-support-values)
        else:
            if float(node.sh_alrt) >= SH_ALRT_THRESHOLD and node.support >= SUPPORT_THRESHOLD:
                color = "Green"
                accepted_clades += 1
            else:
                color = "Red"

            node.img_style["fgcolor"] = color
            node.img_style["size"] = 12
            node.add_face(TextFace(node.sh_alrt, ftype=LABEL_FONT, fsize=LABEL_FONT_SIZE), 0, position="branch-top")
            node.add_face(TextFace(node.support, ftype=LABEL_FONT, fsize=LABEL_FONT_SIZE), 0, position="branch-bottom")
            total_clades += 1

        node.img_style["vt_line_width"] = 2
        node.img_style["hz_line_width"] = 2
    paint_loci(t, loci)
    print(f"{round(accepted_clades/total_clades * 100, 2)}% of clades were accepted (shown with green dots)")


def main():
    # Branch support is stored in node.name
    t = Tree(NEWICK_PATH, format=0)
    ts = TreeStyle()
    ts.mode = 'c' if CIRCULAR_TREE else 'r'
    ts.show_leaf_name = False
    ts.arc_span = 360
    ts.scale = 2000

    make_layouts(t)

    t.render(PDF_PATH, tree_style=ts, w=1000)


if __name__ == "__main__":
    main()
