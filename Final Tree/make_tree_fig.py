#! python3
import sys
print(sys.executable)
from ete3 import Tree, TreeStyle, TextFace, RectFace
from colorir import Palette
import re

NEWICK_PATH = "./final.tree"
# Has to be a pdf since svg is bugged at the momment
PDF_PATH = "./tree.pdf"

ALIGN_LABELS = True
CIRCULAR_TREE = True
ULTRAMETRIC = True
SHOW_DUPS = False
BRANCH_SCALE = 750
LINE_WIDTH = 2
LABEL_FONT_SIZE = 9
LABEL_FONT = "Calibri"
LEGEND_FONT_SIZE = 18
LEGEND_SQ_SIZE = 50

SUPPORT_THRESHOLD = 70

colors = Palette.load(palettes_dir="..")


def make_layouts(t):
    # Creates backgrounds for these specific nodes  
    trim6 = t.get_common_ancestor("XP_037696214.1", "XP_012642053.1")
    trim34 = t.get_common_ancestor("XP_012414743.2", "NP_001192122.1")
    trim22 = t.get_common_ancestor("XP_023501633.1", "XP_016005550.2")
    trim5 = t.get_common_ancestor("XP_004389130.1", "XP_037696204.1")
    trim6.img_style["bgcolor"] = colors.trim6
    trim34.img_style["bgcolor"] = colors.trim34
    trim22.img_style["bgcolor"] = colors.trim22
    trim5.img_style["bgcolor"] = colors.trim5

    cached_sps = t.get_cached_content(store_attr="species")

    accepted_clades = 0
    total_clades = 0
    for node in t.traverse():
        if node.is_root():
            node.img_style["fgcolor"] = "Gray"
        elif node.is_leaf():
            node.img_style["size"] = 0
            if ALIGN_LABELS:
                pos = "aligned"
            else:
                pos = "branch-right"
            sp = re.search(r"[A-Z][a-z]+ [a-z]+", node.species).group()
            sp_face = TextFace(text=sp, fsize=LABEL_FONT_SIZE, ftype=LABEL_FONT, fstyle="italic")
            sp_face.margin_left = 10
            sp_face.margin_right = 10
            sp_face.hz_align = 0
            node.add_face(sp_face, column=0, position=pos)
            id_face = TextFace(text=node.name, fsize=LABEL_FONT_SIZE, ftype=LABEL_FONT, fgcolor=colors.protlabel)
            id_face.margin_right = 10
            node.add_face(id_face, column=1, position=pos)
        else:
            if SHOW_DUPS and node not in [trim6.up, trim6.up.up]:
                for sp in cached_sps[node]:
                    if sp in cached_sps[node.children[0]] and sp in cached_sps[node.children[1]]:
                        node.img_style["vt_line_color"] = colors.red
                        break
            if node.support >= SUPPORT_THRESHOLD:
                support_face = TextFace(int(node.support), ftype=LABEL_FONT, fsize=LABEL_FONT_SIZE)
                support_face.margin_right = 6
                node.add_face(support_face, 1, position="branch-top")
                accepted_clades += 1
            node.img_style["size"] = 0
            total_clades += 1

        node.img_style["vt_line_width"] = LINE_WIDTH
        node.img_style["hz_line_width"] = LINE_WIDTH

    print(f"{round(accepted_clades/total_clades * 100, 2)}% of clades were accepted (bootstrap >= {SUPPORT_THRESHOLD})")


def main():
    t = Tree(NEWICK_PATH, format=0)
    t.ladderize()
    ts = TreeStyle()
    ts.allow_face_overlap = True
    ts.draw_guiding_lines = True
    ts.show_leaf_name = False
    ts.arc_span = 360
    ts.mode = 'c' if CIRCULAR_TREE else 'r'
    ts.scale = BRANCH_SCALE
    ts.show_scale = not ULTRAMETRIC
    if ULTRAMETRIC: t.convert_to_ultrametric()

    for name in ["trim6", "trim34", "trim5", "trim22"]:
        t_face = TextFace(text=name.upper(), fsize=LEGEND_FONT_SIZE, ftype=LABEL_FONT)
        t_face.margin_top = 10
        t_face.margin_right = 10
        t_face.hz_align = 2
        sq_face = RectFace(width=LEGEND_SQ_SIZE, height=LEGEND_SQ_SIZE, fgcolor=colors.black, bgcolor=colors.get_color(name))
        sq_face.margin_top = 10
        sq_face.margin_right = 10
        ts.legend.add_face(t_face, 0)
        ts.legend.add_face(sq_face, 1)

    make_layouts(t)

    t.render(PDF_PATH, tree_style=ts, w=1000)


if __name__ == "__main__":
    main()
