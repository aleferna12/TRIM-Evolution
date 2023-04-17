import json
import re
import pandas as pd
from colorir import RGBGrad, Palette, ColorFormat, HexRGB
from ete3 import Tree, TreeStyle, RectFace, TextFace, StaticItemFace
from ete3.treeview.faces import _label_painter
from PyQt5.QtWidgets import QGraphicsPolygonItem
from PyQt5.QtGui import QPolygonF, QBrush, QColor, QPen
from PyQt5.QtCore import QPointF

P_THRESHOLD = 0.05
MAX_POS_SEL = 8

COLLAPSE = 0
COLLAPSE_WIDTH_CLADE = 4
BRANCH_SCALE = 400
THIN_REC_WIDTH = 2
THICK_REC_WIDTH = 8
SHOW_INFO = False
SHOW_NODES = True
INFO_FONT_SIZE = 6
LABEL_FONT_SIZE = 9
LABEL_FONT = "Calibri"
LEGEND_FONT_SIZE = 18

FILE_SUFFIX = ""

colors = Palette.load(palettes_dir="..", color_format=ColorFormat(HexRGB, include_a=True))


class TriangleItem(QGraphicsPolygonItem):
    def __init__(self, width, height, label="", color=colors.black):

        self.pol = QPolygonF()
        self.pol.append(QPointF(width, 0))
        self.pol.append(QPointF(width, height))
        self.pol.append(QPointF(0, height/2))

        self.label = label
        QGraphicsPolygonItem.__init__(self, self.pol)

        self.setBrush(QBrush(QColor(color)))
        self.setPen(QPen(QColor(color)))

    # Makes up for API changes between QT4 and QT5
    def rect(self):
        return self.boundingRect()

    def paint(self, p, option, widget):
        super(TriangleItem, self).paint(p, option, widget)
        _label_painter(self, p, option, widget)


def collapse_node(node, trim_stats):
    # Collapse node if doesnt have interesting children
    if (
        COLLAPSE
        and not node.is_root()
        and len(node) >= COLLAPSE
        and all(trim_stats[child.name]["Corrected P-value"] > P_THRESHOLD for child in node.traverse())
    ):
        leaves = node.get_leaves()
        dist = node.get_distance(leaves[0])
        polygon_item = TriangleItem(
            dist * BRANCH_SCALE,
            COLLAPSE_WIDTH_CLADE * len(leaves),
            color=colors.protlabel
        )
        node.add_face(StaticItemFace(polygon_item), position="branch-right", column=0)
        node.children = []
        return True
    return False


ts = TreeStyle()
ts.scale = BRANCH_SCALE
ts.show_leaf_name = False
df = pd.read_csv("../cluster_seqs.csv", sep=";", index_col=0)
node_df = pd.DataFrame(columns=["trim", "name", "pos_sel", "aa"])
grad = RGBGrad([colors.negsel, colors.neutsel, colors.possel], use_linear_rgb=True)
pat = re.compile(r"([XN]M_\d+?_\d+)_(\w+?_[a-z]+)")

for trim in [6, 34, 5, 22, 222]:
    with open(f"./aBSREL Results/trim{trim}.json") as file:
        trim_res = json.load(file)
    trim_stats = trim_res["branch attributes"]["0"]
    inp_newick = trim_res["input"]["trees"]["0"] + ";"
    tree = Tree(inp_newick, format=1)
    if trim == 34:
        # Adjust tree dist since trim34 is topology only
        for node in tree.iter_descendants():
            node.dist = 0.05

    tree.convert_to_ultrametric()
    tree.ladderize()
    tree.img_style["size"] = 0
    tree.img_style["vt_line_width"] = THIN_REC_WIDTH

    for node in tree.iter_descendants(is_leaf_fn=lambda n: collapse_node(n, trim_stats)):
        node_stats = trim_stats[node.name]

        if node.is_leaf():
            match = re.search(pat, node.name)
            species = match.group(2).replace("_", " ").lower().capitalize()
            seq_id = match.group(1)[::-1].replace("_", ".", 1)[::-1]
            seq_id = df.index[df["RNAId"] == seq_id].item()
            sp_face = TextFace(text=species, fsize=LABEL_FONT_SIZE, ftype=LABEL_FONT, fstyle="italic")
            sp_face.margin_left = 10
            sp_face.margin_right = 10
            sp_face.hz_align = 0
            id_face = TextFace(text=seq_id, fsize=LABEL_FONT_SIZE, ftype=LABEL_FONT, fgcolor=colors.protlabel)
            id_face.margin_right = 10
            node.add_face(sp_face, column=0, position="aligned")
            node.add_face(id_face, column=1, position="aligned")

        node.img_style["size"] = 0
        p = node_stats["Corrected P-value"]
        if p > P_THRESHOLD:
            height = THIN_REC_WIDTH
        else:
            height = THICK_REC_WIDTH
            pos_sel = list(node_stats["Rate Distributions"][-1])
            pos_sel[0] = "> 8.0" if pos_sel[0] > 8 else str(round(pos_sel[0], 1))
            pos_sel[1] = str(round(pos_sel[1] * 100, 1)) + '%'
            if SHOW_NODES:
                if not node.is_leaf():
                    name = "Br" + str(sum(node_df["name"].str.match(r"Br\d+")) + 1)
                    face = TextFace(name, fsize=INFO_FONT_SIZE)
                    face.margin_right = 3
                    node.add_face(face, column=0, position="branch-top")
                else:
                    name = seq_id  # Previously obtained seq_id
            node_df.loc[len(node_df)] = [f"TRIM{trim}", name] + pos_sel
            if SHOW_INFO:
                face = TextFace(f"{pos_sel[0]}|{pos_sel[1]}", fsize=INFO_FONT_SIZE)
                face.margin_right = 3
                node.add_face(face, column=0, position="branch-bottom")
        node.img_style["hz_line_width"] = height
        node.img_style["vt_line_width"] = THIN_REC_WIDTH

        for i, r_class in enumerate(node_stats["Rate Distributions"]):
            r_class[0] = min(r_class[0], MAX_POS_SEL)
            r_color = grad.perc((r_class[0] / MAX_POS_SEL) ** (1/3))
            width = r_class[1] * node.dist * BRANCH_SCALE
            rect_face = RectFace(width=width, height=height, bgcolor=r_color, fgcolor=colors.transp)
            node.add_face(rect_face, column=i, position="float")

    tree.render(f"Figures/trim{trim}_absrel{FILE_SUFFIX}.pdf", tree_style=ts, w=1000)

node_df.sort_values(by=["trim", "name"], inplace=True)
node_df.to_csv("Figures/absrel_table.csv", sep=';', index=False)