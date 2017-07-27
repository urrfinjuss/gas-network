#!/usr/bin/env python
import gi
import networkx as nx
import matplotlib.pyplot as plt
from numpy import arange, sin, pi
from matplotlib.backends.backend_gtk3agg import FigureCanvasGTK3Agg as FigureCanvas
try:
	import pygraphviz
	from networkx.drawing.nx_agraph import graphviz_layout
except ImportError:
	try:
		import pydotplus
		from networkx.drawing.nx_pydot import graphviz_layout
	except ImportError:
		raise ImportError("This example needs Graphviz and either PyGraphviz or PyDotPlus")

gi.require_version('Gtk', '3.0')
from gi.repository import Gtk

class Handler:
	def onDeleteWindow(self, *args):
		Gtk.main_quit(*args)

	def onButtonPressed(self, button):
		print("Hello World!")



builder = Gtk.Builder()
builder.add_from_file("../glade_widgets/main_window.glade")
builder.connect_signals(Handler())

window = builder.get_object("window1")

# -- do the plot
f = plt.figure(figsize=(5,4), dpi=100)
a = f.add_subplot(111)
plt.axis('off')


# -- do the network graph
G = nx.balanced_tree(3, 5)
pos = graphviz_layout(G, prog='twopi', args='')
nx.draw(G, pos, node_size=20, alpha = 0.5, node_color="blue", with_labels=False)
xlim = a.get_xlim()
ylim = a.get_ylim()

#t = arange(0.0, 3.0, 0.01)
#s = sin(2*pi*t)
#a.plot(t, s)
# -- get the small window
swindow1 = builder.get_object("scrolled_window1")
swindow1.set_border_width(10)

canvas = FigureCanvas(f)
canvas.set_size_request(800,600)
#canvas.show()
swindow1.add_with_viewport(canvas)

window.show_all()
Gtk.main()

