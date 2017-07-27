#!/usr/bin/env python
import gi
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import time

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
builder.add_from_file("../glade_widgets/pipe_window.glade")
builder.connect_signals(Handler())

window = builder.get_object("window1")
sw     = builder.get_object("scr_win")
sw.set_border_width(10)

# -- read pipe data file
name = "pipe_000"
for indx in [0, 100]:
	filename = "../../result/%s/data_%03d.txt" % (name, indx);
	print("Loading:\t%s" % filename)
	raw = np.loadtxt(filename, usecols=(0,1,2));
	nrows, ncols = raw.shape;
	L = raw[nrows-1,0];

	# -- do the plot
	f = plt.figure(figsize=(5,4), dpi=100)
	a1 = f.add_subplot(211)
	a1.plot( raw[:,0], raw[:,1])
	a1.set_xlim(0, L)
	a1.set_ylim(0, 10.0)
	a1.set_ylabel(r"Pressure, MPa")
	a1.set_title('Pipe %3d' % indx)
	a1.axis('on')

	a2 = f.add_subplot(212)
	a2.plot( raw[:,0], raw[:,2])
	a2.set_xlim(0, L)
	a2.set_ylim(-250.0, 250.0)
	a2.set_ylabel("Flux, kg/m^2/s")
	a2.set_xlabel('Distance, km')
	a2.axis('on')

	# -- plot to GTK backend
	canvas = FigureCanvas(f)
	canvas.set_size_request(500,400)
	sw.add_with_viewport(canvas)
	window.show_all()
	time.sleep(1)

Gtk.main()

