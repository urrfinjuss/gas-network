#!/usr/bin/env python

from Tkinter import Tk, Frame, BOTH, Text, Menu, END, Label, Canvas, Button, Toplevel, Message
from ttk import Style
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from scipy.interpolate import interp1d
from numpy import amax, amin

import time
import tkFileDialog 
import subprocess 
import tkMessageBox as box
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib as mpl
import os, os.path
from sys import exit


class Example(Frame):
	def __init__(self, parent, n):
		Frame.__init__(self, parent)
		self.parent = parent;
		self.parent.maxsize(1200,720)
		self.parent.minsize(1200,720)
		self.initUI()
		self.placeWindow(n)


	def initUI(self):
		self.parent.title("Gas Model")
		self.font1 = {'fontname'   : 'Helvetica', 'color'      : 'r', 'fontweight' : 'normal', 'fontsize'   : 14}
		self.font2 = {'fontname'   : 'Helvetica', 'color'      : 'r', 'fontweight' : 'bold', 'fontsize'   : 14}
		self.npts = 64;  
		self.sounds = 377.9683;
		self.skip = 1000.
		self.dt = 5000.*self.skip/(self.npts * self.sounds)
		menubar = Menu(self.parent)
		menubar.config(font=self.font1)
		self.parent.config(menu=menubar)

		fileMenu = Menu(menubar, tearoff=0)
		fileMenu.config(font=self.font1)
		fileMenu.add_command(label="Load", command=self.onLoad)
		fileMenu.add_command(label="Close", command=self.onClose)
		fileMenu.add_command(label="Exit", command=self.quit)
		menubar.add_cascade(label="File", menu=fileMenu)
		
		paramsMenu = Menu(menubar, tearoff=0)
		paramsMenu.config(font=self.font1)
		paramsMenu.add_command(label="Set Parameters", command=self.onSetParams)
		menubar.add_cascade(label="Parameters", menu=paramsMenu)

		title = Label(self.parent, text="Gas Network", fg = "black", font=self.font2)
	        title.grid(sticky="N", pady=4, padx=5, row=0, column=0, columnspan=2)

		self.f = plt.figure(figsize=(5,4), dpi=80)
		self.cmap = plt.cm.jet #ok
		self.a = self.f.add_subplot(111)
		self.a.axis([-1, 1, -1, 1])
		self.a.axis('off')
		self.canvas = FigureCanvasTkAgg(self.f, master=self.parent)

		self.G = nx.Graph()
		nx.draw_networkx(self.G, pos=nx.spring_layout(self.G), ax=self.a)
		self.canvas.show()
		self.canvas.get_tk_widget().config(width=800-20, height=600-100)
		self.canvas.get_tk_widget().grid(sticky="NW", pady=4, padx=5, row=1, column=0, columnspan=6, rowspan=6) 

		self.txt = Text(self.parent, width=54, height=33);
		self.txt.grid(sticky="NW", pady=4, padx=5, row=1, column=8, columnspan=4, rowspan=6)	

		self.log = Text(self.parent, width=112, height=6);
		self.log.grid(sticky="NW", pady=4, padx=5, row=7, column=0, columnspan=6,rowspan=2)

		RunButton = Button(self.parent, text="Run", width=10, command=self.onRun).grid(sticky='SW', pady=4, padx=5, row=7, column=8)
		DrawButton = Button(self.parent, text="Draw", width=10, command=self.onDraw).grid(sticky='SW', pady=4, padx=5, row=7, column=9)
		StopButton = Button(self.parent, text="Stop", width=10, command=self.onStop).grid(sticky='SW', pady=4, padx=5, row=8, column=9)
		QuitButton = Button(self.parent, text="Quit", width=10, command=self.quit ).grid(sticky='SW', pady=4, padx=5, row=7, column=10)
		self.Loaded = False
		self.Animate = False


	def placeWindow(self, n):
		w = 1200
		h = 720
		sw = self.parent.winfo_screenwidth()
		sh = self.parent.winfo_screenheight()
		if (n == 1) :
			x = sw - w;
			y = sh - h;
		else:
			x = 0;
			y = 0;
		self.parent.geometry('%dx%d+%d+%d' % (w, h, x, y))
		print "Calling Tk Window. Native resolution %d x %d" % (sw, sh)
	def onSetParams(self):
		top = Toplevel(self.parent)
		top.title("Model Setup")
		msg = Message(top, text="")
		msg.pack()
		
		applyb = Button(top, text="Apply", command=top.destroy)
		closeb = Button(top, text="Close", command=top.destroy)
		applyb.pack()
		closeb.pack()
		

	def onOpen(self):      
		ftypes = [('Gas Model files', '*.txt'), ('All files', '*')]
		dlg = tkFileDialog.Open(self, filetypes = ftypes)
		fl = dlg.show()
        
        	if fl != '':
			text = self.readFile(fl)
			self.log.insert(END, "Opened file: %s\n" % fl)
			self.txt.insert(END, text)
			

	def onLoad(self):      
		ftypes = [('Gas Model files', '*.txt'), ('All files', '*')]
		dlg = tkFileDialog.Open(self, filetypes = ftypes)
		fl = dlg.show()
		self.Filename = fl
        	if fl != '':
	        	self.Loaded = True;
			self.log.insert(END, "Loaded model from file:\t%s\n" % fl)
			text = self.read_model(fl)
			self.log.insert(END, "Checking for data along pipes ... ")
			for j in range(self.G.number_of_edges()):
				fname = "%s/data_%03d.txt" % (os.getcwd(), j)
				try:
					f = open(fname,"r")
					if (j == self.G.number_of_edges() - 1):
						self.log.insert(END, "OK\n")
					f.close()
				except IOError:
					self.log.insert(END, "\nFailed to open %s" % fname)
			self.refine_graph();
			self.plotGraph()

	def onClose(self):      
		self.txt.delete(1.0, END)
		self.G.clear()
		self.a.clear()
		self.a.axis([-1, 1, -1, 1])
		self.a.axis('off')
		#nx.draw_networkx(self.G, pos=nx.spring_layout(self.G), ax=self.a)
		self.canvas.draw()
		self.log.insert(END, "Closed\n")
		self.Loaded = False;

	def onClose(self):      
		self.txt.delete(1.0, END)
		self.G.clear()
		self.a.clear()
		self.a.axis([-1, 1, -1, 1])
		self.a.axis('off')
		#nx.draw_networkx(self.G, pos=nx.spring_layout(self.G), ax=self.a)
		self.canvas.draw()
		self.log.insert(END, "Closed\n")
		self.Loaded = False;

	def onStop(self):
		self.Animate = False
		self.log.insert(END, "Stopped\n")

	def onRun(self):

		if self.Loaded:
			self.log.insert(END, "Running\n")
			fhin = open("python.cfg","w")			
			fhin.write("# This is an auto-generated input file for gas network simulaion\n")
			fhin.write("# Pipe lengths must be increments of 10 km, npts number of gridpoints per 1 km of pipe\n")
			fhin.write("# Initial data is read from files data_%03d.txt\n# where the files data_000.txt-data_007.txt must be present in \n# the directory in which code is executed\n\n")
			fhin.write("[Configuration]\n")
			nl = len(self.Filename.split('//'));
			self.log.insert(END, "Writing %s" % self.Filename.split('//')[nl-1])
			fhin.write("#name=\t%s\n" % self.Filename.split('//')[nl-1]) #'model9d.txt'
			fhin.write("#nprof=\tdata\n")
			fhin.write("#nincd=\tincidence.txt\n")
			fhin.write("#nodes=\t%d\n" % self.G.number_of_nodes())
			fhin.write("#ncmps=\t%d\n" % self.ncomps)
			fhin.write("#links=\t%d\n" % self.G.number_of_edges())
			fhin.write("#npts=\t%d\n" % self.npts)
			fhin.write("#distr_gamma=\t%d\n" % 1)
			fhin.write("#fix_pressure=\t%d\n" % 1)
			fhin.write("#drawpng=\t0\n")
			fhin.write("#skipnum=\t%d\n\n" % self.skip)
			fhin.write("[Noise]\n")
			fhin.write("#corr_time=\t10\n")
			fhin.write("#corr_dist=\t20\n")
			fhin.write("#amplitude=\t0\n\n")
			fhin.write("# Simulation parameters: sound speed [m/s],\n# dissipation coefficient [1],\n# simulation time in hours.\n\n")
			fhin.write("[Simulation Parameters]\n")
			fhin.write("#sound_speed=\t%.7f\n" % self.sounds)
			fhin.write("#dissip_coef=\t0.01\n")
			fhin.write("#simul_time=\t2.\n\n")
			fhin.write("# Numeration of pipes goes according incidence.txt:\n# take first column, first nonzero entry in the first column\n# is pipe 1, the second nonzero entry in the first column is pipe 2 and so on.\n#\n# All entries must be in SI units, except pipe lengths which must be given in km.\n\n")
			fhin.write("[Pipe Parameters]\n")
			fhin.write("#lengths\n")
			for m in range(self.G.number_of_edges()):
				for i,j in self.G.edges():
					if (self.G[i][j]['key'] == m):
						fhin.write("%d " % self.G[i][j]['length'])
			fhin.write("\n#diameters\n")
			for m in range(self.G.number_of_edges()):
				for i,j in self.G.edges():
					if (self.G[i][j]['key'] == m):
						fhin.write("%.4f " % (self.G[i][j]['diameter']/39.370))
			fhin.write("\n\n# nodes are numerated by the rows of matrix.cfg\n# e.g. first row corresponds to node 1 and etc\n\n")
			fhin.write("[Node Parameters]\n")
			fhin.write("#compressions\n")
			for m in range(self.G.number_of_nodes()):
				fhin.write("1.0 ")
			fhin.write("\n")
			fhin.close()
			A = nx.adjacency_matrix(self.G);
			U = np.triu(A); L = np.tril(A);
			np.savetxt('incidence.txt', L - U, '%.1f')
			subprocess.call(["./simulate.x", "python.cfg"])

		else:
			self.log.insert(END, "No Model Loaded\n")			

	def onDraw(self):
		self.log.insert(END, "Drawing ... ")
		if self.Loaded:
			for j in range(self.G.number_of_edges()):
				path = "%s/pipe_%03d/" % (os.getcwd(),j)
				num_files = len([f for f in os.listdir(path)if os.path.isfile(os.path.join(path, f))])
			self.log.insert(END, "Found %3d frames\n" % num_files)
			n = 1
			self.Animate = True
			while n < num_files+1:
				self.a.cla()
				self.a.axis([-1, 1, -1, 1])
				self.a.axis('off')
				self.a.set_title("Gas Network colored by values of pressure.")
				self.a.text(0.45,0.05, "Time t = %.4f hours" % (1.*n*self.dt/3600.),transform=self.a.transAxes )
				self.updateColors(n)
				self.parent.update()
				n = n+1
				if not self.Animate:
					break
		else:
			self.log.insert(END,"No model Loaded\n")

	def ispresent(self, node1, node2):
		for e in self.G.edges():
			if node1 in e and node2 in e:
				return True
		return False

	def read_model(self, filename):
		f = open(filename, "r")
		self.DG = nx.DiGraph()
		str =  f.readline(); lst =  str.split()
		nv = int(lst[0]); ne = int(lst[1]); self.ncomps = int(lst[2]);
		for k in range(nv):
			lst = f.readline().split()
			R = np.sqrt(float(lst[3])*float(lst[3]) + float(lst[2])*float(lst[2]))
			Arg = 0; #(float(lst[3]) + float(lst[2]))*2.*np.pi/132;
			self.G.add_node(k, xy = (float(lst[3]))/20. - 1., float(lst[2])/70. - 1.) )
		for k in range(ne):
			lst = f.readline().split()
			self.G.add_edge(int(lst[1]), int(lst[2]), diameter=float(lst[3]), 
				length=float(lst[4]), DW=float(lst[5]), key=1)
		node_pos=[]; edge_dir=[];
		for k in range(self.ncomps):
			lst = f.readline().split()
			node_pos.append(int(lst[1]))
			edge_dir.append(int(lst[2]))
			print node_pos[k], edge_dir[k]

		list = self.G.edges() 
		list = sorted(list, key = lambda x: (x[0], x[1]))
		for k in range(ne):
			self.G[list[k][0]][list[k][1]]['key'] = k;
		
		elabels = dict([((u,v,),int(d['key'])) for u,v,d in self.G.edges(data=True)])

		for u,v,d in self.G.edges_iter(data=True):
			for k in range(self.ncomps):
				if self.G[u][v]['key']==edge_dir[k]:
					#print "Compressor here %d" % edge_dir[k]
					#self.DG.add_edge(u,v)
					self.DG.add_node(u)
					#self.DG.add_node(u+128)
					#self.DG.add_edge(u, u+128)
					(x0, y0) = self.G.node[u]['xy']
					(x1, y1) = self.G.node[v]['xy']
					norm = 0.1/np.sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))
					self.DG.node[u]['xy'] = (x0,y0)
					#self.DG.node[u+128]['xy'] = (x0 + norm*(x1-x0), y0 + norm*(y1-y0))
			
		f.close()
		self.log.insert(END, "Number of vertices %d\n" % self.G.number_of_nodes())
		self.log.insert(END, "Number of edges %d\n" % self.G.number_of_edges())
		self.log.insert(END, "Number of compressors %d\n" % self.ncomps)

		f = open(filename, "r")
		text = f.read()
		self.txt.insert(END, text)


	def refine_graph(self):
		self.H = nx.Graph();
		self.H = self.G.copy(); 

		offset = 512;
		self.inc = 512;

		#pdata = []
		#fdata = []
		#xdata = []
		pmax_loc = np.empty([self.H.number_of_edges(), 1], dtype=np.float64)
		pmin_loc = np.empty([self.H.number_of_edges(), 1], dtype=np.float64)
		for j in range(self.H.number_of_edges()):
			fname = "%s/data_%03d.txt" % (os.getcwd(), j)			
			cdata = np.loadtxt(fname).transpose()
			#xdata.append(cdata[0])
			#pdata.append(cdata[1])
			#fdata.append(cdata[2])
			pmax_loc[j] = amax(cdata[1])
			pmin_loc[j] = amin(cdata[1])
		#self.pmax = amax(pmax_loc) 
		#self.pmin = amin(pmin_loc) 

		PsiTokPa = 6.89475729;
		self.pmax = 800*PsiTokPa;
		self.pmin = 300*PsiTokPa;

		#print self.pmax, self.pmin
		for i,j in self.G.edges():
			x1 = self.G.node[i]['xy'][0]
			y1 = self.G.node[i]['xy'][1]

			x2 = self.G.node[j]['xy'][0]
			y2 = self.G.node[j]['xy'][1]
			
			L = self.G[i][j]['length']; #np.sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
			N = int(2*L);


			dx = (x2-x1)/N;  # according to model.txt
			dy = (y2-y1)/N;  # according to model.txt
			dl = 1./N;        
			
			tmp = offset;

			fname = "%s/data_%03d.txt" % (os.getcwd(), self.G[i][j]['key'])
			cdata = np.loadtxt(fname).transpose()

			l = cdata[0]/self.G[i][j]['length']

			# colorbar [-1:1]
			cm_c = 0.5*(self.pmax+self.pmin)
			cm_a = (self.pmax - self.pmin)   
			cm = 2*(cdata[1] - cm_c)/cm_a

			#colorbar [0:1]
			cm_c = self.pmin;
			cm = (cdata[1] - cm_c)/cm_a

			#colorbar [pmin:pmax]
			#cm = cdata[1]

			fcm = interp1d(l, cm)
			
			self.H.add_edge(i,tmp, diameter=self.G[i][j]['diameter'], length=L/N, DW=self.G[i][j]['DW'], key=self.G[i][j]['key'], color= fcm(0.5*dl).astype(float) )#np.interp(0.5*dl, l, cm) )
			self.H.node[tmp]['xy'] =(x1 + dx, y1 + dy)
			for k in range(N-2):
				self.H.add_edge(tmp,tmp+1, diameter=self.G[i][j]['diameter'], length=L/N, DW=self.G[i][j]['DW'], key=self.G[i][j]['key'], color=fcm((k+1.5)*dl).astype(float)) #, color=np.interp(0.5*dl+(k+1)*dl, l, cm))
				self.H.node[tmp+1]['xy'] = (x1 + (k+2)*dx, y1 + (k+2)*dy)
				tmp = tmp + 1;
			self.H.add_edge(tmp,j, diameter=self.G[i][j]['diameter'], length=L/N, DW=self.G[i][j]['DW'], key=self.G[i][j]['key'], color=fcm((N-0.5)*dl).astype(float))#np.interp(l-0.5*dl, l, cm))	
			self.H.remove_edge(i,j)

			str = "color_%02d.txt" % self.G[i][j]['key']
			cols = np.empty([N,2], dtype = np.float64 );
			for j in range(N):
				cols[j][0] = (0.5+j)*dl
				cols[j][1] = fcm((j+0.5)*dl).astype(float)
			np.savetxt(str, cols)
			

			offset = offset + self.inc;

	def plotGraph(self):
		self.posH = nx.get_node_attributes(self.H,'xy')
		self.posG = nx.get_node_attributes(self.G,'xy')
		self.posDG = nx.get_node_attributes(self.DG,'xy')

		self.color=[]
		for i,j in self.H.edges_iter():
			self.color.append(self.H[i][j]['color'])
		self.elabels = dict([((u,v,),int(d['key'])) for u,v,d in self.G.edges(data=True)])
		#for i,j in self.H.edges_iter():
		#	print self.H[i][j]['color']
		#print elabels
		#nx.draw_networkx_nodes(self.G, posG, with_labels=False, node_size = 50, ax=self.a)
		#nx.draw_networkx_edges(self.G, posG, edge_labels=elabels, edge_color='b', with_labels=True, width = 2, ax=self.a)

		nx.draw_networkx_nodes(self.G, self.posG, with_labels=False, node_size = 150, node_color='#A0CBE2', ax=self.a)
		nx.draw_networkx_edges(self.H, self.posH, edge_color=self.color, edge_cmap=self.cmap, width = 5, alpha=1, ax=self.a) # , edge_vmin=0, edge_vmax=1
		nx.draw_networkx_nodes(self.DG, self.posDG, node_size=150, alpha=1, ax=self.a, arrows=True, edge_color='r')
		#nx.draw_networkx_edges(self.DG, self.posDG, alpha=1, width = 2, ax=self.a, arrows=True, edge_color='r')
		#nx.draw_networkx_edges(self.H, posH, edge_color='k', width = 0.5, alpha=1., ax=self.a)
		nx.draw_networkx_edge_labels(self.G, self.posG, edge_labels=self.elabels, font_size= 12, alpha = 0.4, rotate=True)
		self.a.set_title("Gas Network colored by values of pressure.")
		self.a.text(0.45,0.05, "Time t = %.4f hours" % 0,transform=self.a.transAxes )
		self.ax1 = self.f.add_axes([ 0.02, 0.05, 0.04, 0.9])
		norm = mpl.colors.Normalize(vmin=3, vmax=7)
		cb1 = mpl.colorbar.ColorbarBase(self.ax1, cmap=self.cmap, norm=norm, orientation='vertical') #, norm=norm, orientation='horizontal')
		cb1.set_ticks([3, 4, 5, 6, 7])
		cb1.set_label('Pressure in MPa', labelpad=3)
		self.canvas.draw()
		

	def readFile(self, filename):
		f = open(filename, "r")
		text = f.read()
		return text


	def updateColors(self, nt):
		offset = 512;
		for i,j in self.G.edges():

			L = self.G[i][j]['length']; #np.sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
			N = int(2*L);
			dl = 1./N;        
			
			tmp = offset;

			fname = "%s/pipe_%03d/data_%03d.txt" % (os.getcwd(), self.G[i][j]['key'], nt)
			cdata = np.loadtxt(fname).transpose()

			l = cdata[0]/L

			# colorbar [-1:1]
			cm_c = 0.5*(self.pmax+self.pmin)
			cm_a = (self.pmax - self.pmin)   
      			cm = 2*(cdata[1] - cm_c)/cm_a

			#colorbar [0:1]
			cm_c = self.pmin;
      			cm = (cdata[1] - cm_c)/cm_a
			#print cm 
			
			#colorbar [pmin:pmax]
			#cm = cdata[1]

			fcm = interp1d(l, cm)
			
			self.H[i][tmp]['color']= fcm(0.5*dl)
			for k in range(N-2):
				self.H[tmp][tmp+1]['color'] = fcm((k+1.5)*dl)
				tmp = tmp + 1;
			self.H[tmp][j]['color'] = fcm((N-0.5)*dl)

			str = "color_%02d.txt" % self.G[i][j]['key']
			cols = np.empty([N,2], dtype = np.float64 );
			for j in range(N):
				cols[j][0] = (0.5+j)*dl
				cols[j][1] = fcm((j+0.5)*dl).astype(float)
			np.savetxt(str, cols)
			

			offset = offset + self.inc;
		self.color=[]
		for i,j in self.H.edges():
			self.color.append(self.H[i][j]['color'])

		nx.draw_networkx_nodes(self.G, self.posG, with_labels=False, node_size = 150, node_color='#A0CBE2', ax=self.a)
		nx.draw_networkx_edges(self.H, self.posH, edge_color=self.color, edge_cmap=self.cmap, width = 5, alpha=1, ax=self.a, with_labels=False) # , edge_vmin=0, edge_vmax=1
		nx.draw_networkx_nodes(self.DG, self.posDG, node_size=150, alpha=1, ax=self.a, arrows=True, edge_color='r')
		#nx.draw_networkx_edges(self.H, posH, edge_color='k', width = 0.5, alpha=1., ax=self.a)
		#nx.draw_networkx_edge_labels(self.G, posG, edge_labels=self.elabels, font_size= 12, alpha = 1, rotate=True)
		self.canvas.draw()


		

def main():

	root = Tk()
	app = Example(root, 0)
	root.mainloop()

if __name__ == '__main__':
	main()
















