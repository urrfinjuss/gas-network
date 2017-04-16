#!/usr/bin/env python
import numpy as np
#import oct2py as op
import networkx as nx
import matplotlib.pyplot as plt

#class GNet():
#	def __init__(self):
#		self.G = nx.Graph()
#		self.read_model('../model9.txt')	# must be passed as argument
#	def read_model(self,filename):
#		print 'Reading', filename
#		fh = open(filename,"r")		
#		str = fh.readline().split();		
#		self.nVertex = int(str[0]);
#		self.nEdge = int(str[1]);
#		self.nComp = int(str[2]);
#		print str
		

#	def read_model(self, filename):
#		f = open(filename, "r")
#		self.DG = nx.DiGraph()
#		str =  f.readline(); lst =  str.split()
#		nv = int(lst[0]); ne = int(lst[1]); self.ncomps = int(lst[2]);
#		for k in range(nv):
#			lst = f.readline().split()
#			self.G.add_node(k, xy = (float(lst[3]), float(lst[2])) )
#		for k in range(ne):
#			lst = f.readline().split()
#			self.G.add_edge(int(lst[1]), int(lst[2]), diameter=float(lst[3]), 
#				length=float(lst[4]), DW=float(lst[5]), key=1)
#		node_pos=[]; edge_dir=[];
#		for k in range(self.ncomps):
#			lst = f.readline().split()
#			node_pos.append(int(lst[1]))
#			edge_dir.append(int(lst[2]))
#			print node_pos[k], edge_dir[k]
#
#		list = self.G.edges() 
#		list = sorted(list, key = lambda x: (x[0], x[1]))
#		for k in range(ne):
#			self.G[list[k][0]][list[k][1]]['key'] = k;
#		
#		elabels = dict([((u,v,),int(d['key'])) for u,v,d in self.G.edges(data=True)])
#
#		for u,v,d in self.G.edges_iter(data=True):
#			for k in range(self.ncomps):
#				if self.G[u][v]['key']==edge_dir[k]:
#					#print "Compressor here %d" % edge_dir[k]
#					#self.DG.add_edge(u,v)
#					self.DG.add_node(u)
#					#self.DG.add_node(u+128)
#					#self.DG.add_edge(u, u+128)
#					(x0, y0) = self.G.node[u]['xy']
#					(x1, y1) = self.G.node[v]['xy']
#					norm = 0.1/np.sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))
#					self.DG.node[u]['xy'] = (x0,y0)
#					#self.DG.node[u+128]['xy'] = (x0 + norm*(x1-x0), y0 + norm*(y1-y0))
#			
#		f.close()
#		self.log.insert(END, "Number of vertices %d\n" % self.G.number_of_nodes())
#		self.log.insert(END, "Number of edges %d\n" % self.G.number_of_edges())
#		self.log.insert(END, "Number of compressors %d\n" % self.ncomps)
#
#		f = open(filename, "r")
#		text = f.read()
#		self.txt.insert(END, text)


def generate_single(Pl,Pr,L, npts, D, c, alpha, fname):
  #L = [50,80,80,80]*1e3; npts = 64;
  #D = [0.9144,0.9144,0.9144,0.6350];
  #c  =  377.9683; alpha = 0.01;


	dx = np.array(1./npts); 
	N = L*npts + 1;
	x = np.linspace(0,L*1e3, N); 
  	beta = alpha*np.power(c,2)/D;
  	f = np.sqrt(  np.abs(   (np.power(Pl*1e6,2)   -   np.power(Pr*1e6,2)  )    )/beta/(L*1e3)   )
        f = f*np.sign(np.power(Pl*1e6,2)-np.power(Pr*1e6,2));
  	p = np.sqrt(  np.power(Pl*1e6,2) - beta*f*np.abs(f)*x   );
  
	print dx,N,alpha,c,D,beta,f,Pl,Pr

 	fn = open(fname, "w");
 	fn.write('# 1. x, km 2. P, MPa 3. Flux, kg/m^2/s 4. q, kg/s\n# N = %d\n\n' % (N-1))
 	for j in range(x.size):
   		fn.write('%.12e\t%.12e\t%.12e\t%.12e\n' % (1e-3*x[j], p[j]*1e-6, f, 0.));
 	
 	fn.close();

def generate_bc(fname, net, Nc):
	f = open(fname, "r")
	fc = open('../c.txt', "w");
	fd = open('../d.txt', "w");
        fc.write('# 1. time, sec 2. Compression Ratios\n\n');
        fd.write('# 1. time, sec 2. Demands, kg/s\n\n');
  
 	while True:
		line = f.readline()
		if not line: break
		fc.write('%.12e' % float(line.split()[0]))
		fd.write('%.12e' % float(line.split()[0]))
		for ic in range(Nc+1): 
			for u,v,d in GasN.edges(data=True):
				if (GasN[u][v]['c_key'] == ic):
					fc.write('\t%.12e' % float(line.split()[ic+1]))
		fc.write('\n')
		#print GasN.nodes(data=True)[1][1]['slack']
		l = 0;
		for j in GasN.nodes():
			#print GasN.nodes(data=True)[j-1][1]['slack'],j,GasN.number_of_nodes()
			l = l+1
 			if GasN.nodes(data=True)[j][1]['slack']:
				#print 'Slack Node %d' % l
				fd.write('\t%.12e' % float(0.))
			else:
				fd.write('\t%.12e' % float(line.split()[Nc+l]))
		fd.write('\n')
	f.close()
	fc.close()
	fd.close()

filename = '../csv/model9.txt'
iname = '../python.cfg'
mname = 'model9.txt'
f_press = '../csv/model9_ex1_ic.csv'
f_presPl = '../csv/model9_ex1_PPin.csv'
f_presPr = '../csv/model9_ex1_PPout.csv'
f_fluxPl = '../csv/model9_ex1_FFin.csv'
f_fluxPr = '../csv/model9_ex1_FFout.csv'
f_params = '../csv/model9_ex1_parameters.csv'

npts = 16
c  =  377.9683
alpha = 0.01
time = 24.

GasN = nx.Graph()
DG = nx.DiGraph()


f = open(filename,"r")
str = f.readline().split()

Ncompressors = int(str[2]);

fp = open(f_press,"r")
line = fp.readline().split()
for k in range(int(str[0])):
	lst = f.readline().split()
	if (int(lst[8]) == 0):
		GasN.add_node(k, xy = (float(lst[3]), float(lst[2])), P = float(line[k]),  slack=False)
	else:
		GasN.add_node(k, xy = (float(lst[3]), float(lst[2])), P = float(line[k]),  slack=True)
fp.close()

fpPl = open(f_presPl,"r")
ffPl = open(f_fluxPl,"r")
fpPr = open(f_presPr,"r")
ffPr = open(f_fluxPr,"r")

linPl = fpPl.readline().split()
linFl = ffPl.readline().split()
linPr = fpPr.readline().split()
linFr = ffPr.readline().split()
for k in range(int(str[1])):
	lst = f.readline().split()
	GasN.add_edge(int(lst[1]), int(lst[2]), diameter=float(lst[3]), compressor = False, 
	length=int(lst[4]), DW=float(lst[5]), key=k, P_l = float(linPl[k]), c_key = -1,
	P_r = float(linPr[k]), flux = float(linFl[k]))

fpPl.close()
ffPl.close()
fpPr.close()
ffPr.close()

#print GasN.nodes(data=True)
for k in range(int(str[2])):
	lst = f.readline().split()
	for u,v,d in GasN.edges_iter(data=True):
       		if (int(lst[1]) == u) and (int(lst[2])==GasN[u][v]['key']):
			GasN[u][v]['compressor'] = True
			GasN[u][v]['c_key'] = int(lst[0]);

f.close()

for u,v,d in GasN.edges(data=True):
        fname = '../data_%03d.txt' % (GasN[u][v]['key'])
        #print fname, GasN[u][v]['P_l'], GasN[u][v]['P_r'], GasN[u][v]['length'],
	generate_single(GasN[u][v]['P_l'], GasN[u][v]['P_r'], GasN[u][v]['length'], 
	npts, GasN[u][v]['diameter'], c, alpha, fname);


generate_bc(f_params, GasN, Ncompressors)

A = nx.adjacency_matrix(GasN);
U = np.triu(A); L = np.tril(A);
np.savetxt('../incidence.txt', L - U, '%.1f')


fhin = open(iname,"w")			
fhin.write("# This is an auto-generated input file for gas network simulaion\n")
fhin.write("# Pipe lengths must be increments of 10 km, npts number of gridpoints per 1 km of pipe\n")
fhin.write("# Initial data is read from files data_%03d.txt\n# where the files data_000.txt-data_007.txt must be present in \n# the directory in which code is executed\n\n")
fhin.write("[Configuration]\n")
fhin.write("#name=\t%s\n" % mname) #'model9d.txt'
fhin.write("#nprof=\tdata\n")
fhin.write("#nincd=\tincidence.txt\n")
fhin.write("#nodes=\t%d\n" % GasN.number_of_nodes())
fhin.write("#ncmps=\t%d\n" % Ncompressors)
fhin.write("#links=\t%d\n" % GasN.number_of_edges())
fhin.write("#npts=\t%d\n" % npts)
fhin.write("#distr_gamma=\t%d\n" % 1)
fhin.write("#fix_pressure=\t%d\n" % 1)
fhin.write("#drawpng=\t0\n")
fhin.write("#skipnum=\t%d\n\n" % 1000)
fhin.write("[Noise]\n")
fhin.write("#corr_time=\t10\n")
fhin.write("#corr_dist=\t20\n")
fhin.write("#amplitude=\t0\n\n")
fhin.write("# Simulation parameters: sound speed [m/s],\n# dissipation coefficient [1],\n# simulation time in hours.\n\n")
fhin.write("[Simulation Parameters]\n")
fhin.write("#sound_speed=\t%.7f\n" % c)
fhin.write("#dissip_coef=\t%.7f\n" % alpha)
fhin.write("#simul_time=\t%.7f\n\n" % time)
fhin.write("# Numeration of pipes goes according incidence.txt:\n# take first column, first nonzero entry in the first column\n# is pipe 1, the second nonzero entry in the first column is pipe 2 and so on.\n#\n# All entries must be in SI units, except pipe lengths which must be given in km.\n\n")
fhin.write("[Pipe Parameters]\n")
fhin.write("#lengths\n")
for m in range(GasN.number_of_edges()):
	for i,j in GasN.edges():
		if (GasN[i][j]['key'] == m):
			fhin.write("%d " % GasN[i][j]['length'])
fhin.write("\n#diameters\n")
for m in range(GasN.number_of_edges()):
	for i,j in GasN.edges():
		if (GasN[i][j]['key'] == m):
			fhin.write("%.4f " % (GasN[i][j]['diameter']))
fhin.write("\n\n# nodes are numerated by the rows of matrix.cfg\n# e.g. first row corresponds to node 1 and etc\n\n")
fhin.write("[Node Parameters]\n")
fhin.write("#compressions\n")
for m in range(GasN.number_of_nodes()):
	fhin.write("1.0 ")
fhin.write("\n")
fhin.close()
		




