#! /usr/bin/octave -qf

function generate_pressure_node(fname, dt, tmax, pressure)
# generates the set of boundary data at the nodes of the network
# for the fixed pressure nodes.
# input:
#   n         number of network nodes, 
#   tmax      max time, hours
#
# the default time increment is 1.0 second.


	nsteps = tmax*3600./dt;
	fh = fopen(fname,"w");
	fprintf(fh, "# 1. Time, secs 2. Pressure, MPa\n");
	fprintf(fh, "# Auto generated by octave script: generate_bc_nodes.x\n\n");
	freq = 6.*pi/tmax/3600;
	for k = 0:nsteps
	  fprintf(fh, "%.12e\t%.12e\n", k*dt, pressure + 0.0*sin(freq*k*dt));
	end
	fclose(fh);

end

function generate_demand_node(fname, dt, tmax, flow)
# generates the set of boundary data at the nodes of the network
# for the fixed demand nodes.
# input:
#   n         number of network nodes, 
#   tmax      max time, hours
#
# the default time increment is 1.0 second.


	nsteps = tmax*3600./dt;
	fh = fopen(fname,"w");
	if (fh ~= -1)
		fprintf(fh, "# 1. Time, secs 2. Demand, kg/s\n");
		fprintf(fh, "# Auto generated by octave script: generate_bc_nodes.x\n\n");
		for k = 0:nsteps
			if (k <= nsteps/2)
				fprintf(fh, "%.12e\t%.12e\n", k*dt, flow);
			else 
				fprintf(fh, "%.12e\t%.12e\n", k*dt, flow + 4.0); % switch from 0.5 to 4.0
			endif
		end
		fclose(fh);
	else 
		fprintf("Failed to open: %s\n", fname);
		exit(0);
	endif
end

conf = "../../python.cfg";
pr = [6.5e+6, 6.5e+6, 6.5e+6, 6.5e+6];
fl = [1.0e+2, 6.0e+1, 4.0e+1, 1.0e+1];
wd = [0.9144, 0.9144, 0.9144, 0.6350];
Sd = 0.25*pi*wd.^2;
fh = fopen(conf, "r");

if (fh)
  while (~feof(fh))
		line = fgets(fh);
  	[value1, value2] = sscanf(line, "%s\t%s", "C");
		if     (strcmp("name=", value1)) network = value2;
		elseif (strcmp("tmax=", value1)) tmax = str2num(value2);
		endif
	endwhile
  fclose(fh);
else
  fprintf("Could not find base configuration file: %s", conf);
	exit(0);
endif
network = strcat("../../", network);

dt = 1.0;			# dt   -- time step size (secs)

fprintf("Call: generate_bc_nodes.x\n");
fprintf("Time discretization:\t\t%.2f secs\n", dt);
fprintf("Final time:\t\t%12.2f hours\n", tmax);

fh = fopen(network, "r");
if (fh ~= -1) 
	txt = fgetl(fh);
	[nN, nP, nC] = sscanf(txt, "%d\t%d\t%d", "C");
	nType = zeros(nN, 1);
	for j = 1:nN
		txt = fgetl(fh);
		[id, ~, ~, nType(j)] = sscanf(txt, "%d\t%lf\t%lf\t%d", "C");
		fname = sprintf("../bc/nodes/node_%03d.txt", id);
		if (nType(j) == 1) 
			fprintf("Calling pressure init for node %d\n", id);
			generate_pressure_node(fname, dt, tmax, 1.0e-06*pr(j));
		elseif (nType(j) == 0) 
			fprintf("Calling demand init for node %d\n", id);
			generate_demand_node(fname, dt, tmax, 0.);
		else
			fprintf("Unknown Node Type %d.\n", id);
			exit(0);
		endif
	endfor
else
	fprintf("Cannot find network configuration file.\n");
	exit(0);
endif
generate_demand_node("../bc/nodes/node_002.txt", dt, tmax, fl(2)*Sd(2)-fl(4)*Sd(4));
generate_demand_node("../bc/nodes/node_003.txt", dt, tmax, fl(3)*Sd(3)+fl(4)*Sd(4));


#generate_bc_nodes_p(num1, 1.0)
#generate_bc_nodes_f(num2, 1.0)
