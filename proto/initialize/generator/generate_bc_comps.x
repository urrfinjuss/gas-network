#! /usr/bin/octave -qf

function generate_comp(fname, dt, tmax)

	nsteps = tmax*3600./dt;
	fh = fopen(fname,"w");
	fprintf(fh, "# 1. Time, secs 2. Compression Ratio\n");
	fprintf(fh, "# Auto generated by octave script: generate_bc_comps.x\n\n");
	for k = 0:nsteps
	  fprintf(fh, "%.12e\t%.12e\n", k*dt, 1.0);
	end
	fclose(fh);

end

conf = "../../python.cfg";
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
	for j = 1:nC
		txt = fgetl(fh);
		fname = sprintf("../bc/comps/comp_%03d.txt", j-1);
		generate_comp(fname, dt, tmax);
	endfor
else
	fprintf("Cannot find network configuration file.\n");
	exit(0);
endif
