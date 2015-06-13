#include "network.h"
#include <mgl2/mgl_cf.h>

static HMGL gr;
static HMDT x, p, m, q, g;
static HMDT xp, pp, mp;
static double cs;

void mgl_init_draw(network *net) {
  xp = mgl_create_data();
  pp = mgl_create_data();
  mp = mgl_create_data();
  x = mgl_create_data();
  p = mgl_create_data();
  m = mgl_create_data();
  q = mgl_create_data();
  g = mgl_create_data();
  cs = net->c;
}

void mgl_draw_pipe(gpipe_ptr in, network *net, char* fname, char* title) {
  char str1[80], str2[80];
  
  gr = mgl_create_graph(800,600);
  mgl_data_create(x, in->N, 1, 1);
  mgl_data_create(p, in->N, 1, 1);
  mgl_data_create(m, in->N, 1, 1);
  //mgl_data_create(q, in->N, 1, 1);

  for (int k = 0; k < in->N; k++) {
     mgl_data_set_value(x, (in->L)*k/(in->N - 1.), k, 0, 0);
     mgl_data_set_value(p, 0.001*(in->p[k]), k, 0, 0);
     mgl_data_set_value(m, (in->f[k])/cs, k, 0, 0);
     mgl_data_set_value(q, (in->q[k]), k, 0, 0);
  }
  
  mgl_add_legend(gr, title,"");

  mgl_subplot(gr, 2, 1, 0, "");
  mgl_title(gr, "Pressure", "#", -2);
  mgl_set_ranges(gr, 0, in->L, -10, 1200, 0, 0);
  mgl_axis(gr,"xyz","","");
  mgl_plot_xy(gr, x, p, "b1","");
  mgl_legend(gr, 3, "", "");
  mgl_box(gr);

  mgl_subplot(gr, 2, 1, 1, "");
  mgl_title(gr, "Flux", "#", -2);
  mgl_set_ranges(gr, 0, in->L, -1200, 1200, 0, 0);
  mgl_axis(gr,"xyz","","");
  mgl_plot_xy(gr, x, m, "b1","");
  mgl_box(gr);

  mgl_write_png_solid(gr, fname,"");
  mgl_delete_graph(gr);
}





void mgl_draw_pressure(network *net, char* fname, char* title) {
  char str1[80], str2[80], msg[80];
  gpipe_ptr in;
  int l, np = 8;   // must be even
  gr = mgl_create_graph(800,600);

  mgl_set_ranges(gr, 0, 1, -10, 800, 0, 0);
  mgl_axis(gr,"xyz","","");
  mgl_label(gr, 'x', "x/L", 0, "");
  mgl_label(gr, 'y', "Pressure, kPa", 0, "");
  for (int nl = 0; nl < net->nlinks; nl ++) {
    in = net->link[nl];
    mgl_data_create(xp, np, 1, 1);
    mgl_data_create(pp, np, 1, 1);
    mgl_data_create(x, in->N, 1, 1);
    mgl_data_create(p, in->N, 1, 1);

    for (int k = 0; k < in->N; k++) {
       mgl_data_set_value(x, k/(in->N - 1.), k, 0, 0);
       mgl_data_set_value(p, 0.001*(in->p[k]), k, 0, 0);
    }
    l = 0;
    for (int k = floor((in->N)/2)-np/2-1; k < floor((in->N)/2)+np/2; k++) {
       mgl_data_set_value(xp, k/(in->N - 1.), l, 0, 0);
       mgl_data_set_value(pp, 0.001*(in->p[k]), l, 0, 0);
       l++;
    }
    sprintf(msg, "pipe %d", nl);
    mgl_plot_xy(gr, x, p, "b","");
    mgl_text_xy(gr, xp, pp, msg, "k:rC", "");
  }
  mgl_box(gr);
  mgl_write_png_solid(gr, fname,"");
  mgl_delete_graph(gr);
}

void mgl_draw_flux(network *net, char* fname, char* title) {
  char str1[80], str2[80], msg[80];
  gpipe_ptr in;
  int l, np = 8;   // must be even
  gr = mgl_create_graph(800,600);

  mgl_set_ranges(gr, 0, 1, -10, 100, 0, 0);
  mgl_axis(gr,"xyz","","");
  mgl_label(gr, 'x', "x/L", 0, "");
  mgl_label(gr, 'y', "Flux, kg/(m^2 s)", 0, "");
  for (int nl = 0; nl < net->nlinks; nl ++) {
    in = net->link[nl];
    mgl_data_create(xp, np, 1, 1);
    mgl_data_create(mp, np, 1, 1);
    mgl_data_create(x, in->N, 1, 1);
    mgl_data_create(m, in->N, 1, 1);
    for (int k = 0; k < in->N; k++) {
       mgl_data_set_value(x, k/(in->N - 1.), k, 0, 0);
       mgl_data_set_value(m, (in->f[k])/cs, k, 0, 0);
    }
    l = 0;
    for (int k = floor((in->N)/2)-np/2-1; k < floor((in->N)/2)+np/2; k++) {
       mgl_data_set_value(xp, k/(in->N - 1.), l, 0, 0);
       mgl_data_set_value(mp, (in->f[k])/cs, l, 0, 0);
       l++;
    }
    sprintf(msg, "pipe %d", nl);
    mgl_plot_xy(gr, x, m, "b","");
  }
  mgl_box(gr);
  mgl_write_png_solid(gr, fname,"");
  mgl_delete_graph(gr);
}



