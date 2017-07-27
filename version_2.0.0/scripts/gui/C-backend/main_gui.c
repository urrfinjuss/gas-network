#include <stdlib.h>
#include <gtk/gtk.h>
#include <glib/gprintf.h>
 
int main(int argc, char *argv[]) {
	GtkWidget       *window;
	GtkBuilder      *builder; 

	gtk_init(&argc, &argv);
	g_printf("Detected GTK+ %d.%d.%d\n", gtk_major_version, gtk_minor_version, gtk_micro_version);
	printf("Required GTK+ 3.20 or higher\n");
	builder = gtk_builder_new();

	gtk_builder_add_from_file (builder, "main_window.glade", NULL);
	window = GTK_WIDGET(gtk_builder_get_object(builder, "window_main"));
	if (window == NULL) {
		printf("Failed to start: Main Window\n");
		exit(0);
	}
	gtk_builder_connect_signals(builder, NULL);
 
	g_object_unref(builder);
 
	gtk_widget_show(window);                
	gtk_main();
 
	return 0;
}
 
// called when window is closed
void on_window_main_destroy() {
	gtk_main_quit();
}

void on_window_main_quit() {
	gtk_main_quit();
}
