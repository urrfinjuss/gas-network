#include "header.h"

void create_sample_hdf5_file() {
  hid_t		file_id;
  herr_t	status;

  file_id = H5Fcreate("file.h5", H5F_ACC_TRUNC, 
  			H5P_DEFAULT, H5P_DEFAULT);
  status = H5Fclose(file_id);
}
