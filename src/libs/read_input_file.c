#include "read_input_file.h"

// Static variables are now private to this file
static DataRow data[MAX_ROWS];
static int total_rows = 0;

int read_file(const char *filename) {
  FILE *file = fopen(filename, "r");
  if (!file) {
    perror("Error opening file");
    return 1;
  }

  total_rows = 0;
  while (fscanf(file, "%lf %lf %lf %s", &data[total_rows].R, &data[total_rows].C00, 
                &data[total_rows].C6, data[total_rows].output_filename) == 4) {
    total_rows++;
    if (total_rows >= MAX_ROWS) {
      fprintf(stderr, "Exceeded maximum number of rows (%d)\n", MAX_ROWS);
      break;
    }
  }

  fclose(file);
  return 0;
}

DataRow get_row(int row_number) {
  if (row_number < 0 || row_number >= total_rows) {
    fprintf(stderr, "Row number %d is out of bounds (0-%d)\n", row_number, total_rows - 1);
    exit(EXIT_FAILURE);
  }
  return data[row_number];
}

int get_total_rows() {
  return total_rows;
}