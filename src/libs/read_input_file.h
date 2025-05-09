#ifndef READ_INPUT_FILE_H
#define READ_INPUT_FILE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_FILENAME_LEN 256
#define MAX_LINES 1000

typedef struct {
  double R;
  double C00;
  double C6;
  char output_filename[MAX_FILENAME_LEN];
} DataRow;

#define MAX_ROWS MAX_LINES

// Function declarations
int read_file(const char *filename);
DataRow get_row(int row_number);
int get_total_rows();

#endif // READ_INPUT_FILE_H

