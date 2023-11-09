// C/C++ headers
#include <cassert>

// Athena++ headers
#include "utils.hpp"

void read_data_table(char const *fname, double** data, int *rows, int *cols)
{
  FILE *fp = fopen(fname, "r");

  // read dimension of the table
  fscanf(fp, "%d%d", rows, cols);

  NewCArray(data, *rows, *cols);

  char buf[256];
  int irow = 0, icol;
  char *p;
  while (fgets(buf, 256, fp) != NULL) {
    p = std::strtok(buf," ,");
    int icol = 0;
    while (p != NULL) {
      sscanf(p, "%lf", &data[irow][icol]);
      p = std::strtok(NULL," ,");
      assert(icol++ < *cols);
    }
    assert(irow++ < *rows);
  }

  fclose(fp);
}

void ReadDataTable(AthenaArray<Real>& data, std::string fname, char c) {
  // remove comment
  std::string str_file = DecommentFile(fname);

  // read first time to determine dimension
  std::stringstream inp(str_file);
  //std::ifstream inp(fname.c_str(), std::ios::in);
  std::string line;
  std::getline(inp, line);
  int rows = 0, cols = 0;
  if (!line.empty()) {
    rows = 1;
    cols = line[0] == c ? 0 : 1;
    for (int i = 1; i < line.length(); ++i)
      if (line[i-1] == c && line[i] != c) cols++;
  }
  while (std::getline(inp, line)) ++rows;
  rows--;
  //inp.close();

  //str_file = DecommentFile(fname);

  // read second time
  data.NewAthenaArray(rows, cols);
  std::stringstream inp2(str_file);
  //inp.open(fname.c_str(), std::ios::in);

  for (int i = 0; i < rows; ++i)
    for (int j = 0; j < cols; ++j)
      inp2 >> data(i,j);
}

void ReadDataTableForInterp(std::string fname, std::vector<Real>& file_x, std::vector<Real>& file_y, int& n_file, bool enforce_ascending) {
  AthenaArray<Real> file_data;
  ReadDataTable(file_data, fname);

  n_file = file_data.GetDim2();
  file_x.reserve(n_file);
  file_y.reserve(n_file);

  for (int i = 0; i < n_file; ++i)
  {
    file_x[i] = file_data(i,0);
    file_y[i] = file_data(i,1);

    if (enforce_ascending && i > 0 && file_x[i] < file_x[i-1]) {
      std::stringstream msg;
      msg << "###### FATAL ERROR in ReadDataTableForInterp" << std::endl
          << "file \"" << fname<< "\" must be in ascending order." << std::endl;
      ATHENA_ERROR(msg);
    }
  }
}