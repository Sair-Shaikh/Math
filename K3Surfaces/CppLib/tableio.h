#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <set>
#include "constants.h"

std::tuple<unsigned**, unsigned**, unsigned**, unsigned**, unsigned*, unsigned*> generate_all_tables();


void write_table(unsigned* table, int size, std::string fname);
void write_table(int* table, int size, std::string fname);
void write_table(unsigned** table, int size1, int size2, std::string fname);
void write_table(unsigned*** table, int size1, int size2, int size3, std::string fname);

int* read_table(int size, std::string fname, int why);
unsigned* read_table(int size, std::string fname);
unsigned** read_table(int size1, int size2, std::string fname);
unsigned*** read_table(int size1, int size2, int size3, std::string fname);

const std::string dirname = "/Users/sairshaikh/Math/K3Surfaces/Fq_tables/";


/////////////////
// Saving
////////////////

void write_table(unsigned* table, int size, std::string fname){

  // Write data to file
  std::ofstream file;
  file.open(dirname + fname, std::ios_base::binary);
  for (unsigned i = 0; i < size; i++)
    file.write((char*)(&table[i]), sizeof(unsigned));
  
  file.close();
  return;
}

void write_table(int* table, int size, std::string fname){

  // Write data to file
  std::ofstream file;
  file.open(dirname + fname, std::ios_base::binary);
  for (unsigned i = 0; i < size; i++)
    file.write((char*)(&table[i]), sizeof(int));
  
  file.close();
  return;
}

void write_table(unsigned** table, int size1, int size2, std::string fname){

  // Write data to file
  std::ofstream file;
  file.open(dirname + fname, std::ios_base::binary);
  for (unsigned i = 0; i < size1; i++)
    for (unsigned j = 0; j < size2; j++)
      file.write((char*)(&table[i][j]), sizeof(unsigned));
  
  file.close();
  return;
}

void write_table(unsigned*** table, int size1, int size2, int size3, std::string fname){

  // Write data to file
  std::ofstream file;
  file.open(dirname + fname, std::ios_base::binary);
  for (unsigned i = 0; i < size1; i++)
    for (unsigned j = 0; j < size2; j++)
      for (unsigned k = 0; k < size3; k++)
	file.write((char*)(&table[i][j][k]), sizeof(unsigned));
  
  file.close();
  return;
}

/////////////////
// Loading
////////////////

unsigned* read_table(int size, std::string fname){

  // malloc.
  unsigned* table = new unsigned[size];
  for (unsigned i = 0; i < size; i++)
    table[i] = 0;

  // Extract data from file.
  std::ifstream file;
  file.open(dirname + fname, std::ios_base::binary);
  file.read((char*)(table), sizeof(unsigned) * size);
  //for (unsigned i = 0; i < size; i++)
  //  file.read((char*)(&table[i]), sizeof(unsigned));
  
  file.close();
  return table;
}

int* read_table(int size, std::string fname, int why){

  // malloc.
  int* table = new int[size];
  for (int i = 0; i < size; i++)
    table[i] = 0;
  
  // Extract data from file
  std::ifstream file;
  file.open(dirname + fname, std::ios_base::binary);
  file.read((char*)(table), sizeof(int) * size);
  //for (unsigned i = 0; i < size; i++)
  
  file.close();
  return table;
}

unsigned** read_table(int size1, int size2, std::string fname){

  unsigned** table = new unsigned*[size1];
  for (unsigned i = 0; i < size1; i++) {
    table[i] = new unsigned[size2];
    for (unsigned j = 0; j < size2; j++)
      table[i][j] = 0;
  }
  
  // Extract data from file
  std::ifstream file;
  file.open(dirname + fname, std::ios_base::binary);
  for (unsigned i = 0; i < size1; i++)
    file.read((char*)(table[i]), sizeof(unsigned) * size2);
    // for (unsigned j = 0; j < size2; j++)
    //   file.read((char*)(&table[i][j]), sizeof(unsigned));
  
  file.close();
  return table;
}

unsigned*** read_table(int size1, int size2, int size3, std::string fname){

  // Extract data from file.
  std::ifstream file;
  int NN = size1 * size2 * size3;
  unsigned* table_as_array = read_table(NN, fname);
  file.close();
  
  // Set up the pointer table.
  unsigned*** table = new unsigned**[size1];
  int M = 0;
  for (unsigned i = 0; i < size1; i++) {
    table[i] = new unsigned*[size2];
    for (unsigned j = 0; j < size2; j++) {
      table[i][j] = (unsigned *)(table_as_array + M);
      M += size3;
    }
  }
  return table;
}

