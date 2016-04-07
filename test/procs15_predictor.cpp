// procs15_predictor.cpp --- Generates BMRB star file from input pdb-files
// Copyright (C) 2016 Anders Steen Christensen, Anders S. Larsen, Lars Andersen Bratholm
//
// This file is part of Phaistos
//
// Phaistos is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Phaistos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Phaistos.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <string.h>
#include <fstream>

#include <boost/type_traits/is_base_of.hpp>

#include "protein/chain_fb.h"
#include "protein/definitions.h"
#include "protein/iterators/pair_iterator_chaintree.h"

#include "energy/energy.h"
#include "energy/observable.h"
#include "energy/observable_collection.h"
#include "energy/term_procs15.h"

using namespace phaistos;


void run_predictor(std::vector< std::string > &pdbfiles, TermProCS15::Settings &settings) {

     // Create chain from PDB filename
     ChainFB chain(pdbfiles[0],definitions::ALL_ATOMS);

     // Add atoms missing in the pdb structure
     chain.add_atoms(definitions::ALL_PHYSICAL_ATOMS);

     using namespace std;
     // generate dummy star file
     ofstream  myfile;
     myfile.open ("temp.str", ios::out);
     for (unsigned int i = 0; i < chain.size(); i++){
          myfile << (i+1) << " " << (i+1) << " " << chain[i].residue_type << " " << "CA" << " " << "C" << " " << "56.0" << "\n";
     }

     myfile.close();

     settings.star_filename = "temp.str";

     TermProCS15* procs15 = new TermProCS15(&chain, settings); 
     for (std::vector< std::string >::iterator pdb = pdbfiles.begin(); pdb != pdbfiles.end(); ++pdb){
          std::cout << "start prediction " << (*pdb) << std::endl;
          // Create chain from PDB filename
          ChainFB chain_loop((*pdb),definitions::ALL_ATOMS);
          // Add atoms missing in the pdb structure
          chain_loop.add_atoms(definitions::ALL_PHYSICAL_ATOMS);

          //change extension for predictions
          int lastindex_dot = (*pdb).find_last_of(".");
          std::string rawname = (*pdb).substr(0, lastindex_dot);
          std::string filename = rawname + ".procs";

          if (!(settings.debug)) {
               FPtable prediction = (*procs15).predict(chain_loop);
               (*procs15).write_nmr_star_format(chain_loop,filename, prediction);
          } else {
               std::vector<FPtable> prediction = (*procs15).predict_full(chain_loop);
               (*procs15).write_nmr_star_format_full(chain_loop,filename, prediction);
          }
     }
     remove("temp.str");
}

// Parse the command line
void get_command_line(int &argc, char *argv[], 
                      TermProCS15::Settings &settings, std::vector< std::string > &pdbfiles) {

     for (int i = 1; i < argc; i++) {
          std::string arg = argv[i];
          std::string arg_name = "";

          // keyword or argument ?
          if (arg.substr(0, 2) == std::string("--")) {
               arg_name = arg.substr(2, arg.length() - 2);
          } else if (arg.substr(0, 1) == std::string("-")) {
               arg_name = arg.substr(1, arg.length() - 1);
          } else {
               continue;
          }

          if (arg_name == std::string("load-ca")) {
               i++;
               settings.load_ca = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("load-cb")) {
               i++;
               settings.load_cb = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("load-co")) {
               i++;
               settings.load_co = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("load-n")) {
               i++;
               settings.load_n = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("load-h")) {
               i++;
               settings.load_hn = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("load-ha")) {
               i++;
               settings.load_ha = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("details")) {
               i++;
               settings.debug = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("data-folder")) {
               i++;
               settings.procsnumpypath = std::string(argv[i]);
          }
          else if (arg_name == std::string("pdb-files")) {
               while((i+1 < argc) && (std::string(argv[i+1]).substr(std::string(argv[i+1]).find_last_of(".") + 1) == "pdb")) {
                    i++;
                    pdbfiles.push_back(std::string(argv[i]));
               }
          }
     }
}



int main(int argc, char *argv[]) {

     if (argc < 2) {
         std::cerr << "Usage: ./procs15_predictor --data-folder <path to numpy datafiles> --pdb-files file1.pdb file2.pdb ..."
                   << "\nfileX.pdb is PDB structures to predict. Prediction of an atom type can by disabled by"
                   << " including \n--load-Y 0, \nwhere Y can be (ha, ca, cb, co, h, n)"
                   << "\nDetailed output can be printed with --details 1"
                   << "\nMemory usage for predictions with all atom types is approximately 29 GB of RAM." << std::endl;
          exit(1);
     }

     std::vector< std::string > pdbfiles;
     TermProCS15::Settings settings;

     get_command_line(argc, argv, settings, pdbfiles);

     run_predictor(pdbfiles, settings);
}
