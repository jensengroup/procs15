// test_camshift.cpp --- Test of camshift energy class
// Copyright (C) 2011 Anders Steen Christensen
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
          int lastindex = (*pdb).find_last_of("."); 
          std::string rawname = (*pdb).substr(0, lastindex);
          std::string filename = rawname + ".procs";

          FPtable prediction = (*procs15).predict(chain_loop);
          (*procs15).write_nmr_star_format(chain_loop,filename, prediction);
     }
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

          if (arg_name == std::string("load-CA-Chemical-shifts")) {
               i++;
               settings.load_ca = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("load-CB-Chemical-shifts")) {
               i++;
               settings.load_cb = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("load-CO-Chemical-shifts")) {
               i++;
               settings.load_co = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("load-N-Chemical-shifts")) {
               i++;
               settings.load_n = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("load-HN-Chemical-shifts")) {
               i++;
               settings.load_hn = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("load-HA-Chemical-shifts")) {
               i++;
               settings.load_ha = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-ha-rc")) {
               i++;
               settings.include_ha_rc = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-hn-rc")) {
               i++;
               settings.include_hn_rc = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-hn-previous")) {
               i++;
               settings.include_hn_previous_residue_correction = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-ha-previous")) {
               i++;
               settings.include_ha_previous_residue_correction = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-ca-previous")) {
               i++;
               settings.include_ca_previous_residue_correction = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-cb-previous")) {
               i++;
               settings.include_cb_previous_residue_correction = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-co-previous")) {
               i++;
               settings.include_co_previous_residue_correction = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-hn-following")) {
               i++;
               settings.include_hn_following_residue_correction = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-ha-following")) {
               i++;
               settings.include_ha_following_residue_correction = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-ca-following")) {
               i++;
               settings.include_ca_following_residue_correction = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-cb-following")) {
               i++;
               settings.include_cb_following_residue_correction = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-co-following")) {
               i++;
               settings.include_co_following_residue_correction = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-n-following")) {
               i++;
               settings.include_n_following_residue_correction = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-ha-primary-hn-hbond")) {
               i++;
               settings.include_ha_primary_hn_hbond = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-ca-primary-hn-hbond")) {
               i++;
               settings.include_ca_primary_hn_hbond = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-co-primary-hn-hbond")) {
               i++;
               settings.include_co_primary_hn_hbond = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-n-primary-hn-hbond")) {
               i++;
               settings.include_n_primary_hn_hbond = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-hn-primary-ha-hbond")) {
               i++;
               settings.include_hn_primary_ha_hbond = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-co-primary-ha-hbond")) {
               i++;
               settings.include_co_primary_ha_hbond = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-hn-secondary-hn-hbond")) {
               i++;
               settings.include_hn_secondary_hn_hbond = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-ha-secondary-hn-hbond")) {
               i++;
               settings.include_ha_secondary_hn_hbond = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-ca-secondary-hn-hbond")) {
               i++;
               settings.include_ca_secondary_hn_hbond = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-hn-secondary-ha-hbond")) {
               i++;
               settings.include_hn_secondary_ha_hbond = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-ha-secondary-ha-hbond")) {
               i++;
               settings.include_ha_secondary_ha_hbond = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-ca-secondary-ha-hbond")) {
               i++;
               settings.include_ca_secondary_ha_hbond = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-co-secondary-ha-hbond")) {
               i++;
               settings.include_co_secondary_ha_hbond = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("include-n-secondary-ha-hbond")) {
               i++;
               settings.include_n_secondary_ha_hbond = (bool)atoi(argv[i]);
          }
          else if (arg_name == std::string("HN-water-bond-correction")) {
               i++;
               settings.use_water_correction = (bool)atoi(argv[i]);
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

     if (argc < 1) {
         std::cerr << "Usage: ./procs15_predictor < --setting=value > --pdb-files file1.pdb file2.pdb"
                   << "\n(fileX.pdb is PDB structure to predict) (<location of data files> path to numpy datafiles) (<setting=value> is procs15 settings, e.g. --include-ha-rc 0)  "
                   << "\nMemory usage for predictions with all atom types is approximately 29 GB of RAM." << std::endl;
          exit(1);
     }

     std::string numpy_path = argv[1];
     std::vector< std::string > pdbfiles;
     TermProCS15::Settings settings;

     get_command_line(argc, argv, settings, pdbfiles);

     run_predictor(pdbfiles, settings);
}
