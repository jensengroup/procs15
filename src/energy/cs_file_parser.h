// Copyright (C) 2014 by Anders Steen Christensen
//
// This file is part of PHAISTOS
//
// PHAISTOS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// PHAISTOS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with PHAISTOS.  If not, see <http://www.gnu.org/licenses/>.
//
#include <iostream>

#include "lib_definitions.h"

using namespace procs15;

namespace cs_parser {

     //! The NMR STAR format file parser 
     //! \param star_filename Name of NMR STAR datafile
     //! \return value_matrix A list of chemical shifts ordered by residue number
     std::vector< std::vector<FPtype> > value_matrix_from_starfile(const std::string& star_filename, 
                                                                   phaistos::ChainFB& chain) {

          using namespace phaistos;
          using namespace definitions;

          FPtable value_matrix(chain.size(), empty_contribution);
          std::vector<std::string> lines = file_to_string_vector(star_filename);

          for (unsigned int i=0; i<lines.size(); ++i) {
               std::string line = lines[i];
               // Remove white space at beginning and end
               boost::trim(line);
               if (line.size() == 0)
                    continue;
               std::vector<std::string> tokens;
               boost::split(tokens, line, boost::is_any_of(" \t"), boost::token_compress_on);

               int residue_index = boost::lexical_cast<int>(tokens[1]);
               std::string res_type = tokens[2];
               FPtype value = boost::lexical_cast<FPtype>(tokens[5]);
               std::string atom_type = tokens[3];
               // force HA2 on GLY
               if ((atom_type == "HA2") && (res_type == "GLY")){
                   atom_type = "HA";
               }

               if (residue_index > chain.size()) {
                   std::cerr << "ProCS15 WARNING: Data "<< line << " not found in chain" << std::endl;
                   continue;
               } 

               ResidueEnum residue_type_chain = chain[residue_index - 1].residue_type;

               bool contains_bugs = false;
               if ((residue_type_chain == definitions::ALA) && (res_type != "ALA")) contains_bugs = true;
               if ((residue_type_chain == definitions::CYS) && (res_type != "CYS")) contains_bugs = true;
               if ((residue_type_chain == definitions::ASP) && (res_type != "ASP")) contains_bugs = true;
               if ((residue_type_chain == definitions::GLU) && (res_type != "GLU")) contains_bugs = true;
               if ((residue_type_chain == definitions::PHE) && (res_type != "PHE")) contains_bugs = true;
               if ((residue_type_chain == definitions::GLY) && (res_type != "GLY")) contains_bugs = true;
               if ((residue_type_chain == definitions::HIS) && (res_type != "HIS")) contains_bugs = true;
               if ((residue_type_chain == definitions::ILE) && (res_type != "ILE")) contains_bugs = true;
               if ((residue_type_chain == definitions::LYS) && (res_type != "LYS")) contains_bugs = true;
               if ((residue_type_chain == definitions::LEU) && (res_type != "LEU")) contains_bugs = true;
               if ((residue_type_chain == definitions::MET) && (res_type != "MET")) contains_bugs = true;
               if ((residue_type_chain == definitions::ASN) && (res_type != "ASN")) contains_bugs = true;
               if ((residue_type_chain == definitions::PRO) && (res_type != "PRO")) contains_bugs = true;
               if ((residue_type_chain == definitions::GLN) && (res_type != "GLN")) contains_bugs = true;
               if ((residue_type_chain == definitions::ARG) && (res_type != "ARG")) contains_bugs = true;
               if ((residue_type_chain == definitions::SER) && (res_type != "SER")) contains_bugs = true;
               if ((residue_type_chain == definitions::THR) && (res_type != "THR")) contains_bugs = true;
               if ((residue_type_chain == definitions::VAL) && (res_type != "VAL")) contains_bugs = true;
               if ((residue_type_chain == definitions::TRP) && (res_type != "TRP")) contains_bugs = true;
               if ((residue_type_chain == definitions::TYR) && (res_type != "TYR")) contains_bugs = true;

               if (contains_bugs) {
                    std::cerr << "ProCS15 ERROR: CS-data mismatch in " << residue_type_chain << residue_index - 1 << std::endl;
                    std::cerr << "ProCS15 ERROR: Trying to parse line:" << line << std::endl;
                    exit(1);
               }


               if (atom_type == "CA") {
                    value_matrix[residue_index - 1][1] = value;
               } else if (atom_type == "CB") {
                    value_matrix[residue_index - 1][5] = value;
               } else if (atom_type == "C") {
                    value_matrix[residue_index - 1][4] = value;
               } else if (atom_type == "N") {
                    value_matrix[residue_index - 1][3] = value;
               } else if (atom_type == "HA") {
                    value_matrix[residue_index - 1][0] = value;
               } else if (atom_type == "H") {
                    value_matrix[residue_index - 1][2] = value;
               } else {
                      std::cerr << "ProCS15 WARNING - NMR-STAR parser: Skipping unknown atom type " << atom_type << std::endl;
               }
          }

          return value_matrix;
     }

} //End namespace camshift_parser
