// Copyright (C) 2014-2015 by Anders S. Christensen, Anders S. Larsen, Lars A. Bratholm
//
// This file is part of PHAISTOS
//
// PHAISTOS free software: you can redistribute it and/or modify
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


#ifndef TERM_PROCS15
#define TERM_PROCS15

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream> 

#include "term_procs15_base.h"
#include "procs15_backend.h"


using namespace procs15;

namespace phaistos {

class TermProCS15: public TermProCS15Base<TermProCS15> {

public:

     // For convenience, define local EnergyTermCommon
     typedef phaistos::TermProCS15Base<TermProCS15> TermProCS15Base;

     // Use settings from base class
     typedef TermProCS15Base::Settings Settings;

     //!Constructor
     TermProCS15(ChainFB *chain, const Settings &settings=Settings(),
             RandomNumberEngine *random_number_engine = &random_global)
          : TermProCS15Base(chain, "procs15", settings, random_number_engine) {
     }

     //!Copy constructor
     TermProCS15(const TermProCS15 &other,
                 RandomNumberEngine *random_number_engine,
                 int thread_index, ChainFB *chain)
          : TermProCS15Base(other, random_number_engine, thread_index, chain) {

     }

     FPlist get_primary_h_bond_corrections(ResidueIterator<ChainFB> &res1, phaistos::ChainFB& chain) {

          //! PRIMARY H-BOND res1 is the donor 
          FPlist h_bond_corrections = vector_utils::make_vector<FPtype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

          for (ResidueIterator<ChainFB> res2(chain); !(res2).end(); ++res2) {
               // Hydrogen bonds to neighbours are included in tripeptide dihedral scan
               if ( abs(res1->index - res2->index) < 2 ) {
                    continue;
               }
               std::vector<HBondAcceptor> acceptors = get_acceptors(*res2, chain);
               if (abs(res1->index - res2->index) == 2){
                    //remove backbone O
                    acceptors.erase(acceptors.begin()); //TODO Due to standard angles, this should still be included. Maybe subtract correction for standard angles.
               }
               for (std::vector<HBondAcceptor>::iterator acceptor = acceptors.begin(); acceptor != acceptors.end(); acceptor++){
                    // Amide hydrogen donor
                    if (res1->has_atom(H)) {
                         const FPtype r_oh = ((*res1)[H]->position - acceptor->acceptor_oxygen->position).norm();
                         if (r_oh < hn_bond_cutoff) {
                              switch (acceptor->acceptor_type){
                                   default:
                                        break;
                                   case AcceptorAmide: {

                                        const Vector_3D h_pos = (*res1)[H]->position;
                                        const Vector_3D o_pos = acceptor->acceptor_oxygen->position;
                                        const Vector_3D a_2 = acceptor->acceptor_second_atom->position;
                                        const Vector_3D a_3 = acceptor->acceptor_third_atom->position;

                                        const FPtype rho = calc_dihedral(h_pos, o_pos, a_2, a_3);
                                        const FPtype theta = calc_angle(h_pos, o_pos, a_2);

                                        const boost::multi_array<float, 1> data_vector =  get_bb_bb_data_vector(r_oh, theta, rho);


                                        for (unsigned int i = 0; i < 6; i++) {
                                             h_bond_corrections[i] += (FPtype)data_vector[i];
                                        }
                                        break;
                                   }

                                   case AcceptorAlcohol: {

                                        const Vector_3D h_pos = (*res1)[H]->position;
                                        const Vector_3D o_pos = acceptor->acceptor_oxygen->position;
                                        const Vector_3D a_2 = acceptor->acceptor_second_atom->position;
                                        const Vector_3D a_3 = acceptor->acceptor_third_atom->position;

                                        const FPtype rho = calc_dihedral(h_pos, o_pos, a_2, a_3);
                                        const FPtype theta = calc_angle(h_pos, o_pos, a_2);

                                        const boost::multi_array<float, 1> data_vector = get_bb_sc_alcohol_data_vector(r_oh, theta, rho);

                                        for (unsigned int i = 0; i < 6; i++) {
                                             h_bond_corrections[i] += (FPtype)data_vector[i];
                                        }

                                        break;
                                   }
                                   case AcceptorCarboxylate: {

                                        const Vector_3D h_pos = (*res1)[H]->position;
                                        const Vector_3D o_pos = acceptor->acceptor_oxygen->position;
                                        const Vector_3D a_2 = acceptor->acceptor_second_atom->position;
                                        const Vector_3D a_3 = acceptor->acceptor_third_atom->position;

                                        const FPtype rho = calc_dihedral(h_pos, o_pos, a_2, a_3);
                                        const FPtype theta = calc_angle(h_pos, o_pos, a_2);

                                        const boost::multi_array<float, 1> data_vector = get_bb_sc_carboxy_data_vector(r_oh, theta, rho);

                                        for (unsigned int i = 0; i < 6; i++) {
                                             h_bond_corrections[i] += (FPtype)data_vector[i];
                                        }

                                        break;
                                   }
                              } //end switch


                         }
                    }

                    //! Alpha hydrogen donor
                    FPtype r_oha = ((*res1)[get_ha_atom_type(res1->residue_type)]->position - acceptor->acceptor_oxygen->position).norm();
                    if (r_oha < ha_bond_cutoff) { 
                         switch (acceptor->acceptor_type){ //! acceptor type switch
                              default:
                                   break;
                              case AcceptorAmide: {
                                   const Vector_3D h_pos = (*res1)[get_ha_atom_type(res1->residue_type)]->position;
                                   const Vector_3D o_pos = acceptor->acceptor_oxygen->position;
                                   const Vector_3D a_2 = acceptor->acceptor_second_atom->position;
                                   const Vector_3D a_3 = acceptor->acceptor_third_atom->position;

                                   const FPtype rho = calc_dihedral(h_pos, o_pos, a_2, a_3);
                                   const FPtype theta = calc_angle(h_pos, o_pos, a_2);

                                   const boost::multi_array<float, 1> data_vector = get_bb_bb_alphah_oxygen_data_vector(r_oha,theta,rho);

                                   for (unsigned int i = 0; i < 6; i++) {
                                        h_bond_corrections[i+6] += (FPtype)data_vector[i];
                                   }
                              break;
                              }

                              case AcceptorAlcohol: {
                                   const Vector_3D h_pos = (*res1)[get_ha_atom_type(res1->residue_type)]->position;
                                   const Vector_3D o_pos = acceptor->acceptor_oxygen->position;
                                   const Vector_3D a_2 = acceptor->acceptor_second_atom->position;
                                   const Vector_3D a_3 = acceptor->acceptor_third_atom->position;

                                   const FPtype rho = calc_dihedral(h_pos, o_pos, a_2, a_3);
                                   const FPtype theta = calc_angle(h_pos, o_pos, a_2);
                                   const boost::multi_array<float, 1> data_vector = get_bb_bb_alphah_alcohol_data_vector(r_oha,theta,rho);

                                   for (unsigned int i = 0; i < 6; i++) {
                                        h_bond_corrections[i+6] += (FPtype)data_vector[i];
                                   }
                              break;
                              }

                              case AcceptorCarboxylate: {
                                   const Vector_3D h_pos = (*res1)[get_ha_atom_type(res1->residue_type)]->position;
                                   const Vector_3D o_pos = acceptor->acceptor_oxygen->position;
                                   const Vector_3D a_2 = acceptor->acceptor_second_atom->position;
                                   const Vector_3D a_3 = acceptor->acceptor_third_atom->position;

                                   const FPtype rho = calc_dihedral(h_pos, o_pos, a_2, a_3);
                                   const FPtype theta = calc_angle(h_pos, o_pos, a_2);
                                   const boost::multi_array<float, 1> data_vector = get_bb_bb_alphah_carboxy_data_vector(r_oha,theta,rho);

                                   for (unsigned int i = 0; i < 6; i++) {
                                        h_bond_corrections[i+6] += (FPtype)data_vector[i];
                                   }
                              break;
                              }
                         } //! end acceptor type switch

                    }
                    if (res1->residue_type == GLY) {
                         r_oha = ((*res1)[HA3]->position - acceptor->acceptor_oxygen->position).norm();
                         if (r_oha < ha_bond_cutoff) { 
                              switch (acceptor->acceptor_type){ //! acceptor type switch
                                   default:
                                        break;
                                   case AcceptorAmide: {
                                        const Vector_3D h_pos = (*res1)[HA3]->position;
                                        const Vector_3D o_pos = acceptor->acceptor_oxygen->position;
                                        const Vector_3D a_2 = acceptor->acceptor_second_atom->position;
                                        const Vector_3D a_3 = acceptor->acceptor_third_atom->position;

                                        const FPtype rho = calc_dihedral(h_pos, o_pos, a_2, a_3);
                                        const FPtype theta = calc_angle(h_pos, o_pos, a_2);

                                        const boost::multi_array<float, 1> data_vector = get_bb_bb_alphah_oxygen_data_vector(r_oha,theta,rho);


                                        //Don't include effect on HA2 for HA3 donor
                                        for (unsigned int i = 0; i < 5; i++) {
                                             h_bond_corrections[i+6] += (FPtype)data_vector[i];
                                        }
                                   break;
                                   }

                                   case AcceptorAlcohol: {
                                        const Vector_3D h_pos = (*res1)[HA3]->position;
                                        const Vector_3D o_pos = acceptor->acceptor_oxygen->position;
                                        const Vector_3D a_2 = acceptor->acceptor_second_atom->position;
                                        const Vector_3D a_3 = acceptor->acceptor_third_atom->position;

                                        const FPtype rho = calc_dihedral(h_pos, o_pos, a_2, a_3);
                                        const FPtype theta = calc_angle(h_pos, o_pos, a_2);
                                        const boost::multi_array<float, 1> data_vector = get_bb_bb_alphah_alcohol_data_vector(r_oha,theta,rho);

                                        //Don't include effect on HA2 for HA3 donor
                                        for (unsigned int i = 0; i < 5; i++) {
                                             h_bond_corrections[i+6] += (FPtype)data_vector[i];
                                        }
                                   break;
                                   }

                                   case AcceptorCarboxylate: {
                                        const Vector_3D h_pos = (*res1)[HA3]->position;
                                        const Vector_3D o_pos = acceptor->acceptor_oxygen->position;
                                        const Vector_3D a_2 = acceptor->acceptor_second_atom->position;
                                        const Vector_3D a_3 = acceptor->acceptor_third_atom->position;

                                        const FPtype rho = calc_dihedral(h_pos, o_pos, a_2, a_3);
                                        const FPtype theta = calc_angle(h_pos, o_pos, a_2);
                                        const boost::multi_array<float, 1> data_vector = get_bb_bb_alphah_carboxy_data_vector(r_oha,theta,rho);
                                        for (unsigned int i = 0; i < 5; i++) {
                                             h_bond_corrections[i+6] += (FPtype)data_vector[i];
                                        }
                                   break;
                                   }
                              } //! end acceptor type switch

                         }//end if
                    }//end if
               } //end acceptor iterator
          } //end Residue iterator

          return h_bond_corrections;
     }

     FPlist get_secondary_h_bond_corrections(ResidueIterator<ChainFB> &res1, phaistos::ChainFB& chain) {

          FPlist h_bond_corrections = vector_utils::make_vector<FPtype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
          //! SECONDARY H-BOND res1 is the acceptor 
          for (ResidueIterator<ChainFB> res2(chain); !(res2).end(); ++res2) {

               Residue *res3;

               if (res1->terminal_status != CTERM) {
                    res3 = res1->get_neighbour(1);

               } else {

                    break;
               }

               if (abs(res1->index - res2->index) < 2) {
                   continue;
               }

               std::vector<HBondDonor> donors = get_donors(*res2);
               //TODO Due to standard angles, this should still be included. Maybe subtract correction for standard angles.
               if (abs(res1->index - res2->index) == 2){ //!only side chain donor remove HN remove HA 
                    donors.erase(donors.begin()); 
                    if (res2->residue_type != PRO){
                         donors.erase(donors.begin()); 
                    }
                    // Remove HA3 as well
                    if (res2->residue_type == GLY){
                         donors.erase(donors.begin()); 
                    }
               }

               for (std::vector<HBondDonor>::iterator donor = donors.begin(); donor != donors.end(); donor++){
                    switch (donor->donor_type){
                         default: { //! DonorAmide DonorAmmonium DonorGuanidinium DonorImidazolium
                              const FPtype r_oh = ((*res1)[O]->position - donor->donor_hydrogen->position).norm();
                              if (r_oh < hn_bond_cutoff) {
                                   const Vector_3D h_pos = donor->donor_hydrogen->position;
                                   const Vector_3D o_pos = (*res1)[O]->position;
                                   const Vector_3D c_pos = (*res1)[C]->position;
                                   const Vector_3D n_pos = (*res3)[N]->position;

                                   const FPtype rho = calc_dihedral(h_pos, o_pos, c_pos, n_pos);
                                   const FPtype theta = calc_angle(h_pos, o_pos, c_pos);

                                   //data_vector = get_bb_bb_data_vector(r_oh, theta, rho);
                                   const boost::multi_array<float, 1> data_vector = get_bb_bb_data_vector(r_oh, theta, rho);

                                   for (unsigned int i = 0; i < 6; i++) {
                                        h_bond_corrections[i] += (FPtype)(data_vector[i+6]);
                                   }
                              }
                              break;
                         }

                         case DonorAlphaHydrogen: {
                              const FPtype r_oha = ((*res1)[O]->position - donor->donor_hydrogen->position).norm();
                              if (r_oha < ha_bond_cutoff){

                                   const Vector_3D h_pos = donor->donor_hydrogen->position;
                                   const Vector_3D o_pos = (*res1)[O]->position;
                                   const Vector_3D c_pos = (*res1)[C]->position;
                                   const Vector_3D n_pos = (*res3)[N]->position;

                                   const FPtype rho = calc_dihedral(h_pos, o_pos, c_pos, n_pos);
                                   const FPtype theta = calc_angle(h_pos, o_pos, c_pos);

                                   const boost::multi_array<float, 1> data_vector = get_bb_bb_alphah_oxygen_data_vector(r_oha,theta,rho);

                                   for (unsigned int i = 0; i < 6; i++) {
                                        h_bond_corrections[i+6] += (FPtype)(data_vector[i+6]);
                                   }
                              }
                              break;
                         }
                    } //end switch
               } //end donor iterator
          } //end residue iterator
          return h_bond_corrections;
     }

     FPlist get_ring_current_corrections(ResidueIterator<ChainFB> &res1, phaistos::ChainFB& chain) {

          FPlist data_vector = vector_utils::make_vector<FPtype>(0.0, 0.0);

          for (ResidueIterator<ChainFB> res2(chain); !(res2).end(); ++res2) {

               if ( abs(res1->index - res2->index) < 2 ) {
                    continue;
               }

               std::vector<AromaticRing> aromatic_rings = get_aromatic_rings((*res2));

               for (std::vector<AromaticRing>::iterator ring = aromatic_rings.begin();
                    ring != aromatic_rings.end();
                    ++ring) {

                    FPtype dist_squared = (ring->ring_center - (*res1)[get_ha_atom_type(res1->residue_type)]->position).norm_squared();
                    if (dist_squared < rc_cutoff2) {
                         data_vector[0] += calc_aromatic_interaction((*res1)[get_ha_atom_type(res1->residue_type)]->position,
                                                                             (*ring),
                                                                             dist_squared);
                    }

                    if (res1->has_atom(H)) {
                         FPtype dist_squared = (ring->ring_center - (*res1)[H]->position).norm_squared();
                         if (dist_squared < rc_cutoff2) {
                              data_vector[1] += calc_aromatic_interaction((*res1)[H]->position,
                                                                                      (*ring),
                                                                                      dist_squared);
                         }
                    }
               }
          }
          return data_vector;
     }


     FPlist get_shieldings(ResidueIterator<ChainFB> &res1) {

          FPlist data_vector(18, 0.0);
          if ((res1->terminal_status == CTERM) || (res1->terminal_status == NTERM)) {
               return data_vector;
          }

          for (std::vector<unsigned int>::iterator atype = atypevector.begin(); atype != atypevector.end(); ++atype){ //loop over atoms types
               const int em1 = res1->get_neighbour(-1)->residue_type;
               const int e0 = res1->residue_type;
               const int ep1 = res1->get_neighbour(1)->residue_type;


               // check that data for the atom type is loaded, and that it's not CB for GLY or H for PRO.
               if ((loaded_tables[em1*6+*atype]) && (loaded_tables[e0*6+*atype]) && (e0*6+*atype!=31) && (e0*6+*atype!=76) && (loaded_tables[ep1*6+*atype])){

                    //! get correction from previous residue
                    const FPtype r1 = (list_of_tables[em1*6+*atype]->get_values(2,read_phi_psi_index(*res1->get_neighbour(-1)),get_chi(*res1->get_neighbour(-1)))
                                                                            -list_of_tables[0+*atype]->get_st(2));
                    //! get backbone term
                    const FPtype r2 = (list_of_tables[e0*6+*atype]->get_values(1,read_phi_psi_index(*res1),get_chi(*res1)));

                    //! get correction from next residue
                    const FPtype r3 = (list_of_tables[ep1*6+*atype]->get_values(0,read_phi_psi_index(*res1->get_neighbour(1)),get_chi(*res1->get_neighbour(1)))
                                                                            -list_of_tables[0+*atype]->get_st(0));


                    data_vector[*atype] = r1;
                    data_vector[*atype+6] = r2;
                    data_vector[*atype+12] = r3;

               }

          }

          return data_vector;
     }





     FPtable predict(phaistos::ChainFB& chain, MoveInfo *move_info=NULL) {

          using namespace phaistos;
          using namespace definitions;

          FPlist primary_h_bond_corrections;
          FPlist secondary_h_bond_corrections;
          FPlist ring_current_corrections;
          FPlist shieldings;

          FPtable chemical_shifts(chain.size(),empty_contribution);



          for (ResidueIterator<ChainFB> res1(chain); !(res1).end(); ++res1) {

               primary_h_bond_corrections = get_primary_h_bond_corrections(res1, chain);
               if ((settings.load_hn) && (settings.use_water_correction == true) && (std::abs(primary_h_bond_corrections[4]) < 0.000001) && (res1->residue_type != PRO)) {
                     chemical_shifts[res1->index][2] += water_bonding_correction;
               }
               secondary_h_bond_corrections = get_secondary_h_bond_corrections(res1, chain);
               ring_current_corrections = get_ring_current_corrections(res1, chain);
               shieldings = get_shieldings(res1);



               if (res1->terminal_status == NTERM) {
                   if (settings.include_hn_secondary_hn_hbond) chemical_shifts[res1->index+1][2] += secondary_h_bond_corrections[4];
                   if (settings.include_n_secondary_hn_hbond)   chemical_shifts[res1->index+1][3] += secondary_h_bond_corrections[3];

               if (settings.include_hn_secondary_ha_hbond) chemical_shifts[res1->index+1][2] += secondary_h_bond_corrections[4+6];
               if (settings.include_n_secondary_ha_hbond)  chemical_shifts[res1->index+1][3] += secondary_h_bond_corrections[3+6];
                   continue;
               }
               if (res1->terminal_status == CTERM) {
                   if (settings.include_co_primary_hn_hbond) chemical_shifts[res1->index-1][4] += primary_h_bond_corrections[2];
                   continue;
               }
               if (settings.include_ca_previous_residue_correction) chemical_shifts[res1->index][1] -= shieldings[0];
               if (settings.include_cb_previous_residue_correction) chemical_shifts[res1->index][5] -= shieldings[1];
               if (settings.include_co_previous_residue_correction) chemical_shifts[res1->index][4] -= shieldings[2];
               if (settings.include_hn_previous_residue_correction) chemical_shifts[res1->index][2] -= shieldings[4];
               if (settings.include_n_previous_residue_correction) chemical_shifts[res1->index][3]  -= shieldings[3];
               if (settings.include_ha_previous_residue_correction) chemical_shifts[res1->index][0] -= shieldings[5];

               chemical_shifts[res1->index][1] -= shieldings[0+6];
               chemical_shifts[res1->index][5] -= shieldings[1+6];
               chemical_shifts[res1->index][4] -= shieldings[2+6];
               chemical_shifts[res1->index][3] -= shieldings[3+6];
               chemical_shifts[res1->index][2] -= shieldings[4+6];
               chemical_shifts[res1->index][0] -= shieldings[5+6];

               if (settings.include_ca_following_residue_correction) chemical_shifts[res1->index][1] -= shieldings[0+12];
               if (settings.include_cb_following_residue_correction) chemical_shifts[res1->index][5] -= shieldings[1+12];
               if (settings.include_co_following_residue_correction) chemical_shifts[res1->index][4] -= shieldings[2+12];
               if (settings.include_hn_following_residue_correction) chemical_shifts[res1->index][2] -= shieldings[4+12];
               if (settings.include_n_following_residue_correction) chemical_shifts[res1->index][3]  -=  shieldings[3+12];
               if (settings.include_ha_following_residue_correction) chemical_shifts[res1->index][0] -= shieldings[5+12];

               if (settings.include_ca_primary_hn_hbond) chemical_shifts[res1->index][1] += primary_h_bond_corrections[0];
               if (settings.include_cb_primary_hn_hbond) chemical_shifts[res1->index][5] += primary_h_bond_corrections[1];
               if (settings.include_hn_primary_hn_hbond) chemical_shifts[res1->index][2] += primary_h_bond_corrections[4];
               if (settings.include_n_primary_hn_hbond)  chemical_shifts[res1->index][3] += primary_h_bond_corrections[3];
               if (settings.include_co_primary_hn_hbond) chemical_shifts[res1->index-1][4] += primary_h_bond_corrections[2];
               if (settings.include_ha_primary_hn_hbond) chemical_shifts[res1->index][0] += primary_h_bond_corrections[5];

               if (settings.include_ca_primary_ha_hbond) chemical_shifts[res1->index][1] += primary_h_bond_corrections[0+6];
               if (settings.include_cb_primary_ha_hbond) chemical_shifts[res1->index][5] += primary_h_bond_corrections[1+6];
               if (settings.include_co_primary_ha_hbond) chemical_shifts[res1->index][4] += primary_h_bond_corrections[2+6];
               if (settings.include_hn_primary_ha_hbond) chemical_shifts[res1->index][2] += primary_h_bond_corrections[4+6];
               if (settings.include_n_primary_ha_hbond)  chemical_shifts[res1->index][3] += primary_h_bond_corrections[3+6];
               if (settings.include_ha_primary_ha_hbond) chemical_shifts[res1->index][0] += primary_h_bond_corrections[5+6];

               if (settings.include_ca_secondary_hn_hbond) chemical_shifts[res1->index][1] += secondary_h_bond_corrections[0];
               if (settings.include_cb_secondary_hn_hbond) chemical_shifts[res1->index][5] += secondary_h_bond_corrections[1];
               if (settings.include_co_secondary_hn_hbond) chemical_shifts[res1->index][4] += secondary_h_bond_corrections[2];
               if (settings.include_hn_secondary_hn_hbond) chemical_shifts[res1->index+1][2] += secondary_h_bond_corrections[4];
               if (settings.include_n_secondary_hn_hbond)   chemical_shifts[res1->index+1][3] += secondary_h_bond_corrections[3];
               if (settings.include_ha_secondary_hn_hbond) chemical_shifts[res1->index][0] += secondary_h_bond_corrections[5];

               if (settings.include_ca_secondary_ha_hbond) chemical_shifts[res1->index][1] += secondary_h_bond_corrections[0+6];
               if (settings.include_cb_secondary_ha_hbond) chemical_shifts[res1->index][5] += secondary_h_bond_corrections[1+6];
               if (settings.include_co_secondary_ha_hbond) chemical_shifts[res1->index][4] += secondary_h_bond_corrections[2+6];
               if (settings.include_hn_secondary_ha_hbond) chemical_shifts[res1->index+1][2] += secondary_h_bond_corrections[4+6];
               if (settings.include_n_secondary_ha_hbond)  chemical_shifts[res1->index+1][3] += secondary_h_bond_corrections[3+6];
               if (settings.include_ha_secondary_ha_hbond) chemical_shifts[res1->index][0] += secondary_h_bond_corrections[5+6];

               if (settings.include_ha_rc) chemical_shifts[res1->index][0] += ring_current_corrections[0];
               if (settings.include_hn_rc) chemical_shifts[res1->index][2] += ring_current_corrections[1];

               //! zero chemical_shifts_table if only non bonded contribution
               for (unsigned int i = 0; i < 6; i++){
                         if (std::abs(shieldings[i+6]) <= 0.000001) {
                              chemical_shifts[res1->index][transform_id(i)] = 0;
                         }
               }
          //! zero non bonding terms for Proline HN
               if (res1->residue_type == PRO){
                    chemical_shifts[res1->index][2] = 0.0;
               }
          }

          //! zero chemical_shifts_table at N and C terminus
          for (unsigned int o = 0; o < 6; o++){
               chemical_shifts[0][o] = 0;
               chemical_shifts[chemical_shifts.size()-1][o] = 0; 
          }




          return chemical_shifts;
     }

     // Dummy function to avoid duplicate reject() code
     void reject_cache() {}

}; // end class TermProCS15

////! Observable specialization for TermProCS15
//template <>
//class Observable<TermProCS15>: public TermProCS15, public ObservableBase {
//
//public:
//
//     //! Local settings class.
//     const class Settings: public TermProCS15::Settings, public ObservableBase::Settings {
//     public:
//
//          //! Constructor. Defines default values for settings object.
//          Settings(){}
//
//          //! Output operator
//          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
//               o << static_cast<const TermProCS15::Settings>(settings);
//               o << static_cast<const ObservableBase::Settings>(settings);
//               return o;
//          }
//     } settings; //!< Local settings object·
//
//     //! Constructor.
//     //! \param energy_term TermCamshift energy term object
//     //! \param settings Local Settings object
//     //! \param reference_energy_function All observables have a pointer to a reference energy function which they can refer to.
//     Observable(const TermProCS15 &energy_term,
//                const ObservableBase::Settings &settings=ObservableBase::Settings(),
//                Energy<ChainFB> *reference_energy_function=NULL)
//          : TermProCS15(energy_term),
//            settings(dynamic_cast<const Settings&>(settings)) {
//     }
//
//     //! Copy Constructor.
//     //! \param other Source object from which copy is made
//     //! \param thread_index Index indicating in which thread|rank the copy exists
//     //! \param chain Molecule chain
//     Observable(const Observable &other, int thread_index, ChainFB *chain)
//          : TermProCS15(other, random_number_engine, thread_index, chain),
//            settings(other.settings) {
//     }
//
//
//     //! Chemical shift RMSD between two sets of chemical shifts
//     //! \param cs1_this A matrix containing chemical shifts
//     //! \param cs2_this A matrix containing chemical shifts
//     //! \return A vector containing RMSDs
//     FPlist calc_rmsds(const FPtable &cs1_this,
//                                    const FPtable &cs2_this) {
//
//          FPlist chi_sq = empty_contribution;
//          std::vector<int> n_bins = vector_utils::make_vector(0, 0, 0, 0, 0, 0);
//
//          for (unsigned int i = 0; i <  std::min(cs1_this.size(), cs2_this.size()); i++) {
//               for (unsigned int j = 0; j < 6; j++) {
//
//                   if ((std::fabs(cs1_this[i][j]) > 0.0001)
//                    && (std::fabs(cs2_this[i][j]) > 0.0001)) {
//
//                         const FPtype diff  = cs1_this[i][j] - cs2_this[i][j];
//                         chi_sq[j] += diff * diff;
//                         n_bins[j] += 1;
//                    }
//               }
//          }
//
//          FPlist rmsds = empty_contribution;
//
//          for (unsigned int j = 0; j < 6; j++) {
//               rmsds[j] = std::sqrt(chi_sq[j] / n_bins[j]);
//          }
//
//          return rmsds;
//
//     }
//
//
////     FPtable get_full_prediction_error(const FPtable &cs1_this, const FPtable &cs2_this) {
////
////          FPlist cs_vector;
////          FPtable prediction_errors(6,cs_vector);
////
////          for (unsigned int i = 0; i <  std::min(cs1_this.size(), cs2_this.size()); i++) {
////               for (unsigned int j = 0; j < 6; j++) {
////
////                   if ((std::fabs(cs1_this[i][j]) > 0.0001)
////                    && (std::fabs(cs2_this[i][j]) > 0.0001)
////                    && !(std::isnan(cs1_this[i][j]))
////                    && !(std::isnan(cs2_this[i][j]))) {
////
////                         const double diff = cs1_this[i][j] - cs2_this[i][j];
////                         prediction_errors[j].push_back(diff);
////                    }
////                   else {
////                       prediction_errors[j].push_back(0.0);
////                   }
////               }
////          }
////
////          return prediction_errors;
////     }
//
//
//     //! Make observation.
//     virtual std::string observe(MoveInfo *move_info=NULL,
//                                 PHAISTOS_LONG_LONG current_iteration=0,
//                                 bool register_only=false) {
//
//          using namespace procs15;
//          this->list_of_tables[4]->get_st(2);
//          //// Energy to be returned
//          //double energy = this->evaluate();
//
////        //  // Calculate new chemical shifts
////        //  this->predicted_chemical_shifts = this->predict(*(this->chain));
//
//          //// Calculate RMSDs
//          //FPlist rmsds = calc_rmsds(this->predicted_chemical_shifts,
//          //                          this->experimental_chemical_shifts);
//
//          // Output stream
//          std::stringstream s;
//          //s << std::fixed << std::setprecision(5) << energy << ":" << rmsds;
//
////          // Add full prediction error to output
////          if (settings.output_full_prediction_vector) {
////               FPtable prediction_error = get_full_prediction_error(this->predicted_chemical_shifts,
////                                                                    this->experimental_chemical_shifts);
////               s << ":" << prediction_error;
////          }
//
//          return s.str();
//
//     }
//
//}; // End class Observable


} // end namespace phaistos

#endif // TERM_PROCS15
