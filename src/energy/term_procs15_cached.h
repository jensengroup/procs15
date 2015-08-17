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

#ifndef TERM_PROCS15_CACHED
#define TERM_PROCS15_CACHED


#include "term_procs15_base.h"
#include "procs15_backend.h"


using namespace lib_definitions;

namespace phaistos {

class TermProCS15Cached: public TermProCS15Base<TermProCS15Cached>,
                           public ProCS15Backend {
private:

     std::vector<bool> loaded_data;

     struct HbPair {

         Atom* hb_don;
         Atom* hb_acc;
         Atom* hb_acc_second;
         Atom* hb_acc_third;

         unsigned int other_don_index; //sidechain donors or HA donors
         unsigned int other_acc_index; //sidechain acceptors

         AcceptorEnum acceptor_type;
         DonorEnum donor_type;
     };

     struct RingHydrogenPair {

         unsigned int hydrogen_residue_index;
         unsigned int ring_residue_index;
         unsigned int ring_nr;

         bool HAatom;
         bool HNatom;

         Atom* HA_atom;
         Atom* HN_atom; 

     };


     //! Put backbone donor in the vector
     std::vector< Atom* > bb_don_atoms;

     // It was easier and faster to make a special case for HA3 instead of including it in bb_don_atoms
     //! Put HA3 donor in the vector
     std::vector< Atom* > ha3_don_atoms;

     //! Put backbone acceptor in the vector
     std::vector< Atom* > bb_acc_atoms;

     //! Put sidechain donor in the vector
     std::vector< Atom* > sc_don_atoms;

     //! Put sidechain acceptor in the vector
     std::vector< Atom* > sc_acc_atoms;

     std::vector< std::vector<  std::vector< HbPair > > > bbdon_bbacc_pairs;
     std::vector< std::vector<  std::vector< HbPair > > > bbdon_scacc_pairs;
     std::vector< std::vector<  std::vector< HbPair > > > scdon_bbacc_pairs;
     std::vector< std::vector<  std::vector< HbPair > > > ha3don_bbacc_pairs;
     std::vector< std::vector<  std::vector< HbPair > > > ha3don_scacc_pairs;

     FPtable bbdon_bbacc_backup;
     FPtable bbdon_scacc_backup;
     FPtable ha3don_bbacc_backup;
     FPtable ha3don_scacc_backup;
     FPtable scdon_bbacc_backup;

     FPtable bbdon_bbacc;
     FPtable bbdon_scacc;
     FPtable ha3don_bbacc;
     FPtable ha3don_scacc;
     FPtable scdon_bbacc;

     std::vector< std::vector<  std::vector< RingHydrogenPair > > > ring_hydrogen_pairs;

     FPtable chemical_shifts_table;
     FPtable non_bonded_chemical_shifts_table;


     std::vector< unsigned int > amide_bond_count_table;


     void find_hydrogen_ring_pairs(){

          //! find hydrogen-ring pairs for Ring current cached
          for (unsigned int start_index = 0; start_index < chain->size(); start_index++) {

               std::vector<  std::vector< RingHydrogenPair > > ringh_pairs_start_index(0);

               for (unsigned int end_index = 0; end_index < chain->size(); end_index++) {

                    std::vector< RingHydrogenPair > ringh_pairs_end_index(0);

                    for (ResidueIterator<ChainFB> hydrogen_residue(*this->chain); !hydrogen_residue.end(); ++hydrogen_residue){

                         for (ResidueIterator<ChainFB> ring_residue(*this->chain); !ring_residue.end(); ++ring_residue){

                              std::vector<AromaticRing> aromatic_rings = get_aromatic_rings((*ring_residue));

                              //! Hydrogens outside modified region
                              if  ( ((hydrogen_residue->index < start_index) || (hydrogen_residue->index > end_index)) &&
                                    ((   ring_residue->index >= start_index) && (ring_residue->index <= end_index)   ) ){

                                   if (abs(hydrogen_residue->index - ring_residue->index) > 1){

                                        unsigned int ring_nr = 0;

                                        for (std::vector<AromaticRing>::iterator ring = aromatic_rings.begin();
                                             ring != aromatic_rings.end();
                                             ++ring) {

                                             RingHydrogenPair pair;
                                             pair.hydrogen_residue_index = hydrogen_residue->index;
                                             pair.ring_residue_index = ring_residue->index;
                                             pair.ring_nr = ring_nr;

                                             // HA (HA2 for GLY)
                                             if ( hydrogen_residue->has_atom(get_ha_atom_type(hydrogen_residue->residue_type)) ){
                                                  pair.HAatom = true;
                                                  pair.HA_atom = (*hydrogen_residue)[get_ha_atom_type(hydrogen_residue->residue_type)];
                                             }
                                             else{
                                                  pair.HAatom = false;
                                                 }
                                             // H
                                             if ( hydrogen_residue->has_atom(H) ){
                                                  pair.HNatom = true;
                                                  pair.HN_atom = (*hydrogen_residue)[H];
                                             }
                                             else{
                                                  pair.HNatom = false;
                                             }
                                             ringh_pairs_end_index.push_back( pair );
                                             ring_nr += 1;


                                        }
                                   }
                              }
                         }

                         //! Hydrogens inside modified region 
                         if ((hydrogen_residue->index >= start_index) && (hydrogen_residue->index <= end_index)) {

                              for (ResidueIterator<ChainFB> ring_residue2(*this->chain); !ring_residue2.end(); ++ring_residue2){

                                   if (abs(hydrogen_residue->index - ring_residue2->index) > 1){

                                        std::vector<AromaticRing> aromatic_rings2 = get_aromatic_rings((*ring_residue2));
                                        int ring_nr2 = 0;
                                        for (std::vector<AromaticRing>::iterator ring = aromatic_rings2.begin();
                                                         ring != aromatic_rings2.end();
                                                         ++ring) {

                                             RingHydrogenPair pair;
                                             pair.hydrogen_residue_index = hydrogen_residue->index;
                                             pair.ring_residue_index = ring_residue2->index;
                                             pair.ring_nr = ring_nr2;

                                             if ( hydrogen_residue->has_atom(get_ha_atom_type(hydrogen_residue->residue_type)) ){
                                                  pair.HAatom = true;
                                                  pair.HA_atom = (*hydrogen_residue)[get_ha_atom_type(hydrogen_residue->residue_type)];
                                             }
                                             else{
                                                  pair.HAatom = false;
                                             }
                                             if ( hydrogen_residue->has_atom(H) ){
                                                  pair.HNatom = true;
                                                  pair.HN_atom = (*hydrogen_residue)[H];
                                             }
                                             else{
                                                  pair.HNatom = false;
                                             }
                                             ringh_pairs_end_index.push_back( pair );

                                             ring_nr2 += 1;
                                             }
                                   }
                              }
                         }
                    }
              ringh_pairs_start_index.push_back( ringh_pairs_end_index );
              }
          ring_hydrogen_pairs.push_back( ringh_pairs_start_index );
          }
     }




     std::vector< RingHydrogenPair > get_hydrogen_ring_pairs( int start_index, int end_index) {

          std::vector< RingHydrogenPair > modified_pairs;

          for (ResidueIterator<ChainFB> hydrogen_residue(*this->chain); !hydrogen_residue.end(); ++hydrogen_residue){

               for (ResidueIterator<ChainFB> ring_residue(*this->chain); !ring_residue.end(); ++ring_residue){

                    std::vector<AromaticRing> aromatic_rings = get_aromatic_rings((*ring_residue));

                    //! Hydrogens outside modified region
                    if  ( ((hydrogen_residue->index < start_index) || (hydrogen_residue->index > end_index)) &&
                          ((   ring_residue->index >= start_index) && (ring_residue->index <= end_index)   ) ){

                         if (abs(hydrogen_residue->index - ring_residue->index) > 1){

                              int ring_nr = 0;

                              for (std::vector<AromaticRing>::iterator ring = aromatic_rings.begin();
                                   ring != aromatic_rings.end();
                                   ++ring) {

                                   RingHydrogenPair pair;
                                   pair.hydrogen_residue_index = hydrogen_residue->index;
                                   pair.ring_residue_index = ring_residue->index;
                                   pair.ring_nr = ring_nr;

                                   if ( hydrogen_residue->has_atom(get_ha_atom_type(hydrogen_residue->residue_type)) ){
                                        pair.HAatom = true;
                                        pair.HA_atom = (*hydrogen_residue)[get_ha_atom_type(hydrogen_residue->residue_type)]; 
                                   }
                                   else pair.HAatom = false;


                                   if ( hydrogen_residue->has_atom(H)){
                                        pair.HNatom = true;
                                        pair.HN_atom = (*hydrogen_residue)[H];
                                   }
                                   else pair.HNatom = false;
                                   modified_pairs.push_back( pair );
                                   ring_nr += 1;


                              }
                         }
                    }
               }

               //! Hydrogens inside modified region 
               if  ((hydrogen_residue->index >= start_index) && (hydrogen_residue->index <= end_index)) {


                    for (ResidueIterator<ChainFB> ring_residue2(*this->chain); !ring_residue2.end(); ++ring_residue2){

                         if (abs(hydrogen_residue->index - ring_residue2->index) > 1){

                               std::vector<AromaticRing> aromatic_rings2 = get_aromatic_rings((*ring_residue2));
                               int ring_nr2 = 0;
                               for (std::vector<AromaticRing>::iterator ring = aromatic_rings2.begin();
                                                ring != aromatic_rings2.end();
                                                ++ring) {

                                    RingHydrogenPair pair;
                                    pair.hydrogen_residue_index = hydrogen_residue->index;
                                    pair.ring_residue_index = ring_residue2->index;
                                    pair.ring_nr = ring_nr2;

                                    if ( hydrogen_residue->has_atom(get_ha_atom_type(hydrogen_residue->residue_type)) ){
                                         pair.HAatom = true;
                                         pair.HA_atom = (*hydrogen_residue)[get_ha_atom_type(hydrogen_residue->residue_type)]; 
                                    }
                                    else pair.HAatom = false;


                                    if ( hydrogen_residue->has_atom(H) ){
                                         pair.HNatom = true;
                                         pair.HN_atom = (*hydrogen_residue)[H];
                                    }
                                    else pair.HNatom = false;
                                    modified_pairs.push_back( pair );

                                    ring_nr2 += 1;
                               }
                         }
                    }
               }
          }
     return modified_pairs;
     }

     // hbonds between mainchain atoms
     std::vector< HbPair > get_modified_bbbb_hbpair(int start_index, int end_index){

          std::vector< HbPair > modified_pairs;

          for (unsigned int i = 0; i< bb_don_atoms.size(); i++ ) {
               Residue *res_don = bb_don_atoms[i]->residue;

               for(unsigned int j = 0; j< bb_acc_atoms.size(); j++ ) {

                   Residue *res_acc = bb_acc_atoms[j]->residue;

                   if (((res_acc->index >= start_index) &&
                        (res_acc->index <= end_index)) ||
                       ((res_don->index >= start_index) &&
                        (res_don->index <= end_index))) {

                         if(abs(res_acc->index - res_don->index) > 2) {
                              switch (bb_don_atoms[i]->atom_type){
                                   // HN to O hbond
                                   case H: {
                                        HbPair pair;
                                        pair.hb_don = bb_don_atoms[i];
                                        pair.hb_acc = bb_acc_atoms[j];
                                        pair.hb_acc_second = (*res_acc)[C];
                                        pair.hb_acc_third  = (*chain)[res_acc->index+1][N];
                                        pair.acceptor_type = AcceptorAmide;
                                        pair.donor_type    = DonorAmide; 
                                        modified_pairs.push_back(pair);
                                        }
                                        break;
                                   // HA to O hbond
                                   case HA: case HA2: {

                                        HbPair pair;
                                        pair.hb_don = bb_don_atoms[i];
                                        pair.hb_acc = bb_acc_atoms[j];
                                        pair.hb_acc_second = (*res_acc)[C];
                                        pair.hb_acc_third  = (*chain)[res_acc->index+1][N];
                                        pair.acceptor_type = AcceptorAmide;
                                        pair.donor_type    = DonorAlphaHydrogen; 
                                        modified_pairs.push_back(pair);
                                        }
                                        break;

                                   default:
                                        std::cout << "(ProCS15Cached)critical fail in get_modified_bbbb_hbpair" << std::endl;
                                        break;
                              }
                         }
                    }
               }
          }

     return modified_pairs;

     }

     // ha3 - backbone amide hbonds
     std::vector< HbPair > get_modified_ha3bb_hbpair(int start_index, int end_index){

          std::vector< HbPair > modified_pairs;

          for (unsigned int i = 0; i< ha3_don_atoms.size(); i++ ) {
               Residue *res_don = ha3_don_atoms[i]->residue;

               for(unsigned int j = 0; j< bb_acc_atoms.size(); j++ ) {

                   Residue *res_acc = bb_acc_atoms[j]->residue;

                   if (((res_acc->index >= start_index) &&
                        (res_acc->index <= end_index)) ||
                       ((res_don->index >= start_index) &&
                        (res_don->index <= end_index))) {

                         if(abs(res_acc->index - res_don->index) > 2) {
                              HbPair pair;
                              pair.hb_don = ha3_don_atoms[i];
                              pair.hb_acc = bb_acc_atoms[j];
                              pair.hb_acc_second = (*res_acc)[C];
                              pair.hb_acc_third  = (*chain)[res_acc->index+1][N];
                              pair.acceptor_type = AcceptorAmide;
                              pair.donor_type    = DonorAlphaHydrogen; 
                              modified_pairs.push_back(pair);
                         }
                    }
               }
          }

     return modified_pairs;

     }

     // Hbonds between sidechain and mainchain
     std::vector< HbPair > get_modified_bbsc_hbpair(unsigned int start_index, unsigned int end_index){

          unsigned int acc_count;
          std::vector< HbPair > modified_pairs;
          for (unsigned int i = 0; i < bb_don_atoms.size(); i++ ) {

               acc_count = 0;
               Residue *res_don = bb_don_atoms[i]->residue;

               for(unsigned int j = 0; j< sc_acc_atoms.size(); j++ ) {
                   Residue *res_acc = sc_acc_atoms[j]->residue;

                   if (((res_acc->index >= start_index) &&
                        (res_acc->index <= end_index)) ||
                       ((res_don->index >= start_index) &&
                        (res_don->index <= end_index))) {
                         if ( j > 0 ){
                              if ((sc_acc_atoms[j]->residue->index) == sc_acc_atoms[j-1]->residue->index){

                                   acc_count += 1;
                              }
                              else{
                                   acc_count = 0;
                              }
                         }
                         if(abs(res_acc->index - res_don->index) > 1) {
                              switch (bb_don_atoms[i]->atom_type){
                                   case H: //! amide hydrogen donor
                                        {
                                        if ((sc_acc_atoms[j]->atom_type == O)
                                                && (sc_acc_atoms[j]->residue->terminal_status == CTERM)) {

                                             HbPair pair;
                                             pair.hb_don = bb_don_atoms[i];
                                             pair.hb_acc = sc_acc_atoms[j];
                                             pair.hb_acc_second = (*res_acc)[C];
                                             pair.hb_acc_third  = (*res_acc)[OXT];
                                             pair.other_acc_index = acc_count;
                                             pair.acceptor_type = AcceptorCarboxylate;
                                             pair.donor_type    = DonorAmide; 
                                             modified_pairs.push_back(pair);

                                        } else if (sc_acc_atoms[j]->atom_type == OXT) {

                                             HbPair pair;
                                             pair.hb_don = bb_don_atoms[i];
                                             pair.hb_acc = sc_acc_atoms[j];
                                             pair.hb_acc_second = (*res_acc)[C];
                                             pair.hb_acc_third  = (*res_acc)[O];
                                             pair.other_acc_index = acc_count;
                                             pair.acceptor_type = AcceptorCarboxylate;
                                             pair.donor_type    = DonorAmide; 
                                             modified_pairs.push_back(pair);

                                        } else {

                                             HbPair pair;
                                             pair.hb_don = bb_don_atoms[i];
                                             pair.hb_acc = sc_acc_atoms[j];
                                             pair.other_acc_index = acc_count;
                                             pair.donor_type  = DonorAmide;

                                             switch (sc_acc_atoms[j]->residue->residue_type) {

                                                  case ASN: {

                                                       pair.hb_acc_second = (*res_acc)[CG];
                                                       pair.hb_acc_third  = (*res_acc)[ND2];
                                                       pair.acceptor_type = AcceptorAmide;
                                                       modified_pairs.push_back(pair);
                                                       break;
                                                  }

                                                  case GLN: {

                                                       pair.hb_acc_second = (*res_acc)[CD];
                                                       pair.hb_acc_third  = (*res_acc)[NE2];
                                                       pair.acceptor_type = AcceptorAmide;
                                                       modified_pairs.push_back(pair);
                                                       break;
                                                  }

                                                  case GLU: {
                                                       pair.hb_acc_second = (*res_acc)[CD];
                                                       pair.hb_acc_third  = (*res_acc)[CG];
                                                       pair.acceptor_type = AcceptorCarboxylate;
                                                       modified_pairs.push_back(pair);
                                                       break;
                                                  }

                                                  case ASP: {

                                                       pair.hb_acc_second = (*res_acc)[CG];
                                                       pair.hb_acc_third  = (*res_acc)[CB];
                                                       pair.acceptor_type = AcceptorCarboxylate;
                                                       modified_pairs.push_back(pair);
                                                       break;
                                                  }

                                                  case THR: {

                                                       pair.hb_acc_second = (*res_acc)[CB];
                                                       pair.hb_acc_third  = (*res_acc)[HG1];
                                                       pair.acceptor_type = AcceptorAlcohol;
                                                       modified_pairs.push_back(pair);
                                                       break;
                                                  }

                                                  case SER: {

                                                       pair.hb_acc_second = (*res_acc)[CB];
                                                       pair.hb_acc_third  = (*res_acc)[HG];
                                                       pair.acceptor_type = AcceptorAlcohol;
                                                       modified_pairs.push_back(pair);
                                                       break;
                                                  }

                                                  case TYR: {
                                                       pair.hb_acc_second = (*res_acc)[CZ];
                                                       pair.hb_acc_third  = (*res_acc)[HH];
                                                       pair.acceptor_type = AcceptorAlcohol;
                                                       modified_pairs.push_back(pair);
                                                       break;
                                                  }

                                                  default:
                                                       break;

                                             } //end switch residue type
                                        } //end else
                                        break;
                                   } //end case
                                   case HA: case HA2://! Donor alpha hydrogen 
                                        {

                                        if ((sc_acc_atoms[j]->atom_type == O)
                                                && (sc_acc_atoms[j]->residue->terminal_status == CTERM)) {

                                                  HbPair pair;
                                                  pair.hb_don = bb_don_atoms[i];
                                                  pair.hb_acc = sc_acc_atoms[j];
                                                  pair.hb_acc_second = (*res_acc)[C];
                                                  pair.hb_acc_third  = (*res_acc)[OXT];
                                                  pair.other_acc_index = acc_count;
                                                  pair.acceptor_type = AcceptorCarboxylate;
                                                  pair.donor_type    = DonorAlphaHydrogen; 
                                                  modified_pairs.push_back(pair);

                                        } else if (sc_acc_atoms[j]->atom_type == OXT) {

                                                  HbPair pair;
                                                  pair.hb_don = bb_don_atoms[i];
                                                  pair.hb_acc = sc_acc_atoms[j];
                                                  pair.hb_acc_second = (*res_acc)[C];
                                                  pair.hb_acc_third  = (*res_acc)[O];
                                                  pair.other_acc_index = acc_count;
                                                  pair.acceptor_type = AcceptorCarboxylate;
                                                  pair.donor_type    = DonorAlphaHydrogen; 
                                                  modified_pairs.push_back(pair);

                                        } else {

                                                  HbPair pair;
                                                  pair.hb_don = bb_don_atoms[i];
                                                  pair.hb_acc = sc_acc_atoms[j];
                                                  pair.other_acc_index = acc_count;
                                                  pair.donor_type = DonorAlphaHydrogen; 

                                             switch (sc_acc_atoms[j]->residue->residue_type) {

                                                  case ASN: {

                                                       pair.hb_acc_second = (*res_acc)[CG];
                                                       pair.hb_acc_third  = (*res_acc)[ND2];
                                                       pair.acceptor_type = AcceptorAmide;
                                                       modified_pairs.push_back(pair);
                                                       break;
                                                  }

                                                  case GLN: {

                                                       pair.hb_acc_second = (*res_acc)[CD];
                                                       pair.hb_acc_third  = (*res_acc)[NE2];
                                                       pair.acceptor_type = AcceptorAmide;
                                                       modified_pairs.push_back(pair);
                                                       break;
                                                  }

                                                  case GLU: {
                                                       pair.hb_acc_second = (*res_acc)[CD];
                                                       pair.hb_acc_third  = (*res_acc)[CG];
                                                       pair.acceptor_type = AcceptorCarboxylate;
                                                       modified_pairs.push_back(pair);
                                                       break;
                                                  }

                                                  case ASP: {

                                                       pair.hb_acc_second = (*res_acc)[CG];
                                                       pair.hb_acc_third  = (*res_acc)[CB];
                                                       pair.acceptor_type = AcceptorCarboxylate;
                                                       modified_pairs.push_back(pair);
                                                       break;
                                                  }

                                                  case THR: {

                                                       pair.hb_acc_second = (*res_acc)[CB];
                                                       pair.hb_acc_third  = (*res_acc)[HG1];
                                                       pair.acceptor_type = AcceptorAlcohol;
                                                       modified_pairs.push_back(pair);
                                                       break;
                                                  }

                                                  case SER: {

                                                       pair.hb_acc_second = (*res_acc)[CB];
                                                       pair.hb_acc_third  = (*res_acc)[HG];
                                                       pair.acceptor_type = AcceptorAlcohol;
                                                       modified_pairs.push_back(pair);
                                                       break;
                                                  }

                                                  case TYR: {
                                                       pair.hb_acc_second = (*res_acc)[CZ];
                                                       pair.hb_acc_third  = (*res_acc)[HH];
                                                       pair.acceptor_type = AcceptorAlcohol;
                                                       modified_pairs.push_back(pair);
                                                       break;
                                                  }

                                                  default:
                                                       break;

                                             }
                                        }
                                        break;
                                   }

                                   default:
                                        std::cout << "(ProCS15Cached)critical fail in get_modified_bbsc_hbpair" << std::endl;
                                        break;
                              }
                         }
                    }
               }
          }


          return modified_pairs;
     }

     // Hbonds between HA3 donor and sidechain acceptor
     std::vector< HbPair > get_modified_ha3sc_hbpair(unsigned int start_index, unsigned int end_index){

          unsigned int acc_count;
          std::vector< HbPair > modified_pairs;
          for (unsigned int i = 0; i < ha3_don_atoms.size(); i++ ) {

               acc_count = 0;
               Residue *res_don = ha3_don_atoms[i]->residue;

               for(unsigned int j = 0; j< sc_acc_atoms.size(); j++ ) {
                   Residue *res_acc = sc_acc_atoms[j]->residue;

                   if (((res_acc->index >= start_index) &&
                        (res_acc->index <= end_index)) ||
                       ((res_don->index >= start_index) &&
                        (res_don->index <= end_index))) {
                         if ( j > 0 ){
                              if ((sc_acc_atoms[j]->residue->index) == sc_acc_atoms[j-1]->residue->index){

                                   acc_count += 1;
                              }
                              else{
                                   acc_count = 0;
                              }
                         }
                         if(abs(res_acc->index - res_don->index) > 1) {
                              if ((sc_acc_atoms[j]->atom_type == O)
                                      && (sc_acc_atoms[j]->residue->terminal_status == CTERM)) {

                                        HbPair pair;
                                        pair.hb_don = ha3_don_atoms[i];
                                        pair.hb_acc = sc_acc_atoms[j];
                                        pair.hb_acc_second = (*res_acc)[C];
                                        pair.hb_acc_third  = (*res_acc)[OXT];
                                        pair.other_acc_index = acc_count;
                                        pair.acceptor_type = AcceptorCarboxylate;
                                        pair.donor_type    = DonorAlphaHydrogen; 
                                        modified_pairs.push_back(pair);

                              } else if (sc_acc_atoms[j]->atom_type == OXT) {

                                        HbPair pair;
                                        pair.hb_don = ha3_don_atoms[i];
                                        pair.hb_acc = sc_acc_atoms[j];
                                        pair.hb_acc_second = (*res_acc)[C];
                                        pair.hb_acc_third  = (*res_acc)[O];
                                        pair.other_acc_index = acc_count;
                                        pair.acceptor_type = AcceptorCarboxylate;
                                        pair.donor_type    = DonorAlphaHydrogen; 
                                        modified_pairs.push_back(pair);

                              } else {

                                        HbPair pair;
                                        pair.hb_don = ha3_don_atoms[i];
                                        pair.hb_acc = sc_acc_atoms[j];
                                        pair.other_acc_index = acc_count;
                                        pair.donor_type = DonorAlphaHydrogen; 

                                   switch (sc_acc_atoms[j]->residue->residue_type) {

                                        case ASN: {

                                             pair.hb_acc_second = (*res_acc)[CG];
                                             pair.hb_acc_third  = (*res_acc)[ND2];
                                             pair.acceptor_type = AcceptorAmide;
                                             modified_pairs.push_back(pair);
                                             break;
                                        }

                                        case GLN: {

                                             pair.hb_acc_second = (*res_acc)[CD];
                                             pair.hb_acc_third  = (*res_acc)[NE2];
                                             pair.acceptor_type = AcceptorAmide;
                                             modified_pairs.push_back(pair);
                                             break;
                                        }

                                        case GLU: {
                                             pair.hb_acc_second = (*res_acc)[CD];
                                             pair.hb_acc_third  = (*res_acc)[CG];
                                             pair.acceptor_type = AcceptorCarboxylate;
                                             modified_pairs.push_back(pair);
                                             break;
                                        }

                                        case ASP: {

                                             pair.hb_acc_second = (*res_acc)[CG];
                                             pair.hb_acc_third  = (*res_acc)[CB];
                                             pair.acceptor_type = AcceptorCarboxylate;
                                             modified_pairs.push_back(pair);
                                             break;
                                        }

                                        case THR: {

                                             pair.hb_acc_second = (*res_acc)[CB];
                                             pair.hb_acc_third  = (*res_acc)[HG1];
                                             pair.acceptor_type = AcceptorAlcohol;
                                             modified_pairs.push_back(pair);
                                             break;
                                        }

                                        case SER: {

                                             pair.hb_acc_second = (*res_acc)[CB];
                                             pair.hb_acc_third  = (*res_acc)[HG];
                                             pair.acceptor_type = AcceptorAlcohol;
                                             modified_pairs.push_back(pair);
                                             break;
                                        }

                                        case TYR: {
                                             pair.hb_acc_second = (*res_acc)[CZ];
                                             pair.hb_acc_third  = (*res_acc)[HH];
                                             pair.acceptor_type = AcceptorAlcohol;
                                             modified_pairs.push_back(pair);
                                             break;
                                        }

                                        default:
                                             break;

                                   }
                              }
                         }
                    }
               }
          }


          return modified_pairs;
     }


     std::vector< HbPair > get_modified_scbb_hbpair(unsigned int start_index, unsigned int end_index) {

          std::vector< HbPair > modified_pairs;

          unsigned int don_index = 0;

          for (unsigned int i = 0; i< sc_don_atoms.size(); i++ ) {

               if (i > 0){ //! generate sidechain index of donor atom
                    if (sc_don_atoms[i]->residue->index == sc_don_atoms[i-1]->residue->index){

                         don_index += 1;

                    }
                    else {
                         don_index = 0;
                    }
               }

               for(unsigned int j = 0; j< bb_acc_atoms.size(); j++ ) {

                    Residue *res_don = sc_don_atoms[i]->residue;
                    Residue *res_acc = bb_acc_atoms[j]->residue;

                    if (((res_acc->index >= start_index) &&
                         (res_acc->index <= end_index)) ||
                        ((res_don->index >= start_index) &&
                         (res_don->index <= end_index))) {

                         if(abs(res_acc->index - res_don->index) > 1) {

                              HbPair pair;
                              pair.hb_don = sc_don_atoms[i];
                              pair.hb_acc = bb_acc_atoms[j];
                              pair.hb_acc_second = (*res_acc)[C];
                              pair.hb_acc_third  = (*chain)[res_acc->index+1][N];
                              pair.other_don_index = don_index;
                              modified_pairs.push_back(pair);
                         }
                    }
               }
          }

          return modified_pairs; 

     }



     void find_hydrogenbond_atoms(){

          //! Import protein definitions (such as residue names)
          using namespace definitions;

          //bbdon is backbone donors that contribute in primary hbond terms
          //N HA atoms. Including CTERM and NTERM for convenience
          unsigned int num_bbdon = this->chain->size();
          //ha3don is HA3 donors
          unsigned int num_ha3don = 0;
          //bbacc is backbone acceptors that contribute in secondary hbond terms
          //N-1 C=O from all but CTERM
          unsigned int num_bbacc = this->chain->size()-1;
          //scdon is donors that don't contribute in primary hbond terms
          //2 from NTERM H1+H2
          unsigned int num_scdon = 2;
          //scacc is acceptors that don't contribute in secondary hbond terms
          //2 from CTERM O+OXT
          unsigned int num_scacc = 2;

          // Count the number of acceptors and donors
          for (ResidueIterator<ChainFB> it(*this->chain); !it.end(); ++it){

              switch (it->residue_type) {
                  case ASN:
                  case GLN:
                      num_scacc+=1;
                      num_scdon+=2;
                      break;
                  case ASP:
                  case GLU:
                      num_scacc+=2;
                      break;
                  case LYS:
                      num_scdon+=3;
                      break;
                  case ARG:
                      num_scdon+=5;
                      break;
                  case SER: case THR: case TYR:
                      num_scacc+=1;
                      break;
                  case TRP:
                      num_scdon+=1;
                      break;
                  case HIS:
                      if ((*it).has_atom(HD1)) {
                           num_scdon+=1;
                      }
                      if ((*it).has_atom(HE2)) {
                           num_scdon+=1;
                      }
                      break;
                  default:
                      break;
              };

              if (it->terminal_status == NTERM && it->residue_type != PRO){
                   //H3
                   num_scdon += 1;
              }


              else if(it->terminal_status != NTERM && it->residue_type != PRO ){
                   //H
                   num_bbdon += 1;
              }
              if(it->residue_type == GLY){
                   // HA3.
                   num_ha3don += 1;
              }
          }

          bb_don_atoms.resize(num_bbdon);
          ha3_don_atoms.resize(num_ha3don);
          bb_acc_atoms.resize(num_bbacc);
          sc_don_atoms.resize(num_scdon);
          sc_acc_atoms.resize(num_scacc);

          unsigned int ibbdon=0;
          unsigned int iha3don=0;
          unsigned int ibbacc=0;
          unsigned int iscdon=0;
          unsigned int iscacc=0;
          // Put the acceptors and donors in the vectors
          for (ResidueIterator<ChainFB> it(*this->chain); !it.end(); ++it){


               if (it->terminal_status == NTERM) {
                    sc_don_atoms[iscdon++] = (*it)[H1];
                    sc_don_atoms[iscdon++] = (*it)[H2];
                    if (it->residue_type != PRO) {
                         sc_don_atoms[iscdon++] = (*it)[H3];
                    }
               } else if (it->residue_type != PRO) {
                    bb_don_atoms[ibbdon++] = (*it)[H];
               }


               if(it->residue_type == GLY){
                    bb_don_atoms[ibbdon++] = (*it)[HA2];
                    ha3_don_atoms[iha3don++] = (*it)[HA3];
               }
               else{
                    bb_don_atoms[ibbdon++] = (*it)[HA];
               }

               if(it->terminal_status != CTERM){

                    bb_acc_atoms[ibbacc++] = (*it)[O];
               }
               else {
                    sc_acc_atoms[iscacc++] = (*it)[O];
                    sc_acc_atoms[iscacc++] = (*it)[OXT];
               }


               if(it->residue_type == ARG){
                    sc_don_atoms[iscdon++] = (*it)[HE];
                    sc_don_atoms[iscdon++] = (*it)[HH11];
                    sc_don_atoms[iscdon++] = (*it)[HH12];
                    sc_don_atoms[iscdon++] = (*it)[HH22];
                    sc_don_atoms[iscdon++] = (*it)[HH21];
               }
               else if(it->residue_type == LYS){
                    sc_don_atoms[iscdon++] = (*it)[HZ1];
                    sc_don_atoms[iscdon++] = (*it)[HZ2];
                    sc_don_atoms[iscdon++] = (*it)[HZ3];
               }
               else if(it->residue_type == ASN){
                    sc_acc_atoms[iscacc++] = (*it)[OD1];
                    sc_don_atoms[iscdon++] = (*it)[HD21];
                    sc_don_atoms[iscdon++] = (*it)[HD22];
               }
               else if(it->residue_type == GLN){
                    sc_acc_atoms[iscacc++] = (*it)[OE1];
                    sc_don_atoms[iscdon++] = (*it)[HE21];
                    sc_don_atoms[iscdon++] = (*it)[HE22];
               }

               else if(it->residue_type == ASP){
                    sc_acc_atoms[iscacc++] = (*it)[OD1];
                    sc_acc_atoms[iscacc++] = (*it)[OD2];
               }
               else if(it->residue_type == GLU){
                    sc_acc_atoms[iscacc++] = (*it)[OE1];
                    sc_acc_atoms[iscacc++] = (*it)[OE2];
               }

               else if(it->residue_type == SER){
                    sc_acc_atoms[iscacc++] = (*it)[OG];
               }
               else if(it->residue_type == THR){
                    sc_acc_atoms[iscacc++] = (*it)[OG1];
               }

               else if(it->residue_type == TYR){
                    sc_acc_atoms[iscacc++] = (*it)[OH];
               }
               else if(it->residue_type == HIS){
                    if ((*it).has_atom(HD1)) {
                         sc_don_atoms[iscdon++] = (*it)[HD1];
                    }
                    if ((*it).has_atom(HE2)) {
                         sc_don_atoms[iscdon++] = (*it)[HE2];
                    }
               }
               else if(it->residue_type == TRP){
                    sc_don_atoms[iscdon++] = (*it)[HE1];
               }



          }

          if (ibbdon != num_bbdon) {
              std::cout << "error bbdon: " << ibbdon << " " << num_bbdon << std::endl;
              exit(1);
          }
          if (iha3don != num_ha3don) {
              std::cout << "error ha3don: " << iha3don << " " << num_ha3don << std::endl;
              exit(1);
          }
          if (ibbacc != num_bbacc) {
              std::cout << "error bbacc: " << ibbacc << " " << num_bbacc << std::endl;
              exit(1);
          }
          if (iscdon != num_scdon) {
              std::cout << "error scdon: " << iscdon << " " << num_scdon << std::endl;
              exit(1);
          }
          if (iscacc != num_scacc) {
              std::cout << "error scacc: " << iscacc << " " << num_scacc << std::endl;
              exit(1);
          }

          for (unsigned int i = 0; i < bb_don_atoms.size(); i++) {

               FPlist scacc;

               for (unsigned int j = 0; j < sc_acc_atoms.size(); j++) {
                    scacc.push_back(0.0);
               }
               bbdon_scacc.push_back(scacc);

               FPlist bbacc;
               for (unsigned int j = 0; j < bb_acc_atoms.size(); j++) {
                    bbacc.push_back(0.0);
               }
               bbdon_bbacc.push_back(bbacc);
          }

          for (unsigned int i = 0; i < sc_don_atoms.size(); i++) {

               FPlist bbacc;

               for (unsigned int j = 0; j < bb_acc_atoms.size(); j++) {
                    bbacc.push_back(0.0);
               }
               scdon_bbacc.push_back(bbacc);
          }

          for (unsigned int i = 0; i < ha3_don_atoms.size(); i++) {

               FPlist scacc;

               for (unsigned int j = 0; j < sc_acc_atoms.size(); j++) {
                    scacc.push_back(0.0);
               }
               ha3don_scacc.push_back(scacc);

               FPlist bbacc;
               for (unsigned int j = 0; j < bb_acc_atoms.size(); j++) {
                    bbacc.push_back(0.0);
               }
               ha3don_bbacc.push_back(bbacc);
          }

     }



public:

     class Ring_current_cached{ //! used in ring current cached calculations
          public:
                FPlist ring_shifts;
                FPlist ring_shifts_old;

          Ring_current_cached(){

               ring_shifts = empty_contribution;
               ring_shifts_old = empty_contribution;
          };

          void backup(){

                   //HA and HN
                   ring_shifts_old[0] = ring_shifts[0];
                   ring_shifts_old[2] = ring_shifts[2];
          }

          void rollback(){

                   ring_shifts[0] = ring_shifts_old[0];
                   ring_shifts[2] = ring_shifts_old[2];

          }

          void clear(){

                   ring_shifts[0] = 0.0;
                   ring_shifts[2] = 0.0;
          }
     };

     std::vector< std::vector< Ring_current_cached > > init_ringcurrent_chemical_shifts_table(phaistos::ChainFB& chain){

          unsigned int length = chain.size();

          std::vector< std::vector< Ring_current_cached > > out(length);

  
          for (unsigned int j = 0; j < length; j++) {

               std::vector< Ring_current_cached > row;

               for (ResidueIterator<ChainFB> res1(chain); !(res1).end(); ++res1) {

                     row.push_back( Ring_current_cached() );

               }

               out[j] = row;
          }

          return out;
     }

     class Bonded_shifts { //! used in bonded cached calculations

          public:

               FPlist shifts;
               FPlist shifts_old;

          Bonded_shifts(){

               shifts = empty_contribution;
               shifts_old = empty_contribution;
          };

          void backup(){

              for (unsigned int i = 0; i < 6; i++) {
                   shifts_old[i] = shifts[i];
              }
          }

          void rollback(){

              for (unsigned int i = 0; i < 6; i++) {
                   shifts[i] = shifts_old[i];
              }
          }

          void clear(){

              for (unsigned int i = 0; i < 6; i++) {
                   shifts[i] = 0.0;
              }
          }
     };

     std::vector< Bonded_shifts > init_bonded_chemical_shifts_table(phaistos::ChainFB& chain){

          std::vector< Bonded_shifts > out;

          for (ResidueIterator<ChainFB> res1(chain); !(res1).end(); ++res1) {

               out.push_back(Bonded_shifts());
          }

          return out;
     }



     // For convenience, define local EnergyTermCommon
     typedef phaistos::TermProCS15Base<TermProCS15Cached> TermProCS15Base;

     // Use settings from base class
     typedef TermProCS15Base::Settings Settings;


     class AmideWaterInteraction {

     public:

          bool water_interaction;
          bool water_interaction_old;

          AmideWaterInteraction(){

               water_interaction_old = true;
               water_interaction = true;
          }

          void backup(){
               water_interaction_old = water_interaction;
          }
     
          void rollback(){
               water_interaction = water_interaction_old;
          };
     };


     class DonorAcceptorInteraction {

     public:

         FPlist donor_interaction;
         FPlist donor_interaction_prev_res;
         FPlist donor_interaction_next_res;

         FPlist acceptor_interaction;
         FPlist acceptor_interaction_prev_res;
         FPlist acceptor_interaction_next_res;

         FPlist donor_interaction_old;
         FPlist donor_interaction_prev_res_old;
         FPlist donor_interaction_next_res_old;

         FPlist acceptor_interaction_old;
         FPlist acceptor_interaction_prev_res_old;
         FPlist acceptor_interaction_next_res_old;

         bool bound; //!bound or not bound
         bool bound_old;

         void backup_acceptor() {
              for (unsigned int i = 0; i < 6; i++) {
                   acceptor_interaction_old[i]          = acceptor_interaction[i];
                   acceptor_interaction_prev_res_old[i] = acceptor_interaction_prev_res[i];
                   acceptor_interaction_next_res_old[i] = acceptor_interaction_next_res[i];
              }
         }

         void backup_donor() {
              for (unsigned int i = 0; i < 6; i++) {
                   donor_interaction_old[i]          = donor_interaction[i];
                   donor_interaction_prev_res_old[i] = donor_interaction_prev_res[i];
                   donor_interaction_next_res_old[i] = donor_interaction_next_res[i];
                                          bound_old  = bound; 
              }
         }

         void rollback_acceptor() {
              for (unsigned int i = 0; i < 6; i++) {
                   acceptor_interaction[i]          = acceptor_interaction_old[i];
                   acceptor_interaction_prev_res[i] = acceptor_interaction_prev_res_old[i];
                   acceptor_interaction_next_res[i] = acceptor_interaction_next_res_old[i];
              }
         }

         void rollback_donor() {
              for (unsigned int i = 0; i < 6; i++) {
                   donor_interaction[i]          = donor_interaction_old[i];
                   donor_interaction_prev_res[i] = donor_interaction_prev_res_old[i];
                   donor_interaction_next_res[i] = donor_interaction_next_res_old[i];
                                          bound  = bound_old;
              }
         }

         void clear_acceptor() {

              for (unsigned int i = 0; i < 6; i++) {
                   acceptor_interaction[i]          = 0.0;
                   acceptor_interaction_prev_res[i] = 0.0;
                   acceptor_interaction_next_res[i] = 0.0;
              }
         }

         void clear_donor() {
              for (unsigned int i = 0; i < 6; i++) {
                   donor_interaction[i]          = 0.0;
                   donor_interaction_prev_res[i] = 0.0;
                   donor_interaction_next_res[i] = 0.0;
                                           bound = false;
              }
         }

         DonorAcceptorInteraction() {

                   bound = false;
                   bound_old = false;

                   donor_interaction                  = empty_contribution;
                   donor_interaction_prev_res         = empty_contribution;
                   donor_interaction_next_res         = empty_contribution;

                   acceptor_interaction               = empty_contribution;
                   acceptor_interaction_prev_res      = empty_contribution;
                   acceptor_interaction_next_res      = empty_contribution;

                   donor_interaction_old              = empty_contribution;
                   donor_interaction_prev_res_old     = empty_contribution;
                   donor_interaction_next_res_old     = empty_contribution;

                   acceptor_interaction_old           = empty_contribution;
                   acceptor_interaction_prev_res_old  = empty_contribution;
                   acceptor_interaction_next_res_old  = empty_contribution;

          };
     };



     class AmideAmideInteraction: public DonorAcceptorInteraction{
     public:

          AmideAmideInteraction(): DonorAcceptorInteraction() {}
     };


     class AlphaOxygenInteraction: public DonorAcceptorInteraction{
     public:

          AlphaOxygenInteraction(): DonorAcceptorInteraction() {}
     };


     class RingCurrentInteraction{ //! wrapper class for Ring_current_cached 
          public:
               std::vector< Ring_current_cached > sidechain_rings;
               unsigned int side_chain_rings_nr;

               RingCurrentInteraction(ResidueEnum t){

                    switch (t){
                         case TRP:
                              sidechain_rings.push_back( Ring_current_cached() );
                              sidechain_rings.push_back( Ring_current_cached() );
                              side_chain_rings_nr = 2;
                              break;

                         case TYR:
                              sidechain_rings.push_back( Ring_current_cached() );
                              side_chain_rings_nr = 1;
                              break;

                         case HIS:
                              sidechain_rings.push_back( Ring_current_cached() );
                              side_chain_rings_nr = 1;
                              break;

                         case PHE:
                              sidechain_rings.push_back( Ring_current_cached() );
                              side_chain_rings_nr = 1;
                              break;

                         default:
                              sidechain_rings.push_back( Ring_current_cached() );
                              side_chain_rings_nr = 0;
                              break;
                     }
               };

          void backup(unsigned int index){

                sidechain_rings[index].backup();
          }

          void rollback(unsigned int index){

                sidechain_rings[index].rollback();
          }

          void clear(unsigned int index){             
                sidechain_rings[index].clear();
          }
     };

     std::vector < std::vector< RingCurrentInteraction > > init_ringcurrentinteraction_new(phaistos::ChainFB& chain){

          unsigned int length = chain.size();

          std::vector < std::vector< RingCurrentInteraction > > out(length);

          for (unsigned int i = 0; i < length; i++){
               std::vector< RingCurrentInteraction > row;

               for (ResidueIterator<ChainFB> res2(chain); !(res2).end(); ++res2) {
                    row.push_back( RingCurrentInteraction(res2->residue_type) );
               }
               out[i] = row;
          }

          return out;
     }

     class AlphaSidechainInteraction{ //! Wrapper class for AlphaOxygenInteraction

          public:
               std::vector< AlphaOxygenInteraction > interactions;
               unsigned int side_chain_acceptors_nr;

               //res1 is acceptor, res2 is donor
               AlphaSidechainInteraction(ResidueEnum t, bool cterm){

               switch (t) {

                    case ASN: case GLN: case SER: case THR: case TYR:
                         interactions.push_back(AlphaOxygenInteraction() );
                         side_chain_acceptors_nr = 1;
                         if (cterm == true){
                              interactions.push_back( AlphaOxygenInteraction() );
                              interactions.push_back( AlphaOxygenInteraction() );
                              side_chain_acceptors_nr = 3;
                         }
                         break;


                    case ASP: case GLU:
                         interactions.push_back( AlphaOxygenInteraction() );
                         interactions.push_back( AlphaOxygenInteraction() );
                         side_chain_acceptors_nr = 2;
                         if (cterm == true){
                              interactions.push_back( AlphaOxygenInteraction() );
                              interactions.push_back( AlphaOxygenInteraction() );
                              side_chain_acceptors_nr = 4;
                         }
                         break;

                    default:
                         interactions.push_back( AlphaOxygenInteraction() ); //Dummy
                         side_chain_acceptors_nr = 0;
                         if (cterm == true){
                              interactions.push_back( AlphaOxygenInteraction() );
                              side_chain_acceptors_nr = 2;
                         }
                         break;
               }
          };

          void backup_acceptor(int index){
               interactions[index].backup_acceptor();
          }
          void backup_donor(int index){
               interactions[index].backup_donor();
          }
          void clear_acceptor(int index){
               interactions[index].clear_acceptor();
          }
          void clear_donor(int index){
               interactions[index].clear_donor();
          }
          void rollback_acceptor(int index){
               interactions[index].rollback_acceptor();
          }
          void rollback_donor(int index){
               interactions[index].rollback_donor();
          }

     };


     class SidechainAmideInteraction{ //! Wrapper class for AmideAmideInteraction

          public:
               std::vector< AmideAmideInteraction > interactions;
               unsigned int side_chain_donors_nr;
 
               // All donors, excluding HA, not contributing to primary HN bond corrections.
               SidechainAmideInteraction(ResidueEnum t, bool nterm, phaistos::ResidueFB& res){

               switch (t){
                    case LYS:
                         interactions.push_back( AmideAmideInteraction() );
                         interactions.push_back( AmideAmideInteraction() );
                         interactions.push_back( AmideAmideInteraction() );
                         side_chain_donors_nr = 3;
                         if (nterm == true){
                              interactions.push_back( AmideAmideInteraction() );
                              interactions.push_back( AmideAmideInteraction() );
                              interactions.push_back( AmideAmideInteraction() );
                              side_chain_donors_nr = 6;
                         }

                         break;

                    case ARG:
                         interactions.push_back( AmideAmideInteraction() );
                         interactions.push_back( AmideAmideInteraction() );
                         interactions.push_back( AmideAmideInteraction() );
                         interactions.push_back( AmideAmideInteraction() );
                         interactions.push_back( AmideAmideInteraction() );
                         side_chain_donors_nr = 5;
                         if (nterm == true){
                              interactions.push_back( AmideAmideInteraction() );
                              interactions.push_back( AmideAmideInteraction() );
                              interactions.push_back( AmideAmideInteraction() );
                              side_chain_donors_nr = 8;
                         }
                         break;

                    case HIS:
                         side_chain_donors_nr = 0;
                         if (res.has_atom(HD1)) {
                              interactions.push_back( AmideAmideInteraction() );
                              side_chain_donors_nr++;
                         }
                         if (res.has_atom(HE2)) {
                              interactions.push_back( AmideAmideInteraction() );
                              side_chain_donors_nr++;
                         }
                         if (nterm == true){
                              interactions.push_back( AmideAmideInteraction() );
                              interactions.push_back( AmideAmideInteraction() );
                              interactions.push_back( AmideAmideInteraction() );
                              side_chain_donors_nr += 3;
                         }
                         break;

                    case TRP:
                         interactions.push_back( AmideAmideInteraction() );
                         side_chain_donors_nr = 1;
                         if (nterm == true){
                              interactions.push_back( AmideAmideInteraction() );
                              interactions.push_back( AmideAmideInteraction() );
                              interactions.push_back( AmideAmideInteraction() );
                              side_chain_donors_nr = 4;
                         }
                         break;

                    case ASN:
                         interactions.push_back( AmideAmideInteraction() );
                         interactions.push_back( AmideAmideInteraction() );
                         side_chain_donors_nr = 2;
                         if (nterm == true){
                              interactions.push_back( AmideAmideInteraction() );
                              interactions.push_back( AmideAmideInteraction() );
                              interactions.push_back( AmideAmideInteraction() );
                              side_chain_donors_nr = 5;
                         }
                         break;

                    case GLN:
                         interactions.push_back( AmideAmideInteraction() );
                         interactions.push_back( AmideAmideInteraction() );
                         side_chain_donors_nr = 2;
                         if (nterm == true){
                              interactions.push_back( AmideAmideInteraction() );
                              interactions.push_back( AmideAmideInteraction() );
                              interactions.push_back( AmideAmideInteraction() ); 
                              side_chain_donors_nr = 5;
                         }
                         break;
                    case PRO:
                         interactions.push_back( AmideAmideInteraction() ); //not used just for housekeeping(unless == Nterm)
                         side_chain_donors_nr = 0;
                         if (nterm == true){
                              interactions.push_back( AmideAmideInteraction() );
                              side_chain_donors_nr = 2;
                         }
                         break;

                    default:
                         interactions.push_back( AmideAmideInteraction() ); //not used just for housekeeping(unless == Nterm)
                         side_chain_donors_nr = 0;
                         if (nterm == true){
                              interactions.push_back( AmideAmideInteraction() );
                              interactions.push_back( AmideAmideInteraction() ); 
                              side_chain_donors_nr = 3;
                         }
                         break;
               }
          };

          void backup_acceptor(int index){
               interactions[index].backup_acceptor();
          }
          void backup_donor(int index){
               interactions[index].backup_donor();
          }
          void clear_acceptor(int index){
               interactions[index].clear_acceptor();
          }
          void clear_donor(int index){
               interactions[index].clear_donor();
          }
          void rollback_acceptor(int index){
               interactions[index].rollback_acceptor();
          }
          void rollback_donor(int index){
               interactions[index].rollback_donor();
          }

     };

     class AmideSidechainInteraction{ //! Wrapper class for AmideAmideInteraction

          public:
               std::vector< AmideAmideInteraction > interactions;
               int side_chain_acceptors_nr;
 
               AmideSidechainInteraction(ResidueEnum t, bool cterm){

               // All acceptors not contributing to secondary HN hbond corrections.
               switch (t){

                    case ASN:
                         interactions.push_back( AmideAmideInteraction() );
                         side_chain_acceptors_nr = 1;
                         if (cterm == true){
                              interactions.push_back( AmideAmideInteraction() );
                              interactions.push_back( AmideAmideInteraction() );
                              side_chain_acceptors_nr = 3;
                         }
                         break;

                    case GLN:
                         interactions.push_back( AmideAmideInteraction() );
                         side_chain_acceptors_nr = 1;
                         if (cterm == true){
                              interactions.push_back( AmideAmideInteraction() );
                              interactions.push_back( AmideAmideInteraction() );
                              side_chain_acceptors_nr = 3;
                         }
                         break;

                    case ASP:
                         interactions.push_back( AmideAmideInteraction() );
                         interactions.push_back( AmideAmideInteraction() );
                         side_chain_acceptors_nr = 2;
                         if (cterm == true){
                              interactions.push_back( AmideAmideInteraction() );
                              interactions.push_back( AmideAmideInteraction() );
                              side_chain_acceptors_nr = 4;
                         }
                         break;

                    case GLU:
                         interactions.push_back( AmideAmideInteraction() );
                         interactions.push_back( AmideAmideInteraction() );
                         side_chain_acceptors_nr = 2;
                         if (cterm == true){
                              interactions.push_back( AmideAmideInteraction() );
                              interactions.push_back( AmideAmideInteraction() );
                              side_chain_acceptors_nr = 4;
                         }
                         break;

                    case SER:
                         interactions.push_back( AmideAmideInteraction() );
                         side_chain_acceptors_nr = 1;
                         if (cterm == true){
                              interactions.push_back( AmideAmideInteraction() );
                              interactions.push_back( AmideAmideInteraction() );
                              side_chain_acceptors_nr = 3;
                         }
                         break;

                    case THR:
                         interactions.push_back( AmideAmideInteraction() );
                         side_chain_acceptors_nr = 1;
                         if (cterm == true){
                              interactions.push_back( AmideAmideInteraction() );
                              interactions.push_back( AmideAmideInteraction() );
                              side_chain_acceptors_nr = 3;
                         }
                         break;

                    case TYR:
                         interactions.push_back( AmideAmideInteraction() ); 
                         side_chain_acceptors_nr = 1;
                         if (cterm == true){
                              interactions.push_back( AmideAmideInteraction() );
                              interactions.push_back( AmideAmideInteraction() );
                              side_chain_acceptors_nr = 3;
                         }
                         break;

                    default:
                         interactions.push_back( AmideAmideInteraction() );// Not used, only kept for bookkeeping.
                         side_chain_acceptors_nr = 0;
                         if (cterm == true){
                              interactions.push_back( AmideAmideInteraction() );
                              side_chain_acceptors_nr = 2;
                         }
                         break;
                }
          };

          void backup_acceptor(int index){
               interactions[index].backup_acceptor();
          }
          void backup_donor(int index){
               interactions[index].backup_donor();
          }
          void clear_acceptor(int index){
               interactions[index].clear_acceptor();
          }
          void clear_donor(int index){
               interactions[index].clear_donor();
          }
          void rollback_acceptor(int index){
               interactions[index].rollback_acceptor();
          }
          void rollback_donor(int index){
               interactions[index].rollback_donor();
          }

     };

     //! make scacc amide bbdon matrix
     std::vector< std::vector< AmideSidechainInteraction > > init_amide_sidechain_interaction_matrix(phaistos::ChainFB& chain){

          unsigned int length = chain.size();
          std::vector< std::vector< AmideSidechainInteraction > > initialized_amide_sidechain_interaction_matrix(length);

          for (unsigned int j = 0; j < length; j++) {
               std::vector< AmideSidechainInteraction > initialized_amide_sidechain_interation_matrix_row;
               for (ResidueIterator<ChainFB> res1(chain); !(res1).end(); ++res1) {
                    bool cterm = false;
                    if (res1->terminal_status == CTERM){
                         cterm = true;
                    }
                     initialized_amide_sidechain_interation_matrix_row.push_back( AmideSidechainInteraction(res1->residue_type, cterm) );

               }
               initialized_amide_sidechain_interaction_matrix[j] = initialized_amide_sidechain_interation_matrix_row;
          }

          return initialized_amide_sidechain_interaction_matrix;
    }

     //! make scdon bbacc matrix
     std::vector< std::vector< SidechainAmideInteraction > > init_sidechain_amide_interaction_matrix(phaistos::ChainFB& chain){

          unsigned int length = chain.size();
          std::vector< std::vector< SidechainAmideInteraction > > initialized_sidechain_amide_interation_matrix(length);

          for (unsigned int j = 0; j < length; j++) {

               std::vector< SidechainAmideInteraction > initialized_sidechain_amide_interation_matrix_row;

               for (ResidueIterator<ChainFB> res1(chain); !(res1).end(); ++res1) {

                      bool nterm = false;
                      if (res1->terminal_status == NTERM){
                         nterm = true;
                      }
                      initialized_sidechain_amide_interation_matrix_row.push_back( SidechainAmideInteraction(res1->residue_type, nterm, *res1) );
               }

               initialized_sidechain_amide_interation_matrix[j] = initialized_sidechain_amide_interation_matrix_row;
          }

          return initialized_sidechain_amide_interation_matrix;
    }


    //! make bbdonbbacc alpha_oxygen matrix interaction between alpha protons and backbone oxygen
    std::vector< std::vector< AlphaOxygenInteraction > > init_alpha_oxygen_interaction_matrix(phaistos::ChainFB& chain) {

          unsigned int length = chain.size();
          std::vector< std::vector< AlphaOxygenInteraction > > initialized_alpha_oxygen_interation_matrix(length);

          for (unsigned int i = 0; i < length; i++) {
               std::vector< AlphaOxygenInteraction > initialized_alpha_oxygen_interation_matrix_row(length);

               for (unsigned int j = 0; j < length; j++) {

                    initialized_alpha_oxygen_interation_matrix_row[j] = AlphaOxygenInteraction();

               }
               initialized_alpha_oxygen_interation_matrix[i] = initialized_alpha_oxygen_interation_matrix_row;
          }

          return initialized_alpha_oxygen_interation_matrix;
    }

    //! make for interaction between backbone alpha hydrogen and side chain acceptors
    std::vector< std::vector< AlphaSidechainInteraction > > init_alpha_sidechain_interaction_matrix(phaistos::ChainFB& chain){

          unsigned int length = chain.size();
          std::vector< std::vector< AlphaSidechainInteraction > > initialized_alpha_sidechain_interaction_matrix(length);

          for (unsigned int j = 0; j < length; j++) {
               std::vector< AlphaSidechainInteraction > initialized_alpha_sidechain_interation_matrix_row;
               for (ResidueIterator<ChainFB> res1(chain); !(res1).end(); ++res1) {

                    bool cterm = false;
                    if (res1->terminal_status == CTERM){
                         cterm = true;
                    }

                    initialized_alpha_sidechain_interation_matrix_row.push_back( AlphaSidechainInteraction(res1->residue_type, cterm) );

               }
               initialized_alpha_sidechain_interaction_matrix[j] = initialized_alpha_sidechain_interation_matrix_row;
          }
          return initialized_alpha_sidechain_interaction_matrix;

    }


    //! make bbdonbbacc amide_amide matrix
    std::vector< std::vector< AmideAmideInteraction > > init_amide_amide_interaction_matrix(phaistos::ChainFB& chain) {

          unsigned int length = chain.size();
          std::vector< std::vector< AmideAmideInteraction > > initialized_amide_amide_interation_matrix(length);

          for (unsigned int i = 0; i < length; i++) {
               std::vector< AmideAmideInteraction > initialized_amide_amide_interation_matrix_row(length);

               for (unsigned int j = 0; j < length; j++) {

                    initialized_amide_amide_interation_matrix_row[j] = AmideAmideInteraction();

               }
               initialized_amide_amide_interation_matrix[i] = initialized_amide_amide_interation_matrix_row;
          }

          return initialized_amide_amide_interation_matrix;
     }


     struct ModifiedPair {

          int acc_index;
          int don_index;
          unsigned int other_don_index; //! index of side chain donor hydrogen atoms
          unsigned int other_acc_index; //! index of side chain acceptor atoms
          DonorEnum donor_type;
          AcceptorEnum acceptor_type;
     };
     
     struct ModifiedRingHydrogenPair {
         int hydrogen_residue_index;
         int ring_residue_index;
         unsigned int ring_nr;

     };


     std::vector< ModifiedPair > bbdon_bbacc_modified_pairs;
     std::vector< ModifiedPair > ha3don_bbacc_modified_pairs;
     std::vector< ModifiedPair > scdon_bbacc_modified_pairs;
     std::vector< ModifiedPair > bbdon_scacc_modified_pairs;
     std::vector< ModifiedPair > ha3don_scacc_modified_pairs;
     std::vector< RingHydrogenPair > rh_pairs;


     std::vector< ModifiedRingHydrogenPair > ring_hydrogen_modified_pairs;


     std::vector< std::vector< AmideAmideInteraction > >     amide_amide_interaction_matrix_bbbb;
     std::vector< std::vector< SidechainAmideInteraction > > sidehchain_amide_interaction_matrix;
     std::vector< std::vector< AmideSidechainInteraction > > amide_sidechain_interaction_matrix;

     std::vector< std::vector< AlphaOxygenInteraction > > alpha_oxygen_interaction_matrix_bbbb;
     std::vector< std::vector< AlphaOxygenInteraction > > ha3_oxygen_interaction_matrix_bbbb;
     std::vector< std::vector< AlphaSidechainInteraction > > alpha_sidechain_interaction_matrix;
     std::vector< std::vector< AlphaSidechainInteraction > > ha3_sidechain_interaction_matrix;

     std::vector< Bonded_shifts > bonded_interactions;

     std::vector< std::vector< Ring_current_cached > > ring_current_interactions;

     std::vector< std::vector< RingCurrentInteraction > > ring_current_interactions_new;

     FPtable bonded_chemical_shifts_table;
     FPtable ring_current_table;

     std::vector< AmideWaterInteraction > amide_water_interaction;

     unsigned int bonded_start;
     unsigned int bonded_end;

     unsigned int ring_start;
     unsigned int ring_end;

     //Constructor
     TermProCS15Cached(ChainFB *chain,
                  const Settings &settings=Settings(),
                  RandomNumberEngine *random_number_engine = &random_global)
          : TermProCS15Base(chain, "procs15-cached", settings, random_number_engine),
          ProCS15Backend(*chain, settings.load_ca, settings.load_cb, settings.load_co, settings.load_n, settings.load_hn, settings.load_ha, settings.procsnumpypath) {

          using namespace definitions;

          bbdon_bbacc_modified_pairs.clear();
          ha3don_bbacc_modified_pairs.clear();
          scdon_bbacc_modified_pairs.clear();
          bbdon_scacc_modified_pairs.clear();
          ha3don_scacc_modified_pairs.clear();

          find_hydrogenbond_atoms();

          find_hydrogen_ring_pairs();

          ring_current_interactions_new = init_ringcurrentinteraction_new(*chain);
          amide_amide_interaction_matrix_bbbb =  init_amide_amide_interaction_matrix(*chain);
          amide_sidechain_interaction_matrix  =  init_amide_sidechain_interaction_matrix(*chain);
          sidehchain_amide_interaction_matrix =  init_sidechain_amide_interaction_matrix(*chain);

          alpha_oxygen_interaction_matrix_bbbb = init_alpha_oxygen_interaction_matrix(*chain); 
          ha3_oxygen_interaction_matrix_bbbb = init_alpha_oxygen_interaction_matrix(*chain); 
          alpha_sidechain_interaction_matrix = init_alpha_sidechain_interaction_matrix(*chain);
          ha3_sidechain_interaction_matrix = init_alpha_sidechain_interaction_matrix(*chain);

          loaded_data.push_back(settings.load_ca);
          loaded_data.push_back(settings.load_cb);
          loaded_data.push_back(settings.load_co);
          loaded_data.push_back(settings.load_n);
          loaded_data.push_back(settings.load_hn);
          loaded_data.push_back(settings.load_ha);

          //! make non_bonded_chemical_shifts_table 
          for (ResidueIterator<ChainFB> res1(*chain); !(res1).end(); ++res1) {
               non_bonded_chemical_shifts_table.push_back(empty_contribution);
          }

          //! make bonded_chemical_shifts_table 
          for (ResidueIterator<ChainFB> res1(*chain); !(res1).end(); ++res1) {
               bonded_chemical_shifts_table.push_back(empty_contribution);
          }

          //! make ring current chemical shifts table 
          for (ResidueIterator<ChainFB> res1(*chain); !(res1).end(); ++res1) {
               ring_current_table.push_back(empty_contribution);
          }
          
         //! make table that count amide number of bonds
          for (ResidueIterator<ChainFB> res1(*chain); !(res1).end(); ++res1) {
               amide_bond_count_table.push_back( 0 );
          }

          for (ResidueIterator<ChainFB> res1(*chain); !(res1).end(); ++res1) {
               amide_water_interaction.push_back( AmideWaterInteraction() );
          }

          bonded_interactions = init_bonded_chemical_shifts_table(*chain);
          ring_current_interactions = init_ringcurrent_chemical_shifts_table(*chain);
          this->predicted_chemical_shifts = this->predict(*chain);

     }



     //Copy constructor
     TermProCS15Cached(const TermProCS15Cached &other,
                  RandomNumberEngine *random_number_engine,
                  int thread_index, ChainFB *chain)
          : TermProCS15Base(other, random_number_engine, thread_index, chain),
          ProCS15Backend(other, *chain),
          bbdon_bbacc_modified_pairs(other.bbdon_bbacc_modified_pairs),
          ha3don_bbacc_modified_pairs(other.ha3don_bbacc_modified_pairs),
          scdon_bbacc_modified_pairs(other.scdon_bbacc_modified_pairs),
          bbdon_scacc_modified_pairs(other.bbdon_scacc_modified_pairs),
          ha3don_scacc_modified_pairs(other.ha3don_scacc_modified_pairs){

          find_hydrogenbond_atoms();
          find_hydrogen_ring_pairs();

          amide_amide_interaction_matrix_bbbb = init_amide_amide_interaction_matrix(*chain);
          amide_sidechain_interaction_matrix =  init_amide_sidechain_interaction_matrix(*chain);
          sidehchain_amide_interaction_matrix =  init_sidechain_amide_interaction_matrix(*chain);

          alpha_oxygen_interaction_matrix_bbbb = init_alpha_oxygen_interaction_matrix(*chain); 
          ha3_oxygen_interaction_matrix_bbbb = init_alpha_oxygen_interaction_matrix(*chain); 
          alpha_sidechain_interaction_matrix = init_alpha_sidechain_interaction_matrix(*chain);
          ha3_sidechain_interaction_matrix = init_alpha_sidechain_interaction_matrix(*chain);

          ring_current_interactions = init_ringcurrent_chemical_shifts_table(*chain); 
          ring_current_interactions_new = init_ringcurrentinteraction_new(*chain);

          amide_water_interaction = other.amide_water_interaction;

          loaded_data = other.loaded_data;

          for (ResidueIterator<ChainFB> res1(*chain); !(res1).end(); ++res1) {
               non_bonded_chemical_shifts_table.push_back(empty_contribution);
          }

          //! make bonded_chemical_shifts_table and interactions 
          bonded_interactions = init_bonded_chemical_shifts_table(*chain);

          for (ResidueIterator<ChainFB> res1(*chain); !(res1).end(); ++res1) {
               bonded_chemical_shifts_table.push_back(empty_contribution);
          }


          //! make ring current chemical shifts table 
          for (ResidueIterator<ChainFB> res1(*chain); !(res1).end(); ++res1) {
               ring_current_table.push_back(empty_contribution);
          }

          for (ResidueIterator<ChainFB> res1(*chain); !(res1).end(); ++res1) {
               amide_bond_count_table.push_back( 0 );
          }

     }


     //! Check whether any of the moved residues involve potential hydrogen bonds
     inline bool hydrogen_bonding(MoveInfo *move_info){

          // Import protein definitions (such as residue names)
          using namespace definitions;

          for (int i=move_info->modified_angles_start; i< move_info->modified_angles_end; i++) {
               ResidueEnum re = (*this->chain)[i].residue_type;
               if (re == ARG ||
                   re == LYS||
                   re == ASN||
                   re == GLN||
                   re == HIS||
                   re == TRP||
                   re == ASP||
                   re == GLU|| 
                   re == SER||
                   re == THR||
                   re == TYR ) {
                         return true;
               }
          }
          return false;
     }


     FPtable predict(phaistos::ChainFB& chain, MoveInfo *move_info=NULL) {


          int start_index = 0;
          int end_index = chain.size() - 1;

          if (move_info) {
              start_index = move_info->modified_positions_start;
              end_index   = move_info->modified_positions_end - 1;
          }

          //! Backbone donor - backbone acceptor interactions
          if ((!move_info) || (move_info->move_type != definitions::SIDECHAIN)){

               bbdon_bbacc_modified_pairs.clear();
               std::vector< HbPair > pairs = get_modified_bbbb_hbpair(start_index, end_index);
               for (std::vector< HbPair>::iterator pair = pairs.begin(); pair < pairs.end(); pair++) {

                    const unsigned int res_don_index = pair->hb_don->residue->index;
                    const unsigned int res_acc_index = pair->hb_acc->residue->index;



                    switch (pair->donor_type){
                         case DonorAmide: //! Amide proton donor
                         {
                              amide_amide_interaction_matrix_bbbb[res_acc_index][res_don_index].backup_acceptor();
                              amide_amide_interaction_matrix_bbbb[res_don_index][res_acc_index].backup_donor();

                              const FPtype r_oh = (pair->hb_acc->position - pair->hb_don->position).norm();

                              if (r_oh < 3.0) {

                                   const FPtype theta = calc_angle(pair->hb_don->position,
                                                                   pair->hb_acc->position,
                                                                   pair->hb_acc_second->position);

                                   const FPtype rho = calc_dihedral(pair->hb_don->position,
                                                                    pair->hb_acc->position,
                                                                    pair->hb_acc_second->position,
                                                                    pair->hb_acc_third->position);

                                   const boost::multi_array<float, 1> data_vector = get_bb_bb_data_vector(r_oh, theta, rho);
 

                                   //primary H hbond terms.
                                   amide_amide_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction[0]          = (FPtype)data_vector[5]; // HA, n
                                   //amide_amide_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction[1]          = (FPtype)data_vector[0]; // CA, n
                                   amide_amide_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction[2]          = (FPtype)data_vector[4]; // H,  n
                                   amide_amide_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction[3]          = (FPtype)data_vector[3]; // N,  n
                                   //amide_amide_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction_prev_res[4] = (FPtype)data_vector[2]; // C,  n-1

                                   //secondary H hbond terms.
                                   amide_amide_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction[0]             = (FPtype)data_vector[11]; // HA, n
                                   amide_amide_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction[1]             = (FPtype)data_vector[ 6]; // CA, n
                                   amide_amide_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_next_res[2]    = (FPtype)data_vector[ 10]; // H,  n+1
                                   amide_amide_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_next_res[3]    = (FPtype)data_vector[ 9]; // N,  n+1
                                   amide_amide_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction[4]             = (FPtype)data_vector[ 8]; // C,  n

                                   amide_amide_interaction_matrix_bbbb[res_don_index][res_acc_index].bound = true;

                              } else {

                                   amide_amide_interaction_matrix_bbbb[res_acc_index][res_don_index].clear_acceptor();
                                   amide_amide_interaction_matrix_bbbb[res_don_index][res_acc_index].clear_donor();
                                   amide_amide_interaction_matrix_bbbb[res_don_index][res_acc_index].bound = false;

                              }
                              //! add to non bonded table

                              unsigned int chain_length = non_bonded_chemical_shifts_table.size() - 1;

                              for (unsigned int i = 0; i < 5; i++) {

                                   non_bonded_chemical_shifts_table[res_don_index][i] -= amide_amide_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction_old[i];
                                   non_bonded_chemical_shifts_table[res_don_index][i] += amide_amide_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction[i];
                                   if (res_don_index != 0){
                                        non_bonded_chemical_shifts_table[res_don_index-1][i] -= amide_amide_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction_prev_res_old[i];
                                        non_bonded_chemical_shifts_table[res_don_index-1][i] += amide_amide_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction_prev_res[i];
                                   }
                                   if (res_don_index != chain_length){
                                        non_bonded_chemical_shifts_table[res_don_index+1][i] -= amide_amide_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction_next_res_old[i];
                                        non_bonded_chemical_shifts_table[res_don_index+1][i] += amide_amide_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction_next_res[i];
                                   }
                                   non_bonded_chemical_shifts_table[res_acc_index][i] -= amide_amide_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_old[i];
                                   non_bonded_chemical_shifts_table[res_acc_index][i] += amide_amide_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction[i];
                                   if (res_acc_index != 0){
                                        non_bonded_chemical_shifts_table[res_acc_index-1][i] -= amide_amide_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_prev_res_old[i];
                                        non_bonded_chemical_shifts_table[res_acc_index-1][i] += amide_amide_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_prev_res[i];
                                   }
                                   if (res_acc_index != chain_length){
                                        non_bonded_chemical_shifts_table[res_acc_index+1][i] -= amide_amide_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_next_res_old[i];
                                        non_bonded_chemical_shifts_table[res_acc_index+1][i] += amide_amide_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_next_res[i];
                                   }

                              }

                              //! add state of hydrogen bond pair to amide_bond_count_table
                              amide_bond_count_table[res_don_index] -= amide_amide_interaction_matrix_bbbb[res_don_index][res_acc_index].bound_old;
                              amide_bond_count_table[res_don_index] += amide_amide_interaction_matrix_bbbb[res_don_index][res_acc_index].bound;

                         }
                         break;

                         case DonorAlphaHydrogen: //! Alpha hydrogen donor
                         {
                              alpha_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].backup_acceptor();
                              alpha_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].backup_donor();

                              const FPtype r_oha = (pair->hb_acc->position - pair->hb_don->position).norm();

                              if (r_oha < 4.0) {

                                   const FPtype theta = calc_angle(pair->hb_don->position,
                                                                   pair->hb_acc->position,
                                                                   pair->hb_acc_second->position);

                                   const FPtype rho = calc_dihedral(pair->hb_don->position,
                                                                    pair->hb_acc->position,
                                                                    pair->hb_acc_second->position,
                                                                    pair->hb_acc_third->position);

                               const boost::multi_array<float, 1> data_vector = get_bb_bb_alphah_oxygen_data_vector(r_oha,theta,rho);


                                   //primary HA hbond terms.
                                   alpha_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction[0]          = (FPtype)data_vector[5]; // HA, n
                              //     alpha_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction[1]          = (FPtype)data_vector[0]; // CA, n
                                   alpha_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction[2]          = (FPtype)data_vector[4]; // H,  n
                                   //alpha_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction[3]          = (FPtype)data_vector[3]; // N,  n
                                   //alpha_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction[4]          = (FPtype)data_vector[2]; // C,  n
                              //     alpha_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction[5]          = (FPtype)data_vector[1]; // CB, n

                                   //secondary HA hbond terms.
                                   alpha_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction[0]             = (FPtype)data_vector[11]; // HA, n
                                   alpha_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction[1]             = (FPtype)data_vector[ 6]; // CA, n
                                   alpha_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_next_res[2]    = (FPtype)data_vector[ 10]; // H,  n+1
                                   alpha_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_next_res[3]    = (FPtype)data_vector[ 9]; // N,  n+1
                                   alpha_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction[4]             = (FPtype)data_vector[ 8]; // C,  n

                              } else{

                                   alpha_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].clear_acceptor();
                                   alpha_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].clear_donor();
                              }

                              //! add to non bonded table

                              unsigned int chain_length = non_bonded_chemical_shifts_table.size() - 1;

                              for (unsigned int i = 0; i < 6; i++) {

                                   non_bonded_chemical_shifts_table[res_don_index][i] -= alpha_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction_old[i];
                                   non_bonded_chemical_shifts_table[res_don_index][i] += alpha_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction[i];

                                   if (res_don_index != 0){
                                        non_bonded_chemical_shifts_table[res_don_index-1][i] -= alpha_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction_prev_res_old[i];
                                        non_bonded_chemical_shifts_table[res_don_index-1][i] += alpha_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction_prev_res[i];
                                   }
                                   if (res_don_index != chain_length){
                                        non_bonded_chemical_shifts_table[res_don_index+1][i] -= alpha_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction_next_res_old[i];
                                        non_bonded_chemical_shifts_table[res_don_index+1][i] += alpha_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction_next_res[i];
                                   }

                                   if (i != 6){

                                   non_bonded_chemical_shifts_table[res_acc_index][i] -= alpha_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_old[i];
                                   non_bonded_chemical_shifts_table[res_acc_index][i] += alpha_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction[i];
                                   if (res_acc_index != 0){
                                        non_bonded_chemical_shifts_table[res_acc_index-1][i] -= alpha_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_prev_res_old[i];
                                        non_bonded_chemical_shifts_table[res_acc_index-1][i] += alpha_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_prev_res[i];
                                   }
                                   if (res_acc_index != chain_length){
                                        non_bonded_chemical_shifts_table[res_acc_index+1][i] -= alpha_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_next_res_old[i];
                                        non_bonded_chemical_shifts_table[res_acc_index+1][i] += alpha_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_next_res[i];
                                   }

                                   }
                              }

                         break;
                         } //end DonorAlphaHydrogen
 
                    default:
                         break;
                    } //end switch
                    ModifiedPair this_pair;
                    this_pair.don_index = res_don_index;
                    this_pair.acc_index = res_acc_index;
                    this_pair.acceptor_type = pair->acceptor_type; 
                    this_pair.donor_type = pair->donor_type;
                    bbdon_bbacc_modified_pairs.push_back(this_pair);
               }//end for

               //HA3
               ha3don_bbacc_modified_pairs.clear();
               pairs = get_modified_ha3bb_hbpair(start_index, end_index);
               for (std::vector< HbPair>::iterator pair = pairs.begin(); pair < pairs.end(); pair++) {

                    const unsigned int res_don_index = pair->hb_don->residue->index;
                    const unsigned int res_acc_index = pair->hb_acc->residue->index;

                    ha3_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].backup_acceptor();
                    ha3_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].backup_donor();

                    const FPtype r_oha = (pair->hb_acc->position - pair->hb_don->position).norm();

                    if (r_oha < 4.0) {

                         const FPtype theta = calc_angle(pair->hb_don->position,
                                                         pair->hb_acc->position,
                                                         pair->hb_acc_second->position);

                         const FPtype rho = calc_dihedral(pair->hb_don->position,
                                                          pair->hb_acc->position,
                                                          pair->hb_acc_second->position,
                                                          pair->hb_acc_third->position);

                     const boost::multi_array<float, 1> data_vector = get_bb_bb_alphah_oxygen_data_vector(r_oha,theta,rho);


                         //primary HA hbond terms.
                    //     ha3_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction[1]          = (FPtype)data_vector[0]; // CA, n
                         ha3_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction[2]          = (FPtype)data_vector[4]; // H,  n
                         //ha3_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction[3]          = (FPtype)data_vector[3]; // N,  n
                         //ha3_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction[4]          = (FPtype)data_vector[2]; // C,  n
                    //     ha3_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction[5]          = (FPtype)data_vector[1]; // CB, n

                         //secondary HA hbond terms.
                         ha3_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction[0]             = (FPtype)data_vector[11]; // HA, n
                         ha3_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction[1]             = (FPtype)data_vector[ 6]; // CA, n
                         ha3_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_next_res[2]    = (FPtype)data_vector[ 10]; // H,  n+1
                         ha3_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_next_res[3]    = (FPtype)data_vector[ 9]; // N,  n+1
                         ha3_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction[4]             = (FPtype)data_vector[ 8]; // C,  n

                    } else{
                         ha3_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].clear_acceptor();
                         ha3_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].clear_donor();
                    }

                    //! add to non bonded table

                    unsigned int chain_length = non_bonded_chemical_shifts_table.size() - 1;

                    for (unsigned int i = 0; i < 6; i++) {

                         //Don't include HA
                         if (i != 0){
                         non_bonded_chemical_shifts_table[res_don_index][i] -= ha3_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction_old[i];
                         non_bonded_chemical_shifts_table[res_don_index][i] += ha3_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction[i];

                         if (res_don_index != 0){
                              non_bonded_chemical_shifts_table[res_don_index-1][i] -= ha3_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction_prev_res_old[i];
                              non_bonded_chemical_shifts_table[res_don_index-1][i] += ha3_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction_prev_res[i];
                         }
                         if (res_don_index != chain_length){
                              non_bonded_chemical_shifts_table[res_don_index+1][i] -= ha3_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction_next_res_old[i];
                              non_bonded_chemical_shifts_table[res_don_index+1][i] += ha3_oxygen_interaction_matrix_bbbb[res_don_index][res_acc_index].donor_interaction_next_res[i];
                         }
                         }

                         //Don't include CB
                         if (i != 6){

                         non_bonded_chemical_shifts_table[res_acc_index][i] -= ha3_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_old[i];
                         non_bonded_chemical_shifts_table[res_acc_index][i] += ha3_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction[i];
                         if (res_acc_index != 0){
                              non_bonded_chemical_shifts_table[res_acc_index-1][i] -= ha3_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_prev_res_old[i];
                              non_bonded_chemical_shifts_table[res_acc_index-1][i] += ha3_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_prev_res[i];
                         }
                         if (res_acc_index != chain_length){
                              non_bonded_chemical_shifts_table[res_acc_index+1][i] -= ha3_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_next_res_old[i];
                              non_bonded_chemical_shifts_table[res_acc_index+1][i] += ha3_oxygen_interaction_matrix_bbbb[res_acc_index][res_don_index].acceptor_interaction_next_res[i];
                         }

                         }
                    }
                    ModifiedPair this_pair;
                    this_pair.don_index = res_don_index;
                    this_pair.acc_index = res_acc_index;
                    this_pair.acceptor_type = pair->acceptor_type; 
                    this_pair.donor_type = pair->donor_type;
                    bbdon_bbacc_modified_pairs.push_back(this_pair);
               }//end HA3
          }//end if
          else{ //! if only sidechain modification
               bbdon_bbacc_modified_pairs.clear();
               ha3don_bbacc_modified_pairs.clear();
          }


          //! Side chain donor backbone acceptor and backbone donor side chain acceptor
          if ((!move_info) || !(move_info->move_type == definitions::SIDECHAIN && !hydrogen_bonding(move_info))) {

               scdon_bbacc_modified_pairs.clear();
               std::vector< HbPair > pairs1 = get_modified_scbb_hbpair(start_index, end_index);

     
               for (std::vector< HbPair >::iterator pair = pairs1.begin(), pairs1_end = pairs1.end(); //! Side chain donor backbone acceptor
                    pair != pairs1_end; ++pair) {

                    const unsigned int res_don_index = pair->hb_don->residue->index;
                    const unsigned int res_acc_index = pair->hb_acc->residue->index;
                    const unsigned int other_don_index = pair->other_don_index;

                    sidehchain_amide_interaction_matrix[res_acc_index][res_don_index].backup_acceptor(pair->other_don_index);

                    const FPtype r_oh = (pair->hb_acc->position - pair->hb_don->position).norm();

                    if (r_oh < 3.0) {

                         const FPtype theta = calc_angle(pair->hb_don->position,
                                                         pair->hb_acc->position,
                                                         pair->hb_acc_second->position);

                         const FPtype rho = calc_dihedral(pair->hb_don->position,
                                                          pair->hb_acc->position,
                                                          pair->hb_acc_second->position,
                                                          pair->hb_acc_third->position);

                         const boost::multi_array<float, 1> data_vector = get_bb_bb_data_vector(r_oh, theta, rho);



                         //secondary H hbond.
                         sidehchain_amide_interaction_matrix[res_acc_index][res_don_index].interactions[pair->other_don_index].acceptor_interaction[0]          = (FPtype)data_vector[11];// HA, n
                         sidehchain_amide_interaction_matrix[res_acc_index][res_don_index].interactions[pair->other_don_index].acceptor_interaction[1]          = (FPtype)data_vector[ 6];// CA, n
                         sidehchain_amide_interaction_matrix[res_acc_index][res_don_index].interactions[pair->other_don_index].acceptor_interaction_next_res[2] = (FPtype)data_vector[ 10];// Hn  n+1
                         sidehchain_amide_interaction_matrix[res_acc_index][res_don_index].interactions[pair->other_don_index].acceptor_interaction_next_res[3] = (FPtype)data_vector[ 9];// N,  n+1
                         sidehchain_amide_interaction_matrix[res_acc_index][res_don_index].interactions[pair->other_don_index].acceptor_interaction[4]          = (FPtype)data_vector[ 8];// C,  n
                    }
                    else{
                         sidehchain_amide_interaction_matrix[res_acc_index][res_don_index].clear_acceptor(pair->other_don_index);
                    }
                    //! add results to non bonded
                    unsigned int chain_length = non_bonded_chemical_shifts_table.size() - 1;

                    for (unsigned int i = 0; i < 5; i++) {
                         non_bonded_chemical_shifts_table[res_acc_index][i] -= 
                                           sidehchain_amide_interaction_matrix[res_acc_index][res_don_index].interactions[other_don_index].acceptor_interaction_old[i];

                         non_bonded_chemical_shifts_table[res_acc_index][i] += 
                                           sidehchain_amide_interaction_matrix[res_acc_index][res_don_index].interactions[other_don_index].acceptor_interaction[i];

                         if (res_acc_index != 0){
                              non_bonded_chemical_shifts_table[res_acc_index-1][i] -= 
                                                sidehchain_amide_interaction_matrix[res_acc_index][res_don_index].interactions[other_don_index].acceptor_interaction_prev_res_old[i];

                              non_bonded_chemical_shifts_table[res_acc_index-1][i] += 
                                                sidehchain_amide_interaction_matrix[res_acc_index][res_don_index].interactions[other_don_index].acceptor_interaction_prev_res[i];
                         }
                         if (res_acc_index != chain_length){
                              non_bonded_chemical_shifts_table[res_acc_index+1][i] -= 
                                                sidehchain_amide_interaction_matrix[res_acc_index][res_don_index].interactions[other_don_index].acceptor_interaction_next_res_old[i];

                              non_bonded_chemical_shifts_table[res_acc_index+1][i] += 
                                                sidehchain_amide_interaction_matrix[res_acc_index][res_don_index].interactions[other_don_index].acceptor_interaction_next_res[i];
                          }

                    }

                    ModifiedPair this_pair;
                    this_pair.don_index = res_don_index;
                    this_pair.acc_index = res_acc_index;
                    this_pair.other_don_index = pair->other_don_index;
                    scdon_bbacc_modified_pairs.push_back(this_pair);

               }//end for


               bbdon_scacc_modified_pairs.clear();
               std::vector< HbPair > pairs2 = get_modified_bbsc_hbpair(start_index, end_index);
               for (std::vector< HbPair >::iterator pair = pairs2.begin(), pairs2_end = pairs2.end(); //! backbone donor side chain acceptor
                    pair != pairs2_end; ++pair) {

                    const unsigned int res_don_index = pair->hb_don->residue->index;
                    const unsigned int res_acc_index = pair->hb_acc->residue->index;

                    const unsigned int other_acc_index =  pair->other_acc_index;

                    switch (pair->donor_type){ //! donor type switch
                         default:
                              break;
                         case DonorAmide:
                              { //! Amide proton donor
                              amide_sidechain_interaction_matrix[res_don_index][res_acc_index].backup_donor(other_acc_index);
                              const FPtype r_oh = (pair->hb_acc->position - pair->hb_don->position).norm();

                              if (r_oh < 3.0) {

                                   const FPtype theta = calc_angle(pair->hb_don->position,
                                                                   pair->hb_acc->position,
                                                                   pair->hb_acc_second->position);

                                   const FPtype rho = calc_dihedral(pair->hb_don->position,
                                                                    pair->hb_acc->position,
                                                                    pair->hb_acc_second->position,
                                                                    pair->hb_acc_third->position);

                                   amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].bound = true;

                                   switch (pair->acceptor_type){
                                        default:
                                             break;
                                        case AcceptorAmide:
                                             {
                                             const boost::multi_array<float, 1> data_vector = get_bb_bb_data_vector(r_oh, theta, rho);


                                             //primary H hbond.
                                            amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[0]           = (FPtype)data_vector[5]; // HA, n
                                            // amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[1]           = (FPtype)data_vector[0]; // CA, n
                                             amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[2]           = (FPtype)data_vector[4]; // H,  n
                                             amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[3]           = (FPtype)data_vector[3]; // N,  n
                                           //  amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction_prev_res[4]  = (FPtype)data_vector[2]; // C,  n-1
                                             break;
                                        }
                                        case AcceptorAlcohol:
                                             {

                                             const boost::multi_array<float, 1> data_vector = get_bb_sc_alcohol_data_vector(r_oh, theta, rho);
                                           // primary H hbond.

                                             amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[0]           = (FPtype)data_vector[5]; // HA, n
                                           //  amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[1]           = (FPtype)data_vector[0]; // CA,  n
                                             amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[2]           = (FPtype)data_vector[4]; // H,  n
                                             amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[3]           = (FPtype)data_vector[3]; // N,  n
                                           //  amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction_prev_res[4]  = (FPtype)data_vector[2]; // C-1 n
                                             break;
                                        }

                                        case AcceptorCarboxylate:
                                             {
                                             const boost::multi_array<float, 1> data_vector = get_bb_sc_carboxy_data_vector(r_oh, theta, rho);


                                           //primary H hbond
                                             amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[0]           = (FPtype)data_vector[5]; // HA, n
                                           //  amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[1]           = (FPtype)data_vector[0]; // CA,  n
                                             amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[2]           = (FPtype)data_vector[4]; // H,  n
                                             amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[3]           = (FPtype)data_vector[3]; // N,  n
                                           //  amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction_prev_res[4]  = (double)data_vector[2]; // C-1 n
                                             break;
                                        }
                                   }

                              } else {
                                   amide_sidechain_interaction_matrix[res_don_index][res_acc_index].clear_donor(other_acc_index);
                                   amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].bound = false;
                              }

                              int chain_length = non_bonded_chemical_shifts_table.size() - 1;

                              //! add to non bonded
                              for (unsigned int i = 0; i < 5; i++) {
                                   non_bonded_chemical_shifts_table[res_don_index][i] -= 
                                                       amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction_old[i];

                                   non_bonded_chemical_shifts_table[res_don_index][i] += 
                                                       amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[i];

                                   if (res_don_index != 0){
                                        non_bonded_chemical_shifts_table[res_don_index-1][i] -= 
                                                            amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction_prev_res_old[i];

                                        non_bonded_chemical_shifts_table[res_don_index-1][i] +=
                                                             amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction_prev_res[i];
                                   }
                                   if (res_don_index != chain_length){
                                        non_bonded_chemical_shifts_table[res_don_index+1][i] -= 
                                                            amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction_next_res_old[i];

                                        non_bonded_chemical_shifts_table[res_don_index+1][i] += 
                                                            amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction_next_res[i];
                                   }
                              }

                              //! add state of hydrogen bond pair to amide_bond_count_table
                              amide_bond_count_table[res_don_index] -= amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].bound_old;
                              amide_bond_count_table[res_don_index] += amide_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].bound;

                         } //! end amide proton donor
                         break;
                         case DonorAlphaHydrogen:
                              { //! alpha hydrogen donor

                              alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].backup_donor(other_acc_index);
                              const FPtype r_oha = (pair->hb_acc->position - pair->hb_don->position).norm();                    

                              if (r_oha < 4.0) {

                                   const FPtype theta = calc_angle(pair->hb_don->position,
                                                                   pair->hb_acc->position,
                                                                   pair->hb_acc_second->position);

                                   const FPtype rho = calc_dihedral(pair->hb_don->position,
                                                                    pair->hb_acc->position,
                                                                    pair->hb_acc_second->position,
                                                                    pair->hb_acc_third->position);

                                   switch (pair->acceptor_type){
                                        default:
                                             break;
                                        case AcceptorAmide:
                                             {

                                             const boost::multi_array<float, 1> data_vector = get_bb_bb_alphah_oxygen_data_vector(r_oha, theta, rho);


                                //primary HA hbond.
                                             alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[0]           = (FPtype)data_vector[5]; // HA, n
                                   //          alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[1]           = (FPtype)data_vector[0]; // CA, n
                                             alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[2]           = (FPtype)data_vector[4]; // H,  n
//                                             alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[3]           = (FPtype)data_vector[3]; // N,  n
                                   //          alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[4]           = (FPtype)data_vector[2]; // C,  n
                                   //          alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[5]           = (FPtype)data_vector[1]; // CB,  n
                                        }
                                        break;


                                        case AcceptorAlcohol: 
                                             {

                                             const boost::multi_array<float, 1> data_vector = get_bb_bb_alphah_alcohol_data_vector(r_oha, theta, rho);

                                   //primary HA hbond
                                             alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[0]           = (FPtype)data_vector[5]; // HA, n
                                     //        alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[1]           = (FPtype)data_vector[0]; // CA, n
                                             alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[2]           = (FPtype)data_vector[4]; // H,  n
                                     //        alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[3]           = (FPtype)data_vector[3]; // N,  n
                                     //        alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[4]           = (FPtype)data_vector[2]; // C,  n                           
                                     //        alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[5]           = (FPtype)data_vector[1]; // CB,  n                                         

                                        }
                                        break;



                                        case AcceptorCarboxylate:
                                             {

                                             const boost::multi_array<float, 1> data_vector = get_bb_bb_alphah_carboxy_data_vector(r_oha, theta, rho);

                                 //primary HA
                                             alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[0]           = (FPtype)data_vector[5]; // HA, n
                                     //        alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[1]           = (FPtype)data_vector[0]; // CA, n
                                             alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[2]           = (FPtype)data_vector[4]; // H,  n
//                                             alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[3]           = (FPtype)data_vector[3]; // N,  n
                                     //        alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[4]           = (FPtype)data_vector[2]; // C,  n
                                     //        alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[5]           = (double)data_vector[1]; // CB,  n
                                             }
                                             break;

                                   }//end switch

                              } else {
                                   alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].clear_donor(other_acc_index);
                              }

                              int chain_length = non_bonded_chemical_shifts_table.size() - 1;

                              //! add to non bonded
                              for (unsigned int i = 0; i < 6; i++) {     
                                   non_bonded_chemical_shifts_table[res_don_index][i] -= 
                                                       alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction_old[i];

                                   non_bonded_chemical_shifts_table[res_don_index][i] += 
                                                       alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[i];

                                   if (res_don_index != 0){
                                        non_bonded_chemical_shifts_table[res_don_index-1][i] -= 
                                                            alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction_prev_res_old[i];

                                        non_bonded_chemical_shifts_table[res_don_index-1][i] +=
                                                             alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction_prev_res[i];
                                   }
                                   if (res_don_index != chain_length){
                                        non_bonded_chemical_shifts_table[res_don_index+1][i] -= 
                                                            alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction_next_res_old[i];

                                        non_bonded_chemical_shifts_table[res_don_index+1][i] += 
                                                            alpha_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction_next_res[i];
                                   }
                              }

                         } //! end alpha hydrogen donor
                         break;
                    } //! end donor type switch

                    ModifiedPair this_pair;
                    this_pair.don_index = res_don_index;
                    this_pair.acc_index = res_acc_index;
                    this_pair.other_acc_index = other_acc_index;
                    this_pair.donor_type = pair->donor_type;
                    this_pair.acceptor_type = pair->acceptor_type;
                    bbdon_scacc_modified_pairs.push_back(this_pair);

               }//end for

               //HA3 special case
               ha3don_scacc_modified_pairs.clear();
               pairs2 = get_modified_ha3sc_hbpair(start_index, end_index);
               for (std::vector< HbPair >::iterator pair = pairs2.begin(), pairs2_end = pairs2.end(); //! backbone donor side chain acceptor
                    pair != pairs2_end; ++pair) {

                    const unsigned int res_don_index = pair->hb_don->residue->index;
                    const unsigned int res_acc_index = pair->hb_acc->residue->index;

                    const unsigned int other_acc_index =  pair->other_acc_index;
                    ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].backup_donor(other_acc_index);
                    const FPtype r_oha = (pair->hb_acc->position - pair->hb_don->position).norm();

                    if (r_oha < 4.0) {

                         const FPtype theta = calc_angle(pair->hb_don->position,
                                                         pair->hb_acc->position,
                                                         pair->hb_acc_second->position);

                         const FPtype rho = calc_dihedral(pair->hb_don->position,
                                                          pair->hb_acc->position,
                                                          pair->hb_acc_second->position,
                                                          pair->hb_acc_third->position);

                         switch (pair->acceptor_type){
                              default:
                                   break;
                              case AcceptorAmide:
                                   {

                                   const boost::multi_array<float, 1> data_vector = get_bb_bb_alphah_oxygen_data_vector(r_oha, theta, rho);


                                 //primary HA
                         //          ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[1]           = (FPtype)data_vector[0]; // CA, n
                                     ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[2]           = (FPtype)data_vector[4]; // H,  n
                         //          ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[3]           = (FPtype)data_vector[3]; // N,  n
                         //          ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[4]           = (FPtype)data_vector[2]; // C,  n
                         //          ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[5]           = (FPtype)data_vector[1]; // CB,  n
                              }
                              break;


                              case AcceptorAlcohol: 
                                   {

                                   const boost::multi_array<float, 1> data_vector = get_bb_bb_alphah_alcohol_data_vector(r_oha, theta, rho);

                         //primary HA hbond
                           //        ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[1]           = (FPtype)data_vector[0]; // CA, n
                                     ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[2]           = (FPtype)data_vector[4]; // H,  n
                           //        ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[3]           = (FPtype)data_vector[3]; // N,  n
                           //        ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[4]           = (FPtype)data_vector[2]; // C,  n                           
                           //        ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[5]           = (FPtype)data_vector[1]; // CB,  n                                         

                              }
                              break;



                              case AcceptorCarboxylate:
                                   {

                                   const boost::multi_array<float, 1> data_vector = get_bb_bb_alphah_carboxy_data_vector(r_oha, theta, rho);

                       //primary HA
                           //        ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[1]           = (FPtype)data_vector[0]; // CA, n
                                     ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[2]           = (FPtype)data_vector[4]; // H,  n
                           //        ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[3]           = (FPtype)data_vector[3]; // N,  n
                           //        ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[4]           = (FPtype)data_vector[2]; // C,  n
                           //        ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[5]           = (double)data_vector[1]; // CB,  n
                                   }
                                   break;

                         }//end switch

                    } else {
                         ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].clear_donor(other_acc_index);
                    }

                    int chain_length = non_bonded_chemical_shifts_table.size() - 1;

                    //! add to non bonded
                    //Ignore HA
                    for (unsigned int i = 1; i < 6; i++) {
                         non_bonded_chemical_shifts_table[res_don_index][i] -=
                                             ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction_old[i];

                         non_bonded_chemical_shifts_table[res_don_index][i] += 
                                             ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction[i];

                         if (res_don_index != 0){
                              non_bonded_chemical_shifts_table[res_don_index-1][i] -= 
                                                  ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction_prev_res_old[i];

                              non_bonded_chemical_shifts_table[res_don_index-1][i] +=
                                                   ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction_prev_res[i];
                         }
                         if (res_don_index != chain_length){
                              non_bonded_chemical_shifts_table[res_don_index+1][i] -= 
                                                  ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction_next_res_old[i];

                              non_bonded_chemical_shifts_table[res_don_index+1][i] += 
                                                  ha3_sidechain_interaction_matrix[res_don_index][res_acc_index].interactions[other_acc_index].donor_interaction_next_res[i];
                         }
                    }

                    ModifiedPair this_pair;
                    this_pair.don_index = res_don_index;
                    this_pair.acc_index = res_acc_index;
                    this_pair.other_acc_index = other_acc_index;
                    this_pair.donor_type = pair->donor_type;
                    this_pair.acceptor_type = pair->acceptor_type;
                    ha3don_scacc_modified_pairs.push_back(this_pair);

               }//end HA3
          }
          else{
               scdon_bbacc_modified_pairs.clear();
               bbdon_scacc_modified_pairs.clear();
               ha3don_scacc_modified_pairs.clear();
          }



          //! Bonded contributions cached 
          for (ResidueIterator<ChainFB> res1(chain); !(res1).end(); ++res1) {

               if ((res1->index >= (start_index-1)) && (res1->index <= (end_index+1))){

                    bonded_interactions[res1->index].backup(); //! backup prediction

                    if ( (res1->index != 0) && (res1->get_neighbour(1) != NULL)){

                         for (std::vector<unsigned int>::iterator atype = atypevector.begin();atype != atypevector.end(); ++atype){ //loop over atoms types

                          const ResidueEnum em1 = res1->get_neighbour(-1)->residue_type;
                              const ResidueEnum e0 = res1->residue_type;
                              const ResidueEnum ep1 = res1->get_neighbour(1)->residue_type;

                              if ( loaded_tables[em1*6 + *atype] &&  loaded_tables[e0*6 + *atype] && (e0*6 + *atype != 31) && (e0*6 + *atype != 76)  && loaded_tables[ep1*6 + *atype]  ){

                                   const FPtype r1 = (list_of_tables[em1*6 + *atype]->get_values( 2, read_phi_psi_index(*res1->get_neighbour(-1)), get_chi(*res1->get_neighbour(-1)))
                                                                                           -list_of_tables[0 + *atype]->get_st(2) );

                                   const FPtype r2 = (list_of_tables[e0*6 + *atype]->get_values( 1, read_phi_psi_index(*res1), get_chi(*res1)) );


                               const FPtype r3 = (list_of_tables[ep1*6 + *atype]->get_values( 0, read_phi_psi_index(*res1->get_neighbour(1)), get_chi(*res1->get_neighbour(1)))
                                                                                           -list_of_tables[0 + *atype]->get_st(0) );


                                   const FPtype rt =  r1 + r2 + r3;

                                   bonded_interactions[res1->index].shifts[transform_id(*atype)] = -rt;
                              }

                         }
                    }
               }
          }



          if (start_index != 0){
              bonded_start = start_index - 1; 
          }
          else{
               bonded_start = start_index;
          }

          unsigned int size = non_bonded_chemical_shifts_table.size()-1;
          if (end_index != size){
               bonded_end = end_index + 1;
          }
          else{
               bonded_end = end_index;
          }

          //! add bonded contributions
          for (unsigned int j = bonded_start; j <= bonded_end; j++){
               for (unsigned int i = 0; i < 6; i++){
                    bonded_chemical_shifts_table[j][i] = bonded_interactions[j].shifts[i];
               }
          }

          //! zero non bonding terms for Proline HN
          for (unsigned int j = 0; j < non_bonded_chemical_shifts_table.size(); j++){
               if (chain[j].residue_type == PRO){
                    non_bonded_chemical_shifts_table[j][transform_id(4)] = 0.0;
               }
          }

          chemical_shifts_table.clear(); 

          chemical_shifts_table = add_tables(bonded_chemical_shifts_table,
                                          non_bonded_chemical_shifts_table);

          for (ResidueIterator<ChainFB> res1(chain); !(res1).end(); ++res1) {
               if ( amide_bond_count_table[res1->index] == 0 ){
                    chemical_shifts_table[res1->index][transform_id(ProCS15_HN)] += 2.07;
               } 
          }

          rh_pairs.clear();
          rh_pairs = get_hydrogen_ring_pairs(start_index,end_index);

          //! Ring current cached
          for (std::vector< RingHydrogenPair >::iterator pair = rh_pairs.begin();
                    pair != rh_pairs.end(); ++pair) {

               const unsigned int hydrogen = pair->hydrogen_residue_index;
               const unsigned int ring     = pair->ring_residue_index;
               const unsigned int ring_nr  = pair->ring_nr;

               std::vector<AromaticRing> aromatic_rings = get_aromatic_rings( chain[ring] );

               ring_current_interactions_new[hydrogen][ring].backup(ring_nr);

               if (pair->HAatom){
                    FPtype dist_squared = (aromatic_rings[ring_nr].ring_center - pair->HA_atom->position).norm_squared();
                    if (dist_squared < 64.0) {
                         ring_current_interactions_new[hydrogen][ring].sidechain_rings[ring_nr].ring_shifts[0]
                                                           = calc_aromatic_interaction(pair->HA_atom->position, aromatic_rings[ring_nr] , dist_squared);
                    }
                    else{
                         ring_current_interactions_new[hydrogen][ring].sidechain_rings[ring_nr].ring_shifts[0] = 0.0;
                    }
               }
               if (pair->HNatom){
                    FPtype dist_squared = (aromatic_rings[ring_nr].ring_center - pair->HN_atom->position).norm_squared();

                    if (dist_squared < 64.0) {


                         ring_current_interactions_new[hydrogen][ring].sidechain_rings[ring_nr].ring_shifts[2]
                                                           = calc_aromatic_interaction(pair->HN_atom->position, aromatic_rings[ring_nr], dist_squared);

                    }
                    else{

                         ring_current_interactions_new[hydrogen][ring].sidechain_rings[ring_nr].ring_shifts[2] = 0.0;
                    }
               }


               ring_current_table[hydrogen][0] -= ring_current_interactions_new[hydrogen][ring].sidechain_rings[ring_nr].ring_shifts_old[0]; 
               ring_current_table[hydrogen][0] += ring_current_interactions_new[hydrogen][ring].sidechain_rings[ring_nr].ring_shifts[0];

               ring_current_table[hydrogen][2] -= ring_current_interactions_new[hydrogen][ring].sidechain_rings[ring_nr].ring_shifts_old[2]; 
               ring_current_table[hydrogen][2] += ring_current_interactions_new[hydrogen][ring].sidechain_rings[ring_nr].ring_shifts[2];

          }

          //! add ring current to chemical shift tables
          for (unsigned int i = 0; i < chemical_shifts_table.size(); i++){

               chemical_shifts_table[i][2] += ring_current_table[i][2];
               chemical_shifts_table[i][0] += ring_current_table[i][0];

          }

          //! zero chemical_shifts_table at N and C terminus
          for (unsigned int o = 0; o < 6; o++){
               chemical_shifts_table[0][o] = 0;
               chemical_shifts_table[chemical_shifts_table.size()-1][o] = 0; 
          }


          //! zero chemical_shifts_table if only non bonded contribution
          for (unsigned int i = 0; i < 6; i++){
              if (loaded_data[i] == false){
                   for (unsigned int j = 0; j < chemical_shifts_table.size(); j++){
                        chemical_shifts_table[j][transform_id(i)] = 0;
                   }
              }
          }

          return chemical_shifts_table;
     }


     //! Reject last energy evaluation, and roll back changes in sigma, sigma_sigma and chi va
     void reject_cache() {

          //! bonded reject
          for (int j = bonded_start; j <= bonded_end; j++){


                    bonded_chemical_shifts_table[j] = bonded_interactions[j].shifts_old;

               bonded_interactions[j].rollback();
          }


          //! ring current cached reject
          for (std::vector< RingHydrogenPair >::iterator pair = rh_pairs.begin(), pairs1_end = rh_pairs.end();
                    pair != pairs1_end; ++pair) {

                const unsigned int hydrogen = pair->hydrogen_residue_index;
                const unsigned int ring     = pair->ring_residue_index;
                const unsigned int ring_nr  = pair->ring_nr;

                ring_current_table[hydrogen][0] += ring_current_interactions_new[hydrogen][ring].sidechain_rings[ring_nr].ring_shifts_old[0];
                ring_current_table[hydrogen][0] -= ring_current_interactions_new[hydrogen][ring].sidechain_rings[ring_nr].ring_shifts[0];

                ring_current_table[hydrogen][2] += ring_current_interactions_new[hydrogen][ring].sidechain_rings[ring_nr].ring_shifts_old[2];
                ring_current_table[hydrogen][2] -= ring_current_interactions_new[hydrogen][ring].sidechain_rings[ring_nr].ring_shifts[2];

                ring_current_interactions_new[hydrogen][ring].rollback(ring_nr);

          }

          int chain_length = non_bonded_chemical_shifts_table.size() - 1;

          //! backbone donor backbone acceptor reject
          for (std::vector< ModifiedPair >::iterator pair = bbdon_bbacc_modified_pairs.begin(); pair != bbdon_bbacc_modified_pairs.end(); ++pair) {
               switch (pair->donor_type){
                    case DonorAmide:

                         for (unsigned int i = 0; i < 5; i++) {

                              non_bonded_chemical_shifts_table[pair->don_index][i] += amide_amide_interaction_matrix_bbbb[pair->don_index][pair->acc_index].donor_interaction_old[i];
                              non_bonded_chemical_shifts_table[pair->don_index][i] -= amide_amide_interaction_matrix_bbbb[pair->don_index][pair->acc_index].donor_interaction[i];
                              if (pair->don_index != chain_length){
                                   non_bonded_chemical_shifts_table[pair->don_index+1][i] += amide_amide_interaction_matrix_bbbb[pair->don_index][pair->acc_index].donor_interaction_next_res_old[i];
                                   non_bonded_chemical_shifts_table[pair->don_index+1][i] -= amide_amide_interaction_matrix_bbbb[pair->don_index][pair->acc_index].donor_interaction_next_res[i];
                              }
                              if (pair->don_index != 0){ 
                                   non_bonded_chemical_shifts_table[pair->don_index-1][i] += amide_amide_interaction_matrix_bbbb[pair->don_index][pair->acc_index].donor_interaction_prev_res_old[i];
                                   non_bonded_chemical_shifts_table[pair->don_index-1][i] -= amide_amide_interaction_matrix_bbbb[pair->don_index][pair->acc_index].donor_interaction_prev_res[i];
                              }

                              non_bonded_chemical_shifts_table[pair->acc_index][i] += amide_amide_interaction_matrix_bbbb[pair->acc_index][pair->don_index].acceptor_interaction_old[i];
                              non_bonded_chemical_shifts_table[pair->acc_index][i] -= amide_amide_interaction_matrix_bbbb[pair->acc_index][pair->don_index].acceptor_interaction[i]; 
                              if (pair->acc_index != chain_length){
                                   non_bonded_chemical_shifts_table[pair->acc_index+1][i] += amide_amide_interaction_matrix_bbbb[pair->acc_index][pair->don_index].acceptor_interaction_next_res_old[i];
                                   non_bonded_chemical_shifts_table[pair->acc_index+1][i] -= amide_amide_interaction_matrix_bbbb[pair->acc_index][pair->don_index].acceptor_interaction_next_res[i]; 
                              }
                              if (pair->acc_index != 0){
                                   non_bonded_chemical_shifts_table[pair->acc_index-1][i] += amide_amide_interaction_matrix_bbbb[pair->acc_index][pair->don_index].acceptor_interaction_prev_res_old[i];
                                   non_bonded_chemical_shifts_table[pair->acc_index-1][i] -= amide_amide_interaction_matrix_bbbb[pair->acc_index][pair->don_index].acceptor_interaction_prev_res[i];
                              }

                         }


                         amide_bond_count_table[pair->don_index] += amide_amide_interaction_matrix_bbbb[pair->don_index][pair->acc_index].bound_old;
                         amide_bond_count_table[pair->don_index] -= amide_amide_interaction_matrix_bbbb[pair->don_index][pair->acc_index].bound;

                         amide_amide_interaction_matrix_bbbb[pair->don_index][pair->acc_index].rollback_donor();
                         amide_amide_interaction_matrix_bbbb[pair->acc_index][pair->don_index].rollback_acceptor();
                         break;

                    case DonorAlphaHydrogen:

                         for (unsigned int i = 0; i < 6; i++) {

                              non_bonded_chemical_shifts_table[pair->don_index][i] += alpha_oxygen_interaction_matrix_bbbb[pair->don_index][pair->acc_index].donor_interaction_old[i];
                              non_bonded_chemical_shifts_table[pair->don_index][i] -= alpha_oxygen_interaction_matrix_bbbb[pair->don_index][pair->acc_index].donor_interaction[i];
                              if (pair->don_index != chain_length){
                                   non_bonded_chemical_shifts_table[pair->don_index+1][i] += alpha_oxygen_interaction_matrix_bbbb[pair->don_index][pair->acc_index].donor_interaction_next_res_old[i];
                                   non_bonded_chemical_shifts_table[pair->don_index+1][i] -= alpha_oxygen_interaction_matrix_bbbb[pair->don_index][pair->acc_index].donor_interaction_next_res[i];
                              }
                              if (pair->don_index != 0){ 
                                   non_bonded_chemical_shifts_table[pair->don_index-1][i] += alpha_oxygen_interaction_matrix_bbbb[pair->don_index][pair->acc_index].donor_interaction_prev_res_old[i];
                                   non_bonded_chemical_shifts_table[pair->don_index-1][i] -= alpha_oxygen_interaction_matrix_bbbb[pair->don_index][pair->acc_index].donor_interaction_prev_res[i];
                              }

                              non_bonded_chemical_shifts_table[pair->acc_index][i] += alpha_oxygen_interaction_matrix_bbbb[pair->acc_index][pair->don_index].acceptor_interaction_old[i];
                              non_bonded_chemical_shifts_table[pair->acc_index][i] -= alpha_oxygen_interaction_matrix_bbbb[pair->acc_index][pair->don_index].acceptor_interaction[i]; 
                              if (pair->acc_index != chain_length){            
                                   non_bonded_chemical_shifts_table[pair->acc_index+1][i] += alpha_oxygen_interaction_matrix_bbbb[pair->acc_index][pair->don_index].acceptor_interaction_next_res_old[i];
                                   non_bonded_chemical_shifts_table[pair->acc_index+1][i] -= alpha_oxygen_interaction_matrix_bbbb[pair->acc_index][pair->don_index].acceptor_interaction_next_res[i]; 
                              }
                              if (pair->acc_index != 0){
                                   non_bonded_chemical_shifts_table[pair->acc_index-1][i] += alpha_oxygen_interaction_matrix_bbbb[pair->acc_index][pair->don_index].acceptor_interaction_prev_res_old[i];
                                   non_bonded_chemical_shifts_table[pair->acc_index-1][i] -= alpha_oxygen_interaction_matrix_bbbb[pair->acc_index][pair->don_index].acceptor_interaction_prev_res[i];
                              }

                         }
                         alpha_oxygen_interaction_matrix_bbbb[pair->don_index][pair->acc_index].rollback_donor();
                         alpha_oxygen_interaction_matrix_bbbb[pair->acc_index][pair->don_index].rollback_acceptor();
                         break;

                    default:
                         break;
               }
          }

          //! HA3 donor backbone acceptor reject
          for (std::vector< ModifiedPair >::iterator pair = ha3don_bbacc_modified_pairs.begin(); pair != ha3don_bbacc_modified_pairs.end(); ++pair) {
               for (unsigned int i = 0; i < 6; i++) {

                    non_bonded_chemical_shifts_table[pair->don_index][i] += ha3_oxygen_interaction_matrix_bbbb[pair->don_index][pair->acc_index].donor_interaction_old[i];
                    non_bonded_chemical_shifts_table[pair->don_index][i] -= ha3_oxygen_interaction_matrix_bbbb[pair->don_index][pair->acc_index].donor_interaction[i];
                    if (pair->don_index != chain_length){
                         non_bonded_chemical_shifts_table[pair->don_index+1][i] += ha3_oxygen_interaction_matrix_bbbb[pair->don_index][pair->acc_index].donor_interaction_next_res_old[i];
                         non_bonded_chemical_shifts_table[pair->don_index+1][i] -= ha3_oxygen_interaction_matrix_bbbb[pair->don_index][pair->acc_index].donor_interaction_next_res[i];
                    }
                    if (pair->don_index != 0){ 
                         non_bonded_chemical_shifts_table[pair->don_index-1][i] += ha3_oxygen_interaction_matrix_bbbb[pair->don_index][pair->acc_index].donor_interaction_prev_res_old[i];
                         non_bonded_chemical_shifts_table[pair->don_index-1][i] -= ha3_oxygen_interaction_matrix_bbbb[pair->don_index][pair->acc_index].donor_interaction_prev_res[i];
                    }

                    non_bonded_chemical_shifts_table[pair->acc_index][i] += ha3_oxygen_interaction_matrix_bbbb[pair->acc_index][pair->don_index].acceptor_interaction_old[i];
                    non_bonded_chemical_shifts_table[pair->acc_index][i] -= ha3_oxygen_interaction_matrix_bbbb[pair->acc_index][pair->don_index].acceptor_interaction[i]; 
                    if (pair->acc_index != chain_length){
                         non_bonded_chemical_shifts_table[pair->acc_index+1][i] += ha3_oxygen_interaction_matrix_bbbb[pair->acc_index][pair->don_index].acceptor_interaction_next_res_old[i];
                         non_bonded_chemical_shifts_table[pair->acc_index+1][i] -= ha3_oxygen_interaction_matrix_bbbb[pair->acc_index][pair->don_index].acceptor_interaction_next_res[i]; 
                    }
                    if (pair->acc_index != 0){
                         non_bonded_chemical_shifts_table[pair->acc_index-1][i] += ha3_oxygen_interaction_matrix_bbbb[pair->acc_index][pair->don_index].acceptor_interaction_prev_res_old[i];
                         non_bonded_chemical_shifts_table[pair->acc_index-1][i] -= ha3_oxygen_interaction_matrix_bbbb[pair->acc_index][pair->don_index].acceptor_interaction_prev_res[i];
                    }

               }
               ha3_oxygen_interaction_matrix_bbbb[pair->don_index][pair->acc_index].rollback_donor();
               ha3_oxygen_interaction_matrix_bbbb[pair->acc_index][pair->don_index].rollback_acceptor();
          }

          //! sidechain donor backbone acceptor reject
          for (std::vector< ModifiedPair >::iterator pair = scdon_bbacc_modified_pairs.begin(); pair != scdon_bbacc_modified_pairs.end(); ++pair) {

               unsigned int donator_index = pair->other_don_index;

               for (unsigned int i = 0; i < 5; i++) {

                    non_bonded_chemical_shifts_table[pair->acc_index][i] += 
                                        sidehchain_amide_interaction_matrix[pair->acc_index][pair->don_index].interactions[donator_index].acceptor_interaction_old[i];

                    non_bonded_chemical_shifts_table[pair->acc_index][i] -= 
                                        sidehchain_amide_interaction_matrix[pair->acc_index][pair->don_index].interactions[donator_index].acceptor_interaction[i];

                    if (pair->acc_index != chain_length){ 
                         non_bonded_chemical_shifts_table[pair->acc_index+1][i] += 
                                        sidehchain_amide_interaction_matrix[pair->acc_index][pair->don_index].interactions[donator_index].acceptor_interaction_next_res_old[i];

                         non_bonded_chemical_shifts_table[pair->acc_index+1][i] -= 
                                        sidehchain_amide_interaction_matrix[pair->acc_index][pair->don_index].interactions[donator_index].acceptor_interaction_next_res[i];
                    }
                    if (pair->acc_index != 0){
                         non_bonded_chemical_shifts_table[pair->acc_index-1][i] += 
                                        sidehchain_amide_interaction_matrix[pair->acc_index][pair->don_index].interactions[donator_index].acceptor_interaction_prev_res_old[i];

                         non_bonded_chemical_shifts_table[pair->acc_index-1][i] -= 
                                        sidehchain_amide_interaction_matrix[pair->acc_index][pair->don_index].interactions[donator_index].acceptor_interaction_prev_res[i];
                    }

               }
               sidehchain_amide_interaction_matrix[pair->acc_index][pair->don_index].rollback_acceptor(donator_index);

          }

          //! backbone donor sidechain acceptor reject
          for (std::vector< ModifiedPair >::iterator pair = bbdon_scacc_modified_pairs.begin(); pair != bbdon_scacc_modified_pairs.end(); ++pair) {

               unsigned int acceptor_index = pair->other_acc_index;
               switch (pair->donor_type){ //! donor type switch
                    default:
                         break;
                    case DonorAmide:
                         for (unsigned int i = 0; i < 5; i++) {

                              non_bonded_chemical_shifts_table[pair->don_index][i] += 
                                                  amide_sidechain_interaction_matrix[pair->don_index][pair->acc_index].interactions[acceptor_index].donor_interaction_old[i];
                              non_bonded_chemical_shifts_table[pair->don_index][i] -= 
                                                  amide_sidechain_interaction_matrix[pair->don_index][pair->acc_index].interactions[acceptor_index].donor_interaction[i];

                              if (pair->don_index != chain_length){
                                   non_bonded_chemical_shifts_table[pair->don_index+1][i] += 
                                                       amide_sidechain_interaction_matrix[pair->don_index][pair->acc_index].interactions[acceptor_index].donor_interaction_next_res_old[i];

                                   non_bonded_chemical_shifts_table[pair->don_index+1][i] -= 
                                                       amide_sidechain_interaction_matrix[pair->don_index][pair->acc_index].interactions[acceptor_index].donor_interaction_next_res[i];
                              }
 
                              if (pair->don_index != 0){
                                   non_bonded_chemical_shifts_table[pair->don_index-1][i] += 
                                                       amide_sidechain_interaction_matrix[pair->don_index][pair->acc_index].interactions[acceptor_index].donor_interaction_prev_res_old[i];

                                   non_bonded_chemical_shifts_table[pair->don_index-1][i] -= 
                                                       amide_sidechain_interaction_matrix[pair->don_index][pair->acc_index].interactions[acceptor_index].donor_interaction_prev_res[i];
                              }

                         }

                         amide_bond_count_table[pair->don_index] += amide_sidechain_interaction_matrix[pair->don_index][pair->acc_index].interactions[acceptor_index].bound_old;
                         amide_bond_count_table[pair->don_index] -= amide_sidechain_interaction_matrix[pair->don_index][pair->acc_index].interactions[acceptor_index].bound;

                         amide_sidechain_interaction_matrix[pair->don_index][pair->acc_index].rollback_donor(acceptor_index);
                         break;

                    case DonorAlphaHydrogen:
                         for (unsigned int i = 0; i < 6; i++) {

                              non_bonded_chemical_shifts_table[pair->don_index][i] += 
                                                  alpha_sidechain_interaction_matrix[pair->don_index][pair->acc_index].interactions[acceptor_index].donor_interaction_old[i];
                              non_bonded_chemical_shifts_table[pair->don_index][i] -= 
                                                  alpha_sidechain_interaction_matrix[pair->don_index][pair->acc_index].interactions[acceptor_index].donor_interaction[i];

                              if (pair->don_index != chain_length){
                                   non_bonded_chemical_shifts_table[pair->don_index+1][i] += 
                                                       alpha_sidechain_interaction_matrix[pair->don_index][pair->acc_index].interactions[acceptor_index].donor_interaction_next_res_old[i];

                                   non_bonded_chemical_shifts_table[pair->don_index+1][i] -= 
                                                       alpha_sidechain_interaction_matrix[pair->don_index][pair->acc_index].interactions[acceptor_index].donor_interaction_next_res[i];
                              }
 
                              if (pair->don_index != 0){
                                   non_bonded_chemical_shifts_table[pair->don_index-1][i] += 
                                                       alpha_sidechain_interaction_matrix[pair->don_index][pair->acc_index].interactions[acceptor_index].donor_interaction_prev_res_old[i];

                                   non_bonded_chemical_shifts_table[pair->don_index-1][i] -= 
                                                       alpha_sidechain_interaction_matrix[pair->don_index][pair->acc_index].interactions[acceptor_index].donor_interaction_prev_res[i];
                              }
                         }
                         alpha_sidechain_interaction_matrix[pair->don_index][pair->acc_index].rollback_donor(acceptor_index);
                         break;

               } //! donor type end
          }

          //! HA3 donor sidechain acceptor reject
          for (std::vector< ModifiedPair >::iterator pair = ha3don_scacc_modified_pairs.begin(); pair != ha3don_scacc_modified_pairs.end(); ++pair) {

               unsigned int acceptor_index = pair->other_acc_index;
               //skip HA
               for (unsigned int i = 1; i < 6; i++) {

                    non_bonded_chemical_shifts_table[pair->don_index][i] += 
                                        ha3_sidechain_interaction_matrix[pair->don_index][pair->acc_index].interactions[acceptor_index].donor_interaction_old[i];
                    non_bonded_chemical_shifts_table[pair->don_index][i] -= 
                                        ha3_sidechain_interaction_matrix[pair->don_index][pair->acc_index].interactions[acceptor_index].donor_interaction[i];

                    if (pair->don_index != chain_length){
                         non_bonded_chemical_shifts_table[pair->don_index+1][i] += 
                                             ha3_sidechain_interaction_matrix[pair->don_index][pair->acc_index].interactions[acceptor_index].donor_interaction_next_res_old[i];

                         non_bonded_chemical_shifts_table[pair->don_index+1][i] -= 
                                             ha3_sidechain_interaction_matrix[pair->don_index][pair->acc_index].interactions[acceptor_index].donor_interaction_next_res[i];
                    }
 
                    if (pair->don_index != 0){
                         non_bonded_chemical_shifts_table[pair->don_index-1][i] += 
                                             ha3_sidechain_interaction_matrix[pair->don_index][pair->acc_index].interactions[acceptor_index].donor_interaction_prev_res_old[i];

                         non_bonded_chemical_shifts_table[pair->don_index-1][i] -= 
                                             ha3_sidechain_interaction_matrix[pair->don_index][pair->acc_index].interactions[acceptor_index].donor_interaction_prev_res[i];
                    }
               }
               ha3_sidechain_interaction_matrix[pair->don_index][pair->acc_index].rollback_donor(acceptor_index);
          }

     }


     void NMR_STAR_format(phaistos::ChainFB& chain){

          //! print NMR star format

          int atom_count = 1;

          for (unsigned int i = 0; i < chemical_shifts_table.size(); i++){
               for (unsigned int j = 0; j < 6; j++){
                    switch (j){
                         case 0: //!CA
                              if (loaded_tables[0] == true){
                                   std::cout << atom_count << "\t" << (i+1) << "\t" << chain[i].residue_type << "\t" << "CA" << "\t" << "C" << "\t" << chemical_shifts_table[i][transform_id(j)] << std::endl;
                                   atom_count+= 1;
                              }

                              break;
                         case 1: //!CB
                              if (chain[i].residue_type != GLY && loaded_tables[1] == true){
                                   std::cout << atom_count << "\t" << (i+1) << "\t" << chain[i].residue_type << "\t" << "CB" << "\t" << "C" << "\t" << chemical_shifts_table[i][transform_id(j)] << std::endl;
                                   atom_count+= 1;
                              }
                              break;
                         case 2: //!C
                              if (loaded_tables[2] == true){
                                   std::cout << atom_count << "\t" << (i+1) << "\t" << chain[i].residue_type << "\t" << "C" << "\t" << "C" << "\t" << chemical_shifts_table[i][transform_id(j)] << std::endl;
                              }
                              break;
                         case 3: //!N
                              if (loaded_tables[3] == true){
                              std::cout << atom_count << "\t" << (i+1) << "\t" << chain[i].residue_type << "\t" << "N" << "\t" << "N" << "\t" << chemical_shifts_table[i][transform_id(j)] << std::endl;
                              atom_count+= 1;
                              }
                              break;
                         case 4: //!HN
                              if (chain[i].residue_type != PRO && loaded_tables[4] == true){
                              std::cout << atom_count << "\t" << (i+1) << "\t" << chain[i].residue_type << "\t" << "H" << "\t" << "H" << "\t" << chemical_shifts_table[i][transform_id(j)] << std::endl;
                              atom_count+= 1;
                              }
                              break;
                         case 5: //!HA
                              if (loaded_tables[5] == true){
                              std::cout << atom_count << "\t" << (i+1) << "\t" << chain[i].residue_type << "\t" << "HA" << "\t" << "H" << "\t" << chemical_shifts_table[i][transform_id(j)] << std::endl;
                              atom_count+= 1;
                              }
                              break;
                    }
               }
          }
     }


}; // end class TermProCS15Cached

//! Observable specialization for TermProCS15Cached
template <>
class Observable<TermProCS15Cached>: public TermProCS15Cached, public ObservableBase {

public:

     //! Local settings class.
     const class Settings: public TermProCS15Cached::Settings, public ObservableBase::Settings {
     public:

          //! Constructor. Defines default values for settings object.
          Settings(){}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << static_cast<const TermProCS15Cached::Settings>(settings);
               o << static_cast<const ObservableBase::Settings>(settings);
               return o;
          }
     } settings; //!< Local settings object

     //! Constructor.
     //! \param energy_term TermCamshift energy term object
     //! \param settings Local Settings object
     //! \param reference_energy_function All observables have a pointer to a reference energy function which they can refer to.
     Observable(const TermProCS15Cached &energy_term,
                const ObservableBase::Settings &settings=ObservableBase::Settings(),
                Energy<ChainFB> *reference_energy_function=NULL)
          : TermProCS15Cached(energy_term),
            settings(dynamic_cast<const Settings&>(settings)) {
     }

     //! Copy Constructor.
     //! \param other Source object from which copy is made
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     Observable(const Observable &other, int thread_index, ChainFB *chain)
          : TermProCS15Cached(other, random_number_engine, thread_index, chain),
            settings(other.settings) {
     }


     //! Chemical shift RMSD between two sets of chemical shifts
     //! \param cs1_this A matrix containing chemical shifts
     //! \param cs2_this A matrix containing chemical shifts
     //! \return A vector containing RMSDs
     FPlist calc_rmsds(const FPtable &cs1_this,
                                    const FPtable &cs2_this) {

          FPlist chi_sq = empty_contribution;
          std::vector<int> n_bins = vector_utils::make_vector(0, 0, 0, 0, 0, 0);

          for (unsigned int i = 0; i <  std::min(cs1_this.size(), cs2_this.size()); i++) {
               for (unsigned int j = 0; j < 6; j++) {

                   if ((std::fabs(cs1_this[i][j]) > 0.0001)
                    && (std::fabs(cs2_this[i][j]) > 0.0001)) {

                         const FPtype diff  = cs1_this[i][j] - cs2_this[i][j];
                         chi_sq[j] += diff * diff;
                         n_bins[j] += 1;
                    }
               }
          }

          FPlist rmsds = empty_contribution;

          for (unsigned int j = 0; j < 6; j++) {
               rmsds[j] = std::sqrt(chi_sq[j] / n_bins[j]);
          }

          return rmsds;

     }


     FPtable get_full_prediction_error(const FPtable &cs1_this, const FPtable &cs2_this) {

          FPlist cs_vector;
          FPtable prediction_errors(6,cs_vector);

          for (unsigned int i = 0; i <  std::min(cs1_this.size(), cs2_this.size()); i++) {
               for (unsigned int j = 0; j < 6; j++) {

                   if ((std::fabs(cs1_this[i][j]) > 0.0001)
                    && (std::fabs(cs2_this[i][j]) > 0.0001)
                    && !(std::isnan(cs1_this[i][j]))
                    && !(std::isnan(cs2_this[i][j]))) {

                         const double diff = cs1_this[i][j] - cs2_this[i][j];
                         prediction_errors[j].push_back(diff);
                    }
                   else {
                       prediction_errors[j].push_back(0.0);
                   }
               }
          }

          return prediction_errors;
     }


     //! Make observation.
     virtual std::string observe(MoveInfo *move_info=NULL,
                                 PHAISTOS_LONG_LONG current_iteration=0,
                                 bool register_only=false) {

          // Energy to be returned
          double energy = this->evaluate();

          // Calculate new chemical shifts
          this->predicted_chemical_shifts = this->predict(*(this->chain));

          // Calculate RMSDs
          FPlist rmsds = calc_rmsds(this->predicted_chemical_shifts,
                                    this->experimental_chemical_shifts);

          // Output stream
          std::stringstream s;
          s << std::fixed << std::setprecision(5) << energy << ":" << rmsds;

          // Add full prediction error to output
          if (settings.output_full_prediction_vector) {
               FPtable prediction_error = get_full_prediction_error(this->predicted_chemical_shifts,
                                                                    this->experimental_chemical_shifts);
               s << ":" << prediction_error;
          }

          return s.str();

     }

}; // End class Observable

} // end namespace phaistos

#endif // TERM_PROCS15_CACHED

