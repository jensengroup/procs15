// Copyright (C) 2014-2015 by Anders Steen Christensen, Anders S Larsen, Lars Andersen Bratholm
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

#ifndef PROCS15_BACKEND
#define PROCS15_BACKEND

#include <boost/assign/list_of.hpp>

#include "cnpy_phaistos.h"
#include "array_wrapper.h"



namespace phaistos {

using namespace procs15;
using namespace definitions;

class ProCS15Backend {

public:


     //! Vector of pointers to array wrappers that hold numpy data
     std::vector< ArrayWrapper *> list_of_tables;

     //! Which tables (data files) are loaded.
     std::vector<bool> loaded_tables;

     std::vector<unsigned int> atypevector;

     HbondDataWrapper* hbond_array;
     HbondDataWrapper* hbond_array_alcohol;
     HbondDataWrapper* hbond_array_carboxy;
     HbondDataWrapper* halphabond_array_oxygen;
     HbondDataWrapper* halphabond_array_carboxy;
     HbondDataWrapper* halphabond_array_alcohol;



     // Define enum types for hydrogen bond acceptor names
     enum AcceptorEnum {AcceptorAmide = 0,
                       AcceptorAlcohol,
                       AcceptorCarboxylate,
                       NoneBondType,
                       ACCEPTOR_ENUM_SIZE};
     // Define enum types for hydrogen bond donor names
     enum DonorEnum {DonorAmide = 0,
                     DonorAlcohol,
                     DonorAmmonium,
                     DonorGuanidinium,
                     DonorImidazolium,
                     DonorIndole,
                     DonorThiol,
                     DonorAlphaHydrogen,
                     DONOR_ENUM_SIZE};


     enum ProCS15AtomTypes{ ProCS15_CA = 0,
                            ProCS15_CB ,
                            ProCS15_C  ,
                            ProCS15_N  ,
                            ProCS15_HN ,
                            ProCS15_HA };


     // Class that holds the first information on hydrogen bond acceptors
     class HBondAcceptor {
         public:
            Atom* acceptor_oxygen;
            Atom* acceptor_second_atom;
            Atom* acceptor_third_atom;
            unsigned int acceptor_type;
     };


     // Class that holds the first information on hydrogen bond acceptors
     class HBondDonor {
         public:
            Atom* donor_hydrogen;
            Atom* donor_second_atom;
            unsigned int donor_type;
     };




     // Class that holds necessary information about an aromatic ring
     class AromaticRing {
          public:
               Vector_3D ring_center;
               Vector_3D ring_normal_vector;
               unsigned int ring_type;
     };

     // Class that holds variables for a residue
     class ResidueVariable {

          public:
               unsigned int residue_type;
               std::vector<unsigned int> phi_psi_index;
               std::vector<unsigned int> chi_index;
               std::vector<HBondDonor> donors;
               std::vector<HBondAcceptor> acceptors;
               std::vector<AromaticRing> aromatic_rings;

     };



     FPtype calc_aromatic_interaction(Vector_3D h_position,
                                      AromaticRing& ring,
                                      FPtype dist_squared) {

           FPtype dist = sqrt(dist_squared);
           FPtype dist_third_power = dist * dist* dist;

           FPtype theta = calc_angle(h_position, ring.ring_center, ring.ring_normal_vector);
           FPtype cos_theta = cos(theta);
           FPtype three_cos_theta_squared = 3.0 * cos_theta * cos_theta;
           FPtype interaction = (1.0 - three_cos_theta_squared) 
                                * ring_current_intensities[ring.ring_type]
                                * b_factor / dist_third_power;
           return interaction;
     }


    inline unsigned int to_index_with_bounds(FPtype value, unsigned int lower, unsigned int upper) {

        unsigned int converted_value = (unsigned int)value + 0.5;

        if (value < lower) {

             converted_value = lower;

        } else if (value > upper) {

             converted_value = upper;
        }

        return converted_value;

    }


     boost::multi_array<float, 1> get_bb_bb_data_vector(const FPtype &r_oh, const FPtype &theta, const FPtype &rho) {

             const std::vector<unsigned int> indexes = read_hbond_index(r_oh, theta, rho);

                 return hbond_array->fourd_array[indexes[0]][indexes[1]][indexes[2]];
     }

     boost::multi_array<float, 1> get_bb_sc_carboxy_data_vector(const FPtype r_oh, const FPtype theta, const FPtype rho) {

             const std::vector<unsigned int> indexes = read_hbond_index_carboxy(r_oh, theta, rho);
                 return hbond_array_carboxy->fourd_array[indexes[0]][indexes[1]][indexes[2]];
     }

     boost::multi_array<float, 1> get_bb_sc_alcohol_data_vector(const FPtype r_oh, const FPtype theta, const FPtype rho) {

             const std::vector<unsigned int> indexes = read_hbond_index_alcohol(r_oh, theta, rho);
                 return hbond_array_alcohol->fourd_array[indexes[0]][indexes[1]][indexes[2]];
     }


     boost::multi_array<float, 1> get_bb_bb_alphah_oxygen_data_vector(const FPtype r_oh, const FPtype theta, const FPtype rho) {

             const std::vector<unsigned int> indexes = read_halphabond_index_oxygen(r_oh, theta, rho);
                 return halphabond_array_oxygen->fourd_array[indexes[0]][indexes[1]][indexes[2]];
     }

     boost::multi_array<float, 1> get_bb_bb_alphah_carboxy_data_vector(const FPtype r_oh, const FPtype theta, const FPtype rho) {

             const std::vector<unsigned int> indexes = read_halphabond_index_carboxy(r_oh, theta, rho);
                 return halphabond_array_carboxy->fourd_array[indexes[0]][indexes[1]][indexes[2]];
     }

     boost::multi_array<float, 1> get_bb_bb_alphah_alcohol_data_vector(const FPtype r_oh, const FPtype theta, const FPtype rho) {

             const std::vector<unsigned int> indexes = read_halphabond_index_alcohol(r_oh, theta, rho);
                 return halphabond_array_alcohol->fourd_array[indexes[0]][indexes[1]][indexes[2]];
     }




     std::vector<unsigned int> read_hbond_index(FPtype r, FPtype theta, FPtype rho) {

          // Convert theta and rho to 
          theta = theta / M_PI * 180.0;
          rho = rho / M_PI * 180.0;
          const FPtype r_min   = 1.5;
          const FPtype r_delta = 0.01;

          const FPtype theta_min   = 90.0;
          const FPtype theta_delta = 1.0;

          const FPtype rho_min   = -180.0;
          const FPtype rho_delta = 1.0;



          int r_index     = round((r - r_min) / r_delta);

          if (r_index < 0)   r_index = 0;
          if (r_index > 250) r_index = 250;


          int theta_index = round((theta -  theta_min) / theta_delta);

          if (theta_index < 0)   theta_index = 0;
          if (theta_index > 90) theta_index = 90;


          if (rho > 180) rho = 360 - rho;
          int rho_index   = round((rho  - rho_min)/ rho_delta);

          if (rho_index < 0)   rho_index = 0;
          if (rho_index > 361) rho_index = 360;




          return vector_utils::make_vector<unsigned int>(r_index, theta_index, rho_index);
     }

     std::vector<unsigned int> read_hbond_index_carboxy(FPtype r, FPtype theta, FPtype rho) {

          // Convert theta and rho to 
          theta = theta / M_PI * 180.0;
          rho = rho / M_PI * 180.0;
          const FPtype r_min   = 1.5;
          const FPtype r_delta = 0.01;

          const FPtype theta_min   = 90.0;
          const FPtype theta_delta = 1.0;

          const FPtype rho_min   = 0.0;
          const FPtype rho_delta = 1.0;



          int r_index     = round((r - r_min) / r_delta);

          if (r_index < 0)   r_index = 0;
          if (r_index > 250) r_index = 250;

	  
          int theta_index = round((theta -  theta_min) / theta_delta);

          if (theta_index < 0)   theta_index = 0;
          if (theta_index > 90) theta_index = 90;


          rho = rho + 180; // for numpy array index
          int rho_index   = round((rho  - rho_min)/ rho_delta);
          if (rho_index < 0)   rho_index = 0;
          if (rho_index > 361) rho_index = 360;




          return vector_utils::make_vector<unsigned int>(r_index, theta_index, rho_index);
     }


     std::vector<unsigned int> read_hbond_index_alcohol(FPtype r, FPtype theta, FPtype rho) {

          // Convert theta and rho to 
          theta = theta / M_PI * 180.0;
          rho = rho / M_PI * 180.0;
          const FPtype r_min   = 1.5;
          const FPtype r_delta = 0.01;

          const FPtype theta_min   = 90.0;
          const FPtype theta_delta = 1.0;

          const FPtype rho_min   = 0.0;
          const FPtype rho_delta = 1.0;

          int r_index  = round((r - r_min) / r_delta);

          if (r_index < 0)   r_index = 0;
          if (r_index > 250) r_index = 250;

	  
          int theta_index = round((theta -  theta_min) / theta_delta);

          if (theta_index < 0)   theta_index = 0;
          if (theta_index > 90) theta_index = 90;


          rho = rho + 180; // for numpy array index
          int rho_index   = round((rho  - rho_min)/ rho_delta);
          if (rho_index < 0)   rho_index = 0;
          if (rho_index > 361) rho_index = 360;



          return vector_utils::make_vector<unsigned int>(r_index, theta_index, rho_index);
     }

     std::vector<unsigned int> read_halphabond_index_oxygen(FPtype r, FPtype theta, FPtype rho) {

          // Convert theta and rho to 
          theta = theta / M_PI * 180.0;
          rho = rho / M_PI * 180.0;
          const FPtype r_min   = 1.8;
          const FPtype r_delta = 0.01;

          const FPtype theta_min   = 90.0;
          const FPtype theta_delta = 1.0;

          const FPtype rho_min   = 0.0;
          const FPtype rho_delta = 1.0;



          int r_index     = round((r - r_min) / r_delta);

          if (r_index < 0)   r_index = 0;
          if (r_index > 221) r_index = 221;

	  
          int theta_index = round((theta -  theta_min) / theta_delta);

          if (theta_index < 0)   theta_index = 0;
          if (theta_index > 90) theta_index = 90;


          rho = rho + 180; // for numpy array index
          int rho_index   = round((rho  - rho_min)/ rho_delta);
          if (rho_index < 0)   rho_index = 0;
          if (rho_index > 361) rho_index = 360;



          return vector_utils::make_vector<unsigned int>(r_index, theta_index, rho_index);
     }

     std::vector<unsigned int> read_halphabond_index_carboxy(FPtype r, FPtype theta, FPtype rho) {

          // Convert theta and rho to 
          theta = theta / M_PI * 180.0;
          rho = rho / M_PI * 180.0;
          const FPtype r_min   = 1.8;
          const FPtype r_delta = 0.01;

          const FPtype theta_min   = 90.0;
          const FPtype theta_delta = 1.0;

          const FPtype rho_min   = 0.0;
          const FPtype rho_delta = 1.0;



          int r_index     = round((r - r_min) / r_delta);

          if (r_index < 0)   r_index = 0;
          if (r_index > 221) r_index = 221;

	  
          int theta_index = round((theta -  theta_min) / theta_delta);

          if (theta_index < 0)   theta_index = 0;
          if (theta_index > 90) theta_index = 90;


          rho = rho + 180; // for numpy array index
          int rho_index   = round((rho  - rho_min)/ rho_delta);
          if (rho_index < 0)   rho_index = 0;
          if (rho_index > 361) rho_index = 360;




          return vector_utils::make_vector<unsigned int>(r_index, theta_index, rho_index);
     }

     std::vector<unsigned int> read_halphabond_index_alcohol(FPtype r, FPtype theta, FPtype rho) {

          // Convert theta and rho to 
          theta = theta / M_PI * 180.0;
          rho = rho / M_PI * 180.0;
          const FPtype r_min   = 1.8;
          const FPtype r_delta = 0.01;

          const FPtype theta_min   = 90.0;
          const FPtype theta_delta = 1.0;

          const FPtype rho_min   = 0.0;
          const FPtype rho_delta = 1.0;



          int r_index     = round((r - r_min) / r_delta);

          if (r_index < 0)   r_index = 0;
          if (r_index > 221) r_index = 221;

	  
          int theta_index = round((theta -  theta_min) / theta_delta);

          if (theta_index < 0)   theta_index = 0;
          if (theta_index > 90) theta_index = 90;


          rho = rho + 180; // for numpy array index
          int rho_index   = round((rho  - rho_min)/ rho_delta);
          if (rho_index < 0)   rho_index = 0;
          if (rho_index > 361) rho_index = 360;


          return vector_utils::make_vector<unsigned int>(r_index, theta_index, rho_index);
     }

     AtomEnum get_ha_atom_type(ResidueEnum residue_type){

          AtomEnum atomenum;

          switch (residue_type){

               case GLY:
                    atomenum = HA2;
                    return atomenum;
                    break;

               default:
                    atomenum = HA;
                    return atomenum;
                    break;
          }
     }



    std::vector<unsigned int> read_phi_psi_index(phaistos::ResidueFB& res1){
          FPtype phi = res1.get_phi();
          phi = ((phi*rad_to_degrees) > 0.0) ? floor( (phi*rad_to_degrees)  + 0.5) : 360 + ceil((phi*rad_to_degrees)  - 0.5);
          FPtype psi = res1.get_psi();
          psi = ((psi*rad_to_degrees) > 0.0) ? floor( (psi*rad_to_degrees)  + 0.5) : 360 + ceil((psi*rad_to_degrees)  - 0.5);

          return vector_utils::make_vector<unsigned int>(phi, psi);
     }




     std::vector<unsigned int> read_chi_index(phaistos::ResidueFB& res1) {
          //using namespace definitions;

          switch (res1.get_sidechain_dof_values(SIDECHAIN_ATOMS).size()) {
               case 1: {
                    FPtype chi1 = res1.get_sidechain_dof_values(SIDECHAIN_ATOMS)[0];
                    chi1  = ((chi1*rad_to_degrees) > 0.0) ? floor( (chi1*rad_to_degrees)  + 0.5) : 360 + ceil((chi1*rad_to_degrees)  - 0.5);

                    //return boost::assign::list_of<unsigned int>((unsigned int)(chi1));
                    return vector_utils::make_vector<unsigned int>(chi1);
                    break;
               }
               case 2: {
                    FPtype chi1 = res1.get_sidechain_dof_values(SIDECHAIN_ATOMS)[0];
                    chi1  = ((chi1*rad_to_degrees) > 0.0) ? floor( (chi1*rad_to_degrees)  + 0.5) : 360 + ceil((chi1*rad_to_degrees)  - 0.5);


                    FPtype chi2 = res1.get_sidechain_dof_values(SIDECHAIN_ATOMS)[1];
                    chi2  = ((chi2*rad_to_degrees) > 0.0) ? floor( (chi2*rad_to_degrees)  + 0.5) : 360 + ceil((chi2*rad_to_degrees)  - 0.5);


                    return vector_utils::make_vector<unsigned int>(chi1,chi2);
                    //return boost::assign::list_of<unsigned int>((unsigned int)(chi1 ))
                    //                                           ((unsigned int)(chi2 ));
                    break;
               }
               case 3: {
                  FPtype chi1 = res1.get_sidechain_dof_values(SIDECHAIN_ATOMS)[0];
                  chi1  = ((chi1*rad_to_degrees) > 0.0) ? floor( (chi1*rad_to_degrees)  + 0.5) : 360 + ceil((chi1*rad_to_degrees)  - 0.5);


                  FPtype chi2 = res1.get_sidechain_dof_values(SIDECHAIN_ATOMS)[1];
                  chi2  = ((chi2*rad_to_degrees) > 0.0) ? floor( (chi2*rad_to_degrees)  + 0.5) : 360 + ceil((chi2*rad_to_degrees)  - 0.5);

                  FPtype chi3 = res1.get_sidechain_dof_values(SIDECHAIN_ATOMS)[2];
                  chi3  = ((chi3*rad_to_degrees) > 0.0) ? floor( (chi3*rad_to_degrees)  + 0.5) : 360 + ceil((chi3*rad_to_degrees)  - 0.5);

                    return vector_utils::make_vector<unsigned int>(chi1,chi2,chi3);
                  //return boost::assign::list_of<unsigned int>((unsigned int)(chi1 ))
                  //                                             ((unsigned int)(chi2 ))
                  //                                             ((unsigned int)(chi3 ));
                    break;
               }
               case 4: {
                  FPtype chi1 = res1.get_sidechain_dof_values(SIDECHAIN_ATOMS)[0];
                  chi1  = ((chi1*rad_to_degrees) > 0.0) ? floor( (chi1*rad_to_degrees)  + 0.5) : 360 + ceil((chi1*rad_to_degrees)  - 0.5);


                  FPtype chi2 = res1.get_sidechain_dof_values(SIDECHAIN_ATOMS)[1];
                  chi2  = ((chi2*rad_to_degrees) > 0.0) ? floor( (chi2*rad_to_degrees)  + 0.5) : 360 + ceil((chi2*rad_to_degrees)  - 0.5);

                  FPtype chi3 = res1.get_sidechain_dof_values(SIDECHAIN_ATOMS)[2];
                  chi3  = ((chi3*rad_to_degrees) > 0.0) ? floor( (chi3*rad_to_degrees)  + 0.5) : 360 + ceil((chi3*rad_to_degrees)  - 0.5);

                  FPtype chi4 = res1.get_sidechain_dof_values(SIDECHAIN_ATOMS)[3];
                  chi4  = ((chi4*rad_to_degrees) > 0.0) ? floor( (chi4*rad_to_degrees)  + 0.5) : 360 + ceil((chi4*rad_to_degrees)  - 0.5);

                    return vector_utils::make_vector<unsigned int>(chi1,chi2,chi3,chi4);
                    //return boost::assign::list_of<unsigned int>((unsigned int)(chi1 ))
                    //                                           ((unsigned int)(chi2 ))
                    //                                           ((unsigned int)(chi3 ))
                    //                                           ((unsigned int)(chi4 ));
                    break;
               }
               case 5: {
                  FPtype chi1 = res1.get_sidechain_dof_values(SIDECHAIN_ATOMS)[0];
                  chi1  = ((chi1*rad_to_degrees) > 0.0) ? floor( (chi1*rad_to_degrees)  + 0.5) : 360 + ceil((chi1*rad_to_degrees)  - 0.5);


                  FPtype chi2 = res1.get_sidechain_dof_values(SIDECHAIN_ATOMS)[1];
                  chi2  = ((chi2*rad_to_degrees) > 0.0) ? floor( (chi2*rad_to_degrees)  + 0.5) : 360 + ceil((chi2*rad_to_degrees)  - 0.5);

                  FPtype chi3 = res1.get_sidechain_dof_values(SIDECHAIN_ATOMS)[2];
                  chi3  = ((chi3*rad_to_degrees) > 0.0) ? floor( (chi3*rad_to_degrees)  + 0.5) : 360 + ceil((chi3*rad_to_degrees)  - 0.5);

                  FPtype chi4 = res1.get_sidechain_dof_values(SIDECHAIN_ATOMS)[3];
                  chi4  = ((chi4*rad_to_degrees) > 0.0) ? floor( (chi4*rad_to_degrees)  + 0.5) : 360 + ceil((chi4*rad_to_degrees)  - 0.5);
                  
		  FPtype chi5 = res1.get_sidechain_dof_values(SIDECHAIN_ATOMS)[4];
                  chi5  = ((chi5*rad_to_degrees) > 0.0) ? floor( (chi5*rad_to_degrees)  + 0.5) : 360 + ceil((chi5*rad_to_degrees)  - 0.5);


                    return vector_utils::make_vector<unsigned int>(chi1,chi2,chi3,chi4,chi5);
                    //return boost::assign::list_of<unsigned int>((unsigned int)(chi1 ))
                    //                                           ((unsigned int)(chi2 ))
                    //                                           ((unsigned int)(chi3 ))
                    //                                           ((unsigned int)(chi4 ))
                    //                                           ((unsigned int)(chi5 ));
                    break;
               }



               default: {
                    std::vector<unsigned int> empty_vector;
                    return empty_vector;
                    break;
               }
          }
          std::vector<unsigned int> empty_vector;
          return empty_vector;

     }


     std::vector<unsigned int> get_chi(phaistos::ResidueFB& R){ //resize vectors with chi

          //using namespace definitions;

          std::vector<unsigned int> r;

          switch (R.residue_type){

               case ALA:
                    return r;
             	    break;

               case GLY:
                    return r;
                    break;

               case VAL: //! Special case for Val
                    {
                    const Vector_3D n_pos   = (R)[N]->position;
                    const Vector_3D ca_pos  = (R)[CA]->position; 
                    const Vector_3D cb_pos  = (R)[CB]->position;
                    const Vector_3D cg1_pos = (R)[CG1]->position;
                    const Vector_3D cg2_pos = (R)[CG2]->position;

                    const FPtype dihedral_cg1_cg2 = calc_dihedral(cg1_pos, ca_pos, cb_pos, cg2_pos);

                    if ( dihedral_cg1_cg2 > 0){

                         const FPtype chi1_double = calc_dihedral(n_pos,ca_pos,cb_pos,cg1_pos);
                         const unsigned int chi1_int  = ((chi1_double*rad_to_degrees) > 0.0) ? floor( (chi1_double*rad_to_degrees)  + 0.5) : 360 + ceil((chi1_double*rad_to_degrees)  - 0.5);


                         r.push_back(chi1_int);

                         return r;
                    } 
                    else {

                         const FPtype chi1_double = calc_dihedral(n_pos,ca_pos,cb_pos,cg2_pos);
                         const unsigned int chi1_int  = ((chi1_double*rad_to_degrees) > 0.0) ? floor( (chi1_double*rad_to_degrees)  + 0.5) : 360 + ceil((chi1_double*rad_to_degrees)  - 0.5);


                         r.push_back(chi1_int);

                         return r;
                    }

                    break;
                    }

               case THR:
             	    r = read_chi_index(R);
             	    r.resize(2);
             	    return r;
             	    break;

               case SER:
             	    r = read_chi_index(R);
             	    r.resize(1);
             	    return r;
             	    break;

               case CYS:
                    r = read_chi_index(R);
             	    r.resize(1);
             	    return r;
             	    break;

               case ILE:
             	    r = read_chi_index(R);
             	    r.resize(2);
             	    return r;
             	    break;

               case LEU:
                    r = read_chi_index(R);
             	    r.resize(2);
             	    return r;
             	    break;

               case ASN:
                    r = read_chi_index(R);
                    r.resize(2);
                    return r;
                    break;

               case TYR:
                    r = read_chi_index(R);
                    r.resize(2);
                    return r;
                    break;

               case PRO:
                    r = read_chi_index(R);
                    r.resize(0);
                    return r;
                    break;

               case LYS:
                    r = read_chi_index(R);
                    r.resize(4);
                    return r;
                    break;

               case GLN:
                    r = read_chi_index(R);
                    r.resize(3);
                    return r;
                    break;

               case MET:
                    r = read_chi_index(R);
                    r.resize(3);
                    return r;
                    break;

               default:
                    r = read_chi_index(R);
                    return r;
                    break;
          }
     }


     int transform_id(int l){ // Transform atom type id from ca,cb,c,n,h,ha to ha,ca,h,n,c,cb
          switch(l){
               case 0:
	            return 1;
		    break;
	       case 1:
	            return 5;
		    break;
	       case 2:
	            return 4;
	            break;
	       case 3:
	            return 3;
		    break;
	       case 4:
	            return 2;
		    break;
	       case 5:
	            return 0;
		    break;
	       default:
	            std::cout << "ERROR transform id";
		    return 999;
		    break;
          }
     }



    std::vector<HBondDonor> get_donors(phaistos::ResidueFB& res1) {

          //using namespace definitions;

          std::vector<HBondDonor> hbond_donors;

          if (res1.terminal_status == NTERM) {

               HBondDonor donor1;
               HBondDonor donor2;

               donor1.donor_hydrogen    = (res1)[H1];
               donor1.donor_second_atom = (res1)[N];
               donor1.donor_type        = DonorAmmonium;

               donor2.donor_hydrogen    = (res1)[H2];
               donor2.donor_second_atom = (res1)[N];
               donor2.donor_type        = DonorAmmonium;

               hbond_donors.push_back(donor1);
               hbond_donors.push_back(donor2);

               if (res1.residue_type!=PRO) {

                    HBondDonor donor3;

                    donor3.donor_hydrogen    = (res1)[H3];
                    donor3.donor_second_atom = (res1)[N];
                    donor3.donor_type        = DonorAmmonium;

                    hbond_donors.push_back(donor3);
               }

          } else if (res1.has_atom(H)) {

               HBondDonor donor1;

               donor1.donor_hydrogen    = (res1)[H];
               donor1.donor_second_atom = (res1)[N];
               donor1.donor_type        = DonorAmide;

               hbond_donors.push_back(donor1);
          }
          //! get alpha hydrogens
	  if (res1.residue_type != GLY){
 
               HBondDonor donor1;
               donor1.donor_hydrogen    = (res1)[HA];
               donor1.donor_second_atom = (res1)[CA];
               donor1.donor_type = DonorAlphaHydrogen;
               hbond_donors.push_back(donor1);

          }
	  else if (res1.residue_type == GLY){

               HBondDonor donor1;
               donor1.donor_hydrogen    = (res1)[HA2];
               donor1.donor_second_atom = (res1)[CA];
               donor1.donor_type = DonorAlphaHydrogen;

               HBondDonor donor2;
               donor2.donor_hydrogen    = (res1)[HA3];
               donor2.donor_second_atom = (res1)[CA];
               donor2.donor_type = DonorAlphaHydrogen;
  
               hbond_donors.push_back(donor1);
               hbond_donors.push_back(donor2);
          }

          switch (res1.residue_type) {

               case ASN: {
                   HBondDonor donor1, donor2;
                   donor1.donor_hydrogen    = (res1)[HD21];
                   donor1.donor_second_atom = (res1)[ND2];
                   donor1.donor_type        = DonorAmide;

                   donor2.donor_hydrogen    = (res1)[HD22];
                   donor2.donor_second_atom = (res1)[ND2];
                   donor2.donor_type        = DonorAmide;

                   hbond_donors.push_back(donor1);
                   hbond_donors.push_back(donor2);

                   break;
               }

               case GLN: {
                   HBondDonor donor1, donor2;
                   donor1.donor_hydrogen    = (res1)[HE21];
                   donor1.donor_second_atom = (res1)[NE2];
                   donor1.donor_type        = DonorAmide;

                   donor2.donor_hydrogen    = (res1)[HE22];
                   donor2.donor_second_atom = (res1)[NE2];
                   donor2.donor_type        = DonorAmide;

                   hbond_donors.push_back(donor1);
                   hbond_donors.push_back(donor2);

                   break;
               }

               case HIS: {
                   if (res1.has_atom(HD1)) {
                        HBondDonor donor1;
                        donor1.donor_hydrogen    = (res1)[HD1];
                        donor1.donor_second_atom = (res1)[ND1];
                        donor1.donor_type        = DonorIndole;
                        hbond_donors.push_back(donor1);
                   }

                   if (res1.has_atom(HE2)) {
                        HBondDonor donor2;
                        donor2.donor_hydrogen    = (res1)[HE2];
                        donor2.donor_second_atom = (res1)[NE2];
                        donor2.donor_type        = DonorIndole;
                        hbond_donors.push_back(donor2);
                   }

                   break;
               }

               case TRP: {
                   HBondDonor donor1;
                   donor1.donor_hydrogen    = (res1)[HE1];
                   donor1.donor_second_atom = (res1)[NE1];
                   donor1.donor_type        = DonorIndole;
                   hbond_donors.push_back(donor1);
                   break;
               }

               case LYS: {
                   HBondDonor donor1, donor2, donor3;
                   donor1.donor_hydrogen    = (res1)[HZ1];
                   donor1.donor_second_atom = (res1)[NZ];
                   donor1.donor_type        = DonorAmmonium;

                   donor2.donor_hydrogen    = (res1)[HZ2];
                   donor2.donor_second_atom = (res1)[NZ];
                   donor2.donor_type        = DonorAmmonium;

                   donor3.donor_hydrogen    = (res1)[HZ3];
                   donor3.donor_second_atom = (res1)[NZ];
                   donor3.donor_type        = DonorAmmonium;

                   hbond_donors.push_back(donor1);
                   hbond_donors.push_back(donor2);
                   hbond_donors.push_back(donor3);

                   break;
               }

               case ARG: {
                   HBondDonor donor1, donor2, donor3, donor4, donor5;
                   donor1.donor_hydrogen    = (res1)[HH11];
                   donor1.donor_second_atom = (res1)[NH1];
                   donor1.donor_type        = DonorGuanidinium;

                   donor2.donor_hydrogen    = (res1)[HH12];
                   donor2.donor_second_atom = (res1)[NH1];
                   donor2.donor_type        = DonorGuanidinium;

                   donor3.donor_hydrogen    = (res1)[HH21];
                   donor3.donor_second_atom = (res1)[NH2];
                   donor3.donor_type        = DonorGuanidinium;

                   donor4.donor_hydrogen    = (res1)[HH22];
                   donor4.donor_second_atom = (res1)[NH2];
                   donor4.donor_type        = DonorGuanidinium;

                   donor5.donor_hydrogen    = (res1)[HE];
                   donor5.donor_second_atom = (res1)[NE];
                   donor5.donor_type        = DonorGuanidinium;

                   hbond_donors.push_back(donor1);
                   hbond_donors.push_back(donor2);
                   hbond_donors.push_back(donor3);
                   hbond_donors.push_back(donor4);
                   hbond_donors.push_back(donor5);

                   break;
               }

               default: {

                   break;
               }
          }

          return hbond_donors;

     }



     std::vector<HBondAcceptor> get_acceptors(phaistos::ResidueFB& res1, phaistos::ChainFB& chain) {

          //using namespace definitions;

          std::vector<HBondAcceptor> hbond_acceptors;


          if (res1.terminal_status == CTERM) {
               HBondAcceptor acceptor1;
               HBondAcceptor acceptor2;

               acceptor1.acceptor_oxygen      = (res1)[O];
               acceptor1.acceptor_second_atom = (res1)[C];
               acceptor1.acceptor_third_atom  = (res1)[OXT];
               acceptor1.acceptor_type        = AcceptorCarboxylate;

               acceptor2.acceptor_oxygen      = (res1)[OXT];
               acceptor2.acceptor_second_atom = (res1)[C];
               acceptor2.acceptor_third_atom  = (res1)[O];
               acceptor2.acceptor_type        = AcceptorCarboxylate;

               hbond_acceptors.push_back(acceptor1);
               hbond_acceptors.push_back(acceptor2);

          } else {
               HBondAcceptor acceptor1;

               acceptor1.acceptor_oxygen      = (res1)[O];
               acceptor1.acceptor_second_atom = (res1)[C];
               acceptor1.acceptor_third_atom  = chain[res1.index+1][N];
               acceptor1.acceptor_type        = AcceptorAmide;

               hbond_acceptors.push_back(acceptor1);

          }


          switch (res1.residue_type) {

               HBondAcceptor acceptor1;
               HBondAcceptor acceptor2;

               case SER:{
                   acceptor1.acceptor_oxygen      = (res1)[OG];
                   acceptor1.acceptor_second_atom = (res1)[CB];
                   acceptor1.acceptor_third_atom  = (res1)[HG];
                   acceptor1.acceptor_type        = AcceptorAlcohol;

                   hbond_acceptors.push_back(acceptor1);
                   break;}

               case THR:{
                   acceptor1.acceptor_oxygen      = (res1)[OG1];
                   acceptor1.acceptor_second_atom = (res1)[CB];
                   acceptor1.acceptor_third_atom  = (res1)[HG1];
                   acceptor1.acceptor_type        = AcceptorAlcohol;

                   hbond_acceptors.push_back(acceptor1);
                   break;}

               case TYR:{
                   acceptor1.acceptor_oxygen      = (res1)[OH];
                   acceptor1.acceptor_second_atom = (res1)[CZ];
                   acceptor1.acceptor_third_atom  = (res1)[HH];
                   acceptor1.acceptor_type        = AcceptorAlcohol;

                   hbond_acceptors.push_back(acceptor1);
                   break;}

               case ASN:{
                   acceptor1.acceptor_oxygen      = (res1)[OD1];
                   acceptor1.acceptor_second_atom = (res1)[CG];
                   acceptor1.acceptor_third_atom  = (res1)[ND2];
                   acceptor1.acceptor_type        = AcceptorAmide;

                   hbond_acceptors.push_back(acceptor1);
                   break;}

               case GLN:{
                   acceptor1.acceptor_oxygen      = (res1)[OE1];
                   acceptor1.acceptor_second_atom = (res1)[CD];
                   acceptor1.acceptor_third_atom  = (res1)[NE2];
                   acceptor1.acceptor_type        = AcceptorAmide;

                   hbond_acceptors.push_back(acceptor1);
                   break;}

               case ASP:{
                   acceptor1.acceptor_oxygen      = (res1)[OD1];
                   acceptor1.acceptor_second_atom = (res1)[CG];
                   acceptor1.acceptor_third_atom  = (res1)[CB];
                   acceptor1.acceptor_type        = AcceptorCarboxylate;

                   acceptor2.acceptor_oxygen      = (res1)[OD2];
                   acceptor2.acceptor_second_atom = (res1)[CG];
                   acceptor2.acceptor_third_atom  = (res1)[CB];
                   acceptor2.acceptor_type        = AcceptorCarboxylate;

                   hbond_acceptors.push_back(acceptor1);
                   hbond_acceptors.push_back(acceptor2);
                   break;}

               case GLU:{
                   acceptor1.acceptor_oxygen      = (res1)[OE1];
                   acceptor1.acceptor_second_atom = (res1)[CD];
                   acceptor1.acceptor_third_atom  = (res1)[CG];
                   acceptor1.acceptor_type        = AcceptorCarboxylate;

                   acceptor2.acceptor_oxygen      = (res1)[OE2];
                   acceptor2.acceptor_second_atom = (res1)[CD];
                   acceptor2.acceptor_third_atom  = (res1)[CG];
                   acceptor2.acceptor_type        = AcceptorCarboxylate;

                   hbond_acceptors.push_back(acceptor1);
                   hbond_acceptors.push_back(acceptor2);
                   break;}

               default:
                   break;
          }

          return hbond_acceptors;

     }


     std::vector<AromaticRing> get_aromatic_rings(phaistos::ResidueFB& res1) {

          using namespace definitions;

          std::vector<AromaticRing> aromatic_rings;

          switch (res1.residue_type) {
               case PHE: {
                    AromaticRing current_ring;

                    Vector_3D CE1sc, CE2sc, CGsc, CZsc;
                    CE1sc = (res1)[CE1]->position;
                    CE2sc = (res1)[CE2]->position;
                    CGsc  = (res1)[CG]->position;
                    CZsc  = (res1)[CZ]->position;
                    Vector_3D Phe_Centroid = (CGsc + CZsc)*0.5;
                    Vector_3D Phe_Normal = cross_product((CE1sc - CGsc), (CE2sc - CGsc)).normalize();
                    Phe_Normal = Phe_Normal + Phe_Centroid;

                    current_ring.ring_center        = Phe_Centroid;
                    current_ring.ring_normal_vector = Phe_Normal;
                    current_ring.ring_type          = RingBenzene;

                    aromatic_rings.push_back(current_ring);
                    break;
               }

               case TYR: {
                    AromaticRing current_ring;

                    Vector_3D CE1sc, CE2sc, CGsc, CZsc;
                    CE1sc = (res1)[CE1]->position;
                    CE2sc = (res1)[CE2]->position;
                    CGsc  = (res1)[CG]->position;
                    CZsc  = (res1)[CZ]->position;
                    Vector_3D Tyr_Centroid = (CGsc + CZsc)*0.5;
                    Vector_3D Tyr_Normal = cross_product((CE1sc - CGsc), (CE2sc - CGsc)).normalize();
                    Tyr_Normal = Tyr_Normal + Tyr_Centroid;

                    current_ring.ring_center        = Tyr_Centroid;
                    current_ring.ring_normal_vector = Tyr_Normal;
                    current_ring.ring_type          = RingPhenol;

                    aromatic_rings.push_back(current_ring);
                    break;
               }
               case HIS: {
                    AromaticRing current_ring;

                    Vector_3D CE1sc, CGsc, CD2sc;
                    CE1sc = (res1)[CE1]->position;
                    CD2sc = (res1)[CD2]->position;
                    CGsc  = (res1)[CG]->position;
                    Vector_3D His_Middle = (CGsc + CD2sc)*0.5;
                    Vector_3D His_Centroid = His_Middle*0.5527864045 + CE1sc*0.447213955;
                    Vector_3D His_Normal = cross_product((CGsc - CE1sc),(CD2sc - CE1sc)).normalize();
                    His_Normal = His_Normal + His_Centroid;

                    current_ring.ring_center        = His_Centroid;
                    current_ring.ring_normal_vector = His_Normal;
                    current_ring.ring_type          = RingImidazolium;


                    aromatic_rings.push_back(current_ring);
                    break;
               }
               case TRP: {
                    AromaticRing current_ring1, current_ring2;

                    Vector_3D NE1sc, CE2sc, CGsc, CD2sc;
                    NE1sc = (res1)[NE1]->position;
                    CE2sc = (res1)[CE2]->position;
                    CD2sc = (res1)[CD2]->position;
                    CGsc  = (res1)[CG]->position;

                    Vector_3D Trp5_Middle = (CGsc + CD2sc)*0.5;
                    Vector_3D Trp5_Centroid = Trp5_Middle*0.5527864045 + NE1sc*0.447213955;
                    Vector_3D Trp5_Normal = cross_product((CE2sc - CGsc),(NE1sc - CGsc)).normalize();
                    Trp5_Normal = Trp5_Normal + Trp5_Centroid;

                    current_ring1.ring_center        = Trp5_Centroid;
                    current_ring1.ring_normal_vector = Trp5_Normal;
                    current_ring1.ring_type          = RingIndoleSmall;

                    Vector_3D CE3sc, CZ2sc, CH2sc;
                    CE3sc = (res1)[CE3]->position;
                    CZ2sc = (res1)[CZ2]->position;
                    CH2sc = (res1)[CH2]->position;
                    Vector_3D Trp6_Centroid = (CE3sc + CZ2sc)*0.5;
                    Vector_3D Trp6_Normal = cross_product((CH2sc - CE3sc), (CE2sc - CE3sc)).normalize();
                    Trp6_Normal = Trp6_Normal + Trp6_Centroid;

                    current_ring2.ring_center        = Trp6_Centroid;
                    current_ring2.ring_normal_vector = Trp6_Normal;
                    current_ring2.ring_type          = RingIndoleLarge;

                    aromatic_rings.push_back(current_ring1);
                    aromatic_rings.push_back(current_ring2);
                    break;
               }

               default:
                   break;
          }

          return aromatic_rings;

     }


     FPlist add_vectors(FPlist& list1, FPlist& list2) {

          return vector_utils::make_vector((list1[0] + list2[0]),
                                           (list1[1] + list2[1]),
                                           (list1[2] + list2[2]),
                                           (list1[3] + list2[3]),
                                           (list1[4] + list2[4]),
                                           (list1[5] + list2[5]));
     }

     FPlist subtract_vectors(FPlist& list1, FPlist& list2) {

          return vector_utils::make_vector((list1[0] - list2[0]),
                                           (list1[1] - list2[1]),
                                           (list1[2] - list2[2]),
                                           (list1[3] - list2[3]),
                                           (list1[4] - list2[4]),
                                           (list1[5] - list2[5]));
     }



     FPtable add_tables(FPtable& table1, FPtable& table2) {

          assert(table1.size() == table2.size());

          FPtable return_table;

          for (unsigned int i = 0; i < table1.size(); i++) {
               return_table.push_back(add_vectors(table1[i],
                                                  table2[i]));
          }

          return return_table;


     }



     void check_path( bool &load_ca, bool &load_cb, bool &load_co, bool &load_n, bool &load_hn, bool &load_ha, std::string &numpypath ){

          if (!file_exists(numpypath+"delta_hbond_251_91_361_12.npy")){
               std::cerr << "ERROR (ProCS15) File " << numpypath << "delta_hbond_251_91_361_12.npy" << " doesn't exist check data-folder keyword" << std::endl;
               assert(false);
          }

          if (!file_exists(numpypath+"delta_hbond_Carboxylate_251_91_361_6.npy")){          
               std::cerr << "ERROR (ProCS15) File " << numpypath << "delta_hbond_Carboxylate_251_91_361_6.npy" << " doesn't exist check data-folder keyword" << std::endl;
               assert(false);
          }

          if (!file_exists(numpypath+"delta_hbond_Alcohol_251_91_361_6.npy")){          
               std::cerr << "ERROR (ProCS15) File " << numpypath << "delta_hbond_Carboxylate_251_91_361_6.npy" << " doesn't exist check data-folder keyword" << std::endl;
               assert(false);
          }

          if (!file_exists(numpypath+"delta_halphabond_oxygen_221_91_361_12.npy")){          
               std::cerr << "ERROR (ProCS15) File " << numpypath << "delta_halphabond_oxygen_221_91_361_12.npy" << " doesn't exist check data-folder keyword" << std::endl;
               assert(false);
          }

          if (!file_exists(numpypath+"delta_halphabond_carboxy_221_91_361_6.npy")){          
               std::cerr << "ERROR (ProCS15) File " << numpypath << "delta_halphabond_carboxy_221_91_361_6.npy" << " doesn't exist check data-folder keyword" << std::endl;
               assert(false);
          }

          if (!file_exists(numpypath+"delta_halphabond_alcohol_221_91_361_6.npy")){          
               std::cerr << "ERROR (ProCS15) File " << numpypath << "delta_halphabond_alcohol_221_91_361_6.npy" << " doesn't exist check data-folder keyword" << std::endl;
               assert(false);
          }

          std::string procs15_residue_names[20] = {"ALA","CYS","ASP",
                                                   "GLU","PHE","GLY",
                                                   "HIS","ILE","LYS",
                                                   "LEU","MET","ASN",
                                                   "PRO","GLN","ARG",
                                                   "SER","THR","VAL",
                                                   "TRP","TYR"};
   
          if ( load_ca == true ) {          
               for ( unsigned int i = 0; i < 20; i++){
                    if (!file_exists(numpypath+procs15_residue_names[i]+"_ca.npy")){
                         std::cerr << "ERROR (ProCS15) File " << numpypath << procs15_residue_names[i]+"_ca.npy" << " doesn't exist check data-folder keyword" << std::endl;
                    } 
               } 
          }

          if ( load_cb == true ) {           
               for ( int i = 0; i < 20; i++){
                    if (!file_exists(numpypath+procs15_residue_names[i]+"_cb.npy")){
                         std::cerr << "ERROR (ProCS15) File " << numpypath << procs15_residue_names[i]+"_cb.npy" << " doesn't exist check data-folder keyword" << std::endl;
                    } 

               } 
          }

          if ( load_co == true ) { 
               for ( int i = 0; i < 20; i++){
                    if (!file_exists(numpypath+procs15_residue_names[i]+"_co.npy")){
                         std::cerr << "ERROR (ProCS15) File " << numpypath << procs15_residue_names[i]+"_co.npy" << " doesn't exist check data-folder keyword" << std::endl;
                    } 
               } 
          }

          if ( load_n == true ) {           
               for ( int i = 0; i < 20; i++){
                    if (!file_exists(numpypath+procs15_residue_names[i]+"_nh.npy")){
                         std::cerr << "ERROR (ProCS15) File " << numpypath << procs15_residue_names[i]+"_nh.npy" << " doesn't exist check data-folder keyword" << std::endl;
                    } 
               } 
          }

          if ( load_hn == true ) {           
               for ( int i = 0; i < 20; i++){
                    if (!file_exists(numpypath+procs15_residue_names[i]+"_hn.npy")){
                         std::cerr << "ERROR (ProCS15) File " << numpypath << procs15_residue_names[i]+"_hn.npy" << " doesn't exist check data-folder keyword" << std::endl;
                    } 
               } 
          }

          if ( load_ha == true ) {           
               for ( int i = 0; i < 20; i++){
                    if (!file_exists(numpypath+procs15_residue_names[i]+"_ha.npy")){
                         std::cerr << "ERROR (ProCS15) File " << numpypath << procs15_residue_names[i]+"_ha.npy" << " doesn't exist check data-folder keyword" << std::endl;
                    } 
               } 
          }
     } 


     void load_data_linearinterpolation(bool &load_ca, bool &load_cb, bool &load_co, bool &load_n, bool &load_hn, bool &load_ha, std::string &path){


          if (load_ca){
               std::cout << "CA " << std::flush;
               list_of_tables[0] =  new ArrayWrapperStandard(cnpy::npy_load(path+"ALA_ca.npy").get_3d_array()); loaded_tables[0] = true;
               list_of_tables[6] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"CYS_ca.npy").get_4d_array(),gridsize4d); loaded_tables[6] = true;
               list_of_tables[12] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ASP_ca.npy").get_5d_array(),gridsize5d ); loaded_tables[12] = true;
               list_of_tables[18] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"GLU_ca.npy").get_6d_array(),gridsize6d ); loaded_tables[18] = true;
               list_of_tables[24] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"PHE_ca.npy").get_5d_array(),gridsize5d); loaded_tables[24] = true;
               list_of_tables[30] = new ArrayWrapperStandard(cnpy::npy_load(path+"GLY_ca.npy").get_3d_array()); loaded_tables[30] = true;
               list_of_tables[36] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"HIS_ca.npy").get_5d_array(),gridsize5d); loaded_tables[36] = true;
               list_of_tables[42] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ILE_ca.npy").get_5d_array(),gridsize5d); loaded_tables[42] = true;
               list_of_tables[48] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"LYS_ca.npy").get_7d_array(),gridsize7d); loaded_tables[48] = true;
               list_of_tables[54] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"LEU_ca.npy").get_5d_array(),gridsize5d); loaded_tables[54] = true;
               list_of_tables[60] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"MET_ca.npy").get_6d_array(),gridsize6d ); loaded_tables[60] = true;
               list_of_tables[66] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ASN_ca.npy").get_5d_array(),gridsize5d); loaded_tables[66] = true;
               list_of_tables[72] =  new ArrayWrapperStandard(cnpy::npy_load(path+"PRO_ca.npy").get_3d_array()); loaded_tables[72] = true;
               list_of_tables[78] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"GLN_ca.npy").get_6d_array(),gridsize6d ); loaded_tables[78] = true;
               list_of_tables[84] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ARG_ca.npy").get_7d_array(),gridsize7d); loaded_tables[84] = true;
               list_of_tables[90] = new ArrayWrapperInterpolate(cnpy::npy_load(path+"SER_ca.npy").get_4d_array(),gridsize4d); loaded_tables[90] = true;
               list_of_tables[96] = new ArrayWrapperInterpolate(cnpy::npy_load(path+"THR_ca.npy").get_5d_array(),gridsize5d); loaded_tables[96] = true;
               list_of_tables[102] = new ArrayWrapperInterpolate(cnpy::npy_load(path+"VAL_ca.npy").get_4d_array(),gridsize4d); loaded_tables[102] = true;
               list_of_tables[108] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"TRP_ca.npy").get_5d_array(),gridsize5d); loaded_tables[108] = true;
               list_of_tables[114] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"TYR_ca.npy").get_5d_array(),gridsize5d); loaded_tables[114] = true;
          }
          if (load_cb){
               std::cout << "CB " << std::flush;               
               list_of_tables[1] =  new ArrayWrapperStandard(cnpy::npy_load(path+"ALA_cb.npy").get_3d_array());  loaded_tables[1] = true;
               list_of_tables[7] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"CYS_cb.npy").get_4d_array(),gridsize4d); loaded_tables[7] = true;
               list_of_tables[13] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ASP_cb.npy").get_5d_array(),gridsize5d); loaded_tables[13] = true;
               list_of_tables[19] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"GLU_cb.npy").get_6d_array(),gridsize6d); loaded_tables[19] = true;
               list_of_tables[25] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"PHE_cb.npy").get_5d_array(),gridsize5d); loaded_tables[25] = true;
               list_of_tables[31] = new ArrayWrapperStandard(cnpy::npy_load(path+"GLY_cb.npy").get_3d_array()); loaded_tables[31] = true;
               list_of_tables[37] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"HIS_cb.npy").get_5d_array(),gridsize5d); loaded_tables[37] = true;
               list_of_tables[43] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ILE_cb.npy").get_5d_array(),gridsize5d); loaded_tables[43] = true;
               list_of_tables[49] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"LYS_cb.npy").get_7d_array(),gridsize7d); loaded_tables[49] = true;
               list_of_tables[55] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"LEU_cb.npy").get_5d_array(),gridsize5d); loaded_tables[55] = true;
               list_of_tables[61] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"MET_cb.npy").get_6d_array(),gridsize6d); loaded_tables[61] = true;
               list_of_tables[67] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ASN_cb.npy").get_5d_array(),gridsize5d); loaded_tables[67] = true;
               list_of_tables[73] =  new ArrayWrapperStandard(cnpy::npy_load(path+"PRO_cb.npy").get_3d_array()); loaded_tables[73] = true;
               list_of_tables[79] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"GLN_cb.npy").get_6d_array(),gridsize6d); loaded_tables[79] = true;
               list_of_tables[85] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ARG_cb.npy").get_7d_array(),gridsize7d); loaded_tables[85] = true;
               list_of_tables[91] = new ArrayWrapperInterpolate(cnpy::npy_load(path+"SER_cb.npy").get_4d_array(),gridsize4d); loaded_tables[91] = true;
               list_of_tables[97] = new ArrayWrapperInterpolate(cnpy::npy_load(path+"THR_cb.npy").get_5d_array(),gridsize5d); loaded_tables[97] = true;
               list_of_tables[103] = new ArrayWrapperInterpolate(cnpy::npy_load(path+"VAL_cb.npy").get_4d_array(),gridsize4d); loaded_tables[103] = true;
               list_of_tables[109] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"TRP_cb.npy").get_5d_array(),gridsize5d); loaded_tables[109] = true;
               list_of_tables[115] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"TYR_cb.npy").get_5d_array(),gridsize5d); loaded_tables[115] = true;
          }

          if (load_co){
               std::cout << "CO " << std::flush;
               list_of_tables[2] =  new ArrayWrapperStandard(cnpy::npy_load(path+"ALA_co.npy").get_3d_array());  loaded_tables[2] = true;
               list_of_tables[8] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"CYS_co.npy").get_4d_array(),gridsize4d); loaded_tables[8] = true;
               list_of_tables[14] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ASP_co.npy").get_5d_array(),gridsize5d); loaded_tables[14] = true;
               list_of_tables[20] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"GLU_co.npy").get_6d_array(),gridsize6d); loaded_tables[20] = true;
               list_of_tables[26] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"PHE_co.npy").get_5d_array(),gridsize5d); loaded_tables[26] = true;
               list_of_tables[32] = new ArrayWrapperStandard(cnpy::npy_load(path+"GLY_co.npy").get_3d_array()); loaded_tables[32] = true;
               list_of_tables[38] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"HIS_co.npy").get_5d_array(),gridsize5d); loaded_tables[38] = true;
               list_of_tables[44] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ILE_co.npy").get_5d_array(),gridsize5d); loaded_tables[44] = true;
               list_of_tables[50] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"LYS_co.npy").get_7d_array(),gridsize7d); loaded_tables[50] = true;
               list_of_tables[56] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"LEU_co.npy").get_5d_array(),gridsize5d); loaded_tables[56] = true;
               list_of_tables[62] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"MET_co.npy").get_6d_array(),gridsize6d); loaded_tables[62] = true;
               list_of_tables[68] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ASN_co.npy").get_5d_array(),gridsize5d); loaded_tables[68] = true;
               list_of_tables[74] =  new ArrayWrapperStandard(cnpy::npy_load(path+"PRO_co.npy").get_3d_array()); loaded_tables[74] = true;
               list_of_tables[80] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"GLN_co.npy").get_6d_array(),gridsize6d); loaded_tables[80] = true;
               list_of_tables[86] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ARG_co.npy").get_7d_array(),gridsize7d); loaded_tables[86] = true;
               list_of_tables[92] = new ArrayWrapperInterpolate(cnpy::npy_load(path+"SER_co.npy").get_4d_array(),gridsize4d); loaded_tables[92] = true;
               list_of_tables[98] = new ArrayWrapperInterpolate(cnpy::npy_load(path+"THR_co.npy").get_5d_array(),gridsize5d); loaded_tables[98] = true;
               list_of_tables[104] = new ArrayWrapperInterpolate(cnpy::npy_load(path+"VAL_co.npy").get_4d_array(),gridsize4d); loaded_tables[104] = true;
               list_of_tables[110] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"TRP_co.npy").get_5d_array(),gridsize5d); loaded_tables[110] = true;
               list_of_tables[116] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"TYR_co.npy").get_5d_array(),gridsize5d); loaded_tables[116] = true;
          }
          if (load_n){
               std::cout << "N " << std::flush;
               list_of_tables[3] =  new ArrayWrapperStandard(cnpy::npy_load(path+"ALA_nh.npy").get_3d_array());  loaded_tables[3] = true;
               list_of_tables[9] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"CYS_nh.npy").get_4d_array(),gridsize4d); loaded_tables[9] = true;
               list_of_tables[15] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ASP_nh.npy").get_5d_array(),gridsize5d); loaded_tables[15] = true;
               list_of_tables[21] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"GLU_nh.npy").get_6d_array(),gridsize6d); loaded_tables[21] = true;
               list_of_tables[27] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"PHE_nh.npy").get_5d_array(),gridsize5d); loaded_tables[27] = true;
               list_of_tables[33] = new ArrayWrapperStandard(cnpy::npy_load(path+"GLY_nh.npy").get_3d_array()); loaded_tables[33] = true;
               list_of_tables[39] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"HIS_nh.npy").get_5d_array(),gridsize5d); loaded_tables[39] = true;
               list_of_tables[45] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ILE_nh.npy").get_5d_array(),gridsize5d); loaded_tables[45] = true;
               list_of_tables[51] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"LYS_nh.npy").get_7d_array(),gridsize7d); loaded_tables[51] = true;
               list_of_tables[57] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"LEU_nh.npy").get_5d_array(),gridsize5d); loaded_tables[57] = true;
               list_of_tables[63] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"MET_nh.npy").get_6d_array(),gridsize6d); loaded_tables[63] = true;
               list_of_tables[69] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ASN_nh.npy").get_5d_array(),gridsize5d); loaded_tables[69] = true;
               list_of_tables[75] =  new ArrayWrapperStandard(cnpy::npy_load(path+"PRO_nh.npy").get_3d_array()); loaded_tables[75] = true;
               list_of_tables[81] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"GLN_nh.npy").get_6d_array(),gridsize6d); loaded_tables[81] = true;
               list_of_tables[87] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ARG_nh.npy").get_7d_array(),gridsize7d); loaded_tables[87] = true;
               list_of_tables[93] = new ArrayWrapperInterpolate(cnpy::npy_load(path+"SER_nh.npy").get_4d_array(),gridsize4d); loaded_tables[93] = true;
               list_of_tables[99] = new ArrayWrapperInterpolate(cnpy::npy_load(path+"THR_nh.npy").get_5d_array(),gridsize5d); loaded_tables[99] = true;
               list_of_tables[105] = new ArrayWrapperInterpolate(cnpy::npy_load(path+"VAL_nh.npy").get_4d_array(),gridsize4d); loaded_tables[105] = true;
               list_of_tables[111] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"TRP_nh.npy").get_5d_array(),gridsize5d); loaded_tables[111] = true;
               list_of_tables[117] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"TYR_nh.npy").get_5d_array(),gridsize5d); loaded_tables[117] = true;
          }	
          if (load_hn){
               std::cout << "HN " << std::flush;
               list_of_tables[4] =  new ArrayWrapperStandard(cnpy::npy_load(path+"ALA_hn.npy").get_3d_array());  loaded_tables[4] = true;
               list_of_tables[10] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"CYS_hn.npy").get_4d_array(),gridsize4d); loaded_tables[10] = true;
               list_of_tables[16] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ASP_hn.npy").get_5d_array(),gridsize5d); loaded_tables[16] = true;
               list_of_tables[22] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"GLU_hn.npy").get_6d_array(),gridsize6d); loaded_tables[22] = true;
               list_of_tables[28] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"PHE_hn.npy").get_5d_array(),gridsize5d); loaded_tables[28] = true;
               list_of_tables[34] = new ArrayWrapperStandard(cnpy::npy_load(path+"GLY_hn.npy").get_3d_array()); loaded_tables[34] = true;
               list_of_tables[40] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"HIS_hn.npy").get_5d_array(),gridsize5d); loaded_tables[40] = true;
               list_of_tables[46] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ILE_hn.npy").get_5d_array(),gridsize5d); loaded_tables[46] = true;
               list_of_tables[52] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"LYS_hn.npy").get_7d_array(),gridsize7d); loaded_tables[52] = true;
               list_of_tables[58] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"LEU_hn.npy").get_5d_array(),gridsize5d); loaded_tables[58] = true;
               list_of_tables[64] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"MET_hn.npy").get_6d_array(),gridsize6d); loaded_tables[64] = true;
               list_of_tables[70] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ASN_hn.npy").get_5d_array(),gridsize5d); loaded_tables[70] = true;
               list_of_tables[76] =  new ArrayWrapperStandard(cnpy::npy_load(path+"PRO_hn.npy").get_3d_array()); loaded_tables[76] = true;
               list_of_tables[82] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"GLN_hn.npy").get_6d_array(),gridsize6d); loaded_tables[82] = true;
               list_of_tables[88] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ARG_hn.npy").get_7d_array(),gridsize7d); loaded_tables[88] = true;
               list_of_tables[94] = new ArrayWrapperInterpolate(cnpy::npy_load(path+"SER_hn.npy").get_4d_array(),gridsize4d); loaded_tables[94] = true;
               list_of_tables[100] = new ArrayWrapperInterpolate(cnpy::npy_load(path+"THR_hn.npy").get_5d_array(),gridsize5d); loaded_tables[100] = true;
               list_of_tables[106] = new ArrayWrapperInterpolate(cnpy::npy_load(path+"VAL_hn.npy").get_4d_array(),gridsize4d); loaded_tables[106] = true;
               list_of_tables[112] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"TRP_hn.npy").get_5d_array(),gridsize5d); loaded_tables[112] = true;
               list_of_tables[118] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"TYR_hn.npy").get_5d_array(),gridsize5d); loaded_tables[118] = true;
          }

          if (load_ha){
               std::cout << "HA " << std::flush;
               list_of_tables[5] =  new ArrayWrapperStandard(cnpy::npy_load(path+"ALA_ha.npy").get_3d_array());  loaded_tables[5] = true;
               list_of_tables[11] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"CYS_ha.npy").get_4d_array(),gridsize4d); loaded_tables[11] = true;
               list_of_tables[17] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ASP_ha.npy").get_5d_array(),gridsize5d); loaded_tables[17] = true;
               list_of_tables[23] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"GLU_ha.npy").get_6d_array(),gridsize6d); loaded_tables[23] = true;
               list_of_tables[29] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"PHE_ha.npy").get_5d_array(),gridsize5d); loaded_tables[29] = true;
               list_of_tables[35] = new ArrayWrapperStandard(cnpy::npy_load(path+"GLY_ha.npy").get_3d_array()); loaded_tables[35] = true;
               list_of_tables[41] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"HIS_ha.npy").get_5d_array(),gridsize5d); loaded_tables[41] = true;
               list_of_tables[47] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ILE_ha.npy").get_5d_array(),gridsize5d); loaded_tables[47] = true;
               list_of_tables[53] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"LYS_ha.npy").get_7d_array(),gridsize7d); loaded_tables[53] = true;
               list_of_tables[59] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"LEU_ha.npy").get_5d_array(),gridsize5d); loaded_tables[59] = true;
               list_of_tables[65] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"MET_ha.npy").get_6d_array(),gridsize6d); loaded_tables[65] = true;
               list_of_tables[71] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ASN_ha.npy").get_5d_array(),gridsize5d); loaded_tables[71] = true;
               list_of_tables[77] =  new ArrayWrapperStandard(cnpy::npy_load(path+"PRO_ha.npy").get_3d_array()); loaded_tables[77] = true;
               list_of_tables[83] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"GLN_ha.npy").get_6d_array(),gridsize6d); loaded_tables[83] = true;
               list_of_tables[89] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"ARG_ha.npy").get_7d_array(),gridsize7d); loaded_tables[89] = true;
               list_of_tables[95] = new ArrayWrapperInterpolate(cnpy::npy_load(path+"SER_ha.npy").get_4d_array(),gridsize4d); loaded_tables[95] = true;
               list_of_tables[101] = new ArrayWrapperInterpolate(cnpy::npy_load(path+"THR_ha.npy").get_5d_array(),gridsize5d); loaded_tables[101] = true;
               list_of_tables[107] = new ArrayWrapperInterpolate(cnpy::npy_load(path+"VAL_ha.npy").get_4d_array(),gridsize4d); loaded_tables[107] = true;
               list_of_tables[113] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"TRP_ha.npy").get_5d_array(),gridsize5d); loaded_tables[113] = true;
               list_of_tables[119] =  new ArrayWrapperInterpolate(cnpy::npy_load(path+"TYR_ha.npy").get_5d_array(),gridsize5d); loaded_tables[119] = true;
         }


     }


     //! Constructor
     ProCS15Backend (phaistos::ChainFB &chain, bool load_ca, bool load_cb, bool load_co, bool load_n, bool load_hn, bool load_ha, std::string procsnumpypath) {

          using namespace procs15;
	 
          std::cout << "Loading ProCS15 data files" << std::endl;

          check_path( load_ca, load_cb, load_co, load_n, load_hn, load_ha, procsnumpypath );       

          std::cout << "Loading Hydrogen bond data " << std::endl;

          std::cout << procsnumpypath << "delta_hbond_251_91_361_12.npy"
                       << std::endl;

          hbond_array = new HbondDataWrapper(cnpy::npy_load(procsnumpypath+"delta_hbond_251_91_361_12.npy").get_4d_array()); 


          std::cout << procsnumpypath << "delta_hbond_Carboxylate_251_91_361_6.npy"
                       << std::endl;
          hbond_array_carboxy = new HbondDataWrapper(cnpy::npy_load(procsnumpypath+"delta_hbond_Carboxylate_251_91_361_6.npy").get_4d_array());


          std::cout << procsnumpypath << "delta_hbond_Alcohol_251_91_361_6.npy"
                       << std::endl;
          hbond_array_alcohol = new HbondDataWrapper(cnpy::npy_load(procsnumpypath+"delta_hbond_Alcohol_251_91_361_6.npy").get_4d_array());


          std::cout << procsnumpypath << "delta_halphabond_oxygen_221_91_361_12.npy"
                       << std::endl;
          halphabond_array_oxygen = new HbondDataWrapper(cnpy::npy_load(procsnumpypath+"delta_halphabond_oxygen_221_91_361_12.npy").get_4d_array());


          std::cout << procsnumpypath << "delta_halphabond_carboxy_221_91_361_6.npy"
                       << std::endl;
          halphabond_array_carboxy = new HbondDataWrapper(cnpy::npy_load(procsnumpypath+"delta_halphabond_carboxy_221_91_361_6.npy").get_4d_array());

          std::cout << procsnumpypath << "delta_halphabond_alcohol_221_91_361_6.npy"
                       << std::endl;
          halphabond_array_alcohol = new HbondDataWrapper(cnpy::npy_load(procsnumpypath+"delta_halphabond_alcohol_221_91_361_6.npy").get_4d_array());

          list_of_tables.reserve(6*20);
          for (unsigned i=0; i<120; i++) loaded_tables.push_back(false);

	      std::cout << "Loading Chemical shift data: "; 
          load_data_linearinterpolation( load_ca, load_cb, load_co, load_n, load_hn, load_ha, procsnumpypath );
 
          std::cout << std::endl;

          if (load_ca) { 
              atypevector.push_back(ProCS15_CA);
          }
          if (load_cb) { 
              atypevector.push_back(ProCS15_CB);
          }
          if (load_co) { 
              atypevector.push_back(ProCS15_C);
          }
          if (load_n)  { 
              atypevector.push_back(ProCS15_N);
          }
          if (load_hn) { 
              atypevector.push_back(ProCS15_HN);
          }
          if (load_ha) { 
              atypevector.push_back(ProCS15_HA);
          }
     }

     //! Copy Constructor
     ProCS15Backend (const ProCS15Backend &other, phaistos::ChainFB& chain) {

          hbond_array = other.hbond_array;
          hbond_array_carboxy = other.hbond_array_carboxy;
          hbond_array_alcohol = other.hbond_array_alcohol;
          halphabond_array_oxygen = other.halphabond_array_oxygen;
          halphabond_array_carboxy = other.halphabond_array_carboxy;
          halphabond_array_alcohol = other.halphabond_array_alcohol;

          list_of_tables.reserve(6*20);
          for (unsigned i=0; i < 120; i++)(list_of_tables[i] = other.list_of_tables[i]);

          loaded_tables = other.loaded_tables;
          atypevector = other.atypevector;

     }


}; // End class ProCS15Backend



}; // End namespace phaistos




#endif
