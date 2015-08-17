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

#ifndef TERM_PROCS15_BASE
#define TERM_PROCS15_BASE


#include <boost/random/uniform_smallint.hpp>
#include <boost/random/uniform_int.hpp>

#include "energy/energy_term.h"
#include "protein/chain_fb.h"
#include "energy/energy.h"


#include "cs_file_parser.h"
#include<boost/range/numeric.hpp>

using namespace lib_definitions;

namespace phaistos {

//! Class for the chemical shift energy term
template <typename DERIVED_CLASS>
class TermProCS15Base: public EnergyTermCommon<DERIVED_CLASS, ChainFB> {

private:

     //! For convenience, define local EnergyTermCommon.
     typedef phaistos::EnergyTermCommon<DERIVED_CLASS, ChainFB> EnergyTermCommon;

public:

     //! Pointer to the global random number engine.
     RandomNumberEngine *random_number_engine;

     //! Local setting class
     const class Settings: public EnergyTerm<ChainFB>::Settings {
     public:

         //! An NMR STAR formatted chemical shifts file.
         std::string star_filename;

         // Sample weights or use fixed default values.
         bool sample_weights;

         // How to model intercept
         std::string intercept_model;

         // How to model slope
         std::string slope_model;

         //! Which chemical shift data files should be loaded.
         bool load_ca;
         bool load_cb;
         bool load_co;
         bool load_n;
         bool load_hn;
         bool load_ha;
         bool include_ha_rc;
         bool include_hn_rc;
         bool include_hn_previous_residue_correction;
         bool include_ha_previous_residue_correction;
         bool include_ca_previous_residue_correction;
         bool include_cb_previous_residue_correction;
         bool include_co_previous_residue_correction;
         bool include_n_previous_residue_correction; //=1
         bool include_hn_following_residue_correction;
         bool include_ha_following_residue_correction;
         bool include_ca_following_residue_correction;
         bool include_cb_following_residue_correction;
         bool include_co_following_residue_correction;
         bool include_n_following_residue_correction;
         bool include_hn_primary_hn_hbond; //=1
         bool include_ha_primary_hn_hbond;
         bool include_ca_primary_hn_hbond;
         bool include_cb_primary_hn_hbond;
         bool include_co_primary_hn_hbond;
         bool  include_n_primary_hn_hbond;
         bool include_hn_primary_ha_hbond;
         bool include_ha_primary_ha_hbond; //=1
         bool include_ca_primary_ha_hbond; //=0
         bool include_cb_primary_ha_hbond; //=0
         bool include_co_primary_ha_hbond;
         bool  include_n_primary_ha_hbond; //=0
         bool include_hn_secondary_hn_hbond;
         bool include_ha_secondary_hn_hbond;
         bool include_ca_secondary_hn_hbond;
         bool include_cb_secondary_hn_hbond;
         bool include_co_secondary_hn_hbond; //=1
         bool  include_n_secondary_hn_hbond; //=1
         bool include_hn_secondary_ha_hbond;
         bool include_ha_secondary_ha_hbond;
         bool include_ca_secondary_ha_hbond;
         bool include_cb_secondary_ha_hbond;
         bool include_co_secondary_ha_hbond;
         bool  include_n_secondary_ha_hbond;


         bool output_full_prediction_vector;

         //! Add correction to amide proton chemical shifts when not forming hydrogen bond.
         bool use_water_correction;

         //! Path to procs15 numpy files
         std::string procsnumpypath;


         //! Constructor and a reasonable default values.
         Settings(std::string star_filename="",
                  bool sample_weights = false,
                  std::string intercept_model = "MLE",
                  std::string slope_model = "MLE",
                  bool load_ca = true,
                  bool load_cb = true,
                  bool load_co = true,
                  bool load_n  = true,
                  bool load_hn = true,
                  bool load_ha = true,
                  bool include_ha_rc = true,
                  bool include_hn_rc = true,
                  bool include_hn_previous_residue_correction = true,
                  bool include_ha_previous_residue_correction = true,
                  bool include_ca_previous_residue_correction = true,
                  bool include_cb_previous_residue_correction = true,
                  bool include_co_previous_residue_correction = true,
                  bool  include_n_previous_residue_correction = true,
                  bool include_hn_following_residue_correction = true,
                  bool include_ha_following_residue_correction = true,
                  bool include_ca_following_residue_correction = true,
                  bool include_cb_following_residue_correction = true,
                  bool include_co_following_residue_correction = true,
                  bool  include_n_following_residue_correction = true,
                  bool include_hn_primary_hn_hbond = true,
                  bool include_ha_primary_hn_hbond = true,
                  bool include_ca_primary_hn_hbond = false,
                  bool include_cb_primary_hn_hbond = false,
                  bool include_co_primary_hn_hbond = false,
                  bool  include_n_primary_hn_hbond = true,
                  bool include_hn_primary_ha_hbond = true,
                  bool include_ha_primary_ha_hbond = true,
                  bool include_ca_primary_ha_hbond = false,
                  bool include_cb_primary_ha_hbond = false,
                  bool include_co_primary_ha_hbond = false,
                  bool  include_n_primary_ha_hbond = false,
                  bool include_hn_secondary_hn_hbond = true,
                  bool include_ha_secondary_hn_hbond = true,
                  bool include_ca_secondary_hn_hbond = true,
                  bool include_cb_secondary_hn_hbond = false,
                  bool include_co_secondary_hn_hbond = true,
                  bool  include_n_secondary_hn_hbond = true,
                  bool include_hn_secondary_ha_hbond = true,
                  bool include_ha_secondary_ha_hbond = true,
                  bool include_ca_secondary_ha_hbond = true,
                  bool include_cb_secondary_ha_hbond = false,
                  bool include_co_secondary_ha_hbond = true,
                  bool  include_n_secondary_ha_hbond = true,
                  bool output_full_prediction_vector = false,
                  bool use_water_correction = true,
                  std::string procsnumpypath = "./" )
              : star_filename(star_filename),
                sample_weights(sample_weights),
                intercept_model(intercept_model),
                slope_model(slope_model),
                load_ca(load_ca),
                load_cb(load_cb),
                load_co(load_co),
                load_n(load_n),
                load_hn(load_hn),
                load_ha(load_ha),
                include_ha_rc(include_ha_rc),
                include_hn_rc(include_hn_rc),
                include_hn_previous_residue_correction(include_hn_previous_residue_correction),
                include_ha_previous_residue_correction(include_ha_previous_residue_correction),
                include_ca_previous_residue_correction(include_ca_previous_residue_correction),
                include_cb_previous_residue_correction(include_cb_previous_residue_correction),
                include_co_previous_residue_correction(include_co_previous_residue_correction),
                include_n_previous_residue_correction( include_n_previous_residue_correction),
                include_hn_following_residue_correction(include_hn_following_residue_correction),
                include_ha_following_residue_correction(include_ha_following_residue_correction),
                include_ca_following_residue_correction(include_ca_following_residue_correction),
                include_cb_following_residue_correction(include_cb_following_residue_correction),
                include_co_following_residue_correction(include_co_following_residue_correction),
                 include_n_following_residue_correction( include_n_following_residue_correction),
                include_hn_primary_hn_hbond(  include_hn_primary_hn_hbond),
                include_ha_primary_hn_hbond(  include_ha_primary_hn_hbond),
                include_ca_primary_hn_hbond(  include_ca_primary_hn_hbond),
                include_cb_primary_hn_hbond(  include_cb_primary_hn_hbond),
                include_co_primary_hn_hbond(  include_co_primary_hn_hbond),
                 include_n_primary_hn_hbond(   include_n_primary_hn_hbond),
                include_hn_primary_ha_hbond(  include_hn_primary_ha_hbond),
                include_ha_primary_ha_hbond(  include_ha_primary_ha_hbond),
                include_ca_primary_ha_hbond(  include_ca_primary_ha_hbond),
                include_cb_primary_ha_hbond(  include_cb_primary_ha_hbond),
                include_co_primary_ha_hbond(  include_co_primary_ha_hbond),
                 include_n_primary_ha_hbond(   include_n_primary_ha_hbond),
              include_hn_secondary_hn_hbond(include_hn_secondary_hn_hbond),
              include_ha_secondary_hn_hbond(include_ha_secondary_hn_hbond),
              include_ca_secondary_hn_hbond(include_ca_secondary_hn_hbond),
              include_cb_secondary_hn_hbond(include_cb_secondary_hn_hbond),
              include_co_secondary_hn_hbond(include_co_secondary_hn_hbond),
               include_n_secondary_hn_hbond( include_n_secondary_hn_hbond),
              include_hn_secondary_ha_hbond(include_hn_secondary_ha_hbond),
              include_ha_secondary_ha_hbond(include_ha_secondary_ha_hbond),
              include_ca_secondary_ha_hbond(include_ca_secondary_ha_hbond),
              include_cb_secondary_ha_hbond(include_cb_secondary_ha_hbond),
              include_co_secondary_ha_hbond(include_co_secondary_ha_hbond),
               include_n_secondary_ha_hbond( include_n_secondary_ha_hbond),
                output_full_prediction_vector(output_full_prediction_vector),
                use_water_correction(use_water_correction),
                procsnumpypath(procsnumpypath){}


         //! Output operator.
         friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
              o << "star-filename:"    << settings.star_filename << "\n";
              o << "sample-weights:"   << settings.sample_weights << "\n";
              o << "intercept-model:"  << settings.intercept_model << "\n";
              o << "slope_model:"      << settings.slope_model << "\n";
              o << "load_cb:"      << settings.load_cb << "\n";
              o << "load_ca:"      << settings.load_ca << "\n";
              o << "load_co:"      << settings.load_co << "\n";
              o << "load_n:"      << settings.load_n << "\n";
              o << "load_hn:"      << settings.load_hn << "\n";
              o << "load_ha:"      << settings.load_ha << "\n";
              o << "HN-water-bond-correction:"   << settings.use_water_correction << "\n";
              o << static_cast<const typename EnergyTerm<ChainFB>::Settings>(settings);
              return o;
         }
     } settings;

     // If last move was none move or not.
     bool none_move;

     //! Weights in the energy term in the current MC step.
     FPlist weights;

     //! Weights in the energy term in the previous MC step.
     FPlist weights_previous;

     //! Parameters for scaling the predicted chemical shieldings.
     FPlist slopes;
     FPlist slopes_previous;
     FPlist intercepts;
     FPlist intercepts_previous;

     //! Counter for printing at a specified frequency.
     unsigned int n_print;

     //! A table containing the chemical shift from the user-defined input file
     FPtable experimental_chemical_shifts;

     //! Holds the predicted chemical shifts.
     FPtable predicted_chemical_shifts;
     FPtable predicted_chemical_shifts_previous;

     //! Local copy of thread id.
     int thread_id;

     //! Returns samples according to the Normal distribution.·
     //! \param m Mean value.
     //! \param std standard deviation
     FPtype rnom(FPtype m, FPtype std, RandomNumberEngine *rne) {
          boost::normal_distribution<> nom(m, std);
          boost::variate_generator<RandomNumberEngine&, boost::normal_distribution<> > norml( *rne, nom );
          return norml();
     }

     //! Returns a random integer from {min, min+1, ... max-1, max}.
     //! \param min Lowest number in the range.
     //! \param max Highest number in the range.
     //! \return A random number.
     int rand_int(const int min, const int max, RandomNumberEngine *rne) {
          boost::uniform_smallint<> distribution(min, max);
          boost::variate_generator<RandomNumberEngine&, boost::uniform_smallint<> > generator(*rne, distribution);
          return generator();
     }

     //! Runs the ProCS15 predictor on a protein structure. This is a virtual
     //! function, and the implementation is found in the derived classes, TermProCS15 and TermProCS15Cached.
     //! \param chain The chain object
     //! \return A Matrix containing six chemical shifts for each residue
     virtual FPtable predict(phaistos::ChainFB& chain, MoveInfo *move_info=NULL)=0;


     //! Constructor
     TermProCS15Base(ChainFB *chain, std::string name, const Settings &settings=Settings(),
                 RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, name, settings, random_number_engine),
            random_number_engine(random_number_engine),
            settings(settings) {
          using namespace lib_definitions;


          ///! Check that filename exists
          if (!file_exists(settings.star_filename)) {
               std::cerr << "ERROR (ProCS15): File \"" << settings.star_filename << "\" not found. Exiting\n";
               exit(1);
          }

          // HA, CA, H, N, C, CB
          this->weights_previous = 
          this->weights = vector_utils::make_vector<FPtype>(0.40, 1.65, 0.60, 4.15, 1.70, 2.00);

          this->slopes_previous = 
          this->slopes = vector_utils::make_vector<FPtype>(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);

          this->intercepts_previous = 
          this->intercepts = vector_utils::make_vector<FPtype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

          //! Read-in chemical shifts and store in this table.
          this->experimental_chemical_shifts = cs_parser::value_matrix_from_starfile(settings.star_filename, *chain);

          this->n_print = 10000;

          this->none_move = false;

     }

     //! Copy constructor
     TermProCS15Base(const TermProCS15Base &other,
             RandomNumberEngine *random_number_engine,
             int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            random_number_engine(other.random_number_engine),
            settings(other.settings),
            none_move(other.none_move),
            weights(other.weights),
            weights_previous(other.weights_previous),
            slopes(other.slopes),
            slopes_previous(other.slopes_previous),
            intercepts(other.intercepts),
            intercepts_previous(other.intercepts_previous),
            n_print(other.n_print),
            experimental_chemical_shifts(other.experimental_chemical_shifts),
            predicted_chemical_shifts(other.predicted_chemical_shifts),
            predicted_chemical_shifts_previous(other.predicted_chemical_shifts_previous) {

            // Give each thread a defined random_number_engine (so runs can be replicated).
            this->random_number_engine = random_number_engine;

            // Don't copy the thread id (as these are obviously different for each thread).
            this->thread_id = thread_index;
 
     }


     // print boost multiarray
     template <typename Array>
     void print(std::ostream& os, const Array& A)
     {
           typename Array::const_iterator i;
             os << "[";
               for (i = A.begin(); i != A.end(); ++i) {
                       print(os, *i);
                           if (boost::next(i) != A.end())
                                     os << ',';
                             }
                 os << "]" << std::endl;;
     }
     void print(std::ostream& os, const double& x)
     {
           os << x;
     }
     void print(std::ostream& os, const float& x)
     {
           os << x;
     }



     double get_normal_marginalized_slope_marginalized_intercept(FPtable &predicted_cs, FPtable &experimental_cs) {
          double energy = 0.0;
          double cs1;
          double cs2;
          for (unsigned int j=0; j < 6; j++) {
               double chi_sq = 0.0;
               unsigned int n = 0;
               double xmean = 0.0;
               double ymean = 0.0;
               for (unsigned int i=1;i < this->chain->size()-1;i++) {
                   cs1 = predicted_cs[i][j];
                   cs2 = experimental_cs[i][j];
                   if ((std::abs(cs1) > 0.000001) &&
                       (std::abs(cs2) > 0.000001)) {

                       xmean += cs1;
                       ymean += cs2;
                       n++;
                   }
               }

               xmean /= (double)n;
               ymean /= (double)n;

               double b1 = 0.0;
               double a = 0.0;
               double b0;
               for (unsigned int i=1;i < this->chain->size()-1;i++) {
                   cs1 = predicted_cs[i][j];
                   cs2 = experimental_cs[i][j];
                   if ((std::abs(cs1) > 0.000001) &&
                       (std::abs(cs2) > 0.000001)) {

                       b1 += (cs1 - xmean) * (cs2 - ymean);
                       a += (cs1 - xmean) * (cs1 - xmean);
                   }
               }
               b1 /= a;
               b0 = ymean-b1*xmean;
               for (unsigned int i=1;i < this->chain->size()-1;i++) {
                   cs1 = predicted_cs[i][j];
                   cs2 = experimental_cs[i][j];
                   if ((std::abs(cs1) > 0.000001) &&
                       (std::abs(cs2) > 0.000001)) {

                       chi_sq += (cs2-(b0+b1*cs1))*(cs2-(b0+b1*cs1));
                   }
               }
               if (n > 0) energy += 0.5 * n * std::log(chi_sq);
          }
          return energy;
     }

     double get_energy(FPtable &predicted_cs, FPtable &experimental_cs) {
          return get_normal_marginalized_slope_marginalized_intercept(predicted_cs, experimental_cs);
     }

     //! Update single weight.
     void sample_parameter(FPtype &parameter, FPtype stdev) {

          FPtype delta = rnom(0.0,stdev, this->random_number_engine);

          parameter *= std::exp(delta);
     }

     void update_parameter() {

          // Get random atom type
          unsigned int i = rand_int(0,5, this->random_number_engine);

          if (settings.sample_weights) {
              this->sample_parameter(this->weights[i], 0.15);
          }
          if (settings.intercept_model == "sample") {
              this->sample_parameter(this->intercepts[i], 0.15);
          }
          if (settings.slope_model == "sample") {
              this->sample_parameter(this->slopes[i], 0.05);
          }
     }




     //! Calculate the energy and updates the precicted chemical shift table or parameters for the energy term.
     //! \param move_info A MoveInfo object
     //! \return The energy associated with the chemical shift prediction.
     double evaluate(MoveInfo *move_info=NULL) {

         // Set none_move to false as default
          this->none_move = false; 

          // Check if current move was a none-move
          if (move_info) {
               if (move_info->modified_angles.empty() == true) {

                    this->none_move = true;
                    this->update_parameter();

               }
          }
          this->predicted_chemical_shifts = predict(*(this->chain), move_info);

          return get_energy(this->predicted_chemical_shifts,this->experimental_chemical_shifts);

     }

     //! Accept last energy evaluation. Save predictions and model parameters.
     void accept() {

          // If the move was a none-move, energy term parameters
          if (this->none_move) {

               this->weights_previous = this->weights;
               this->intercepts_previous = this->intercepts;
               this->slopes_previous = this->slopes;

          // Else it was a physical move, and backup chemical shifts
          } else {

               this->predicted_chemical_shifts_previous = this->predicted_chemical_shifts;
          }

          print_parameters();

     }

     //! Virtual function to avoid duplicate reject() function in cached and uncached terms
     virtual void reject_cache()=0;

     //! Reject last energy evaluation. Restore predictions and model parameters.
     void reject() {

          // Reject cache
          reject_cache();


          if (this->none_move) {

               this->weights = this->weights_previous;
               this->intercepts = this->intercepts_previous;
               this->slopes = this->slopes_previous;

          } else {

               this->predicted_chemical_shifts = this->predicted_chemical_shifts_previous;
          }

          print_parameters();

     }

}; // End class TermProCS15Base



} // end namespace phaistos
#endif
