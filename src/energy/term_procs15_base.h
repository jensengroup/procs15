// Copyright (C) 2014-2014 by Anders S. Christensen, Anders S. Larsen, Lars A. Bratholm
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
#include <iomanip>
#include <stdlib.h>

#include "procs15_backend.h"

namespace phaistos {
using namespace procs15;

//! Class for the chemical shift energy term
template <typename DERIVED_CLASS>
class TermProCS15Base: public EnergyTermCommon<DERIVED_CLASS, ChainFB>,
                       public ProCS15Backend {

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

         // Probability density function used in energy_term
         std::string energy_type;

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

         //! Add correction to amide proton chemical shifts when not forming hydrogen bond.
         bool use_water_correction;

         //! Path to procs15 numpy files
         std::string procsnumpypath;


         //! Constructor and a reasonable default values.
         Settings(std::string star_filename="",
                  std::string energy_type = "marginalized",
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
                  bool use_water_correction = true,
                  std::string procsnumpypath = "./" )
              : star_filename(star_filename),
                energy_type(energy_type),
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
                use_water_correction(use_water_correction),
                procsnumpypath(procsnumpypath){}


         //! Output operator.
         friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
              o << "star-filename:"    << settings.star_filename << "\n";
              o << "energy-type:"   << settings.energy_type << "\n";
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

     //! Holds memberships in mixture model
     std::vector<bool> member_previous;
     std::vector<bool> member;

     //! Local copy of thread id.
     int thread_id;

     //! Ringcurrent cutoff distance
     double rc_cutoff;

     //! HA hbond distance (max 4)
     double ha_hbond_cutoff;

     //! HN hbond distance (max 4)
     double hn_hbond_cutoff;

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
            ProCS15Backend(*chain, settings.load_ca, settings.load_cb, settings.load_co, settings.load_n, settings.load_hn, settings.load_ha, settings.procsnumpypath),
            random_number_engine(random_number_engine),
            settings(settings) {
          using namespace procs15;


          ///! Check that filename exists
          if (!file_exists(settings.star_filename)) {
               std::cerr << "ERROR (ProCS15): File \"" << settings.star_filename << "\" not found. Exiting\n";
               exit(1);
          }

          // HA, CA, H, N, C, CB (twice in case of mixture model)
          this->weights_previous = 
          this->weights = vector_utils::make_vector<FPtype>(0.36, 1.83, 0.61, 3.78, 1.66, 1.97, 0.20, 0.80, 0.30, 2.00, 0.90, 1.00);

          this->slopes_previous = 
          this->slopes = vector_utils::make_vector<FPtype>(0.467, 1.02, 0.154, 0.463, 0.385, 0.939);

          this->intercepts_previous = 
          this->intercepts = vector_utils::make_vector<FPtype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

          // When slopes and intercepts are fixed, the energy term is parametrised slightly different:
          if ((settings.energy_type == "fix-slope-intercept") || (settings.energy_type == "fix-all")) {
               this->intercepts_previous = 
               this->intercepts = vector_utils::make_vector<FPtype>(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
          }

          //! Read-in chemical shifts and store in this table.
          this->experimental_chemical_shifts = cs_parser::value_matrix_from_starfile(settings.star_filename, *chain);

          this->n_print = 100000;

          if (settings.energy_type == "mixture") {
               for (unsigned int i=0; i < this->chain->size(); i++) {
                    this->member.push_back((bool)rand_int(0, 1, random_number_engine));
               }
          } else {
               this->member.resize(this->chain->size(), true);
          }
          this->member_previous = this->member;

          this->none_move = false;

          // Don't use CYS CB shifts, if more than two CYS residues are in the chain.
          unsigned int count_cys = 0;
          for (unsigned int i = 0; i < this->experimental_chemical_shifts.size(); i++) {
               if ((*(this->chain))[i].residue_type == definitions::CYS) {
                   count_cys++;
               }
          }
          if (count_cys > 1) {
               for (unsigned int i = 0; i < this->experimental_chemical_shifts.size(); i++) {

                    if ((*(this->chain))[i].residue_type == definitions::CYS){
                         this->experimental_chemical_shifts[i][5] = 0.0;
                    }
               }
          }

          rc_cutoff = 8;
          ha_hbond_cutoff = 4;
          hn_hbond_cutoff = 3;

     }

     //! Copy constructor
     TermProCS15Base(const TermProCS15Base &other,
             RandomNumberEngine *random_number_engine,
             int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            ProCS15Backend(other, *chain),
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
            predicted_chemical_shifts_previous(other.predicted_chemical_shifts_previous),
            rc_cutoff(other.rc_cutoff),
            ha_hbond_cutoff(other.ha_hbond_cutoff),
            hn_hbond_cutoff(other.hn_hbond_cutoff),
            member(other.member),
            member_previous(other.member_previous) {

            // Give each thread a defined random_number_engine (so runs can be replicated).
            this->random_number_engine = random_number_engine;

            // Don't copy the thread id (as these are obviously different for each thread).
            this->thread_id = thread_index;
 
     }

     double get_normal_mle_slope_mle_intercept_mle_sigma(FPtable &predicted_cs, FPtable &experimental_cs) {
          double energy = 0.0;
          double cs1;
          double cs2;
          for (unsigned int j=0; j < 6; j++) {
               unsigned int n = 0;
               double xsum = 0.0;
               double ysum = 0.0;
               double xysum = 0.0;
               double x2sum = 0.0;
               double y2sum = 0.0;
               for (unsigned int i=1;i < this->chain->size()-1;i++) {
                   cs1 = predicted_cs[i][j];
                   cs2 = experimental_cs[i][j];
                   if ((std::abs(cs1) > 0.000001) &&
                       (std::abs(cs2) > 0.000001)) {

                       xsum += cs1;
                       ysum += cs2;
                       xysum += cs1*cs2;
                       x2sum += cs1*cs1;
                       y2sum += cs2*cs2;
                       n++;
                   }
               }

               if (n > 0) {
                    double chi_sq = y2sum + ( xysum * (2*xsum * ysum - n* xysum) - x2sum * ysum * ysum ) / (n * x2sum - xsum*xsum);
                    energy += 0.5*(n+1)*std::log(chi_sq);
               }
          }
          return energy;
     }


     double get_normal_marg_slope_marg_intercept_sampled_sigma(FPtable &predicted_cs, FPtable &experimental_cs) {
          double energy = 0.0;
          double cs1;
          double cs2;
          for (unsigned int j=0; j < 6; j++) {
               unsigned int n = 0;
               double xsum = 0.0;
               double ysum = 0.0;
               double xysum = 0.0;
               double x2sum = 0.0;
               double y2sum = 0.0;
               for (unsigned int i=1;i < this->chain->size()-1;i++) {
                   cs1 = predicted_cs[i][j];
                   cs2 = experimental_cs[i][j];
                   if ((std::abs(cs1) > 0.000001) &&
                       (std::abs(cs2) > 0.000001)) {

                       xsum += cs1;
                       ysum += cs2;
                       xysum += cs1*cs2;
                       x2sum += cs1*cs1;
                       y2sum += cs2*cs2;
                       n++;
                   }
               }

               if (n > 0) {
                    // update the slopes vector with the maximum likelihood estimates
                    // Note that the slopes are still marginalized in the energy
                    this->slopes[j] = (n * xysum - xsum * ysum) / (n * x2sum - xsum*xsum);
                    double chi_sq = (n * y2sum - ysum * ysum - this->slopes[j] * (n * xysum - xsum * ysum)) / n;
                    energy += (n-1)*std::log(this->weights[j]) + 0.5 * std::log(n * x2sum - xsum*xsum)
                           + 0.5*chi_sq/(this->weights[j]*this->weights[j]);
                    // noninformative prior slope and offset
                    energy += 1.5 * std::log(this->slopes[j]*this->slopes[j] + 1);
               }
          }
          return energy;
     }

     double get_normal_marg_slope_marg_intercept_marg_sigma(FPtable &predicted_cs, FPtable &experimental_cs) {
          double energy = 0.0;
          double cs1;
          double cs2;
          for (unsigned int j=0; j < 6; j++) {
               unsigned int n = 0;
               double xsum = 0.0;
               double ysum = 0.0;
               double xysum = 0.0;
               double xsum2 = 0.0;
               double ysum2 = 0.0;
               for (unsigned int i=1;i < this->chain->size()-1;i++) {
                   cs1 = predicted_cs[i][j];
                   cs2 = experimental_cs[i][j];
                   if ((std::abs(cs1) > 0.000001) &&
                       (std::abs(cs2) > 0.000001)) {

                       xsum += cs1;
                       ysum += cs2;
                       xysum += cs1*cs2;
                       xsum2 += cs1*cs1;
                       ysum2 += cs2*cs2;
                       n++;
                   }
               }

               if (n > 0) {
                    double chi_sq = ysum2 + ( xysum * (2*xsum * ysum - n* xysum) - xsum2 * ysum * ysum ) / (n * xsum2 - xsum*xsum);
                    energy += 0.5*(n-2)*std::log(chi_sq) + 0.5 * std::log(n * xsum2 - xsum*xsum);
               }
          }
          return energy;
     }

     double get_normal_sampled_slope_sampled_intercept_sampled_sigma(FPtable &predicted_cs, FPtable &experimental_cs) {
          double energy = 0.0;
          double cs1;
          double cs2;
          for (unsigned int j=0; j < 6; j++) {
               unsigned int n = 0;
               double xsum = 0.0;
               double ysum = 0.0;
               double xysum = 0.0;
               double xsum2 = 0.0;
               double ysum2 = 0.0;
               for (unsigned int i=1;i < this->chain->size()-1;i++) {
                   cs1 = predicted_cs[i][j];
                   cs2 = experimental_cs[i][j];
                   if ((std::abs(cs1) > 0.000001) &&
                       (std::abs(cs2) > 0.000001)) {

                       xsum += cs1;
                       ysum += cs2;
                       xysum += cs1*cs2;
                       xsum2 += cs1*cs1;
                       ysum2 += cs2*cs2;
                       n++;
                   }
               }

               if (n > 0) {
                    // the intercept is sampled as the deviation from the MLE intercept, to reduce the correlation between slope and intercept.
                    double chi_sq = n*this->intercepts[j]*this->intercepts[j] + this->slopes[j]*this->slopes[j]*(xsum2-xsum*xsum/n) + 2*this->slopes[j] * (xsum*ysum/n-xysum) + ysum2 - ysum*ysum/n;
                    energy += (n+1)*std::log(this->weights[j]) + 1.5*std::log(1+this->slopes[j]*this->slopes[j]) + 0.5*chi_sq/(this->weights[j]*this->weights[j]);
               }
          }
          return energy;
     }

     double get_normal_sampled_slope_marg_intercept_sampled_sigma(FPtable &predicted_cs, FPtable &experimental_cs) {
          double energy = 0.0;
          double cs1;
          double cs2;
          for (unsigned int j=0; j < 6; j++) {
               unsigned int n = 0;
               double xsum = 0.0;
               double ysum = 0.0;
               double xysum = 0.0;
               double x2sum = 0.0;
               double y2sum = 0.0;
               for (unsigned int i=1;i < this->chain->size()-1;i++) {
                   cs1 = predicted_cs[i][j];
                   cs2 = experimental_cs[i][j];
                   if ((std::abs(cs1) > 0.000001) &&
                       (std::abs(cs2) > 0.000001)) {

                       xsum += cs1;
                       ysum += cs2;
                       xysum += cs1*cs2;
                       x2sum += cs1*cs1;
                       y2sum += cs2*cs2;
                       n++;
                   }
               }

               if (n > 0) {
                    double chi_sq = this->slopes[j]*this->slopes[j]*(x2sum-xsum*xsum/n) + 2*this->slopes[j] * (xsum*ysum/n-xysum) + y2sum - ysum*ysum/n;
                    energy += n*std::log(this->weights[j]) + 1.5*std::log(1+this->slopes[j]*this->slopes[j]) + 0.5*chi_sq/(this->weights[j]*this->weights[j]);
               }
          }
          return energy;
     }

     double get_jeffreys_prior_sigma(FPtable &predicted_cs, FPtable &experimental_cs) {
          double energy = 0.0;
          double cs1;
          double cs2;
          for (unsigned int j = 0; j < 6; j++) {
               // check that at least one prediction and experimental value exist for given type
               // could be done in a prettier fashion
               for (unsigned int i=1;i < this->chain->size()-1;i++) {
                   cs1 = predicted_cs[i][j];
                   cs2 = experimental_cs[i][j];
                   if ((std::abs(cs1) > 0.000001) &&
                       (std::abs(cs2) > 0.000001)) {
                       energy += std::log(this->weights[j]);
                       break;
                   }
               }
          }
          return energy;
     }

     // beta prime prior on sigma^2. Using the transformation p(sigma) = p(sigma^2)*2*sigma
     double get_beta_prime_prior_sigma() {
          double energy = 0.0;
          // alpha and beta fitted from vasco dataset, removing gross outliers.
          std::vector<double> alphas = vector_utils::make_vector<double>(9.48, 38.3, 12.8, 188., 72.5, 59.6);
          std::vector<double> betas =  vector_utils::make_vector<double>(69.9, 11.7, 34.8, 13.8, 27.3, 15.8);
          for (unsigned int i = 0; i < 6; i++) {
               double sigma2 = this->weights[i] * this->weights[i];
               energy += (0.5 - alphas[i]) * std::log(sigma2) + (alphas[i] + betas[i]) * std::log(1+sigma2);
          }
          return energy;
     }

     // lognormal prior on sigma. 95% confidence interval between shiftx2 values and values from charmm minimized structures
     double get_lognormal_prior_sigma() {
          double energy = 0.0;
          std::vector<double> mus = vector_utils::make_vector<double>(-1.104,0.184,-0.714,1.146,0.172,0.347);
          //std::vector<double> sigmas =  vector_utils::make_vector<double>(0.0956,0.177,0.103,0.160,0.183,0.177); //95%
          std::vector<double> sigmas =  vector_utils::make_vector<double>(0.0727,0.135,0.0787,0.121,0.139,0.135);
          for (unsigned int i = 0; i < 6; i++) {
               energy += 0.5 * (std::log(sigmas[i]) - mus[i]) * (std::log(sigmas[i]) - mus[i]) / (this->weights[i] * this->weights[i])
                       + std::log(this->weights[i]);
          }
          return energy;
     }

     // normal prior on the slope
     double get_normal_prior_slope() {
          double energy = 0.0;
          // slopes fitted from vasco dataset, removing gross outliers.
          std::vector<double> mus = vector_utils::make_vector<double>(0.457, 1.02, 0.170, 0.453, 0.379, 0.939);
          std::vector<double> sigmas = vector_utils::make_vector<double>(0.119, 0.065, 0.107, 0.104, 0.124, 0.019);
          for (unsigned int i = 0; i < 6; i++) {
               energy += 0.5 * (slopes[i] - mus[i]) * (slopes[i] - mus[i]) / (sigmas[i]*sigmas[i]);
          }
          return energy;
     }

     double get_energy(FPtable &predicted_cs, FPtable &experimental_cs) {
          if (settings.energy_type == "marginalized") {
               double energy = get_normal_marg_slope_marg_intercept_marg_sigma(predicted_cs, experimental_cs);
          }
          if (settings.energy_type == "mle") return get_normal_mle_slope_mle_intercept_mle_sigma(predicted_cs, experimental_cs);
          if (settings.energy_type == "sample-all") {
               double energy = get_normal_sampled_slope_sampled_intercept_sampled_sigma(predicted_cs, experimental_cs);
               energy += get_beta_prime_prior_sigma();
               return energy;
          }
          if (settings.energy_type == "marginalized-intercept") {
               double energy = get_normal_sampled_slope_marg_intercept_sampled_sigma(predicted_cs, experimental_cs);
               energy += get_beta_prime_prior_sigma();
               return energy;
          }
          if (settings.energy_type == "sample-weights") {
               double energy = get_normal_marg_slope_marg_intercept_sampled_sigma(predicted_cs, experimental_cs);
               energy += get_beta_prime_prior_sigma();
               return energy;
          }
          return 0.0;
     }

     //! Update single weight.
     void sample_parameter(FPtype &parameter, FPtype stdev) {

          FPtype delta = rnom(0.0,stdev, this->random_number_engine);

          parameter *= std::exp(delta);
     }

     //! Update single weight.
     void sample_parameter_add(FPtype &parameter, FPtype stdev) {

          FPtype delta = rnom(0.0,stdev, this->random_number_engine);

          parameter += delta;
     }

     void update_parameter() {
          if (settings.energy_type == "sample-all") {
               // Get random parameter
               unsigned int i = rand_int(0,2, this->random_number_engine);
               // Get random atom type
               unsigned int j = rand_int(0,5, this->random_number_engine);
               if (i == 0) this->sample_parameter(this->weights[j], 0.30);
               else if (i == 1) this->sample_parameter(this->slopes[j], 0.30);
               else if (i == 2) this->sample_parameter_add(this->intercepts[j], 0.50);
          }
          else if (settings.energy_type == "marginalized-intercept") {
               // Get random parameter
               unsigned int i = rand_int(0,1, this->random_number_engine);
               // Get random atom type
               unsigned int j = rand_int(0,5, this->random_number_engine);
               if (i == 0) this->sample_parameter(this->weights[j], 0.30);
               else if (i == 1) this->sample_parameter(this->slopes[j], 0.30);
          }
          else if (settings.energy_type == "sample-weights") {
               // Get random atom type
               unsigned int j = rand_int(0,5, this->random_number_engine);
               this->sample_parameter(this->weights[j], 0.30);
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
          if (!this->none_move) {
               this->predicted_chemical_shifts = predict(*(this->chain), move_info);
          }

          return get_energy(this->predicted_chemical_shifts,this->experimental_chemical_shifts);

     }


     //! Print info on sampled parameters.
     void print_parameters() {
          if (this->n_print > 50000) {
               double cs1;
               double cs2;
               double chi_sq_first;
               double chi_sq_second;
               FPlist rmsd_first = empty_contribution;
               FPlist slope_first = empty_contribution;
               FPlist rmsd_second = empty_contribution;
               // Holds which residues that makes a contribution to the chemical shift energy.
               std::vector<unsigned int> contributions(this->chain->size(),0);
               for (unsigned int j=0; j < 6; j++) {
                    unsigned int n_first = 0;
                    double xsum_first = 0.0;
                    double ysum_first = 0.0;
                    double xysum_first = 0.0;
                    double xsum2_first = 0.0;
                    double ysum2_first = 0.0;
                    unsigned int n_second = 0;
                    double xsum_second = 0.0;
                    double ysum_second = 0.0;
                    double xysum_second = 0.0;
                    double xsum2_second = 0.0;
                    double ysum2_second = 0.0;
                    for (unsigned int i=1;i < this->chain->size()-1;i++) {
                        cs1 = this->predicted_chemical_shifts[i][j];
                        cs2 = this->experimental_chemical_shifts[i][j];
                        if ((std::abs(cs1) > 0.000001) &&
                            (std::abs(cs2) > 0.000001)) {
                            contributions[i] = 1;
                            if (this->member[i]) {
                                 xsum_first += cs1;
                                 ysum_first += cs2;
                                 xysum_first += cs1*cs2;
                                 xsum2_first += cs1*cs1;
                                 ysum2_first += cs2*cs2;
                                 n_first++;
                            } else {
                                 xsum_second += cs1;
                                 ysum_second += cs2;
                                 xysum_second += cs1*cs2;
                                 xsum2_second += cs1*cs1;
                                 ysum2_second += cs2*cs2;
                                 n_second++;
                            }
                        }
                    }
                    // Use MLE expression
                    if (n_first > 2) {
                         chi_sq_first = (n_first * ysum2_first - ysum_first * ysum_first - this->slopes[j] * (n_first * xysum_first - xsum_first * ysum_first)) / n_first;
                         //chi_sq_first = ysum2_first + ( xysum_first * (2*xsum_first * ysum_first - n_first* xysum_first) - xsum2_first * ysum_first * ysum_first ) / (n_first * xsum2_first - xsum_first*xsum_first);
                         rmsd_first[j] = std::sqrt(chi_sq_first/(n_first-2));
                         //slope_first[j] = (n_first * xysum_first - xsum_first*ysum_first)/
                         //                 (n_first*xsum2_first - xsum_first*xsum_first);
                    }
                    if (n_second > 0) {
                         chi_sq_second = ysum2_second + ( xysum_second * (2*xsum_second * ysum_second - n_second* xysum_second) - xsum2_second * ysum_second * ysum_second ) / (n_second * xsum2_second - xsum_second*xsum_second);
                         rmsd_second[j] = std::sqrt(chi_sq_second/(n_second-2));
                    }
               }
               FPlist weights1;
               FPlist weights2;
               for (unsigned int i=0; i < 6; i++) {
                    weights1.push_back(this->weights[i]);
                    weights2.push_back(std::sqrt(this->weights[i+6]*this->weights[i+6]+this->weights[i]*this->weights[i]));
               }
               if (settings.energy_type == "sample-all") {

                    std::cout << std::setprecision(3)
                              << "\n# PROCS15 -- THREAD #" << this->thread_id
                              << " -- WEIGHTS: " << weights1
                              << " -- RMSD: " << rmsd_first
                              << std::endl;
                    std::cout << std::setprecision(3)
                              << "\t\t -- SLOPES: " << this->slopes
                              << " -- INTERCEPTS: " << this->intercepts
                              << std::endl;
               }
               else if (settings.energy_type == "marginalized-intercept") {

                    std::cout << std::setprecision(3)
                              << "\n# PROCS15 -- THREAD #" << this->thread_id
                              << " -- WEIGHTS: " << weights1
                              << " -- RMSD: " << rmsd_first
                              << std::endl;
                    std::cout << std::setprecision(3)
                              << "\t\t -- SLOPES: " << this->slopes
                              << std::endl;
               }
               else if (settings.energy_type == "sample-weights") {

                    std::cout << std::setprecision(3)
                              << "\n# PROCS15 -- THREAD #" << this->thread_id
                              << " -- WEIGHTS: " << weights1
                              << " -- RMSD: " << rmsd_first
                              << "\t\t -- SLOPES: " << this->slopes
                              << std::endl;
               }

               else {
                    std::cout << std::setprecision(3)
                              << "\n# PROCS15 -- THREAD #" << this->thread_id
                              << " -- RMSD: " << rmsd_first
                              << std::endl;
               }

               this->n_print = 0;
          } else {
               this->n_print++;
          }
     }

     void write_nmr_star_format(phaistos::ChainFB& chain, std::string filename, FPtable &predicted_cs){

          int atom_count = 1;

          std::ofstream myfile;
          myfile.open(filename.c_str());
          myfile << std::fixed;
          myfile << std::setprecision(2);

          for (unsigned int i = 0; i < predicted_cs.size(); i++){
               for (unsigned int j = 0; j < 6; j++){
                    switch (j){
                         case 1: //!CA
                              if (loaded_tables[0] == true){
                                   myfile << atom_count << "\t" << (i+1) << "\t" << chain[i].residue_type << "\t" << "CA" << "\t" << "C" << "\t" << predicted_cs[i][1] << "\n";
                                   atom_count+= 1;
                              }

                              break;
                         case 5: //!CB
                              if (chain[i].residue_type != GLY && loaded_tables[1] == true){
                                   myfile << atom_count << "\t" << (i+1) << "\t" << chain[i].residue_type << "\t" << "CB" << "\t" << "C" << "\t" << predicted_cs[i][5] << "\n";
                                   atom_count+= 1;
                              }
                              break;
                         case 4: //!C
                              if (loaded_tables[2] == true){
                                   myfile << atom_count << "\t" << (i+1) << "\t" << chain[i].residue_type << "\t" << "C" << "\t" << "C" << "\t" << predicted_cs[i][4] << "\n";
                                   atom_count+= 1;
                              }
                              break;
                         case 3: //!N
                              if (loaded_tables[3] == true){
                                   myfile << atom_count << "\t" << (i+1) << "\t" << chain[i].residue_type << "\t" << "N" << "\t" << "N" << "\t" << predicted_cs[i][3] << "\n";
                                   atom_count+= 1;
                              }
                              break;
                         case 2: //!HN
                              if (chain[i].residue_type != PRO && loaded_tables[4] == true){
                                   myfile << atom_count << "\t" << (i+1) << "\t" << chain[i].residue_type << "\t" << "H" << "\t" << "H" << "\t" << predicted_cs[i][2] << "\n";
                                   atom_count+= 1;
                              }
                              break;
                         case 0: //!HA
                              if (loaded_tables[5] == true){
                                   if (chain[i].residue_type == GLY) {
                                        myfile << atom_count << "\t" << (i+1) << "\t" << chain[i].residue_type << "\t" << "HA2" << "\t" << "H" << "\t" << predicted_cs[i][0] << "\n";
                                   } else {
                                        myfile << atom_count << "\t" << (i+1) << "\t" << chain[i].residue_type << "\t" << "HA" << "\t" << "H" << "\t" << predicted_cs[i][0] << "\n";
                                   }

                                   atom_count++;
                              }
                              break;
                    }
               }
          }
          myfile.close();
     }

     void write_nmr_star_format_full(phaistos::ChainFB& chain, std::string filename, std::vector<FPtable> &contributions){

          int atom_count = 1;

          std::ofstream myfile;
          myfile.open(filename.c_str());
          myfile << std::fixed;
          myfile << std::setprecision(2);
          myfile << "\t\t\t\t\tSUM \tHB1 \tHaB1\tHB2 \tHaB2\tRC  \tpre \tshld\tfol \twat\n";

          for (unsigned int i = 0; i < contributions.size(); i++){
               for (unsigned int j = 0; j < 6; j++){
                    switch (j){
                         case 1: //!CA
                              if (loaded_tables[0] == true){
                                   myfile << atom_count << "\t" << (i+1) << "\t" << chain[i].residue_type << "\t" << "CA" << "\t" << "C" << "\t" << contributions[i].back()[j];
                                   for (unsigned int k = 0; k < contributions[i].size()-1; k++) {
                                        myfile << "\t" << contributions[i][k][j];
                                   }
                                   myfile << "\n";
                                   atom_count+= 1;
                              }

                              break;
                         case 5: //!CB
                              if (chain[i].residue_type != GLY && loaded_tables[1] == true){
                                   myfile << atom_count << "\t" << (i+1) << "\t" << chain[i].residue_type << "\t" << "CB" << "\t" << "C" << "\t" << contributions[i].back()[j];
                                   for (unsigned int k = 0; k < contributions[i].size()-1; k++) {
                                        myfile << "\t" << contributions[i][k][j];
                                   }
                                   myfile << "\n";
                                   atom_count+= 1;
                              }
                              break;
                         case 4: //!C
                              if (loaded_tables[2] == true){
                                   myfile << atom_count << "\t" << (i+1) << "\t" << chain[i].residue_type << "\t" << "C" << "\t" << "C" << "\t" << contributions[i].back()[j];
                                   for (unsigned int k = 0; k < contributions[i].size()-1; k++) {
                                        myfile << "\t" << contributions[i][k][j];
                                   }
                                   atom_count+= 1;
                                   myfile << "\n";
                              }
                              break;
                         case 3: //!N
                              if (loaded_tables[3] == true){
                                   myfile << atom_count << "\t" << (i+1) << "\t" << chain[i].residue_type << "\t" << "N" << "\t" << "N" << "\t" << contributions[i].back()[j];
                                   for (unsigned int k = 0; k < contributions[i].size()-1; k++) {
                                        myfile << "\t" << contributions[i][k][j];
                                   }
                                   atom_count+= 1;
                                   myfile << "\n";
                              }
                              break;
                         case 2: //!HN
                              if (chain[i].residue_type != PRO && loaded_tables[4] == true){
                                   myfile << atom_count << "\t" << (i+1) << "\t" << chain[i].residue_type << "\t" << "H" << "\t" << "H" << "\t" << contributions[i].back()[j];
                                   for (unsigned int k = 0; k < contributions[i].size()-1; k++) {
                                        myfile << "\t" << contributions[i][k][j];
                                   }
                                   atom_count+= 1;
                                   myfile << "\n";
                              }
                              break;
                         case 0: //!HA
                              if (loaded_tables[5] == true){
                                   if (chain[i].residue_type == GLY) {
                                        myfile << atom_count << "\t" << (i+1) << "\t" << chain[i].residue_type << "\t" << "HA2" << "\t" << "H" << "\t" << contributions[i].back()[j];
                                   } else {
                                        myfile << atom_count << "\t" << (i+1) << "\t" << chain[i].residue_type << "\t" << "HA" << "\t" << "H" << "\t" << contributions[i].back()[j];
                                   }
                                   for (unsigned int k = 0; k < contributions[i].size()-1; k++) {
                                        myfile << "\t" << contributions[i][k][j];
                                   }
                                   myfile << "\n";

                                   atom_count++;
                              }
                              break;
                    }
               }
          }
          myfile.close();
     }

     //! Accept last energy evaluation. Save predictions and model parameters.
     void accept() {

          // If the move was a none-move, energy term parameters
          if (this->none_move) {

               this->weights_previous = this->weights;
               this->intercepts_previous = this->intercepts;
               this->slopes_previous = this->slopes;
               this->member_previous = this->member;

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
               this->member = this->member_previous;

          } else {

               this->predicted_chemical_shifts = this->predicted_chemical_shifts_previous;
          }

          print_parameters();

     }

}; // End class TermProCS15Base



} // end namespace phaistos
#endif
