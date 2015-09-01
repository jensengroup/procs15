namespace module_procs15 {

// Module: energy term initialization
template <typename SETTINGS_MODIFIER>
struct EnergyOptions {

     // Constructor - general case: do nothing
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   CHAIN_TYPE *chain,
                   DBN_TYPE *dbn) {
     }

     // Constructor - ChainFB specific case
     template <typename DBN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   ChainFB *chain,
                   DBN_TYPE *dbn) {

          // Import namespace for make_vector
          using namespace boost::fusion;

          // Define defaults for the different modes
          ModeDefinitions mode_definitions(target, chain);


          // ProCS15 energy term
          for (int counter = occurrences[prefix+"-procs15"]; counter > 0; counter--) {

               typedef TermProCS15 EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));
               //boost::shared_ptr<Settings> settings(new Settings());

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "ProCS15 chemical shift energy (" + prefix + ")",
                         prefix+"-procs15", settings,
                         make_vector(
                              make_vector(std::string("star-filename"),
                                          std::string("Chemical shifts file"),
                                          &settings->star_filename),
                              make_vector(std::string("energy-type"),
                                          std::string("Which probability model to use in the energy term (marginalized (default), mle, sample-all, marginalized-intercept, sample-weights, mixture"),
                                          &settings->energy_type),
                              make_vector(std::string("load-CA-Chemical-shifts"),
                          std::string("Whether to load and do predictions for alpha carbon chemical shifts (=true)"),
                                          &settings->load_ca),
                              make_vector(std::string("load-CB-Chemical-shifts"),
                          std::string("Whether to load and do predictions for beta carbon chemical shifts (=true)"),
                                          &settings->load_cb),
                              make_vector(std::string("load-CO-Chemical-shifts"),
                          std::string("Whether to load and do predictions for carbonyl carbon chemical shifts (=true)"),
                                          &settings->load_co),
                              make_vector(std::string("load-N-Chemical-shifts"),
                          std::string("Whether to load and do predictions for amide chemical shifts (=true)"),
                                          &settings->load_n),
                              make_vector(std::string("load-HN-Chemical-shifts"),
                          std::string("Whether to load and do predictions for amide proton chemical shifts (=true)"),
                                          &settings->load_hn),
                              make_vector(std::string("load-HA-Chemical-shifts"),
                          std::string("Whether to load and do predictions for alpha hydrogen chemical shifts (=true)"),
                                          &settings->load_ha),
                              make_vector(std::string("include-ha-rc"),
                          std::string(""),
                                          &settings->include_ha_rc),
                              make_vector(std::string("include-hn-rc"),
                          std::string(""),
                                          &settings->include_hn_rc),
                              make_vector(std::string("include-hn-previous"),
                          std::string(""),
                                          &settings->include_hn_previous_residue_correction),
                              make_vector(std::string("include-ha-previous"),
                          std::string(""),
                                          &settings->include_ha_previous_residue_correction),
                              make_vector(std::string("include-ca-previous"),
                          std::string(""),
                                          &settings->include_ca_previous_residue_correction),
                              make_vector(std::string("include-cb-previous"),
                          std::string(""),
                                          &settings->include_cb_previous_residue_correction),
                              make_vector(std::string("include-co-previous"),
                          std::string(""),
                                          &settings->include_co_previous_residue_correction),
                          //    make_vector(std::string("include-n-previous"),
                          //std::string(""),
                          //                &settings->include_n_previous_residue_correction),
                              make_vector(std::string("include-hn-following"),
                          std::string(""),
                                          &settings->include_hn_following_residue_correction),
                              make_vector(std::string("include-ha-following"),
                          std::string(""),
                                          &settings->include_ha_following_residue_correction),
                              make_vector(std::string("include-ca-following"),
                          std::string(""),
                                          &settings->include_ca_following_residue_correction),
                              make_vector(std::string("include-cb-following"),
                          std::string(""),
                                          &settings->include_cb_following_residue_correction),
                              make_vector(std::string("include-co-following"),
                          std::string(""),
                                          &settings->include_co_following_residue_correction),
                              make_vector(std::string("include-n-following"),
                          std::string(""),
                                          &settings->include_n_following_residue_correction),
                          //    make_vector(std::string("include-hn-primary-hn-hbond"),
                          //std::string(""),
                          //                &settings->include_hn_primary_hn_hbond),
                              make_vector(std::string("include-ha-primary-hn-hbond"),
                          std::string(""),
                                          &settings->include_ha_primary_hn_hbond),
                              make_vector(std::string("include-ca-primary-hn-hbond"),
                          std::string(""),
                                          &settings->include_ca_primary_hn_hbond),
                              make_vector(std::string("include-co-primary-hn-hbond"),
                          std::string(""),
                                          &settings->include_co_primary_hn_hbond),
                              make_vector(std::string("include-n-primary-hn-hbond"),
                          std::string(""),
                                          &settings->include_n_primary_hn_hbond),
                              make_vector(std::string("include-hn-primary-ha-hbond"),
                          std::string(""),
                                          &settings->include_hn_primary_ha_hbond),
                          //    make_vector(std::string("include-ha-primary-ha-hbond"),
                          //std::string(""),
                          //                &settings->include_ha_primary_ha_hbond),
                          //    make_vector(std::string("include-ca-primary-ha-hbond"),
                          //std::string(""),
                          //                &settings->include_ca_primary_ha_hbond),
                          //    make_vector(std::string("include-cb-primary-ha-hbond"),
                          //std::string(""),
                          //                &settings->include_cb_primary_ha_hbond),
                              make_vector(std::string("include-co-primary-ha-hbond"),
                          std::string(""),
                                          &settings->include_co_primary_ha_hbond),
                          //    make_vector(std::string("include-n-primary-ha-hbond"),
                          //std::string(""),
                          //                &settings->include_n_primary_ha_hbond),
                              make_vector(std::string("include-hn-secondary-hn-hbond"),
                          std::string(""),
                                          &settings->include_hn_secondary_hn_hbond),
                              make_vector(std::string("include-ha-secondary-hn-hbond"),
                          std::string(""),
                                          &settings->include_ha_secondary_hn_hbond),
                              make_vector(std::string("include-ca-secondary-hn-hbond"),
                          std::string(""),
                                          &settings->include_ca_secondary_hn_hbond),
                          //    make_vector(std::string("include-co-secondary-hn-hbond"),
                          //std::string(""),
                          //                &settings->include_co_secondary_hn_hbond),
                          //    make_vector(std::string("include-n-secondary-hn-hbond"),
                          //std::string(""),
                          //                &settings->include_n_secondary_hn_hbond),
                              make_vector(std::string("include-hn-secondary-ha-hbond"),
                          std::string(""),
                                          &settings->include_hn_secondary_ha_hbond),
                              make_vector(std::string("include-ha-secondary-ha-hbond"),
                          std::string(""),
                                          &settings->include_ha_secondary_ha_hbond),
                              make_vector(std::string("include-ca-secondary-ha-hbond"),
                          std::string(""),
                                          &settings->include_ca_secondary_ha_hbond),
                              make_vector(std::string("include-co-secondary-ha-hbond"),
                          std::string(""),
                                          &settings->include_co_secondary_ha_hbond),
                              make_vector(std::string("include-n-secondary-ha-hbond"),
                          std::string(""),
                                          &settings->include_n_secondary_ha_hbond),
                              make_vector(std::string("output-full-prediction-vector"),
                          std::string("(Observable only) Output all the differences between predicted and experimental chemical shifts (=false)"),
                                          &settings->output_full_prediction_vector),
                              make_vector(std::string("HN-water-bond-correction"),
                                          std::string("Whether to use the correction for surface HN bonding with water (=true)"),
                                          &settings->use_water_correction),
                              make_vector(std::string("data-folder"),
                                          std::string("Path to folder containing procs15 QM data files"),
                                          &settings->procsnumpypath)
                              )),
                    super_group, counter==1);
          }

          // ProCS15-cached energy term
          for (int counter = occurrences[prefix+"-procs15-cached"]; counter > 0; counter--) {

               typedef TermProCS15Cached EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "ProCS15 (cached) chemical shift energy (" + prefix + ")",
                         prefix+"-procs15-cached", settings,
                         make_vector(
                              make_vector(std::string("star-filename"),
                                          std::string("Chemical shifts file"),
                                          &settings->star_filename),
                              make_vector(std::string("energy-type"),
                                          std::string("Which probability model to use in the energy term (marginalized (default), mle, sample-all, marginalized-intercept, sample-weights, mixture"),
                                          &settings->energy_type),
                              make_vector(std::string("load-CA-Chemical-shifts"),
                          std::string("Whether to load and do predictions for alpha carbon chemical shifts (=true)"),
                                          &settings->load_ca),
                              make_vector(std::string("load-CB-Chemical-shifts"),
                          std::string("Whether to load and do predictions for beta carbon chemical shifts (=true)"),
                                          &settings->load_cb),
                              make_vector(std::string("load-CO-Chemical-shifts"),
                          std::string("Whether to load and do predictions for carbonyl carbon chemical shifts (=true)"),
                                          &settings->load_co),
                              make_vector(std::string("load-N-Chemical-shifts"),
                          std::string("Whether to load and do predictions for amide chemical shifts (=true)"),
                                          &settings->load_n),
                              make_vector(std::string("load-HN-Chemical-shifts"),
                          std::string("Whether to load and do predictions for amide proton chemical shifts (=true)"),
                                          &settings->load_hn),
                              make_vector(std::string("load-HA-Chemical-shifts"),
                          std::string("Whether to load and do predictions for alpha hydrogen chemical shifts (=true)"),
                                          &settings->load_ha),
                              make_vector(std::string("data-folder"),
                                          std::string("Path to folder containing procs15 QM data files"),
                                          &settings->procsnumpypath)
                              )),
                    super_group, counter==1);
          }

    }
};

}
