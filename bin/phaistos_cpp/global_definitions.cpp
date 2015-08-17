namespace module_procs15 {

// Module: energy term initialization
struct EnergyInitialization {


     // Constructor - general case: do nothing
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     EnergyInitialization(const Options &options, CHAIN_TYPE *chain, DBN_TYPE *dbn,
                                       Energy<CHAIN_TYPE> *energy, std::vector<RandomNumberEngine *> *random_number_generators, 
                                       std::string prefix="") {
     }

     // Constructor - template specific case
     template <typename DBN_TYPE>
     EnergyInitialization(const Options &options, ChainFB *chain, DBN_TYPE *dbn,
                                       Energy<ChainFB> *energy, std::vector<RandomNumberEngine *> *random_number_generators, 
                                       std::string prefix="") {

          Options::OptionValue option;

          // ProCS15 energy term
          option = options[prefix+"-procs15"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef typename TermProCS15::Settings Settings;

               // Add energy term
               energy->add_term(new TermProCS15(chain,
                                                 options.get_settings<Settings>(option,i)));
          }

          // ProCS15-cached energy term
          option = options[prefix+"-procs15-cached"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef typename TermProCS15Cached::Settings Settings;

               // Add energy term
               energy->add_term(new TermProCS15Cached(chain,
                                                 options.get_settings<Settings>(option,i)));
          }


     }

};


}
