// Copyright (c) 2013-2015 by Anders Steen Christensen, Lars Andersen Bratholm
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef LIB_DEFINITIONS
#define LIB_DEFINITIONS

#include <vector>
#include <boost/assign/list_of.hpp>

namespace procs15 {


     // Define single precision or double precision,
     // may have impact on speed.
     typedef double FPtype;
     typedef std::vector<FPtype> FPlist;
     typedef std::vector<FPlist> FPtable;

     const FPtype rad_to_degrees = 57.2957795131;

     enum RingEnum {RingBenzene = 0,
                    RingPhenol,
                    RingImidazolium,
                    RingIndoleSmall,
                    RingIndoleLarge,
                    RING_ENUM_SIZE};

     const std::vector<FPtype> ring_current_intensities = boost::assign::list_of(1.00)  // Benzene
                                                                                (0.81)  // Phenol
                                                                                (0.69)  // Imidazolium
                                                                                (0.57)  // IndoleSmall
                                                                                (1.02); // IndoleLarge
     const FPtype b_factor = 30.42;

     const FPtype hn_bond_cutoff = 3.0;

     const FPtype ha_bond_cutoff = 4.0;

     const FPtype rc_cutoff2 = 64.0;

     const FPtype water_bonding_correction = 2.07;

     const unsigned int gridsize4d = 5;
     const unsigned int gridsize5d = 20;
     const unsigned int gridsize6d = 20;
     const unsigned int gridsize7d = 20;


     FPlist empty_contribution(6, 0.0);
};

#endif
