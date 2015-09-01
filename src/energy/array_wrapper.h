//  Copyright (C) 2014-2015 by Anders S. Larsen
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


#ifndef ARRAY_WRAPPER
#define ARRAY_WRAPPER


#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>

#include <boost/lexical_cast.hpp>
#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include "cnpy_phaistos.h"

#include "lib_definitions.h"

namespace procs15 {

class HbondDataWrapper { //! wraps data for hydrogen bonding
     public: 
          boost::multi_array<float, 4> fourd_array;

          HbondDataWrapper(const boost::multi_array<float, 4>&  array_to_be_wrapped) {
               fourd_array.resize(boost::extents[array_to_be_wrapped.shape()[0]]
                                               [array_to_be_wrapped.shape()[1]]
                                               [array_to_be_wrapped.shape()[2]]
                                               [array_to_be_wrapped.shape()[3]]);
               fourd_array = array_to_be_wrapped;
          }
};



class ArrayWrapper { 
     public: 
          boost::multi_array<float, 1> oned_array;
          boost::multi_array<float, 2> twod_array;
          boost::multi_array<float, 3> threed_array;
          boost::multi_array<float, 4> fourd_array;

          virtual FPtype get_values(unsigned int a)=0;

          virtual FPtype get_values(unsigned int a,
                                   unsigned int b)=0;

          virtual FPtype get_values(unsigned int a,
                                   unsigned int b,
	         		               unsigned int c)=0; 


          virtual FPtype get_values(unsigned int a,
                                   unsigned int b,
                                   unsigned int c,
                                   unsigned int d)=0;

          virtual FPtype get_st(int i)=0;

          virtual void set_tension(FPtype newtension)=0;

          virtual void set_bias(FPtype newbias)=0;

          unsigned int numpy(int a){ //convert 0:360 to -180:180 format of numpy file for use in array
	        if (a>= 180){
      	            a = a - 180;
                }
                else{
                    a = a + 180;
                }
                return (unsigned int) a;
          } 

          virtual FPtype get_values(const int i, const std::vector<unsigned int>& d, const std::vector<unsigned int>& chi )=0; 
};


class ArrayWrapperStandard:public ArrayWrapper { 
     public: 
          boost::multi_array<float, 1> oned_array;
          boost::multi_array<float, 2> twod_array;
          boost::multi_array<float, 3> threed_array;
          boost::multi_array<float, 4> fourd_array;



          FPtype get_values(unsigned int a) {
               return (FPtype)oned_array[a];
          }

          FPtype get_values(unsigned int a,
                            unsigned int b) {
               return (FPtype)twod_array[a][b];
          }

          FPtype get_values(unsigned int a,
                            unsigned int b,
                            unsigned int c) {
               return (FPtype)threed_array[a][b][c];
          }

          FPtype get_values(unsigned int a,
                            unsigned int b,
                            unsigned int c,
                            unsigned int d) {
               return (FPtype)fourd_array[a][b][c][d];
          } 


          void set_tension(FPtype newtension){

          }

          void set_bias(FPtype newbias){

          }


          ArrayWrapperStandard(const boost::multi_array<float, 1>&  array_to_be_wrapped) {
               oned_array.resize(boost::extents[array_to_be_wrapped.shape()[0]]);
               oned_array = array_to_be_wrapped;
          }

          ArrayWrapperStandard(const boost::multi_array<float, 2>& array_to_be_wrapped) {
               twod_array.resize(boost::extents[array_to_be_wrapped.shape()[0]]
                                              [array_to_be_wrapped.shape()[1]]);
               twod_array = array_to_be_wrapped;
          }

          ArrayWrapperStandard(const boost::multi_array<float, 3>&  array_to_be_wrapped) {
               threed_array.resize(boost::extents[array_to_be_wrapped.shape()[0]]
                                                [array_to_be_wrapped.shape()[1]]
                                                [array_to_be_wrapped.shape()[2]]);
               threed_array = array_to_be_wrapped;
          }

          ArrayWrapperStandard(const boost::multi_array<float, 4>&  array_to_be_wrapped) {
               fourd_array.resize(boost::extents[array_to_be_wrapped.shape()[0]]
                                               [array_to_be_wrapped.shape()[1]]
                                               [array_to_be_wrapped.shape()[2]]
                                               [array_to_be_wrapped.shape()[3]]);
              fourd_array = array_to_be_wrapped;
          }
          FPtype get_values(const int i, const std::vector<unsigned int>& d, const std::vector<unsigned int>& chi ){ // i = residue index, d = dihidral , chi = chi      
               switch(chi.size()){
                    case 0:
                    return (FPtype)threed_array[i][numpy(d[0])][numpy(d[1])];
                    break;
               case 1:        
                    return (FPtype)fourd_array[i][numpy(chi[0])][numpy(d[0])][numpy(d[1])];
                    break;
               default:
                    std::cerr << "ERROR wrong arraysize" << std::endl;
                    return 0.0;
                    break;
               }
          }

          FPtype get_st(int i){
               return (FPtype)threed_array[i][(numpy(+360-120))][numpy(140)]; 
          }
          ~ArrayWrapperStandard(){
               std::cerr << "Destructor called unexpectedly" << std::endl; 
          }
};




class ArrayWrapperInterpolate:public ArrayWrapper{
     private:
          unsigned int step;
     public:
          boost::multi_array<float, 1> oned_array;
          boost::multi_array<float, 2> twod_array;
          boost::multi_array<float, 3> threed_array;
          boost::multi_array<float, 4> fourd_array;
          boost::multi_array<float, 5> fived_array;
          boost::multi_array<float, 6> sixd_array;
          boost::multi_array<float, 7> sevend_array;


          ArrayWrapperInterpolate(const boost::multi_array<float, 1>&  array_to_be_wrapped, unsigned int g) {
               oned_array.resize(boost::extents[array_to_be_wrapped.shape()[0]]);
               oned_array = array_to_be_wrapped;
	       step = g;
          }

          ArrayWrapperInterpolate(const boost::multi_array<float, 2>& array_to_be_wrapped, unsigned int g) {
               twod_array.resize(boost::extents[array_to_be_wrapped.shape()[0]]
                                               [array_to_be_wrapped.shape()[1]]);
               twod_array = array_to_be_wrapped;
	       step = g;
          }

          ArrayWrapperInterpolate(const boost::multi_array<float, 3>&  array_to_be_wrapped, unsigned int g) {
               threed_array.resize(boost::extents[array_to_be_wrapped.shape()[0]]
                                                [array_to_be_wrapped.shape()[1]]
                                                [array_to_be_wrapped.shape()[2]]);
               threed_array = array_to_be_wrapped;
 	       step = g;
          }

          ArrayWrapperInterpolate(const boost::multi_array<float, 4>&  array_to_be_wrapped, unsigned int g) {
               fourd_array.resize(boost::extents[array_to_be_wrapped.shape()[0]]
                                               [array_to_be_wrapped.shape()[1]]
                                               [array_to_be_wrapped.shape()[2]]
                                               [array_to_be_wrapped.shape()[3]]);
               fourd_array = array_to_be_wrapped;
	       step = g;
          }
          ArrayWrapperInterpolate(const boost::multi_array<float, 5>&  array_to_be_wrapped, unsigned int g) {
               fived_array.resize(boost::extents[array_to_be_wrapped.shape()[0]]
                                               [array_to_be_wrapped.shape()[1]]
                                               [array_to_be_wrapped.shape()[2]]
                                               [array_to_be_wrapped.shape()[3]]
	         			      [array_to_be_wrapped.shape()[4]]);
               fived_array = array_to_be_wrapped;
	       step = g;
          }
          ArrayWrapperInterpolate(const boost::multi_array<float, 6>&  array_to_be_wrapped, unsigned int g) {
               sixd_array.resize(boost::extents[array_to_be_wrapped.shape()[0]]
                                               [array_to_be_wrapped.shape()[1]]
                                               [array_to_be_wrapped.shape()[2]]
                                               [array_to_be_wrapped.shape()[3]]
                                               [array_to_be_wrapped.shape()[4]]
                                               [array_to_be_wrapped.shape()[5]]);

               sixd_array = array_to_be_wrapped;
	       step = g;
          }
          ArrayWrapperInterpolate(const boost::multi_array<float, 7>&  array_to_be_wrapped, unsigned int g) {
               sevend_array.resize(boost::extents[array_to_be_wrapped.shape()[0]]
                                               [array_to_be_wrapped.shape()[1]]
                                               [array_to_be_wrapped.shape()[2]]
                                               [array_to_be_wrapped.shape()[3]]
                                               [array_to_be_wrapped.shape()[4]]
                                               [array_to_be_wrapped.shape()[5]]
                                               [array_to_be_wrapped.shape()[6]]);

               sevend_array = array_to_be_wrapped;
	       step = g;
          }


          void set_tension(FPtype newtension){

          }

          void set_bias(FPtype newbias){

          }

          FPtype get_values(unsigned int a) {
               return (FPtype)oned_array[a];
          }

          FPtype get_values(unsigned int a,
                            unsigned int b) {
               return (FPtype)twod_array[a][b];
          }

          FPtype get_values(unsigned int a,
                            unsigned int b,
                            unsigned int c) {
               return (FPtype)threed_array[a][b][c];
          }

          FPtype get_values(unsigned int a,
                            unsigned int b,
                            unsigned int c,
                            unsigned int d) {
               return (FPtype)fourd_array[a][b][c][d];
          }
          FPtype get_values(unsigned int a,
                            unsigned int b,
                            unsigned int c,
                            unsigned int d,
	         	            unsigned int e) {
               return (FPtype)fived_array[a][b][c][d][e];
          }
          FPtype get_values(unsigned int a,
                            unsigned int b,
                            unsigned int c,
                            unsigned int d,
                            unsigned int e,
	         	            unsigned int f) {
              return (FPtype)sixd_array[a][b][c][d][e][f];
          }
          FPtype get_values(unsigned int a,
                            unsigned int b,
                            unsigned int c,
                            unsigned int d,
                            unsigned int e,
                            unsigned int f,
 	         	            unsigned int g){
               return (FPtype)sevend_array[a][b][c][d][e][f][g];
          }

          unsigned int ru(int p){    //get upperbound for use in array index 
               return (unsigned int)((numpy(p)+step)/step);
          }
          unsigned int  rl(int p){
               return (unsigned int)(numpy(p)/step); //get lowerbound for use in array index 
          }
	  FPtype upperbound(int p){
	       return (FPtype)((numpy(p)+step)/step)*step; //get upperbbound
	  }
          FPtype lowerbound(int p){
               return (FPtype)(numpy(p)/step)*step; //get lowerbound
          }
	  unsigned int nearest(int a){
	       if ( (numpy(a)-lowerbound(a))>(step/2) ){
	            a = upperbound(a)/step; 
	       }
	       else{
	            a = lowerbound(a)/step;
               }
	       return  a;
          }
          FPtype get_st(int i){         
               return (FPtype)threed_array[i][3][16];
	  }
	  FPtype interpolate4d(int i, FPtype x, FPtype y, FPtype z, FPtype t, int phi, int psi, int grid[][2]){
               FPtype c0; FPtype c1; FPtype c00; FPtype c10; FPtype c01; FPtype c11;
               c00 = (sevend_array[i][phi][psi][grid[0][0]][grid[1][0]][grid[2][0]][grid[3][0]]*(1-x)+sevend_array[i][phi][psi][grid[0][1]][grid[1][0]][grid[2][0]][grid[3][0]]*x)*(1-y)+
                     (sevend_array[i][phi][psi][grid[0][0]][grid[1][1]][grid[2][0]][grid[3][0]]*(1-x)+sevend_array[i][phi][psi][grid[0][1]][grid[1][1]][grid[2][0]][grid[3][0]]*x)*y;


               c10 = (sevend_array[i][phi][psi][grid[0][0]][grid[1][0]][grid[2][1]][grid[3][0]]*(1-x)+sevend_array[i][phi][psi][grid[0][1]][grid[1][0]][grid[2][1]][grid[3][0]]*x)*(1-y)+
                     (sevend_array[i][phi][psi][grid[0][0]][grid[1][1]][grid[2][1]][grid[3][0]]*(1-x)+sevend_array[i][phi][psi][grid[0][1]][grid[1][1]][grid[2][1]][grid[3][0]]*x)*y;

               c01 = (sevend_array[i][phi][psi][grid[0][0]][grid[1][0]][grid[2][0]][grid[3][1]]*(1-x)+sevend_array[i][phi][psi][grid[0][1]][grid[1][0]][grid[2][0]][grid[3][1]]*x)*(1-y)+
                     (sevend_array[i][phi][psi][grid[0][0]][grid[1][1]][grid[2][0]][grid[3][1]]*(1-x)+sevend_array[i][phi][psi][grid[0][1]][grid[1][1]][grid[2][0]][grid[3][1]]*x)*y;


               c11 = (sevend_array[i][phi][psi][grid[0][0]][grid[1][0]][grid[2][1]][grid[3][1]]*(1-x)+sevend_array[i][phi][psi][grid[0][1]][grid[1][0]][grid[2][1]][grid[3][1]]*x)*(1-y)+
                     (sevend_array[i][phi][psi][grid[0][0]][grid[1][1]][grid[2][1]][grid[3][1]]*(1-x)+sevend_array[i][phi][psi][grid[0][1]][grid[1][1]][grid[2][1]][grid[3][1]]*x)*y;

	       c0  = c00*(1-z)+c10*z;	
	       c1  = c01*(1-z)+c11*z;

	       return c0*(1-t)+c1*t; //interpolate t and return


          }
	  FPtype interpolate3d(int i, FPtype x, FPtype y, FPtype z, int phi, int psi, int grid[][2]){
                    FPtype c0; FPtype c1;
		    c0 =
		    ( ( sixd_array[i][phi][psi][grid[0][0]][grid[1][0]][grid[2][0]]*(1-x)+sixd_array[i][phi][psi][grid[0][1]][grid[1][0]][grid[2][0]]*x )*(1-y)+
		    ( sixd_array[i][phi][psi][grid[0][0]][grid[1][1]][grid[2][0]]*(1-x)+sixd_array[i][phi][psi][grid[0][1]][grid[1][1]][grid[2][0]]*x )*y);

		    c1 =
		    ( ( sixd_array[i][phi][psi][grid[0][0]][grid[1][0]][grid[2][1]]*(1-x)+sixd_array[i][phi][psi][grid[0][1]][grid[1][0]][grid[2][1]]*x )*(1-y)+
		    ( sixd_array[i][phi][psi][grid[0][0]][grid[1][1]][grid[2][1]]*(1-x)+sixd_array[i][phi][psi][grid[0][1]][grid[1][1]][grid[2][1]]*x )*y);
		   
		    return c0*(1-z)+c1*z;

          }


          FPtype get_values(const int i, const std::vector<unsigned int>& d, const std::vector<unsigned int>& chi ){ // i = residue index, d = dihidral ,chi = chi

               FPtype x; FPtype y; FPtype z; FPtype t; FPtype c0; FPtype c1; FPtype c00; FPtype c10; FPtype c01; FPtype c11; int grid[4][2];
	       FPtype phix; FPtype psiy;

               switch(chi.size()){
                    case 0:
                         x = ( numpy(d[0])-lowerbound(d[0]) ) / ( upperbound(d[0]) -lowerbound(d[0]) );
                         y = ( numpy(d[1])-lowerbound(d[1]) ) / ( upperbound(d[1]) -lowerbound(d[1]) );

                         return ( (1-x)*(1-y)*threed_array[i][rl(d[0])][rl(d[1])]+
                                         x*(1-y)*threed_array[i][ru(d[0])][rl(d[1])]+
                                         x*y*threed_array[i][ru(d[0])][ru(d[1])]+
                                         (1-x)*y*threed_array[i][rl(d[0])][ru(d[1])]);
                    break;
                    case 1:
                         x = ( numpy(d[0])-lowerbound(d[0]) ) / ( upperbound(d[0]) -lowerbound(d[0]) );
                         y = ( numpy(d[1])-lowerbound(d[1]) ) / ( upperbound(d[1]) -lowerbound(d[1]) );
                         z = ( numpy(chi[0])-lowerbound(chi[0]) ) / ( upperbound(chi[0]) -lowerbound(chi[0]) );

                         c0 =
                         ( ( fourd_array[i][rl(chi[0])][rl(d[0])][rl(d[1])]*(1-x)+fourd_array[i][rl(chi[0])][ru(d[0])][rl(d[1])]*x )*(1-y)+
                         ( fourd_array[i][rl(chi[0])][rl(d[0])][ru(d[1])]*(1-x)+fourd_array[i][rl(chi[0])][ru(d[0])][ru(d[1])]*x )*y);
                         c1 =
                         ( ( fourd_array[i][ru(chi[0])][rl(d[0])][rl(d[1])]*(1-x)+fourd_array[i][ru(chi[0])][ru(d[0])][rl(d[1])]*x )*(1-y)+
                         ( fourd_array[i][ru(chi[0])][rl(d[0])][ru(d[1])]*(1-x)+fourd_array[i][ru(chi[0])][ru(d[0])][ru(d[1])]*x )*y);

                         return c0*(1-z)+c1*z;
                         break;
                    case 2:
                         x = ( numpy(d[0])-lowerbound(d[0]) ) / ( upperbound(d[0])-lowerbound(d[0]) ); //phi
                         y = ( numpy(d[1])-lowerbound(d[1]) ) / ( upperbound(d[1])-lowerbound(d[1]) ); //psi
                         z = ( numpy(chi[0])-lowerbound(chi[0]) ) / ( upperbound(chi[0]) -lowerbound(chi[0]) ); //chi1
                         t = ( numpy(chi[1])-lowerbound(chi[1]) ) / ( upperbound(chi[1]) -lowerbound(chi[1]) ); //chi2

		         grid[0][0] = rl(d[0]); grid[0][1] = ru(d[0]); // x
	                 grid[1][0] = rl(d[1]); grid[1][1] = ru(d[1]); // y
	                 grid[2][0] = rl(chi[0]); grid[2][1] = ru(chi[0]); // z
	                 grid[3][0] = rl(chi[1]); grid[3][1] = ru(chi[1]); // t

                         c00 = (fived_array[i][grid[0][0]][grid[1][0]][grid[2][0]][grid[3][0]]*(1-x)+fived_array[i][grid[0][1]][grid[1][0]][grid[2][0]][grid[3][0]]*x)*(1-y)+
                               (fived_array[i][grid[0][0]][grid[1][1]][grid[2][0]][grid[3][0]]*(1-x)+fived_array[i][grid[0][1]][grid[1][1]][grid[2][0]][grid[3][0]]*x)*y;


                         c10 = (fived_array[i][grid[0][0]][grid[1][0]][grid[2][1]][grid[3][0]]*(1-x)+fived_array[i][grid[0][1]][grid[1][0]][grid[2][1]][grid[3][0]]*x)*(1-y)+
                               (fived_array[i][grid[0][0]][grid[1][1]][grid[2][1]][grid[3][0]]*(1-x)+fived_array[i][grid[0][1]][grid[1][1]][grid[2][1]][grid[3][0]]*x)*y;
 
                         c01 = (fived_array[i][grid[0][0]][grid[1][0]][grid[2][0]][grid[3][1]]*(1-x)+fived_array[i][grid[0][1]][grid[1][0]][grid[2][0]][grid[3][1]]*x)*(1-y)+
                               (fived_array[i][grid[0][0]][grid[1][1]][grid[2][0]][grid[3][1]]*(1-x)+fived_array[i][grid[0][1]][grid[1][1]][grid[2][0]][grid[3][1]]*x)*y;


                         c11 = (fived_array[i][grid[0][0]][grid[1][0]][grid[2][1]][grid[3][1]]*(1-x)+fived_array[i][grid[0][1]][grid[1][0]][grid[2][1]][grid[3][1]]*x)*(1-y)+
                               (fived_array[i][grid[0][0]][grid[1][1]][grid[2][1]][grid[3][1]]*(1-x)+fived_array[i][grid[0][1]][grid[1][1]][grid[2][1]][grid[3][1]]*x)*y;
		         //interpolate z
		         c0  = c00*(1-z)+c10*z;	
		         c1  = c01*(1-z)+c11*z;
                         return c0*(1-t)+c1*t; //interpolate t and return
		         break;
	            case 3: // first interpolate chi angles then phi,psi
                         phix = ( numpy(d[0])-lowerbound(d[0]) ) / ( upperbound(d[0]) -lowerbound(d[0]) ); //phi
                         psiy = ( numpy(d[1])-lowerbound(d[1]) ) / ( upperbound(d[1]) -lowerbound(d[1]) ); //psi
                         x = ( numpy(chi[0])-lowerbound(chi[0]) ) / ( upperbound(chi[0])-lowerbound(chi[0]) ); //chi1
                         y = ( numpy(chi[1])-lowerbound(chi[1]) ) / ( upperbound(chi[1])-lowerbound(chi[1]) ); //chi2
                         z = ( numpy(chi[2])-lowerbound(chi[2]) ) / ( upperbound(chi[2]) -lowerbound(chi[2]) ); //chi3
                         grid[0][0] = rl(chi[0]); grid[0][1] = ru(chi[0]); // x
                         grid[1][0] = rl(chi[1]); grid[1][1] = ru(chi[1]); // y
                         grid[2][0] = rl(chi[2]); grid[2][1] = ru(chi[2]); // z


                         return ( (1-phix)*(1-psiy)*interpolate3d(i,x,y,z,rl(d[0]),rl(d[1]),grid)+
                                         phix*(1-psiy)*interpolate3d(i,x,y,z,ru(d[0]),rl(d[1]),grid)+
                                         phix*psiy*interpolate3d(i,x,y,z,ru(d[0]),ru(d[1]),grid)+
                                        (1-phix)*psiy*interpolate3d(i,x,y,z,rl(d[0]),ru(d[1]),grid) );

		         break;
	
                    case 4: // first interpolate chi angles then phi,psi
                         phix = ( numpy(d[0])-lowerbound(d[0]) ) / ( upperbound(d[0]) -lowerbound(d[0]) ); //phi
                         psiy = ( numpy(d[1])-lowerbound(d[1]) ) / ( upperbound(d[1]) -lowerbound(d[1]) ); //psi
                         x = ( numpy(chi[0])-lowerbound(chi[0]) ) / ( upperbound(chi[0])-lowerbound(chi[0]) ); //chi1
                         y = ( numpy(chi[1])-lowerbound(chi[1]) ) / ( upperbound(chi[1])-lowerbound(chi[1]) ); //chi2
                         z = ( numpy(chi[2])-lowerbound(chi[2]) ) / ( upperbound(chi[2]) -lowerbound(chi[2]) ); //chi3
                         t = ( numpy(chi[3])-lowerbound(chi[3]) ) / ( upperbound(chi[3]) -lowerbound(chi[3]) ); //chi4
                         grid[0][0] = rl(chi[0]); grid[0][1] = ru(chi[0]); // x
                         grid[1][0] = rl(chi[1]); grid[1][1] = ru(chi[1]); // y
                         grid[2][0] = rl(chi[2]); grid[2][1] = ru(chi[2]); // z
                         grid[3][0] = rl(chi[3]); grid[3][1] = ru(chi[3]); // t

                         return ( (1-phix)*(1-psiy)*interpolate4d(i,x,y,z,t,rl(d[0]),rl(d[1]),grid)+
                                   phix*(1-psiy)*interpolate4d(i,x,y,z,t,ru(d[0]),rl(d[1]),grid)+
                                   phix*psiy*interpolate4d(i,x,y,z,t,ru(d[0]),ru(d[1]),grid)+
                                  (1-phix)*psiy*interpolate4d(i,x,y,z,t,rl(d[0]),ru(d[1]),grid) );

                         break;



                    default:
                         std::cerr << "ERROR (ProCS15) wrong arraysize" << std::endl;
                         return 0.0;
                         break;
               }
          }
}; //end class ArrayWrapperInterpolate
}; //end namespace
#endif
