// Original work: Copyright (C) 2011  Carl Rogers
// Modifications: Copyright (C) 2013 Anders S. Christensen
//
// Originally released under MIT License
// License at http://www.opensource.org/licenses/mit-license.php
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


#ifndef LIBCNPY_H_
#define LIBCNPY_H_

#include<algorithm>
#include<vector>
#include<cstdio>
#include<iostream>
#include<cassert>
#include "boost/multi_array.hpp"

namespace cnpy {

class NpyArray {
    public:
        char* raw_data;
        double* loaded_data;
        std::vector<unsigned int> shape;
        unsigned int word_size;
        void destruct() {delete[] raw_data;}


    double get_value(unsigned int a) {
        return loaded_data[a];
    }


    double get_value(unsigned int a, unsigned int b) {
        unsigned int index = a * shape[1]
                           + b;
        return loaded_data[index];
    }


    double get_value(unsigned int a, unsigned int b, unsigned int c) {
        unsigned int index = a * shape[1] * shape[2] 
                           + b * shape[2]
                           + c;
        return loaded_data[index];

    }


    double get_value(unsigned int a, unsigned int b, unsigned int c, unsigned int d) {
        unsigned int index = a * shape[1] * shape[2] * shape[3]
                           + b * shape[2] * shape[3]
                           + c * shape[3]
                           + d;

        return loaded_data[index];
    }

    double get_value(unsigned int a, unsigned int b, unsigned int c, unsigned int d, unsigned int e) {
        unsigned int index = a * shape[1] * shape[2] * shape[3] * shape[4]
                           + b * shape[2] * shape[3] * shape[4]
                           + c * shape[3] * shape[4]
                           + d * shape[4]
			   + e;

        return loaded_data[index];
    }
    double get_value(unsigned int a, unsigned int b, unsigned int c, unsigned int d, unsigned int e, unsigned int f) {
        unsigned int index = a * shape[1] * shape[2] * shape[3] * shape[4] * shape[5]
                           + b * shape[2] * shape[3] * shape[4] * shape[5]
                           + c * shape[3] * shape[4] * shape[5]
                           + d * shape[4] * shape[5]
                           + e * shape[5]
			   + f;	
        return loaded_data[index];
    }
    double get_value(unsigned int a, unsigned int b, unsigned int c, unsigned int d, unsigned int e, unsigned int f, unsigned int g){
        unsigned int index = a * shape[1] * shape[2] * shape[3] * shape[4] * shape[5] * shape[6]
                           + b * shape[2] * shape[3] * shape[4] * shape[5] * shape[6]
                           + c * shape[3] * shape[4] * shape[5] * shape[6]
                           + d * shape[4] * shape[5] * shape[6]
                           + e * shape[5] * shape[6]
                           + f * shape[6]
			   + g;	 
        return loaded_data[index];
    }






    boost::multi_array<float, 1> get_1d_array() {
        assert(shape.size() == 1);
        boost::multi_array<float,1> arr(boost::extents[shape[0]]);
        for (unsigned int i = 0; i < shape[0]; i++) {
            arr[i] = (float)get_value(i);
        }
          return arr;
    }

    boost::multi_array<float, 2> get_2d_array() {
        assert(shape.size() == 2);
        boost::multi_array<float,2> arr(boost::extents[shape[0]][shape[1]]);
        for (unsigned int i = 0; i < shape[0]; i++) {
            for (unsigned int j = 0; j < shape[1]; j++) {
                arr[i][j] = (float)get_value(i, j);
            }
        }
          return arr;
    }


    boost::multi_array<float, 3> get_3d_array() {
        assert(shape.size() == 3);
        boost::multi_array<float,3> arr(boost::extents[shape[0]][shape[1]][shape[2]]);
        for (unsigned int i = 0; i < shape[0]; i++) {
            for (unsigned int j = 0; j < shape[1]; j++) {
                for (unsigned int k = 0; k < shape[2]; k++) {
                    arr[i][j][k] = (float)get_value(i, j, k);
                }
            }
        }
          return arr;
    }

    boost::multi_array<float, 4> get_4d_array() {
        assert(shape.size() == 4);
        boost::multi_array<float,4> arr(boost::extents[shape[0]][shape[1]][shape[2]][shape[3]]);
        for (unsigned int i = 0; i < shape[0]; i++) {
            for (unsigned int j = 0; j < shape[1]; j++) {
                for (unsigned int k = 0; k < shape[2]; k++) {
                    for (unsigned int l = 0; l < shape[3]; l++) {
                        arr[i][j][k][l] = (float)get_value(i, j, k, l);
                    }
                }
            }
        }
        return arr;
    }

    boost::multi_array<float, 5> get_5d_array() {
        assert(shape.size() == 5);
        boost::multi_array<float,5> arr(boost::extents[shape[0]][shape[1]][shape[2]][shape[3]][shape[4]]);
        for (unsigned int i = 0; i < shape[0]; i++) {
            for (unsigned int j = 0; j < shape[1]; j++) {
                for (unsigned int k = 0; k < shape[2]; k++) {
                    for (unsigned int l = 0; l < shape[3]; l++) {
	                for (unsigned int m = 0; m < shape[4]; m++) {
                            arr[i][j][k][l][m] = (float)get_value(i, j, k, l, m);
			}
                    }
                }
            }
        }
        return arr;
    }

    boost::multi_array<float, 6> get_6d_array() {
        assert(shape.size() == 6);
        boost::multi_array<float,6> arr(boost::extents[shape[0]][shape[1]][shape[2]][shape[3]][shape[4]][shape[5]]);
        for (unsigned int i = 0; i < shape[0]; i++) {
            for (unsigned int j = 0; j < shape[1]; j++) {
                for (unsigned int k = 0; k < shape[2]; k++) {
                    for (unsigned int l = 0; l < shape[3]; l++) {
                        for (unsigned int m = 0; m < shape[4]; m++) {
			    for (unsigned int o = 0; o < shape[5]; o++) {
                                 arr[i][j][k][l][m][o] = (float)get_value(i, j, k, l, m , o);
			    }
                        }
                    }
                }
            }
        }
        return arr;
    }

    boost::multi_array<float, 7> get_7d_array() {
        assert(shape.size() == 7);
        boost::multi_array<float,7> arr(boost::extents[shape[0]][shape[1]][shape[2]][shape[3]][shape[4]][shape[5]][shape[6]]);
        for (unsigned int i = 0; i < shape[0]; i++) {
            for (unsigned int j = 0; j < shape[1]; j++) {
                for (unsigned int k = 0; k < shape[2]; k++) {
                    for (unsigned int l = 0; l < shape[3]; l++) {
                        for (unsigned int m = 0; m < shape[4]; m++) {
                            for (unsigned int o = 0; o < shape[5]; o++) {
				for (unsigned int p = 0; p < shape[6]; p++) {
                                 arr[i][j][k][l][m][o][p] = (float)get_value(i, j, k, l, m , o, p);
				}
                            }
                        }
                    }
                }
            }
        }
        return arr;
    }





};


    



    void parse_npy_header(FILE* fp, unsigned int& word_size, unsigned int*& shape, unsigned int& ndims) {   
        char buffer[256]; 
        
        // Original code line:
        // fread(buffer,sizeof(char),11,fp);

        // Work around for warnd unknow return type:
        size_t dummy = fread(buffer,sizeof(char),11,fp);
        (void)dummy;

        std::string header = fgets(buffer,256,fp); 
        assert(header[header.size()-1] == '\n'); 
     
        int loc1, loc2; 
     
        // fortran order 
        // loc1 = header.find("fortran_order")+16; 
        // bool fortran_order = (header.substr(loc1,5) == "True" ? true : false); 
        // assert(!fortran_order); 
     

        //shape 
        loc1 = header.find("("); 
        loc2 = header.find(")"); 
        std::string str_shape = header.substr(loc1+1,loc2-loc1-1); 
        if(str_shape[str_shape.size()-1] == ',') ndims = 1; 
        else ndims = std::count(str_shape.begin(),str_shape.end(),',')+1; 
        shape = new unsigned int[ndims]; 
        for(unsigned int i = 0;i < ndims;i++) { 
            loc1 = str_shape.find(","); 
            shape[i] = atoi(str_shape.substr(0,loc1).c_str()); 
            str_shape = str_shape.substr(loc1+1); 
        } 
     
        //endian, word size, data type 
        //byte order code | stands for not applicable.  
        //not sure when this applies except for byte array 
        loc1 = header.find("descr")+9; 
        // bool littleEndian = (header[loc1] == '<' || header[loc1] == '|' ? true : false); 
        // assert(littleEndian); 
     
        // char type = header[loc1+1]; 
        //assert(type == map_type(T)); 
     
        std::string str_ws = header.substr(loc1+2); 
        loc2 = str_ws.find("'"); 
        word_size = atoi(str_ws.substr(0,loc2).c_str()); 
    } 


    NpyArray load_the_npy_file(FILE* fp) {
        unsigned int* shape;
        unsigned int ndims, word_size;
        cnpy::parse_npy_header(fp,word_size,shape,ndims);
        unsigned long long size = 1; //long long so no overflow when multiplying by word_size
        for(unsigned int i = 0;i < ndims;i++) size *= shape[i];
    
        cnpy::NpyArray arr;
        arr.word_size = word_size;
        arr.shape = std::vector<unsigned int>(shape,shape+ndims);
        delete[] shape;
        arr.raw_data = new char[size*word_size];
        int nread = fread(arr.raw_data,word_size,size,fp);
        return arr;
        (void)nread;
    }


    NpyArray npy_load(std::string fname) {
        FILE* fp = fopen(fname.c_str(), "rb");
    
        if(!fp) {
            printf("npy_load: Error! Unable to open file %s!\n",fname.c_str());
            abort();
        }
    
        NpyArray arr = load_the_npy_file(fp);
        arr.loaded_data = reinterpret_cast<double*>(arr.raw_data);
    
        fclose(fp);
        return arr;
    }


}

#endif
