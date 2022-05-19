/* Orignal author: Bastien Cazaux
 * Modifications: Michel Henrichs
 */

#include "pbwt.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <vector>

#include <iomanip>

#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <algorithm>

#include <tuple>
#include <list>

#include <unistd.h>

using namespace std;

typedef long long int int_all;

// g++  -Ofast -std=c++11 haplotype-pbwt-lite3.cpp -o haplotype-pbwt-lite3 && /usr/bin/time  -f "\n\tTime %E\n\tMemory %M Kbytes (%t Kbytes)" ./haplotype-pbwt-lite3 -f Input/thousandgenomes_chr22.bm -b 1000000  > Output/chr22_1000000

void see(int_all* a,int_all n){
    for (int_all i = 0;i<n;i++)
        std::cerr << a[i];
    std::cerr << endl;
}

class LR_file_hap {
private:
    string filename;
    int_all nb_line;
    int_all nb_column;
    int_all file_size;
    int_all current_line;

    ifstream input_stream;
    istreambuf_iterator<char> input_it;
    int_all* current_char;

public:
    int_all get_number_line(void){
        return nb_line;
    }

    int_all get_number_column(void){
        return nb_column;
    }

    bool is_end(void){
        return current_line > nb_line;
    }

    void close(void){
        input_stream.close();
    }

    int_all* next(void){
        if (is_end())
            return NULL;
        else
        {
            for (int_all i = 0; i < nb_column; i++)
            {
                // -'0' to turn ascii char 0/1 to numbers 0/1
                current_char[i] = *input_it++ - '0';
            }
            input_it++; // skip line feeds
            current_line++;
            return current_char;
        }
    }

    LR_file_hap(string filename,int_all buffer_size,bool visible = false){

        if (visible){
            cerr << "Setup file reading." << endl;
        }

        current_line = 0;
        nb_line = 0;
        nb_column = 0;
        file_size = 0;

        int BufferSize_big = 4096;
        int BufferSize = buffer_size;

        ifstream str2;
        char _buffer[BufferSize_big];
        str2.rdbuf()->pubsetbuf(_buffer, BufferSize_big);
        str2.open (filename.c_str(), ifstream::in);

        istreambuf_iterator<char> istart (str2);
        istreambuf_iterator<char> iend;

        // read file to find out dimensions (lines and columns)
        while (istart != iend)
        {
           if (*istart == '\n') nb_line++;
           file_size++;
           istart++;
        }
        str2.close();
        nb_column = (int_all)(file_size/nb_line-1);
        if (visible)
            cerr << "File has " << nb_line << " lines, " << nb_column << " columns." << endl;
        
        // set up input steam iterator for reading
        char* input_buffer = new char[BufferSize];
        input_stream.rdbuf()->pubsetbuf(input_buffer, BufferSize);
        input_stream.open(filename.c_str(),ifstream::in);
        input_stream.seekg(0, ios::beg);
        input_it = input_stream;

        current_char = new int_all[nb_column];

        if (visible)
            cerr << "Done with input setup." << endl;
    }
};

// Class to simulate two arrays
class Dtype {
private:
    // We take the first array for choice = true
    bool choice;
    int_all* first_array;
    int_all* second_array;
    int_all size;

public:

    Dtype(int_all len){
        choice = true;
        size = len;
        first_array = new int_all[size];
        second_array = new int_all[size];

    }

    void change_choice(void){
        choice = !choice;
    }

    int_all get_index(int_all i,bool principal){
        if (choice){
            if (principal){
                return first_array[i];
            }
            else{
                return second_array[i];
            }
        }
        else{
            if (!principal){
                return first_array[i];
            }
            else{
                return second_array[i];
            }
        }
    }

    void set_index(int_all i,int_all val,bool principal){
        if (choice){
            if (principal){
                first_array[i] = val;
            }
            else{
                second_array[i] = val;
            }
        }
        else{
            if (!principal){
                first_array[i] = val;
            }
            else{
                second_array[i] = val;
            }
        }
    }
};


class Pbwtlite {
private:
    // k is the step of the PBWT
    int_all position;

    // a in a array corresponding to the permutation in rev lex order
    Dtype* a;
    Dtype* d;
    // norm position of the lcp
    Dtype* e;
    // number of line (max size)
    int_all size;
    // size of the alphabet
    int_all size_sigma;
    // size of different element of e
    int_all size_r;
    // store position of real element of e
    int_all* s;
    // occurence of element of s
    int_all* t;

    // Aux
    int_all* C;
    int_all* P;


public:
    Pbwtlite(int_all len,int_all len_alph){
        position = 0;

        size = len;
        size_sigma = len_alph;

        a = new Dtype(size);
        d = new Dtype(size);

        C = new int_all[size_sigma+1];
        P = new int_all[size_sigma];

        for(int_all i=0; i<size; ++i){
            a->set_index(i,i,true);
            d->set_index(i,0,true);
            a->set_index(i,0,false);
            d->set_index(i,0,false);
        }

    }

    void change_choice(void){
        a->change_choice();
        d->change_choice();
        position++;
    }

    int_all get_index_a(int_all i){
        return a->get_index(i,true);
    }
    int_all get_index_d(int_all i){
        return d->get_index(i,true);
    }

    int_all get_position(void){
        return position;
    }

    int_all get_size(void){
        return size;
    }

    int_all maxd(int_all j, int_all i){
        if (j != i){
            d->set_index(j,max(d->get_index(j,false),maxd(a->get_index(j,false),i)),false);
            a->set_index(j,i+1,false);
        }
        return d->get_index(j,false);
    }

    int_all update_interval(tuple<int_all,int_all,int_all>* ar, int_all size_ar, int_all* auxLa, int_all* auxLb, int_all* auxPa, int_all* auxPb, int_all val_min_size){
        int_all j = 0;
        int_all j_auxL = 0;
        int_all j_auxP = 0;
        for (int i = 0;i<size-1;i++){
            if (d->get_index(i,true) > d->get_index(i+1,true)){
                auxLa[j_auxL] = i;
                auxLb[j_auxL++] = d->get_index(i,true);
                auxPa[j_auxP] = i;
                auxPb[j_auxP++] = d->get_index(i,true);
            }
            if (d->get_index(i,true) < d->get_index(i+1,true)){
                auxPa[j_auxP] = i;
                auxPb[j_auxP++] = d->get_index(i,true);

                while (j_auxL > 0 && auxPb[j_auxP-1] < d->get_index(i+1,true)){
                    if ((position-auxPb[j_auxP-1])*(i-auxLa[j_auxL-1]+1) >= val_min_size ){
                        ar[j++] = make_tuple(auxLa[j_auxL-1],i,position-auxPb[j_auxP-1]);
                    }
                    j_auxP--;

                    if (auxLb[j_auxL-1] <= auxPb[j_auxP-1] && auxLb[j_auxL-1] <= d->get_index(i+1,true)){
                        j_auxL--;
                    }
                }
            }
        }
        auxPa[j_auxP] = size-1;
        auxPb[j_auxP++] = d->get_index(size-1,true);
        while (j_auxL > 0 && auxPb[j_auxP-1] < position){
            if ((position-auxPb[j_auxP-1])*(size-1-auxLa[j_auxL-1]+1) >= val_min_size){
                ar[j++] = make_tuple(auxLa[j_auxL-1],size-1,position-auxPb[j_auxP -1]);
            }
            j_auxP--;

            if (j_auxP > 0 && auxLb[j_auxL-1] <= auxPb[j_auxP-1]){
                j_auxL--;
            }
        }
        return j;
    }

    friend ostream& operator<<(ostream& os, const Pbwtlite& P){
        os << "PBWT (" << P.position << ")\n\t\033[01;33ma\t";
        for(int_all i=0; i<P.size; ++i){
            os << P.a->get_index(i,true) << "\t";
        }
        os << "\033[0m\n\t\033[01;34md\t";
        for(int_all i=0; i<P.size; ++i){
            os << P.d->get_index(i,true) << "\t";
        }
        os << "\033[0m\n";
        return os;
    }

    void next(int_all* r){
        change_choice();

        fill(C,C+size_sigma+1,0);
        fill(P,P+size_sigma,-1);

        int_all b;

        for(int_all i=0; i<size; ++i){ C[r[i]+1]++; }
        for(int_all i=1; i<size_sigma+1; ++i){ C[i] += C[i-1]; }
        for(int_all i=0; i<size; ++i){
            b = r[a->get_index(i,false)];
            a->set_index(C[b],a->get_index(i,false),true);
            a->set_index(i,i+1,false);
            if (P[b] == -1){
                d->set_index(C[b],position,true);
            }
            else{
                d->set_index(C[b],maxd(P[b]+1,i),true);
            }
            P[b] = i;
            C[b]++;
        }
    }
};

int_all build_max_interval_ar(tuple<int_all,int_all,int_all>* ar, int_all ar_size, tuple<int_all,int_all,int_all>* ar_max, int_all* ni, Pbwtlite Q){
    int_all ni_a[Q.get_size()+1];
    int_all j = 0;
    ni_a[0] = 0;
    for (int_all i = 0 ; i < Q.get_size(); i++){
        ni_a[i+1] = ni[Q.get_index_a(i)];
    }
    for (int_all i = 1 ; i < Q.get_size()+1; i++){
        ni_a[i] = ni_a[i]+ni_a[i-1];
    }

    for (int_all i = 0 ; i < ar_size; i++){
        if (ni_a[get<1>(ar[i])+1]- ni_a[get<0>(ar[i])] != 0 && ni_a[get<1>(ar[i])+1]- ni_a[get<0>(ar[i])] != get<1>(ar[i])-get<0>(ar[i])+1){
            ar_max[j++] = make_tuple(get<0>(ar[i]),get<1>(ar[i]),get<2>(ar[i]));
        }
    }
    return j;
}

int_all add_block_file(tuple<int_all,int_all,int_all>* ar, int_all size, Pbwtlite Q, ostream& file,
    int_all min_allele_count = 0, int_all min_seq_count = 0){
    int_all reported_blocks = 0;
    for (int_all i = 0;i < size;i++){
        int_all weight = (get<2>(ar[i]))*(get<1>(ar[i]) - get<0>(ar[i])+1);
        int_all block_start = (Q.get_position()-get<2>(ar[i]));
        int_all block_end = Q.get_position()-1;
        int_all allele_count = block_end - block_start + 1;
        int_all witness = Q.get_index_a(get<0>(ar[i]));
        int_all seq_count = get<1>(ar[i]) - get<0>(ar[i])+1;
        if (allele_count >= min_allele_count && seq_count >= min_seq_count)
        {
            file << weight << " (" << block_start << "-" << block_end << "," << witness << ":" << seq_count << ")\n";
            reported_blocks++;
        }
    }
    return reported_blocks;
}

void see_pourcentage(float p,int_all nb_block,int max_size){
    int nb = (int)(p*(max_size));
    cerr << "\r";
    for (int i = 0 ; i< nb-1; i++){
        cerr << "|";
    }
    cerr << ">";
    for (int i = nb ; i< max_size; i++){
        cerr << " ";
    }
    cerr << " " << fixed << setprecision(1) << p*100 << " %\tNb blocks:\t" << nb_block <<flush;
}


// haplotype-pbwt file_name minimal_block_height minimal_block_width
void usage() {
    cerr << "Usage haplotype-pbwt-lite" << endl;
    cerr << "\t\t-f matrix_filename : a file of haplotypes" << endl;

    cerr << "\t\t[-b minimal_block_size=2] : smallest size to output (area)" << endl;
    cerr << "\t\t[-n buffer_size=8192] : size of the buffer (input files)" << endl;
    exit(0);
}

int_all calc_haploblocks(string block_file, string out_path,
    int_all min_allele_count, int_all min_seq_count,
    int_all minimal_block_size, int_all buffer_size)
{
    if (minimal_block_size < 2) {
        cerr << "minimal_block_size is too small, setting to 2" << endl;
        minimal_block_size = 2;
    }
    if (buffer_size < 2) {
        cerr << "buffer_size is too small" << endl;
        exit(0);
    }
    
    //cerr << "Parameters:" << endl;
    //cerr << "\tbuffer size:\t\t" << buffer_size << endl;
    //cerr << "\tminimal block size:\t" << minimal_block_size << endl;

   LR_file_hap file (block_file, buffer_size);
   int nb_column = file.get_number_column();
   int nb_line = file.get_number_line();

   //cerr << "nb of lines (= nb of SNPs):\t" << nb_line << endl;
   //cerr << "nb of columns (=nb of chromosomes):\t" << nb_column << endl;

   int_all* ni = file.next();
   
   // cerr << "First line" << endl;
   // see(ni,nb_column);
   
   int_all j = 0;
   Pbwtlite P (nb_column,2);

   // New
   tuple<int_all,int_all,int_all>* ar_inter = new tuple<int_all,int_all,int_all>[nb_column];
   int_all size_ar_inter;
   tuple<int_all,int_all,int_all>* ar_inter_max = new tuple<int_all,int_all,int_all>[nb_column];
   int_all size_ar_inter_max;

   // Aux
   int_all* aux1a = new int_all[nb_column];
   int_all* aux1b = new int_all[nb_column];
   int_all* aux2a = new int_all[nb_column];
   int_all* aux2b = new int_all[nb_column];

   int_all nb_blocks = 0;
   
   std::ofstream out_file;
   out_file.open(out_path, std::ofstream::out | std::ofstream::trunc);
   // write bm matrix dimension because that info is lost otherwise
   out_file << '#' << nb_column << "," << nb_line << '\n';
   
    while (!file.is_end()){
        
        /*if (j % (max(1,int(nb_line/100))) == 0){
            see_pourcentage(float(j)/float(nb_line),nb_blocks,80);
        }*/
        j++;

        P.next(ni);
        size_ar_inter = P.update_interval(ar_inter,size_ar_inter,aux1a,aux1b,aux2a,aux2b,minimal_block_size);
        
        ni = file.next();

        if (!file.is_end()){
            size_ar_inter_max = build_max_interval_ar(ar_inter, size_ar_inter, ar_inter_max, ni, P);
            nb_blocks += add_block_file(ar_inter_max, size_ar_inter_max, P, out_file, min_allele_count, min_seq_count);
        }
        else{
            nb_blocks += add_block_file(ar_inter, size_ar_inter, P, out_file);
        }
    }

   //see_pourcentage(1,nb_blocks,80);

   cerr << "\nNb of blocks:\t" << nb_blocks << endl;
   return nb_blocks;
}
