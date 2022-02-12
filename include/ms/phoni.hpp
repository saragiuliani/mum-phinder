/* ms_pointers - Computes the matching statistics pointers from BWT and Thresholds
    Copyright (C) 2020 Massimiliano Rossi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file phoni.hppint_vector<>(r, 0, log_n);
   \brief phoni.hpp Computes the matching statistics pointers from BWT.
   \authors Sara Giuliani, Giuseppe Romana
   \date 07/2020
*/

#ifndef _MS_POINTERS_HH
#define _MS_POINTERS_HH


/** FLAGS **/
// #define MEASURE_TIME 1  //measure the time for LCE and backward search?
//#define NAIVE_LCE_SCHEDULE 1 //stupidly execute two LCEs without heurstics
#ifndef NAIVE_LCE_SCHEDULE //apply a heuristic
#define SORT_BY_DISTANCE_HEURISTIC 1 // apply a heuristic to compute the LCE with the closer BWT position first
#endif

#ifndef DCHECK_HPP
#define DCHECK_HPP
#include <string>
#include <sstream>
#include <stdexcept>

// # GR libs
#include <numeric>
#include <algorithm>

#ifndef DCHECK
#ifdef NDEBUG
#define ON_DEBUG(x)
#define DCHECK_(x,y,z)
#define DCHECK(x)
#define DCHECK_EQ(x, y)
#define DCHECK_NE(x, y)
#define DCHECK_LE(x, y)
#define DCHECK_LT(x, y)
#define DCHECK_GE(x, y)
#define DCHECK_GT(x, y)
#else//NDEBUG
#define ON_DEBUG(x) x
#define DCHECK_(x,y,z) \
  if (!(x)) throw std::runtime_error(std::string(" in file ") + __FILE__ + ':' + std::to_string(__LINE__) + (" the check failed: " #x) + ", we got " + std::to_string(y) + " vs " + std::to_string(z))
#define DCHECK(x) \
  if (!(x)) throw std::runtime_error(std::string(" in file ") + __FILE__ + ':' + std::to_string(__LINE__) + (" the check failed: " #x))
//#define DCHECK_EQ(x, y)   if (!(x==y)) throw std::runtime_error(std::string(" in file ") + __FILE__ + ':' + std::to_string(__LINE__) + (" the check failed: " #x) + ", we got " + std::to_string(y) + " vs " + std::to_string(z))
#define DCHECK_EQ(x, y) DCHECK_((x) == (y), x,y)
#define DCHECK_NE(x, y) DCHECK_((x) != (y), x,y)
#define DCHECK_LE(x, y) DCHECK_((x) <= (y), x,y)
#define DCHECK_LT(x, y) DCHECK_((x) < (y) ,x,y)
#define DCHECK_GE(x, y) DCHECK_((x) >= (y),x,y )
#define DCHECK_GT(x, y) DCHECK_((x) > (y) ,x,y)
#endif //NDEBUG
#endif //DCHECK
#endif /* DCHECK_HPP */


#include <common.hpp>

#include <malloc_count.h>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

#include <r_index.hpp>

#include<ms_rle_string.hpp>


#include "PlainSlp.hpp"
#include "PoSlp.hpp"
#include "ShapedSlp_Status.hpp"
#include "ShapedSlp.hpp"
#include "ShapedSlpV2.hpp"
#include "SelfShapedSlp.hpp"
#include "SelfShapedSlpV2.hpp"
#include "DirectAccessibleGammaCode.hpp"
#include "IncBitLenCode.hpp"
#include "FixedBitLenCode.hpp"
#include "SelectType.hpp"
#include "VlcVec.hpp"


#include <iostream>
#include <stdexcept>


#ifdef MEASURE_TIME
struct Stopwatch {
	std::chrono::high_resolution_clock::time_point m_val = std::chrono::high_resolution_clock::now();

	void reset() {
		m_val = std::chrono::high_resolution_clock::now();
	}
	double seconds() const {
		return std::chrono::duration<double, std::ratio<1>>(std::chrono::high_resolution_clock::now() - m_val).count();
	}

};
#endif//MEASURE_TIME

using var_t = uint32_t;
using Fblc = FixedBitLenCode<>;
using SelSd = SelectSdvec<>;
using SelMcl = SelectMcl<>;
using DagcSd = DirectAccessibleGammaCode<SelSd>;
using DagcMcl = DirectAccessibleGammaCode<SelMcl>;
using Vlc64 = VlcVec<sdsl::coder::elias_delta, 64>;
using Vlc128 = VlcVec<sdsl::coder::elias_delta, 128>;


template <
      class SlpT = SelfShapedSlp<var_t, DagcSd, DagcSd, SelSd>,
      class sparse_bv_type = ri::sparse_sd_vector,
      class rle_string_t = ms_rle_string_sd
          >
class ms_pointers : ri::r_index<sparse_bv_type, rle_string_t>
{
public:

    unsigned char seq_nt4_table[256] = {
        0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

    // std::vector<size_t> thresholds;
    SlpT slp;

    // std::vector<ulint> samples_start;
    int_vector<> samples_start; // SA sampled where the run begins +1 (analouguosly samples_last)
    int_vector<> lcp_start;
    int_vector<> lcp_last;
    //std::vector<size_t> lcp_start;
    //std::vector<size_t> lcp_last;    // int_vector<> samples_end;
    // std::vector<ulint> samples_last;

    // static const uchar TERMINATOR = 1;
    // bool sais = true;
    // /*
    //  * sparse RLBWT: r (log sigma + (1+epsilon) * log (n/r)) (1+o(1)) bits
    //  */
    // //F column of the BWT (vector of 256 elements)
    // std::vector<ulint> F;
    // //L column of the BWT, run-length compressed
    // rle_string_t bwt;
    // ulint terminator_position = 0;
    // ulint r = 0; //number of BWT runs

    typedef size_t size_type;

    ms_pointers()
        : ri::r_index<sparse_bv_type, rle_string_t>()
        {}

    void build(const std::string& filename)
    {
        verbose("Building the r-index from BWT");

        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

        std::string bwt_fname = filename + ".bwt";

        verbose("RLE encoding BWT and computing SA samples");

        if (true)
        {
            std::string bwt_heads_fname = bwt_fname + ".heads";
            std::ifstream ifs_heads(bwt_heads_fname);
            std::string bwt_len_fname = bwt_fname + ".len";
            std::ifstream ifs_len(bwt_len_fname);
            this->bwt = rle_string_t(ifs_heads, ifs_len);

            ifs_heads.seekg(0);
            ifs_len.seekg(0);
            this->build_F_(ifs_heads, ifs_len);
        }
        else
        {
            std::ifstream ifs(bwt_fname);
            this->bwt = rle_string_t(ifs);

            ifs.seekg(0);
            this->build_F(ifs);
        }



        this->r = this->bwt.number_of_runs();
        ri::ulint n = this->bwt.size();
        // int log_r = bitsize(uint64_t(this->r));
        int log_n = bitsize(uint64_t(this->bwt.size()));

        verbose("Number of BWT equal-letter runs: r = " , this->r);
        verbose("Rate n/r = " , double(this->bwt.size()) / this->r);
        verbose("log2(r) = ......................" , log2(double(this->r)));
        verbose("log2(n/r) = " , log2(double(this->bwt.size()) / this->r));



        read_samples(filename + ".ssa", this->r, n, samples_start);
        read_samples(filename + ".esa", this->r, n, this->samples_last);
        read_lcp(filename + ".lcp", this->r, this->lcp_start, this->lcp_last, n);
        verbose("lcp_start ", lcp_start.size(), "lcp_last ", lcp_last.size(), "samples_start ", samples_start.size());


        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();
        verbose("R-index construction complete là");
        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

//load_grammar(filename);

        verbose("text length: ", slp.getLen());
        verbose("bwt length: ", this->bwt.size());
    }

    void read_lcp(std::string filename, ulint r, int_vector<> &lcp_start,  int_vector<> &lcp_last, size_t bound = std::numeric_limits<size_t>::max()){
  		// Get the size in bits of the max lcp
      int bit_size = sizeof(size_t) * 8;
      if (bound > 0)      
        bit_size = bitsize(bound);

      struct stat filestat;
      FILE *fd;

      if ((fd = fopen(filename.c_str(), "r")) == nullptr)
          error("open() file " + filename + " failed");

      int fn = fileno(fd);
      if (fstat(fn, &filestat) < 0)
          error("stat() file " + filename + " failed");

      if (filestat.st_size % SSABYTES != 0)
          error("invalid file " + filename);

	    // Create the vectors
      lcp_start = int_vector<>(r, 0, bit_size);
      lcp_last = int_vector<>(r, 0, bit_size);

  		// Read the vector
      uint64_t lcp_value = 0;
  		uint64_t prev_lcp_value = 0;
      size_t i = 1; //index on the lcp file
  		int k = 0; //index on lcp_start and last;
  		int prev_run = 0; //index of the last run read
      bool new_run = true;
  		fread((char *)&lcp_value, SSABYTES, 1, fd);
  		assert(lcp_value == 0);
      while (fread((char *)&lcp_value, SSABYTES, 1, fd))
      {
        const size_t curr_run = this->bwt.run_of_position(i++);
        lcp_value = lcp_value <= bound ? lcp_value : bound;
        if ( curr_run != prev_run)
        {
          if (not new_run) lcp_last[k] = prev_lcp_value;
          new_run = true;
          ++k;
          lcp_start[k] = lcp_last[k] = 0;
        } 
        else if ( new_run )
        {
          lcp_start[k] = lcp_value <= bound ? lcp_value : bound;
          new_run = false;
        }
        prev_lcp_value = lcp_value;
        prev_run = curr_run;
      }
      if (not new_run) lcp_last[k] = prev_lcp_value;
      assert(k == (r - 1));
      fclose(fd);

      sdsl::util::bit_compress(lcp_start);
      sdsl::util::bit_compress(lcp_last);
  	}

    void load_grammar(const std::string& filename) {
        {
            verbose("Load Grammar");
            std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

            ifstream fs(filename + ".slp");
            slp.load(fs);

            std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();
            verbose("Memory peak: ", malloc_count_peak());
            verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
        }
        DCHECK_EQ(slp.getLen()+1, this->bwt.size());
    }



    void read_samples(std::string filename, ulint r, ulint n, int_vector<> &samples)
    {
        int log_n = bitsize(uint64_t(n));
        struct stat filestat;
        FILE *fd;

        if ((fd = fopen(filename.c_str(), "r")) == nullptr)
            error("open() file " + filename + " failed");

        int fn = fileno(fd);
        if (fstat(fn, &filestat) < 0)
            error("stat() file " + filename + " failed");

        if (filestat.st_size % SSABYTES != 0)
            error("invilid file " + filename);

        size_t length = filestat.st_size / (2*SSABYTES);
        //Check that the length of the file is 2*r elements of 5 bytes
        assert(length == r);

        // Create the vector
        samples = int_vector<>(r, 0, log_n);

        // Read the vector
        uint64_t left = 0;
        uint64_t right = 0;
        size_t i = 0;
        while (fread((char *)&left, SSABYTES, 1, fd) && fread((char *)&right, SSABYTES, 1,fd))
        {
            //verbose("read_sample pre",right);
            ulint val = (right ? right - 1 : n - 1);
            //verbose("read_sample post",val);
            assert(bitsize(uint64_t(val)) <= log_n);
            //verbose("read_sample",val);
            samples[i++] = val;
        }

        fclose(fd);
    }

    void write_int(ostream& os, const size_t& i){
      os.write(reinterpret_cast<const char*>(&i), sizeof(size_t));
    }

// ms_ref	9 - 3 - - - 3 - -
//		 	0 1 2 3 4 5 6 7 8
// p     	a c g t c g a g c
// mum   	t f t f f f t f f
// mum_idx  0 2 6
// idx		0 1 2  -> sort w.r.t. pos in T ->

// t ............|------*********----------|..................................................
// p ...|------*********----------|....|------*********----------|...*********...
    //int cmpfunction(const void *a, const void *b){ return ms_references[a] - ms_references[b]; }

    void remove_duplicatesMUMs(std::vector<std::tuple<size_t, size_t, size_t> >& mums){
      verbose("Removing MUMs duplicates in the pattern!");

      size_t j = 0;

      if (mums.size() > 1){
        std::sort(mums.begin(), mums.end());
        size_t p = std::get<0>(mums[0]);
        size_t l = std::get<1>(mums[0]);
        size_t p_tmp;
        size_t l_tmp;
        size_t last_i = 0;
        bool unique = true;
        for (size_t i = 1; i < mums.size(); i++){
          if(i%10000 == 0) std::cout << "\rProcessed: " << i << "/" << mums.size() <<" MUMs.";
          p_tmp = std::get<0>(mums[i]);
          l_tmp = std::get<1>(mums[i]);

          if (p == p_tmp){
            if(l == l_tmp){
              unique = false;
            }else{
              if (l < l_tmp){
                l = l_tmp;
                last_i = i;
                unique = true;
              }
            }
          }else{
            if(l < l_tmp + (p_tmp -p)){
              if (unique)
                mums[j++] = mums[last_i];
              p = p_tmp;
              l = l_tmp;
              last_i = i;
              unique = true;
            }
          }
        }

        if(unique)             
          mums[j++] = mums[last_i];
        std::cout << "\r";
        verbose("Processed: ", mums.size(), "/", mums.size(), " MUMs.");

      }else{
        if (mums.size() == 1)
          mums[j++] = mums[0];
      }
      mums.resize(j);
      mums.shrink_to_fit();
    }

    void print_MUMs_to_file(ofstream& mums_file, const std::vector<std::tuple<size_t, size_t, size_t> >& mums)
    {
      for( size_t i = 0; i < mums.size(); ++i)
      {
        const size_t& ref = std::get<0>(mums[i]);
        const size_t& len = std::get<1>(mums[i]);
        const size_t& query = std::get<2>(mums[i]);

        mums_file <<std::setw(8)<<ref + 1 << std::setw(10) << query + 1 <<  std::setw(10) << len << std::endl;
      }

    }

    // Computes the matching statistics pointers for the given pattern
    //std::pair<std::vector<size_t>, std::vector<size_t>>
    // size_t query_mum(const std::string& patternfile, const std::string& len_filename, const std::string& ref_filename, const std::string& mum_pos_t_filename, const std::string& mum_pos_p_filename, const std::string& mum_len_filename, bool mumref = false) {
    size_t query_mum(const std::string& patternfile, const std::string& outfile, bool mumref = false) {
      const char* p;
      size_t m;
      map_file(patternfile.c_str(), p, m);

      auto pattern_at = [&] (const size_t pos) {
        return p[pos];
      };

      auto res = _query_mum(pattern_at, m, mumref);

      ofstream outstream(outfile, std::ios::binary);
      print_MUMs_to_file(outstream, res);
      outstream.close();

      // int_vector<>& ms_references = std::get<0>(res);
      // int_vector<>& ms_lengths = std::get<1>(res);
      // std::vector<bool>& mum = std::get<2>(res);


      // ofstream len_file(len_filename, std::ios::binary);
      // ofstream ref_file(ref_filename, std::ios::binary);
      // ofstream mum_pos_t_file(mum_pos_t_filename, std::ios::out);
      // ofstream mum_pos_p_file(mum_pos_p_filename, std::ios::out);
      // ofstream mum_len_file(mum_len_filename, std::ios::out);

      // //mum = remove_duplicates(mum);
      // for (int i = m-1; i >= 0; i--){
      //   // verbose("m:", m, "i:", i, "ms_ref:", ms_references[i], "ms_len:", ms_lengths[i], "ms_twice:", ms_twice[i]);
      //   write_int(ref_file, ms_references[i]);
      //   write_int(len_file, ms_lengths[i]);
      // }



      // //mum_pos_t_file << pattern.first << endl;
      // for (size_t i = 0; i < m; i++){
      //   if (mum[i]){
      //     mum_pos_t_file << ms_references[i] << " ";
      //   }
      // }
      // mum_pos_t_file << endl;

      // //mum_len_file << pattern.first << endl;
      // for(size_t i = 0; i < m; i++){
      //   if (mum[i])
      //     mum_len_file << ms_lengths[i] << " ";
      // }
      // mum_len_file << endl;

      // for(size_t i = 0; i < m; i++){
      //   if (mum[i])
      //     mum_pos_p_file << i << " ";
      // }
      // mum_len_file << endl;

      return m;
    }

    template <typename string_t>
    std::vector<std::tuple<size_t, size_t, size_t>> query_mum(const string_t &pattern, const size_t m, bool mumref = false)
    {
        auto pattern_at = [&] (const size_t pos) {
          return pattern[pos];
        };

        return _query_mum(pattern_at, m, mumref);
    }


    template <typename Func>
    inline std::vector<std::tuple<size_t, size_t, size_t>> _query_mum(Func pattern_at, const size_t m, bool mumref = false, bool acgt_only = true)
    {
        ON_DEBUG(verbose("DEBUG MODE ON!"));
        const size_t n = slp.getLen();

        size_t last_len = 0;
        size_t last_ref = 0;

        //TODO: we could allocate the file here and store the numbers *backwards* !
        // int_vector<> ms_references(m,0); //! stores at (m-i)-th entry the text position of the match with pattern[m-i..] when computing the matching statistics of pattern[m-i-1..]
        // int_vector<> ms_lengths(m,1);
        // int_vector<> ms_p_lengths(m,0);
        // int_vector<> ms_s_lengths(m,0);
        // int_vector<> ms_twice(m,0);
        // std::vector<bool> mum(m,false);
        size_t ms_ref = 0; //! stores at (m-i)-th entry the text position of the match with pattern[m-i..] when computing the matching statistics of pattern[m-i-1..]
        size_t ms_len = 0;
        size_t ms_p_len = 0;
        size_t ms_s_len = 0;
        size_t ms_twice = 0;
        std::vector<std::tuple<size_t, size_t, size_t>> mums;

        //!todo: we need here a while look in case that we want to support suffixes of the pattern that are not part of the text
        DCHECK_GT(this->bwt.number_of_letter(pattern_at(m-1)), 0);

        //! Start with the last character of the pattern

        // pos = first occurrence in the bwt of the last char of the pattern
        size_t pos = this->bwt_size() - 1;
        ms_ref = this->get_last_run_sample();

        // auto pos = this->bwt.select(0, pattern_at(m-1));

        // {
        //     // run_of_j = i when char in pos is in the i^th run of the bwt
        //     const ri::ulint run_of_j = this->bwt.run_of_position(pos);
        //     // first match: position in the reference of the first character in the i^th run
        //     ms_ref = samples_start[run_of_j];

        //     // # GR
        //     // **********************
        //     // Algo 1, ll. 2-5

        //     ms_p_len = 0;

        //     if (this->bwt.number_of_letter(pattern_at(m-1)) > 1){
        //       ms_s_len = 1;
        //       ms_twice = 1;
        //     }
        //     else{ // pattern_at(m-1) occurs once in the ref => ms_successor and ms_twice do not exists
        //       ms_s_len = 0;
        //       ms_twice = 0;
        //     }

        //     DCHECK_EQ(slp.charAt(ms_ref),  pattern_at(m-1));
        // }
        //verbose("ALGO1pos = ", pos);

        // update of pos following LF: try to match the pattern to the reference backward
        // pos = LF(pos, pattern_at(m-1)); // # Algo 1, l. 7
        //verbose("ALGO1 LF pos = ", pos);
        #ifdef MEASURE_TIME
        double time_lce = 0;
        double time_backwardstep = 0;
        size_t count_lce_total = 0;
        size_t count_lce_skips = 0;
        #endif

        // extend the match
        // for (size_t i = 1; i < m; i++) {
        for (size_t i = 0; i < m; i++) {
            if(i%10000 == 0) std::cout << "\rProcessed: " << i << "/" << m <<" characters.";
            last_len = ms_len;
            last_ref = ms_ref;
            const size_t last_p_lengths = ms_p_len;
            const size_t last_s_lengths = ms_s_len;
            const size_t last_twice = ms_twice;
            const auto c = pattern_at(m - i - 1);
            // DCHECK_EQ(ms_lengths[m-i], last_len); // N/A
            // ON_DEBUG(DCHECK_EQ(ms_references[m-i], last_ref)); // N/A


            const size_t number_of_runs_of_c = this->bwt.number_of_letter(c); // !!!!!!!!!!!! number_of_runs_of_c == NUMBER OF OCCURRENCES OF CHARARCTER C IN REFERENCE
            //verbose("N runs of", c, " = ", number_of_runs_of_c);
            if((number_of_runs_of_c == 0 )|| (acgt_only and seq_nt4_table[c] > 3)) {
                ms_len = 0;
                pos = this->bwt_size() - 1;
                ms_ref = this->get_last_run_sample();
                ms_p_len = 0;
                ms_s_len = 0;
                ms_twice = 0;
            }
            else if (pos < this->bwt.size() && this->bwt[pos] == c) { // # GR Algo 2: MMatch
                ms_len = last_len+1;

                ms_ref = last_ref-1;

                // # GR
                // **********************
                // *****  Algo 2  *******
                // **********************
                const ri::ulint rank = this->bwt.rank(pos, c);

                ri::ulint run0, run1;
                size_t sa0, sa1;
                if (pos > 0 and this->bwt[pos - 1] == c){
                  ms_p_len = ms_p_len + 1;
                }
                else{
                  if(rank > 0) {
                    sa0 = this->bwt.select(rank-1, c);
                    run0 = this->bwt.run_of_position(sa0);
                    DCHECK_LT(sa0, pos);
                    ms_p_len = std::min(static_cast<unsigned long>(ms_p_len)+1, lceToRBounded(slp, (this->samples_last[run0])%n, ms_ref, ms_p_len+1));
                  }
                  else{
                    ms_p_len = 0;
                  }
                }
                if ((pos < n - 1) and this->bwt[pos + 1] == c){
                  ms_s_len = ms_s_len + 1;
                }
                else{
                  if(rank < number_of_runs_of_c - 1) { //Q: number of letters of c? Is it not supposed to be < number_of_runs_of_c)
                    sa1 = this->bwt.select(rank+1, c);
                    run1 = this->bwt.run_of_position(sa1);
                    DCHECK_GT(sa1, pos);
                    ms_s_len = std::min(static_cast<unsigned long>(ms_s_len)+1, lceToRBounded(slp, (this->samples_start[run1])%n, ms_ref, ms_s_len+1));
                  }
                  else{
                    ms_s_len = 0;
                  }
                }
                ms_twice = std::max(ms_p_len, ms_s_len);

                // # GR
                // **********************
                // *****  END 2  ********
                // **********************
            }
            else {
              //verbose("ELSE 3, i = ", i);
              // # GR
              // **********************
              // *****  Algo 3  *******
              // **********************
              const ri::ulint rank = this->bwt.rank(pos, c);
              ri::ulint run0, run1;
              uint64_t lcp0, lcp1; // # GR
              size_t sa0, sa1;
              // # GR ll. 4-10
              //if( i > 0 && i < m-1){verbose("MS = ", ms_lengths[m - i]);}

              if(rank > 0) {
                sa0 = this->bwt.select(rank-1, c);
                //verbose("1", c, "---->", m-i-1);
                //verbose("SA0 = ", sa0, "pos = ", pos);
                DCHECK_LT(sa0, pos);
                run0 = this->bwt.run_of_position(sa0);
                if (pos > 0 and sa0 == pos - 1){
                  //verbose("2", "---->", m-i-1);
                  lcp0 = ms_p_len + 1;
                }
                else{
                  //verbose("3", "---->", m-i-1);
                  //verbose("MS_len = ", ms_lengths[m - i], "LCE = ", lceToRBounded(slp, (this->samples_last[run0]+1)%n, ms_references[m - i], ms_lengths[m - i]));
                  const auto lce = lceToRBounded(slp, (this->samples_last[run0]+1)%n, ms_ref, ms_len);
                  lcp0 = std::min(static_cast<unsigned long>(ms_len), lce ) + 1;
                  //verbose("LCP0 3 = ", lcp0);
                }
              }
              else{
                //verbose("4", "---->", m-i-1);
                lcp0 = 0;
              }
              // # GR ll. 11-17
              if(rank < number_of_runs_of_c) {
                //verbose("5", "---->", m-i-1);
                sa1 = this->bwt.select(rank, c);
                DCHECK_GT(sa1, pos);
                run1 = this->bwt.run_of_position(sa1);
                if (pos < n-1 and sa1 == pos + 1){
                  //verbose("6", "---->", m-i-1);
                  lcp1 = ms_s_len + 1;
                }
                else{
                  //verbose("7", "---->", m-i-1);
                  //verbose("MS_len = ", ms_len, "LCE = ", lceToRBounded(slp, (this->samples_start[run1]+1)%n, ms_references[m - i], ms_len));
                  const size_t t1 = (this->samples_start[run1]+1)%n;
                  lcp1 = std::min(static_cast<unsigned long>(ms_len), lceToRBounded(slp, t1, ms_ref, ms_len)) + 1;
                  //verbose("LCP1 7 = ", lcp1);
                }
              }
              else{
                lcp1 = 0;
              }
              // # GR ll. 18-25
              //verbose("LCP0 = ", lcp0, "LCP1 = ", lcp1);
              if (lcp0 <= lcp1){
                // l. 19
                //verbose("8", "---->", m-i-1, "run1 =", run1, "sa1 =", sa1, this->samples_start[run1], this->samples_start[run1+1]);
                ms_ref  = this->samples_start[run1]; // --------------> PROBLEMA
                ms_len = lcp1;
                //write_ref(this->samples_start[run1]); // --------------> PROBLEMA
                //write_len(lcp1);
                ms_p_len = lcp0;
                if (rank < number_of_runs_of_c - 1) {
                  //verbose("9", "---->", m-i-1);
                  auto s_pos = this->bwt.select(rank + 1, c);
                  // ll. 21-22
                  if(s_pos == sa1 + 1){
                  //verbose("10", "---->", m-i-1);
                    //int_vector<0> abbreviated_ms_length = ms_lengths[m - i] & INT_MAX;
                    const size_t lcpS = this->lcp_start[this->bwt.run_of_position(s_pos)];
                    ms_s_len = std::min(static_cast<unsigned long>(ms_len - 1), static_cast<unsigned long>(lcpS)) +1;
                  }
                  // ll. 23-24
                  else{
                    ms_s_len = std::min(static_cast<unsigned long>(ms_len), lceToRBounded(slp,
                                                  this->samples_start[run1],	// sa[sa1]
                                                  this->samples_start[this->bwt.run_of_position(s_pos)],	// sa
                                                  ms_len));
                  }
                }
                // case comment l. 20
                else{
                  //7verbose("12", "---->", m-i-1);
                  ms_s_len = 0;
                }
                // l. 25
                pos = sa1;
                //verbose("SA1 new pos = ", pos);
              }
              // # GR ll. 26-33
              else{
                // l. 27
                //verbose("13", "---->", m-i-1, "run0 =", run0, "sa0 =", sa0, this->samples_last[run0], this->samples_last[run0+1]);

                ms_ref  = this->samples_last[run0]; 
                ms_len = lcp0;
                ms_s_len = lcp1;
                if (rank > 1) {
                  auto p_pos = this->bwt.select(rank - 2, c);
                  // ll. 29-30
                  if(p_pos == sa0 - 1){
                    ms_p_len = std::min(static_cast<unsigned long>(ms_len - 1), static_cast<unsigned long>(lcp_last[this->bwt.run_of_position(p_pos)]))+1;
                  }
                  // ll. 31-32
                  else{
                    ms_p_len = std::min(static_cast<unsigned long>(ms_len - 1), lceToRBounded(slp,
                                                  (this->samples_last[this->bwt.run_of_position(sa0)]+1)%n,		// sa[sa0]
                                                  (this->samples_last[this->bwt.run_of_position(p_pos)]+1)%n,	// sa[p_pos]
                                                  ms_len - 1)) +1;

                  }
                }
                // case comment l. 28
                else{
                  ms_p_len = 0;
                }
                //l. 33
                pos = sa0;
              }
              ms_twice = std::max(ms_p_len, ms_s_len);

              // # GR
              // **********************
              // *****  END 3  ********
              // **********************
            }



            #ifdef MEASURE_TIME
            Stopwatch s;
            #endif
            pos = LF(pos, c); //! Perform one backward step
            //verbose("LF new pos = ", pos);
            #ifdef MEASURE_TIME
            time_backwardstep += s.seconds();
            #endif


            if (m - i - 1 == 0 && ms_len > ms_twice) 
              mums.push_back( std::make_tuple(ms_ref, ms_len, 0) );
            else if (ms_len <= last_len && last_len > last_twice)
              mums.push_back( std::make_tuple(last_ref, last_len, m - i) ); 

            

        }
        std::cout << "\r";
        verbose("Processed: ", m, "/", m, " characters.");
        verbose("Found ", mums.size(), " MUMs candidates.");

    
      if( not mumref )
        remove_duplicatesMUMs(mums);

      #ifdef MEASURE_TIME
      cout << "Time backwardsearch: " << time_backwardstep << std::endl;
      cout << "Time lce: " << time_lce << std::endl;
      cout << "Count lce: " << count_lce_total << std::endl;
      cout << "Count lce skips: " << count_lce_skips << std::endl;
      #endif
      return mums;
      // return std::make_tuple(ms_references,ms_lengths,mum);
    }
    /*
     * \param i position in the BWT
     * \param c character
     * \return lexicographic rank of cw in bwt
     */
    ulint LF(ri::ulint i, ri::uchar c)
    {
        // //if character does not appear in the text, return empty pair
        // if ((c == 255 and this->F[c] == this->bwt_size()) || this->F[c] >= this->F[c + 1])
        //     return {1, 0};
        //number of c before the interval
        ri::ulint c_before = this->bwt.rank(i, c);
        // number of c inside the interval rn
        ri::ulint l = this->F[c] + c_before;
        //verbose("LF i = ", i, "LF res = ", l);
        return l;
    }

    /* serialize the structure to the ostream
     * \param out     the ostream
     */
    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") // const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        out.write((char *)&this->terminator_position, sizeof(this->terminator_position));
        written_bytes += sizeof(this->terminator_position);
        written_bytes += my_serialize(this->F, out, child, "F");
        written_bytes += this->bwt.serialize(out);
        written_bytes += this->samples_last.serialize(out);

        written_bytes += samples_start.serialize(out, child, "samples_start");
        written_bytes += lcp_start.serialize(out, child, "lcp_start");
        written_bytes += lcp_last.serialize(out, child, "lcp_last");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;

    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream &in, const std::string& filename)
    {
        in.read((char *)&this->terminator_position, sizeof(this->terminator_position));
        my_load(this->F, in);
        this->bwt.load(in);
        this->r = this->bwt.number_of_runs();
        this->samples_last.load(in);
        this->samples_start.load(in);
        this->lcp_start.load(in);
        this->lcp_last.load(in);

        load_grammar(filename);
    }

    // // From r-index
    // ulint get_last_run_sample()
    // {
    //     return (samples_last[r - 1] + 1) % bwt.size();
    // }

    protected :

    vector<ulint> build_F_(std::ifstream &heads, std::ifstream &lengths)
    {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);

        this->F = vector<ulint>(256, 0);
        int c;
        {
          ulint i = 0;
          while ((c = heads.get()) != EOF)
          {
            size_t length;
            lengths.read((char *)&length, 5);
            if (c > TERMINATOR)
              this->F[c] += length;
            else
            {
              this->F[TERMINATOR] += length;
              this->terminator_position = i;
            }
            i++;
          }
        }
        for (ulint i = 255; i > 0; --i)
            this->F[i] = this->F[i - 1];
        this->F[0] = 0;
        for (ulint i = 1; i < 256; ++i)
            this->F[i] += this->F[i - 1];
        return this->F;
    }
    };

#endif /* end of include guard: _MS_POINTERS_HH */
