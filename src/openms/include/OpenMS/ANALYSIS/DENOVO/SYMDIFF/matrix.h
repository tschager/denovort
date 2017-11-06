// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Thomas Tschager $
// $Authors: Yves Frank, Simon Rösch, Thomas Tschager, and Valentin Venzin $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_DENOVO_SYMDIFF_MATRIX_H
#define OPENMS_ANALYSIS_DENOVO_SYMDIFF_MATRIX_H

/**
  @brief This file contains implementations of the dynamic programming tables for the linear, the
  position dependent and the neighborhood based models.

  @note To save space, the assumption is made, that scores are integers and in the interval 
  [-1'000'000, 15'777'215].
*/

#include <string>
#include <iostream>
#include <set>
#include <assert.h>
#include "retention_time.h"
#include "amino_acids.h"

/**
  @brief Binary search that returns the position where @em val was found rather than only a bool if it was
  found.
*/
template<class Iter, class T>
Iter binary_find(Iter begin, Iter end, T val)
{
  // Finds the lower bound in at most log(last - first) + 1 comparisons
  Iter i = lower_bound(begin, end, val);
  if (i != end && !(val < *i)) 
    return i; // found
  else 
    return end; // not found
}

// Used in @em transform() to return all keys of a @em std:map
struct RetrieveKey
{
  template <typename T>
  typename T::first_type operator()(T keyValuePair) const
  {
      return keyValuePair.first;
  }
};

/** 
  @brief Struct containing (a,b,time,score) that is stored in the dynamic programming table for the position
  dependent retention time model. The size of the struct is 8 bytes.

  @note Score must be in interval [-1'000'000, 15'777'215], time must be an integer
*/
unsigned int CHAR_MASK = 0xff000000;
unsigned int SCORE_MASK = 0x00ffffff;
unsigned int SCORE_OFFSET = 1000000;
struct data
{
  int t;
  // First byte corresponds to (a,b) (key from mapping (a,b)->char), remaining three bytes are used for the score
  int key_score;

  // Constructor
  data (int score, int time, char key) : t(time), key_score(((key << 24) & CHAR_MASK) | ((score + SCORE_OFFSET) & SCORE_MASK)) {}

  // Getter and setter
  char get_key() 
  { 
    return (key_score & CHAR_MASK) >> 24; 
  }
  int get_score() 
  { 
    return (key_score & SCORE_MASK) - SCORE_OFFSET; 
  };
  void set_score(int score, char key) 
  { 
    key_score = ((key << 24) & CHAR_MASK) | ((score + SCORE_OFFSET) & SCORE_MASK); 
  }
  
  // Operator, order is according to retention time @em t
  bool operator <(const data &o) 
  { 
    return t < o.t; 
  }
};

/**
  @brief Abstract matrix class implements a dynamic programming table.
  @note Parameters @em a and @em b may not always be necessary; then arbitrary values may be supplied.
*/
class Matrix 
{
  public:

    /**
      @brief Stores all occurring distinct retention time values at Q[v][e] (for all a,b) in @em rt_dim

      @param v Vertex
      @param e Edge
      @param rt_dim Pointer to list that is to be filled with distinct retention time values stored at Q[v][e]
    */
    virtual void get_rt_dim(int v, int e, std::vector<int> &rt_dim) = 0;

    /**
      @brief Stores all scores of Q[v][e] in @em rt_a_pair.
      
      @param v Vertex
      @param e Edge
      @param rt_a_pair Vector of pairs where scores are to be stored
    */
    virtual void get_rt_a_pairs(int v, int e, std::vector<std::pair <int, int> > &rt_a_pair) = 0;
    
    /**
      @brief Get score: Q[v][e][a][b][t]. If position is empty return NEG_INF.

      @param v Vertex
      @param e Edge
      @param a Position of last character of the prefix string
      @param b Position of last character of the suffix string
      @param t Retention time of both prefix and suffix
      @return Score Q[v][e][a][b][t] if present, NEG_INF otherwise
    */
    virtual int get(int v, int e, int a, int b, int t) = 0;
    
    /**
      @brief Set score: Q[v][e][a][b][t] = value. If ti already exists, update existing entry to @em value.

      @param v Vertex
      @param e Edge
      @param a Position of last character of the prefix string
      @param b Position of last character of the suffix string
      @param t Retention time of both prefix and suffix
      @param value Score to be stored
    */
    virtual void set(int v, int e, int a, int b, int t, int value) = 0;
    
    /**
      @param v Vertex
      @param e Edge
      @return Whether Q[v][e][a][b][t] is empty for all retention times t, a and b
    */
    virtual bool is_empty(int v, int e) = 0;

    /**
      @brief Delete all entries of table and free memory.
    */
    virtual void clear() = 0;
};

/** 
  @brief Implementation of matrix class for linear retention time prediction functions. The table has the form
  Q[vertex][edge][time].

  @note Parameters @em a and @em b are not needed and are thus ignored by the implementation
  @note For information about inherited functions refer to parent class @em Matrix
*/
class MatrixLin : public Matrix 
{
  std::vector<std::vector<std::map<int, int> > > Q;
  RetentionTime rt;
  
  public:
    
    void get_rt_dim(int v, int e, std::vector<int> &rt_dim) 
    {
      transform(Q[v][e].begin(), Q[v][e].end(), back_inserter(rt_dim), RetrieveKey());
    }

    // Not used in @em MatrixLin
    void get_rt_a_pairs (int v, int e, std::vector<std::pair <int,int> > &rt_a_pair) {};
      
    int get(int v, int e, int a, int b, int t)
    {
      if (!Q[v][e].count(t))
        return NEG_INF;
      else 
        return Q[v][e][t];
    }
  
    void set(int v, int e, int a, int b, int t, int value)
    {
      Q[v][e][t] = value;
    }
  
    bool is_empty(int v, int e) 
    { 
      return Q[v][e].size() == 0; 
    }

    void clear()
    {
      Q.clear();
      Q.shrink_to_fit();
    }

    /**
      @brief Constructor

      @param vertexCount Number of vertices
      @param edge_count Number of edges
      @param ret Retention time of peptide in question
    */
    MatrixLin (int vertexCount, int edgeCount, RetentionTime ret) 
    {
      rt = ret;
      Q = std::vector<std::vector<std::map<int,int> > >(vertexCount,std::vector<std::map<int, int> >(edgeCount, std::map<int, int>()));
    }
};

/**
  @brief Implementation of matrix class for position dependent retention time prediction functions. The table
  has the form Q[vertex][edge][a][b][time] where a and b are the positions of the last character of the prefix
  and the suffix path (both starting at 'the left') respectively.
  
  At each position Q[v][e] an initially empty list of @em data structs is stored. Recall that a @em data
  struct stores (a, b, time, score). Using a list rather than an array is much more memory efficient since most
  entries would be empty in practice.

  The list at Q[v][e] preserves an ascending order of @em data structs with respect to retention time. Thus,
  binary search is used to find an all entries with time t within the list at Q[v][e]. Then linear search
  is performed among all entries with the same t but distinct (a, b) (usually only a few entries).

  @note For information about inherited functions refer to parent class @em Matrix
*/
class MatrixPos : public Matrix 
{
  
  // Two dimensional table Q[v][e] that stores a list of @em data structs at each position
  std::vector<std::vector<std::vector<data> > > Q;
  std::map<std::pair<int, int>, char> key_map; // Used to store (a, b) keys
  RetentionTime rt;
    
  private:

    /**
      @brief Maps parameter pair @em a, @em b to a 1 byte representation (data struct allows only 8 bytes).

      @param a Position of last character of the prefix string
      @param b Position of last character of the suffix string
      @return key generated for @em a and @em b
    */
    char get_key(int a, int b) 
    { 
      return key_map[std::mp(a, b)];
    }

  public:

    void get_rt_dim(int v, int e, std::vector<int> &rt_dim)
    {
      std::set<int> ts;
      for (int i = 0; i < Q[v][e].size(); i++) 
      {
        ts.insert(Q[v][e][i].t);
      }
      std::copy(ts.begin(), ts.end(), back_inserter(rt_dim)); 
    }

    // Not used in @em MatrixPos
    void get_rt_a_pairs(int v, int e, std::vector<std::pair <int,int> > &rt_a_pair) {}

    int get(int v, int e, int a, int b, int t) 
    {

      int aa = a; int bb = b;
      char key = get_key(aa, bb);

      // Case 1 - no entry with time @em t
      if (!Q[v][e].size() || t > Q[v][e][Q[v][e].size() - 1].t)
      {
        return NEG_INF;
      }

      // Case 2 - Entry with time @em t at end of vector (common enough to be a separate case)
      std::vector<data>::iterator it;
      if (Q[v][e][Q[v][e].size() - 1].t == t) 
      {
        it = Q[v][e].end();
        it--;
        do {
          if (it->get_key() == key) 
          {
            return it->get_score();
          }
          it--;
        } while (it->t == t && it > Q[v][e].begin());
        return NEG_INF;
      }

      // Case 3 - General case
      it = binary_find(Q[v][e].begin(), Q[v][e].end(), data(0, t, key)); // Binary search for entries with @em t
      if (it == Q[v][e].end())
        return NEG_INF; // No entry with time @em t
      
      // Entry with time @em t, now search for entry with correct @em key
      else 
      {
        while (it->t == t && it > Q[v][e].begin())
        {
          it--;
        }
        if (it-> t != t) 
          it++;
        while (it->t == t) 
        {
          if (it->get_key() == key) 
            return it->get_score();
          it++;
        }
        return NEG_INF;
      }
    }

    void set(int v, int e, int a, int b, int t, int value)
    {
      int aa = a; int bb = b;
      char key = get_key(aa, bb);
      std::vector<data>::iterator it;
    
      // Case 1 - no entry has time @em t
      if (!Q[v][e].size() || Q[v][e][Q[v][e].size() - 1].t < t)
        Q[v][e].push_back(data(value, t, key));
    
      // Case 2 - entry with @em t is at last position
      else if (Q[v][e][Q[v][e].size() - 1].t == t) 
      {
        // Check whether one of the entries with time @em t has also the same @em key
        it = Q[v][e].end();
        it--;
        bool found = false;
        do {
          if (it->get_key() == key) {
            it->set_score(value, key); // Found entry with @em key -> update score @em value
            found = true;
            break;
          }
          it--;
        } while (it->t == t && it >= Q[v][e].begin());

        // No entry with @em key -> add new @em data struct
        if (!found) 
          Q[v][e].push_back(data(value,t,key));

      // Case 3 - Unknown whether there is a entry with time @em t
      } else 
      {
        // Binary search for structs with time @em t
        it = binary_find(Q[v][e].begin(),Q[v][e].end(),data(0, t, key));
      
        // No entry found -> insert new value at appropriate location
        if (it == Q[v][e].end()) 
        {
          it--;
          while (it->t > t && it > Q[v][e].begin())
          {
            it--;
          }
          if (it->t < t) 
            it++;
          Q[v][e].insert(it,data(value, t, key));	

        // Entry found with time @em t
        } else 
        {
          // check if one entry has also the same @em key, if so, just update the score, else add new entry
          while (it->t == t && it > Q[v][e].begin()) 
          {
            it--;
          }
          if (it-> t != t) 
            it++;
          bool found = false;
          while (it->t == t) 
          {
            if (it->get_key() == key) 
            {
              it->set_score(value, key); // Update found entry
              found = true;
              break;
            }
            it++;
          }
          if (!found) 
            Q[v][e].insert(it, data(value, t, key)); // Add new entry at appropriate location in list
        }
      }
    }

    bool is_empty(int v, int e) 
    { 
      return Q[v][e].size() == 0; 
    }

    void clear()
    {
      Q.clear();
      Q.shrink_to_fit();
    }

    /**
      @brief Constructor

      @param k Number of vertices
      @param edge_count Number of edges
      @param gamma_l Number of position dependent RT coefficients in the prefix
      @param gamma_r Number of position dependent RT coefficients in the suffix
      @param ret Retention time of peptide in question
    */
    MatrixPos(int k, int edge_count, int gamma_l, int gamma_r, RetentionTime ret) 
    {

      // Initialize @key_map
      char key = 0;		
      for (int a = 0; a <= gamma_l; a++) 
      {
        for (int b = 0; b <= gamma_r; b++) 
        {
          key_map[std::mp(a, b)] = key++;
        }
      }

      rt = ret;    
      Q = std::vector<std::vector<std::vector<data> > > (k, std::vector<std::vector<data> > (edge_count, std::vector<data> ()));
    }	
};

/**
  @brief MatrixNei implements the 4 dimensional dynamic programming table Q[vertex][edge][last_char][time]. 
  
  At each position Q[v][e] an initially empty list of @em data structs is stored. Recall that a @em data
  struct stores (a, b, time, score). In this case we leave the slot for b empty. Using a list rather than an 
  array is much more memory efficient since most entries would be empty in practice.

  The list at Q[v][e] preserves an ascending order of @em data structs with respect to retention time. Thus,
  binary search is used to find an all entries with time t within the list at Q[v][e]. Then linear search
  is perforemd among all entries with the same t but distinct a (usually only a few entries).

  @note For information about inherited functions refer to parent class @em Matrix
*/
class MatrixNei : public Matrix 
{
  std::vector<std::vector<std::vector<data> > > Q; 
  RetentionTime rt;
  
  public:

    void get_rt_dim(int v, int e, std::vector<int> &rt_dim)
    {
      if (Q[v][e].size() < 1) 
        return;
    
      rt_dim.push_back(Q[v][e][0].t);
      for (int i = 0; i < Q[v][e].size(); i++)
      {
        if (Q[v][e][i].t == rt_dim.back())
          continue;
        else 
          rt_dim.push_back(Q[v][e][i].t);
      }  
    }

    void get_rt_a_pairs(int v, int e, std::vector<std::pair <int, int> > &rt_a_pair) 
    {
      for (int i = 0; i < Q[v][e].size(); i++)
      {
        rt_a_pair.push_back(std::make_pair(AA_char2int.at(Q[v][e][i].get_key()), Q[v][e][i].t));
      }
    }
      
    int get(int v, int e, int a, int b, int t) 
    {
      char key = AA_int2char.at(a); // Key used for the @em data struct - does not depend on @em b

      // Case 1 - no entry with time @em t
      if (!Q[v][e].size() || t > Q[v][e][Q[v][e].size() - 1].t)
        return NEG_INF;
    
      // Case 2 - Entry with time @em t at end of vector (common enough to be a separate case)
      std::vector<data>::iterator it;
      if (Q[v][e][Q[v][e].size() - 1].t == t)
      {
        it = Q[v][e].end();
        it--;
        do
        {
          if (it->get_key() == key) 
            return it->get_score();
          it--;
        } while (it->t == t && it >= Q[v][e].begin());
        return NEG_INF;
      }

      // Case 3 - General case
      it = binary_find(Q[v][e].begin(),Q[v][e].end(),data(0, t, key)); // Binary search for entries with @em t
      if (it == Q[v][e].end())
        return NEG_INF; // No entry with time @em t

      // Entry with time @em t, now search for entry with correct @em key
      else
      {
        while (it->t == t && it > Q[v][e].begin()) 
        {
          it--;
        }
        if (it-> t != t)
          it++;
        while (it->t == t) 
        {
          if (it->get_key() == key) 
            return it->get_score();
          it++;
        }
        return NEG_INF;
      }
    }

    void set(int v, int e, int a, int b, int t, int value)
    {
      char key = AA_int2char.at(a);

      std::vector<data>::iterator it;

      // Case 1 - no entry has time @em t
      if (!Q[v][e].size() || Q[v][e][Q[v][e].size() - 1].t < t)
        Q[v][e].push_back(data(value, t, key));
    
      // Case 2 - entry with @em t is at last position
      else if (Q[v][e][Q[v][e].size() - 1].t == t) 
      {
        it = Q[v][e].end(); 
        it--;
        bool found = false;
        do 
        {
          if (it->get_key() == key) 
          {
            it->set_score(value,key);
            found = true;
            break;
          }
          it--;
        } while (it->t == t && it > Q[v][e].begin());
      
        if (!found) 
          Q[v][e].push_back(data(value, t, key));
    
      // Case 3 - Unknown whether there is a entry with time @em t
      } else 
      { 
        // Binary search for structs with time @em t
        it = binary_find(Q[v][e].begin(), Q[v][e].end(),data(0, t, key));
        if (it == Q[v][e].end()) 
        {
          // No entry found with time @em t, insert new value at appropriate location
          it--;
          while (it->t > t && it > Q[v][e].begin()) 
          {
            it--;
          }
          if (it->t < t) 
            it++;
          Q[v][e].insert(it, data(value, t, key));	
      
        // Entry found with time @em t
        } else
        {
          // Check if one entry has also the same @em key, if so, just update the score, else add new entry
          while (it->t == t && it > Q[v][e].begin()) 
          {
            it--;
          }
          if (it-> t != t)
            it++;
          bool found = false;
          while (it->t == t)
          {
            if (it->get_key() == key) 
            {
              it->set_score(value, key); // Update entry
              found = true;
              break;
            }
            it++;
          }
          if (!found) 
            Q[v][e].insert(it,data(value, t, key));	// Add new entry
        }
      }
    }

    bool is_empty(int v, int e) 
    { 
      return Q[v][e].size() == 0; 
    }

    void clear()
    {
      Q.clear();
      Q.shrink_to_fit();
    }

    /**
      @brief Constructor

      @param k Number of vertices
      @param edge_count Number of edges
      @param ret Retention time of peptide in question
    */
    MatrixNei(int k, int edge_count, RetentionTime ret) 
    {
      rt = ret;
      Q = std::vector<std::vector<std::vector<data> > > (k, std::vector<std::vector<data> > (edge_count, std::vector<data> ()));
    } 
};

#endif