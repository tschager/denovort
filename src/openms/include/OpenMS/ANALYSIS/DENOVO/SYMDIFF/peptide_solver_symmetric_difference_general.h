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

#ifndef OPENMS_ANALYSIS_DENOVO_SYMDIFF_PEPTIDESOLVERSYMMETRICDIFFERENCEGENERAL_H
#define OPENMS_ANALYSIS_DENOVO_SYMDIFF_PEPTIDESOLVERSYMMETRICDIFFERENCEGENERAL_H

#include "symdiffscore_multyions.h"
#include "config.h"
#include <iostream>
#include <cmath>
#include <set>

#include "retention_time.h"
#include <cfloat>
#include <algorithm>
#include "matrix.h"
#include <queue>
#include <stack>

/**
  @brief returns true iff the absolute difference between @em a and @em b is less than 10 * DBL_EPSILON.

  @param a First double argument
  @param b Second double argument
*/
bool almost_equal (double a, double b) 
{ 
  return (a == b || (abs(a - b) < 10 * DBL_EPSILON)); 
}

/**
  @brief This method is the core of our symmetric difference algorithm.

  1. From a given set of peaks the matrix spectrum graph is computed. 

  2.2. In a second step, we use dynamic programming with a table @em matrix_l to compute scores of optimal
  path pairs at each (vertex, edge, time) triple, where vertex and edge are the ends of the suffix and the 
  prefix path (or vice versa).

  2.3. Then, we compute a dynamic programming table @em matrix_r where we start at all end-entries of the of
  matrix_l and compute the scores in backwards order. This has the advantage that all path pairs that could
  not be completed to a full path in the matrix spectrum graph do not have an entry in @em matrix_r whereas
  they did in @em matrix_l.

  3. In a last step, we perform a DFS on @em matrix_r to find all peptides with a score at least @em rho times
  as good as the best score.

  @note We may refer to edges and/or vertices of the prefix path as being "on the left side" of the graph
  and edges and/or vertices of the suffix path as being "on the right side" of the graph.

  @param peaks Spectrum containing measurements
  @param last_letter Given last character of peptide
  @param rt Object modeling the retention time of the peptide in question
  @param rho Fraction indicating how high the relative (to the maximal) score of a sequence must be to 
    be included
  @param result Pointer to list where candidate peptides are to be stored in
  @param conf Arguments for algorithm
*/
void find_all_peptides(
  std::vector<peak>& peaks, char last_letter, RetentionTime rt,
  double rho, std::vector<peptide>& result, config conf)
{
  LOG("find_all_peptides start.\n");
  scoreCacheClear();

  int k = peaks.size();
  int edge_count = 0;

  // For all vertices of the graph store their outgoing edges
  std::vector< std::vector<simple_edge> > matrix_edges =
    std::vector< std::vector<simple_edge> >(2 * k, std::vector<simple_edge>() );
  // For all vertices of the graph store their incoming edges
  std::vector< std::vector<simple_edge> > matrix_in_edges = 
    std::vector< std::vector<simple_edge> >(2 * k, std::vector<simple_edge>() );
  std::vector< simple_edge > edges; // List of all edges of the graph

  simple_edge end_edge;
  end_edge.from = 1;
  end_edge.to = 1;
  end_edge.label = "-";
  end_edge.id = edge_count++;
  edges.push_back(end_edge);

  /**
    1. Compute the edges of the spectrum graph
    1.1 Edges on the "left side of the graph" 
      (Starting vertex @em x over all x < k/2, end vertex @em m over all x < m.)
  */
  for (int x = 0; x < k / 2; x++)
  {
    for (int m = k - 1; m > x; m--)
    {
      // List of labels matching the mass m - x
      std::vector<std::string> label = get_label(peaks[m].mz - peaks[x].mz, conf);
      if (label.size() > 0)
      {
        // Create an edge for each label in @em label
        for (int i = 0; i < label.size(); i++)
        {
          std::string l = label[i]; // Single label
          if (l.size() > 1)
          {
            // We have a missing peak (the labels contains more than one character) and have to check if 
            // there is a measured peak at @em peaks[x].mz + @em get_aminoacid_mass(l[1])
            double missing = peaks[x].mz+get_aminoacid_mass(l[1]);
            int l = x; 
            int r = m;
            while (l + 1 < r) //binary search
            {
              int next = (l + r) / 2;
              if (peaks[next].mz < missing)
                l = next;
              else
                r = next;
            }
            // If there is no peak with the missing mass, ignore the edge
            if (fabs(peaks[l].mz - missing) < conf.EPS || fabs(peaks[r].mz - missing) < conf.EPS)
              continue;
          }
          simple_edge e;
          e.from = 2 * x;
          e.to = 2 * m;
          e.label = l;
          e.id = edge_count++;
          matrix_edges[e.from].push_back(e);
          matrix_in_edges[e.to].push_back(e);
          edges.push_back(e);
        }
      }
    }
  }

  /**
    1.2 Edges on the "right side of the graph"
      (Starting vertex @em y over all y > k/2, end vertex over all m < y. Imagine path-pairs where both
      start at 0.)
  */  
  for (int y  = k - 1; y > k / 2; y--)
  {
    for (int m = 0; m < y; m++)
    {
      std::vector<std::string> label = get_label(peaks[y].mz - peaks[m].mz, conf);
      if (label.size() > 0)
      {
        for (int i = 0; i < label.size(); i++)
        {
          std::string l = label[i];
          // Restrict last edge to end with last_letter
          if (y != k - 1 || l[0] == last_letter || (l.size() > 1 && removeBrackets(l).front() == last_letter))
          {
            if (l.size() > 1)
            {
              double missing = peaks[m].mz + get_aminoacid_mass(l[1]);
              int l = m;
              int r = y;
              while (l + 1 < r)
              {
                int next = (l + r) / 2;
                if (peaks[next].mz < missing)
                  l = next;
                else
                  r = next;
              }
              if (fabs(peaks[l].mz - missing) < conf.EPS || fabs(peaks[r].mz - missing) < conf.EPS)
                continue;
            }
            simple_edge e;
            e.from = 2 * (k - 1 - y) + 1;
            e.to = 2 * (k - 1 - m) + 1;
            e.label = l;
            e.id = edge_count++;
            matrix_edges[e.from].push_back(e);
            matrix_in_edges[e.to].push_back(e);
            edges.push_back(e);
          }
        }
      }
    }
  }

  /**
    2. Fill dynamic programming tables
    2.1 Setup
  */

  std::pair<int, int> gamma = rt.get_gamma();
  int gamma_l = gamma.first; // Number of position dependent RT coefficients in the prefix
  int gamma_r = gamma.second; // Number of position dependent RT coefficients in the suffix

  // Initialize matrix based on the retention time model that is used
  Matrix *ml, *mr;
  if (rt.isPosModel())
  {
    ml = new MatrixPos(k, edge_count, gamma_l, gamma_r, rt);
    mr = new MatrixPos(k, edge_count, gamma_l, gamma_r, rt);
  }
  else if (rt.isNeiModel()) 
  {
    ml = new MatrixNei(k, edge_count, rt);
    mr = new MatrixNei(k, edge_count, rt);
    matrix_in_edges[1].push_back(end_edge);
  }
  else
  {
    gamma_l = 0;
    gamma_r = 0;
    ml = new MatrixLin(k, edge_count, rt);
    mr = new MatrixLin(k, edge_count, rt);
  }
  Matrix& matrix_l = *(ml);
  Matrix& matrix_r = *(mr);
  
  // Represents all feasible RT values for a solution
  FeasibleTimes feasible_end = FeasibleTimes(0, rt, conf);
  int feas_rt = feasible_end.get_nr_feasible_rts();

  // To store indices from where we can start backtracking
  int max_score = NEG_INF;	
  std::map<int, std::vector<int> > max_scores;	
  std::map<int, std::vector<int> > max_xs;
  std::map<int, std::vector<int> > max_es;
  std::map<int, std::vector<int> > max_as;
  std::map<int, std::vector<int> > max_bs;
  std::map<int, std::vector<int> > max_ts;

  matrix_l.set(0, 0, 0, 0, 0, conf.START_SCORE);
  int nr_terminals = 0; 
  
  /**
    2.2 Position dependent and linear model
    2.2.1 Fill DP @em matrix_l with scores
  */
  if (rt.isPosModel() || rt.isLinModel())
  {
    for (int x = 0; x < k; x++)
    {
      for (int e = 0; e < edge_count; e++)
      {
        std::vector<int> rt_dim;
        matrix_l.get_rt_dim(x, e, rt_dim);
        for (int a = 0; a <= gamma_l; a++) 
        {
          for (int b = 0; b <= gamma_r; b++) 
          {
            for (std::vector<int>::iterator it = rt_dim.begin(); it != rt_dim.end(); ++it) 
            {
              int t = *it;
              if (matrix_l.get(x, e, a, b, t) == NEG_INF)
                continue;
              
              // Check if end-vertex last edge is terminal (== k - 1)
              else if (edges[e].to / 2 + x / 2 == k - 1)
              {
                // Check whether retention time at terminal vertex is feasible
                if (feasible_end.is_feasible(t)) 
                {
                  int new_score = matrix_l.get(x, e, a, b, t);
                  if (new_score > max_score) 
                     max_score = new_score;
                  if (!max_scores.count(t) || max_scores[t][0] <= new_score) 
                  {
                    // Check whether there are already entries at (x, e, a, b, t) with @em new_score ...
                    if (!max_scores.count(t)) 
                    {
                      // ... no -> add new entry
                      std::vector<int> scores; scores.push_back(new_score);	
                      std::vector<int> xs; xs.push_back(x);
                      std::vector<int> es; es.push_back(e);
                      std::vector<int> as; as.push_back(a);
                      std::vector<int> bs;	bs.push_back(b);
                      std::vector<int> ts; ts.push_back(t);
                      max_scores[t] = scores;
                      max_xs[t] = xs; max_es[t] = es; 
                      max_as[t] = as; max_bs[t] = bs; 
                      max_ts[t] = ts;
                    }
                    else if (max_scores[t][0] < new_score) 
                    {
                      // ... yes but @em new_score is better -> clear list and and add new entry
                      max_scores[t].clear();
                      max_xs[t].clear(); max_es[t].clear(); 
                      max_as[t].clear(); max_bs[t].clear(); 
                      max_ts[t].clear();	
                      max_scores[t].push_back(new_score);
                      max_xs[t].push_back(x); max_es[t].push_back(e); 
                      max_as[t].push_back(a); max_bs[t].push_back(b); 
                      max_ts[t].push_back(t);	
                    }
                    else 
                    {
                      // ... yes but the existing ones have the same score, append new entry to list
                      max_scores[t].push_back(new_score);
                      max_xs[t].push_back(x); max_es[t].push_back(e); 
                      max_as[t].push_back(a); max_bs[t].push_back(b); 
                      max_ts[t].push_back(t);
                    }
                  }
                  nr_terminals++;
                  matrix_r.set(x, e, a, b, t, 0); // Mark entry for backtracking
                }
              } 
              // Vertex is not terminal (< k - 1)
              else 
              {
                for (int i = 0; i < matrix_edges[x].size(); i++) 
                {
                  // Compute new a, b, t for current next edge (x, p)
                  simple_edge out_e = matrix_edges[x][i];
                  int length = out_e.label.size();
                  if (length > 1) 
                    length -= 2; // Subtract '(' and ')'
                  int aa, bb, tt;
                  if (out_e.to % 2) 
                  {
                    // Extension on the "right side"
                    aa = a; 
                    bb = std::min(gamma_r, b + length);
                    tt = t + get_seq_rt(out_e.label, INT_MAX, b);
                  } 
                  else 
                  {
                    // Extension on the "left side"
                    bb = b; 
                    aa = std::min(gamma_l, a + length);
                    tt = t + get_seq_rt(out_e.label, a, INT_MAX);
                  }

                  // Skip if the two paths overlap
                  if (edges[e].to / 2 + out_e.to / 2 > k - 1) 
                    continue;
                  
                  // Check whether the prefix or the suffix path will be longer
                  // Store score in corresponding entry (save edge for the longer of the two)
                  else if (out_e.to > edges[e].to) 
                  {
                    // @em out_e is on longer path
                    int val = std::max(matrix_l.get(edges[e].to, out_e.id, aa, bb, tt),
                      matrix_l.get(x, e, a, b, t) + score(peaks, out_e, edges[e], conf));
                    matrix_l.set(edges[e].to, out_e.id, aa, bb, tt, val);
                  } 
                  else 
                  {
                    // @em edges[e] is on longer path
                    int val = std::max(matrix_l.get(out_e.to, edges[e].id, aa, bb, tt),
                      matrix_l.get(x, e, a, b, t) + score(peaks, out_e, edges[e], conf));
                    matrix_l.set(out_e.to,edges[e].id, aa, bb, tt, val);
                  }
                }
              }
            }
          }
        }
      }
    }

    /**
      2.2.2 Fill DP @em matrix_r with scores in backwards order
      
      @brief We start at entries that correspond to feasible paths in the matrix spectrum graph and compute
      the scores in backwards order. This has the advantage (over @em matrix_l) that there will not be any 
      entries in @em matrix_r that cannot be extended to a feasible path on the matrix spectrum graph.
    */
    for (int x = k - 1; x >= 0; x--)
    {
      for (int e = 0; e < edge_count; e++)
      {
        // If @em matrix_l contains no score at entry (x, e), we can skip this entry
        if (matrix_l.is_empty(x, e))
          continue;
        
        std::vector<int> rt_dim;
        matrix_r.get_rt_dim(x, e, rt_dim);
        for (int a = 0; a <= gamma_l; a++) 
        {
          for (int b = 0; b <= gamma_r; b++) 
          {
            for (std::vector<int>::iterator it = rt_dim.begin(); it != rt_dim.end(); ++it) 
            {
              int t = *it;
              if (matrix_r.get(x, e, a, b, t) > NEG_INF) 
              {
                for (int i = 0; i < matrix_in_edges[x].size(); i++)
                {
                  simple_edge in_e = matrix_in_edges[x][i];
                  simple_edge back_e, next_e;
                  
                  // Check if @em e or @em edges[e] was used for the last extension
                  // Recall that always the shorter of the suffix and the prefix path is extended
                  if (in_e.from < edges[e].from) 
                  {
                    // e edge was the last
                    back_e = edges[e]; 
                    next_e = in_e;
                  }
                  else
                  {
                    // x edge was the last
                    back_e = in_e; 
                    next_e = edges[e];
                  }

                  int length = back_e.label.size();
                  if (length > 1) 
                    length -= 2; // Subtract '(' and ')'
                  int aa, bb, tt;
                  
                  // Check whether the extension was on the "right side"
                  if (back_e.to % 2) 
                  {
                    // Need to try all b, b - 1, b - |back_e|, (if b < gamma_r, b - |back_e| suffices)
                    int start_length = length;
                    if (b == gamma_r) 
                      start_length = 0;
                    for (int j = start_length; j <= length; j++) 
                    {
                      bb = b - j;
                      if (bb < 0) 
                        break;
                      tt = t - get_seq_rt(back_e.label, INT_MAX, bb);

                      // Check whether this is an actual path to follow
                      if (matrix_l.get(back_e.from, next_e.id, a, bb, tt) == NEG_INF) 
                        continue;	
                      int val = std::max(matrix_r.get(back_e.from, next_e.id, a, bb, tt),
                        matrix_r.get(x, e, a, b, t) + score(peaks, back_e, next_e, conf));
                      matrix_r.set(back_e.from, next_e.id, a, bb, tt, val);                  
                    }
                  }
                  // Extension was on the "left side"
                  else 
                  {
                    int start_length = length;
                    // Need to try all a, a - 1, a - |back_e|, (if a < gamma_l, a - |back_e| suffices)
                    if (a == gamma_l) 
                      start_length = 0;
                    for (int j = start_length; j <= length; j++) 
                    {
                      aa = a - j;
                      if (aa < 0) 
                        break;
                      tt = t - get_seq_rt(back_e.label, aa, INT_MAX);

                      // Check whether this is an actual path to follow
                      if (matrix_l.get(back_e.from, next_e.id, aa, b, tt) == NEG_INF) 
                        continue;

                      int val = std::max(matrix_r.get(back_e.from, next_e.id, aa, b, tt),
                        matrix_r.get(x, e, a, b, t) + score(peaks, back_e, next_e, conf));
                      matrix_r.set(back_e.from, next_e.id, aa, b, tt, val);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  /**
    2.3 Neighborhood based model
    2.3.1 Fill matrix_l with scores
  */
  else if (rt.isNeiModel())
  {
    std::string edgelabel;
    for (int x = 0; x < k; x++)
    {
      for (int e = 0; e < edge_count; e++)
      {
        if ( e != 0 && (!(edges[e].from < x) || !(x < edges[e].to)) )
          continue;

        std::vector<int> rt_dim;
        matrix_l.get_rt_dim(x,e,rt_dim);
        for (int a = 0; a <= 20; a++) 
        {
          for (std::vector<int>::iterator it = rt_dim.begin(); it != rt_dim.end(); ++it) 
          {
            int t = *it;
            if (matrix_l.get(x, e, a, 0, t) == NEG_INF)
              continue;
            else if (edges[e].to / 2 + x / 2 == k - 1)
            {
              // Found a terminal vertex
              std::string lastChar;
              if (edges[e].to % 2)
                lastChar = std::string(1, AA_int2char.at(a)) + removeBrackets(edges[e].label).back();
              else
                lastChar = removeBrackets(edges[e].label).back() + std::string(1,AA_int2char.at(a));
              
              int tfinal = t + get_seq_rt_nei(lastChar);
              if (feasible_end.is_feasible(tfinal))
              {
                // Terminal vertex is within the feasible RT interval around measured time
                int new_score = matrix_l.get(x, e, a, 0, t);
                int aa;
                if (edges[e].to % 2)
                  aa = AA_char2int.at(removeBrackets(edges[e].label).front());
                else
                  aa = AA_char2int.at(removeBrackets(edges[e].label).back());
                
                if (new_score > max_score)
                  max_score = new_score;
                if (!max_scores.count(tfinal) || max_scores[tfinal][0] <= new_score) 
                {
                  // Check whether there are already entries at (x, e, a, b, t) with @em new_score ...
                  if (!max_scores.count(tfinal)) 
                  {
                    // ... no -> add new entry
                    std::vector<int> scores; scores.push_back(new_score);	
                    std::vector<int> xs; xs.push_back(x);
                    std::vector<int> es; es.push_back(e);
                    std::vector<int> as; as.push_back(aa);
                    std::vector<int> bs;	bs.push_back(0);
                    std::vector<int> ts; ts.push_back(tfinal);
                    max_scores[tfinal] = scores;
                    max_xs[tfinal] = xs; max_es[tfinal] = es; 
                    max_as[tfinal] = as; max_bs[tfinal] = bs; 
                    max_ts[tfinal] = ts;
                  }
                  else if (max_scores[tfinal][0] < new_score) 
                  {
                    // ... yes but @em new_score is better -> clear list and and add new entry
                    max_scores[tfinal].clear();
                    max_xs[tfinal].clear(); max_es[tfinal].clear(); 
                    max_as[tfinal].clear(); max_bs[tfinal].clear(); 
                    max_ts[tfinal].clear();	
                    max_scores[tfinal].push_back(new_score);
                    max_xs[tfinal].push_back(x); max_es[tfinal].push_back(e); 
                    max_as[tfinal].push_back(aa); max_bs[tfinal].push_back(0); 
                    max_ts[tfinal].push_back(t);	
                  }
                  else 
                  {           
                    // ... yes but the existing ones have the same score, append new entry to list
                    max_scores[tfinal].push_back(new_score);
                    max_xs[tfinal].push_back(x); max_es[tfinal].push_back(e); 
                    max_as[tfinal].push_back(aa); max_bs[tfinal].push_back(0); 
                    max_ts[tfinal].push_back(tfinal);
                  }
                }
                nr_terminals++;
                matrix_r.set(x, e, aa, 0, tfinal, 0); // Set for backtracking
              }
            }
            // Vertex is not terminal 
            else 
            {
              for (int i = 0; i < matrix_edges[x].size(); i++) 
              {
                // Compute new a and t for current next edge (x, p)
                simple_edge out_e = matrix_edges[x][i];
                int length = out_e.label.size();
                if (length > 1) 
                  length -= 2; // Subtract '(' and ')'
                int aa,bb,tt;          
                
                if (out_e.to % 2) 
                {
                  edgelabel = out_e.label;
                  edgelabel = removeBrackets(edgelabel);
                  // Reverse label because its on the "right side" but our path pairs started both at 0
                  std::reverse(edgelabel.begin(), edgelabel.end());
                  
                  // Extension on the "right side side" of the graph
                  if (out_e.to > edges[e].to) 
                    // out_e is on longer path
                    aa = AA_char2int.at(removeBrackets(edges[e].label).back());
                  else
                    // edges[e] is still on longer path
                    aa = AA_char2int.at(removeBrackets(out_e.label).back());
                  edgelabel = edgelabel + std::string(1,AA_int2char.at(a));
                } 
                else 
                {
                  edgelabel = out_e.label;
                  edgelabel = removeBrackets(edgelabel);
                  
                  // Extension on the left side
                  if (out_e.to > edges[e].to) 
                    // @em out_e is on longer path
                    aa = AA_char2int.at(removeBrackets(edges[e].label).back());
                  else
                    // @em edges[e] is still on longer path
                    aa = AA_char2int.at(removeBrackets(out_e.label).back());
                  edgelabel = std::string(1, AA_int2char.at(a)) + edgelabel;
                }
                tt = t + get_seq_rt_nei(edgelabel);  
                
                // Skip if the two paths overlap
                if (edges[e].to / 2 + out_e.to / 2 > k - 1) 
                  continue;

                // Check whether the prefix or the suffix path will be longer
                // Store score in corresponding entry (save edge for the longer of the two)
                else if (out_e.to > edges[e].to) 
                {
                  // @em out_e is on longer path
                  int val = std::max(
                    matrix_l.get(edges[e].to, out_e.id, aa, 0, tt), 
                    matrix_l.get(x, e, a, 0, t) + score(peaks, out_e, edges[e], conf)
                  );
                  matrix_l.set(edges[e].to, out_e.id, aa, 0, tt, val);
                } 
                else 
                {
                  // @em edges[e] is still on longer path
                  int val = std::max(
                    matrix_l.get(out_e.to, edges[e].id, aa, 0, tt), 
                    matrix_l.get(x, e, a, 0, t) + score(peaks, out_e, edges[e], conf)
                  );
                  matrix_l.set(out_e.to, edges[e].id, aa, 0, tt, val);
                }
              }
            }
          }
        }
      }
    }
    
    /**
      2.3.2 Fill DP @em matrix_r with scores in backwards order
      
      @brief We start at entries that correspond to feasible paths in the matrix spectrum graph and compute
      the scores in backwards order. This has the advantage (over @em matrix_l) that there will not be any 
      entries in @em matrix_r that cannot be extended to a feasible path on the matrix spectrum graph.
    */
    std::vector< std::pair<int,int> > rt_a_pairs;
    for (int x = k - 1; x >= 0; x--)
    {
      for (int e = 0; e < edge_count; e++)
      {
        if (e != 0 && (!(edges[e].from < x) || !(x < edges[e].to)))
          continue;
        
        // If @em matrix_l contains no score at entry (x, e), we can skip this entry
        if (matrix_l.is_empty(x,e))
          continue;
        
        rt_a_pairs.clear();
        matrix_r.get_rt_a_pairs(x,e,rt_a_pairs);
        
        int a = 0;
        int t = 0;
        for (auto pair : rt_a_pairs)
        {
          a = pair.first;
          t = pair.second;
            
          if (matrix_r.get(x, e, a, 0, t) > NEG_INF) 
          {
            for (int i = 0; i < matrix_in_edges[x].size(); i++)
            {
              simple_edge in_e = matrix_in_edges[x][i];                
              simple_edge back_e, next_e;
              int aa, bb, tt;
              std::string edgelabel = removeBrackets(in_e.label);
              if (in_e.to % 2)
                // Reverse label because its on the "right side" but our path pairs started both at 0
                std::reverse(edgelabel.begin(), edgelabel.end());
              if (in_e.from < edges[e].from) 
              {
                // @em edges[e] edge was the last
                back_e = edges[e]; 
                next_e = in_e;
                  
                if (back_e.id != 0 && edges[e].to % 2)
                  edgelabel = edgelabel + std::string(1, AA_int2char.at(a));
                else
                  edgelabel = std::string(1, AA_int2char.at(a)) + edgelabel;
              }
              else
              {
                // x edge was the last
                back_e = in_e; 
                next_e = edges[e];
                if (back_e.id != 0 && edges[e].to % 2)
                  edgelabel = edgelabel + std::string(1, AA_int2char.at(a));
                else
                  edgelabel = std::string(1,AA_int2char.at(a)) + edgelabel;
              }
              // If it is a terminal vertex, use the whole label of both edges instead of only the last
              // vertex of back_e
              if (edges[e].to / 2 + x / 2 == k - 1) 
              {
                std::string l;
                if (next_e.to % 2) 
                {
                  l = removeBrackets(next_e.label);
                  std::reverse(l.begin(), l.end());
                  edgelabel = removeBrackets(back_e.label) + l;
                }
                else
                {
                  l = removeBrackets(back_e.label);
                  std::reverse(l.begin(), l.end());
                  edgelabel = removeBrackets(next_e.label) + l;
                }
              }
              tt = t - get_seq_rt_nei(edgelabel);
              aa = AA_char2int.at(removeBrackets(back_e.label).front());                 
              if (back_e.id != 0)
              {
                int val = std::max(matrix_r.get(back_e.from, next_e.id, aa, 0, tt),
                  matrix_r.get(x, e, a, 0, t) + score(peaks, back_e, next_e, conf));
                matrix_r.set(back_e.from, next_e.id, aa, 0, tt, val);
              }
              else
              {
                aa = AA_char2int.at(removeBrackets(next_e.label).front());
                int val = std::max(matrix_r.get(0, 0, aa, 0, tt),
                  matrix_r.get(x, e, a, 0, t) + score(peaks, edges[e], in_e, conf));
                matrix_r.set(0, 0 , aa, 0, tt, val);
              }
            }
          }
        }
      }
    }
  }
  matrix_r.set(0, 0, 0, 0, 0, max_score);
  
  /**
    3. Find all paths with score >= rho*max_score, for this we use DFS
    3.1 Setup - preparations for DFS in @em matrix_r
  */
  std::vector< int > edge_counter = std::vector< int >(k, 0);
  int x = 0;
  int e = 0;
  int a = 0; 
  int b = 0; 
  int t = 0;
  std::vector<simple_edge> edge_stack;
  std::vector<int> scores;
  std::vector<int> e_stack;
  std::vector<int> x_stack;
  std::vector<int> a_stack;
  std::vector<int> b_stack;
  std::vector<int> t_stack;
  scores.push_back(conf.START_SCORE);
  
  /**
    3.2 Linear and position dependent model
    3.2.1 Find all candidate peptides with DFS in @em matrix_r
  */
  if (rt.isPosModel() || rt.isLinModel())
  {
    while (true)
    {
      int r_value = matrix_r.get(x, e, a, b, t);
      // Check if at least one solution is reachable from this position (x, e, a, b, t)
      if (scores.back() + r_value >= rho * max_score && edge_counter[x] < matrix_edges[x].size()) 
      {
        if (edges[e].to / 2 + matrix_edges[x][edge_counter[x]].to / 2 > k - 1)
        {
          edge_counter[x]++;
          continue; // The two paths overlap, abort
        }
        edge_stack.push_back(matrix_edges[x][edge_counter[x]]);
        edge_counter[x]++;
        scores.push_back(scores.back() + score(peaks, edge_stack.back(), edges[e], conf));
        
        // Store current location
        e_stack.push_back(e);
        x_stack.push_back(x);
        a_stack.push_back(a);
        b_stack.push_back(b);
        t_stack.push_back(t);

        int length = edge_stack.back().label.size();
        if (length > 1) // Ignore '(' and ')'
          length -= 2;
        if (edge_stack.back().to < edges[e].to)
          // Short path of the path pair stays short
          x = edge_stack.back().to;
        else
        {
          // Switch - short path of path pair becomes the long one
          x = edges[e].to;
          e = edge_stack.back().id;
        }
        edge_counter[x] = 0;

        if (edge_stack.back().to % 2) 
        {
          // Extension was on the "right side"
          t = t + get_seq_rt(edge_stack.back().label, INT_MAX, b);
          b = std::min(gamma_r, b + length);
        }
        else
        {
          // Extension was on the "left side"
          t = t + get_seq_rt(edge_stack.back().label, a, INT_MAX);
          a = std::min(gamma_l, a + length);	
        }

        // Check whether we have found candidate peptide
        if (edges[e].to / 2 + x / 2 == k - 1 && 
          scores.back() >= rho * max_score && feasible_end.is_feasible(t))
        {
          // Terminal vertex found - write all edge labels currently in the @em edge_stack to @em p
          std::vector<std::string> p;
          for (int i = 0; i < edge_stack.size(); i++)
          {
            simple_edge e = edge_stack[i];
            if (e.from % 2 == 0)
              p.push_back(e.label); 
          }
          for (int i = edge_stack.size() - 1; i >= 0; i--)
          {
            simple_edge e = edge_stack[i];
            if (e.from % 2 != 0)
            {
              std::string l = e.label;
              // Reverse label because its on the "right side" but our path pairs started both at 0
              std::reverse(l.begin(), l.end());
              p.push_back(l);
            }
          }
          peptide pp = std::mp(std::mp(scores.back(), t), p);
          result.push_back(pp);
        }
      } 
      else 
      {
        // No solution reachable from current position
        if (edge_stack.size() == 0)
          break; // We are done
        else
        {
          // Restore old location in order to continue with DFS
          x = x_stack.back(); x_stack.pop_back();
          e = e_stack.back(); e_stack.pop_back();
          a = a_stack.back(); a_stack.pop_back();
          b = b_stack.back(); b_stack.pop_back();
          t = t_stack.back(); t_stack.pop_back();
          scores.pop_back();
          edge_stack.pop_back();
        }
      }
    }
  }

  /**
    3.3 Neighborhood based model
    3.3.1 Find all candidate peptides with DFS in @em matrix_r

    @brief The retention time of double edges has to be stored partially in offset stack; if an edge label 
    has size > 1, the retention time coefficient of all but the first character are stored in the
    corresponding offset as soon as a new edge is added to the same path, the offset is added again in this
    way, the rt of matrix_l and matrix_r are in sync t_stack = matrix_l.time - offsetL - offsetR.
  */
  else if (rt.isNeiModel())
  {
    // Store current location
    e_stack.push_back(e);
    x_stack.push_back(x);
    a_stack.push_back(a);
    t_stack.push_back(t);
    edge_stack.push_back(edges[0]);

    int curscore = scores.back();

    std::vector<int> offset_l;
    std::vector<int> offset_r;
    offset_l.push_back(0);
    offset_r.push_back(0);
    int curoffset_l = 0;
    int curoffset_r = 0;
    
    std::string edgelabel;
    int aa;
    int lastChar;
    int tt;
    while (true)
    {
      if (edge_counter[x] >= matrix_edges[x].size()) 
      {
        x_stack.pop_back(); 
        e_stack.pop_back(); 
        a_stack.pop_back(); 
        t_stack.pop_back(); 
        scores.pop_back();
        offset_l.pop_back();
        offset_r.pop_back();
        edge_stack.pop_back();
        
        // No vertices left
        if (x_stack.size() == 0)
          break;
        
        // Go to parent position
        x = x_stack.back();
        e = e_stack.back();
        a = a_stack.back();
        t = t_stack.back();
        curscore = scores.back();
        curoffset_l = offset_l.back();
        curoffset_r = offset_r.back();
      }
      else
      {
        // Check next edge and increment @em edge_counter
        simple_edge& out_e = matrix_edges[x][edge_counter[x]];
        edge_counter[x]++;
        
        // If the two paths overlap, abort
        if (edges[e].to / 2 + out_e.to / 2 > k - 1)
          continue;

        // Get label and compute @em tt and @em aa
        edgelabel = out_e.label;
        edgelabel = removeBrackets(edgelabel);
        if (out_e.to % 2)
          std::reverse(edgelabel.begin(), edgelabel.end());        
        
        // @em out_e is the last edge connecting the two paths
        if (out_e.to / 2 + edges[e].to / 2 == k - 1) 
          lastChar = AA_char2int.at(edgelabel.front());
        else
        {
          if (out_e.to % 2)
            lastChar = AA_char2int.at(edgelabel.back());
          else
            lastChar = AA_char2int.at(edgelabel.front());
        }
        
        // Check whether the extension was on the "right side" of the graph
        if (out_e.to % 2) 
        {
          if (out_e.to > edges[e].to) 
            // @em out_e is on longer path
            aa = AA_char2int.at(removeBrackets(edges[e].label).back());
          else 
            // @em edges[e] is on longer path            
            aa = AA_char2int.at(removeBrackets(out_e.label).back());
          edgelabel = edgelabel + std::string(1, AA_int2char.at(a)) ;
        } 
        // Extension was on "the left side" of the graph
        else 
        {
          if (out_e.to > edges[e].to) 
            // @em out_e is on longer path
            aa = AA_char2int.at(removeBrackets(edges[e].label).back());
          else 
            // @em edges[e] is on longer path
            aa= AA_char2int.at(removeBrackets(out_e.label).back());
          edgelabel = std::string(1, AA_int2char.at(a)) + edgelabel;
        }
        
        curoffset_l = offset_l.back();
        curoffset_r = offset_r.back();
        tt = t;
        if (out_e.to % 2)
        {
          if (edgelabel.size() <= 2) 
          {
            tt += get_seq_rt_nei(edgelabel) + curoffset_r;
            curoffset_r = 0;
          }
          else 
          {
            tt += get_seq_rt_nei(edgelabel.substr(edgelabel.size()-2)) + curoffset_r;
            curoffset_r = get_seq_rt_nei(edgelabel.substr(0, edgelabel.size() - 1));
          }
          
        }
        else
        {
          if (edgelabel.size() <= 2) 
          {
            tt += get_seq_rt_nei(edgelabel) + curoffset_l;
            curoffset_l = 0;
          }
          else 
          {
            tt += get_seq_rt_nei(edgelabel.substr(0,2)) + curoffset_l;
            curoffset_l = get_seq_rt_nei(edgelabel.substr(1));
          }
        }
        // If @em out_e is the last edge, add all offsets to @em tt
        if (x / 2 + edges[e].to / 2 == k - 1)
          tt += curoffset_l + curoffset_r;
                
        // Check score in @em matrix_r
        int r_value = matrix_r.get(x, e, lastChar, 0, tt);
        
        // Check if we found a new feasible path in the spectrum graph
        if (curscore + r_value >= rho * max_score)
        {
          // Push score of new position and new edge to stack
          curscore += score(peaks, out_e, edges[e], conf);
          
          // Go to new position
          if (out_e.to < edges[e].to)
            // Short path of the path pair stays short
            x = out_e.to;
          else
          {
            // Switch - short path of path pair becomes the long one
            x = edges[e].to;
            e = out_e.id;
          }
          a = aa;
          t = tt;
          edge_counter[x] = 0;

          // Store new position location
          e_stack.push_back(e);
          x_stack.push_back(x);
          a_stack.push_back(a);
          t_stack.push_back(t);
          scores.push_back(curscore);
          offset_l.push_back(curoffset_l);
          offset_r.push_back(curoffset_r);
          edge_stack.push_back(out_e);
          
          // Check whether we found a terminal vertex          
          if (edges[e].to / 2 + x / 2 == k - 1)
          {
            edgelabel = "";
            // Check whether the extension was on "the right side"
            if (edges[e].to % 2)
              edgelabel = std::string(1, AA_int2char.at(a)) + removeBrackets(edges[e].label).back();
            else 
              edgelabel = removeBrackets(edges[e].label).back() + std::string(1, AA_int2char.at(a));
            
            int tfinal = t + get_seq_rt_nei(edgelabel) + curoffset_l + curoffset_r;
            
            // Check whether the retention time is feasible
            if (feasible_end.is_feasible(tfinal))
            {
              std::vector<std::string> p;
              for (int i = 0; i < e_stack.size(); i++)
              {
                simple_edge e = edge_stack.at(i);
                if (e.from % 2 == 0)
                  p.push_back(e.label);
              }
              for (int i = e_stack.size() - 1; i >= 0; i--)
              {
                simple_edge e = edge_stack.at(i);
                if (e.from % 2 != 0 && e.id != 0)
                {
                  std::string l = e.label;
                  // Reverse label because its on the "right side" but our path pairs started both at 0
                  std::reverse(l.begin(),l.end());
                  p.push_back(l);
                }
              }
              peptide pp = std::mp(std::mp(curscore, tfinal), p);
              result.push_back(pp);
            }
          }
        }
      }
    }
  }
  matrix_l.clear();
  matrix_r.clear();
}

#endif
