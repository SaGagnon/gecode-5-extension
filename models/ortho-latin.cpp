/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*
 *  Main authors:
 *     Christian Schulte <schulte@gecode.org>
 *
 *  Copyright:
 *     Christian Schulte, 2004
 *
 *  Last modified:
 *     $Date: 2015-09-11 16:29:45 +0200 (Fri, 11 Sep 2015) $ by $Author: schulte $
 *     $Revision: 14672 $
 *
 *  This file is part of Gecode, the generic constraint
 *  development environment:
 *     http://www.gecode.org
 *
 *  Permission is hereby granted, free of charge, to any person obtaining
 *  a copy of this software and associated documentation files (the
 *  "Software"), to deal in the Software without restriction, including
 *  without limitation the rights to use, copy, modify, merge, publish,
 *  distribute, sublicense, and/or sell copies of the Software, and to
 *  permit persons to whom the Software is furnished to do so, subject to
 *  the following conditions:
 *
 *  The above copyright notice and this permission notice shall be
 *  included in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#include <gecode/driver.hh>
#include <gecode/int.hh>

#include "cbs.hpp"

using namespace Gecode;


//template<class View>
//ViewArray<View>& combine(const Home& home, const ViewArray<View>& a,
//                         const ViewArray<View>& b) {
//  ViewArray<View> c(a.size() + b.size());
//  for (int i=0; i<a.size(); i++)
//    c[i] = a[i];
//  for (int i=0; i<b.size(); i++)
//    c[i+a.size()] = b[i];
//
//  return std::move(c);
//}

using HashMap = SharedHashMap<VarId, VarId>;


using CBSBrancherInt = CBSBrancher<Int::IntView>;

class OLBrancher : public CBSBrancherInt {
protected:
  const int n;
  HashMap x1_to_z;
  HashMap x2_to_z;
public:
  OLBrancher(const Home& home, ViewArray<Int::IntView>& x,
             HashMap& x1_to_z0, HashMap& x2_to_z0, int n0)
    : CBSBrancherInt(home, x), n(n0), x1_to_z(x1_to_z0), x2_to_z(x2_to_z0) {}

  OLBrancher(Space& home, bool share, OLBrancher& b)
    : CBSBrancherInt(home, share, b), n(b.n) {
    x1_to_z.update(home,share,b.x1_to_z);
    x1_to_z.update(home,share,b.x2_to_z);
  }
  Candidate getChoice(Space& home) override {
      struct Best {
        Val val;
        Dens dens;
      };

    using ZVarId = VarId;
    __gnu_cxx::hash_map<ZVarId, std::pair<Best,Best>> aggr;

    for_every_log_entry([&](PropId prop_id, SlnCnt slnCnt,
                            VarId var_id, Val val, SlnCnt dens) {
      assert(!x1_to_z.isIn(var_id) || !x2_to_z.isIn(var_id));
      if (x1_to_z.isIn(var_id)) {
        auto zVarId = x1_to_z[var_id];

        if (aggr.find(zVarId) == aggr.end()) {
          auto b = Best{-1, 0};
          aggr[zVarId] = std::make_pair(b, b);
        }

        if (aggr[zVarId].first.dens < dens)
          aggr[zVarId].first = Best{val,dens};
      } else if (x2_to_z.isIn(var_id)) {
        auto zVarId = x2_to_z[var_id];

        if (aggr.find(zVarId) == aggr.end()) {
          auto b = Best{-1, 0};
          aggr[zVarId] = std::make_pair(b, b);
        }

        if (aggr[zVarId].second.dens < dens)
          aggr[zVarId].second = Best{val,dens};
      }
    });

    using E = decltype(aggr)::value_type;
    auto best =
      *std::max_element(aggr.begin(), aggr.end(), [](const E& a,
                                                     const E& b) {

        return a.second.first.dens + a.second.second.dens <
          b.second.first.dens + b.second.second.dens;
    });

    auto zVarIdx = best.first;
    auto v1 = best.second.first.val;
    auto v2 = best.second.second.val;

    auto zVal = (v2-1)*n + v1-1;
    return Candidate{varpos[zVarIdx], zVal};
  }
  Brancher* copy(Space& home, bool share) override  {
    return new (home) OLBrancher(home,share,*this);
  }
  static void post(Home home, ViewArray<Int::IntView>& x,
                   HashMap& x1_to_z, HashMap& x2_to_z, int n) {
    (void) new (home) OLBrancher(home, x, x1_to_z, x2_to_z, n);
  }
};

void olbrancher(Home home, const IntVarArgs& z, const IntVarArgs& x1,
                const IntVarArgs& x2, int n) {
  auto s = (int)z.size() + x1.size() + x2.size();
  ViewArray<Int::IntView> y(home, s);

  {
    int i = 0;
    for (int j=0; j< z.size(); j++) y[i++] =  z[j];
    for (int j=0; j<x1.size(); j++) y[i++] = x1[j];
    for (int j=0; j<x2.size(); j++) y[i++] = x2[j];
  }


  HashMap x1_to_z;
  x1_to_z.init();
  HashMap x2_to_z;
  x2_to_z.init();

  for (int i=0; i<z.size(); i++) {
    x1_to_z[x1[i].varimp()->id()] = z[i].varimp()->id();
    x2_to_z[x2[i].varimp()->id()] = z[i].varimp()->id();
  }

  OLBrancher::post(home, y, x1_to_z, x2_to_z, n);
}

/**
 * \brief %Example: Orthogonal latin squares
 *
 * \ingroup Example
 */
class OrthoLatinSquare : public Script {
protected:
  /// Size of squares
  const int n;
  /// Fields of first square
  IntVarArray x1;
  /// Fields of second square
  IntVarArray x2;

public:
  enum {
    BRANCH_NONE,
    BRANCH_CBS_MAX_SD,
    BRANCH_CBS_A_AVG_SD
  };
  /// Access field at position \a i and \a j in first square
  IntVar& y1(int i, int j) {
    return x1[i*n+j];
  }
  /// Access field at position \a i and \a j in first square
  const IntVar& y1(int i, int j) const {
    return x1[i*n+j];
  }
  /// Access field at position \a i and \a j in second square
  IntVar& y2(int i, int j) {
    return x2[i*n+j];
  }
  /// Access field at position \a i and \a j in second square
  const IntVar& y2(int i, int j) const {
    return x2[i*n+j];
  }

  /// Actual model
  OrthoLatinSquare(const SizeOptions& opt)
    : Script(opt),
      n(opt.size()),
      x1(*this,n*n,1,n), x2(*this,n*n,1,n) {
    const int nn = n*n;
    IntVarArgs z(*this,nn,0,n*n-1);

    distinct(*this, z, opt.ipl());
    // Connect
    {
      IntArgs mod(n*n);
      IntArgs div(n*n);
      for (int i=0; i<n; i++)
        for (int j=0; j<n; j++) {
          mod[i*n+j] = j+1;
          div[i*n+j] = i+1;
        }
      for (int i = nn; i--; ) {
        element(*this, div, z[i], x2[i]);
        element(*this, mod, z[i], x1[i]);
      }
    }

    // Rows
    for (int i = n; i--; ) {
      IntVarArgs ry(n);
      for (int j = n; j--; )
        ry[j] = y1(i,j);
      distinct(*this, ry, opt.ipl());
      for (int j = n; j--; )
        ry[j] = y2(i,j);
      distinct(*this, ry, opt.ipl());
    }
    for (int j = n; j--; ) {
      IntVarArgs cy(n);
      for (int i = n; i--; )
        cy[i] = y1(i,j);
      distinct(*this, cy, opt.ipl());
      for (int i = n; i--; )
        cy[i] = y2(i,j);
      distinct(*this, cy, opt.ipl());
    }

    for (int i = 1; i<n; i++) {
      IntVarArgs ry1(n);
      IntVarArgs ry2(n);
      for (int j = n; j--; ) {
        ry1[j] = y1(i-1,j);
        ry2[j] = y2(i,j);
      }
      rel(*this, ry1, IRT_GQ, ry2);
    }

    if (opt.branching() == BRANCH_CBS_MAX_SD) {
      olbrancher(*this,z,x1,x2,n);
    } else if (opt.branching() == BRANCH_CBS_A_AVG_SD) {
      cbsbranch(*this, z, CBSBranchingHeuristic::A_AVG_SD);
    }

    branch(*this, z, INT_VAR_SIZE_MIN(), INT_VAL_SPLIT_MIN());
  }

  /// Constructor for cloning \a s
  OrthoLatinSquare(bool share, OrthoLatinSquare& s)
    : Script(share,s), n(s.n) {
      x1.update(*this, share, s.x1);
      x2.update(*this, share, s.x2);
  }

  /// Copy during cloning
  virtual Space*
  copy(bool share) {
    return new OrthoLatinSquare(share,*this);
  }
  /// Print solution
  virtual void
  print(std::ostream& os) const {
    for (int i = 0; i<n; i++) {
      os << "\t";
      for (int j = 0; j<n; j++) {
        os.width(2);
        os << y1(i,j) << "  ";
      }
      os << std::endl;
    }
    os << std::endl;
    for (int i = 0; i<n; i++) {
      os << "\t";
      for (int j = 0; j<n; j++) {
        os.width(2);
        os << y2(i,j) << "  ";
      }
      os << std::endl;
    }
    os << std::endl;
  }

};

/**
 * \brief Main function
 * \relates OrthoLatinSquare
 */
int
main(int argc, char* argv[]) {
  SizeOptions opt("OrthoLatinSquare");
  opt.size(7);
  opt.ipl(IPL_DOM);


  opt.branching(OrthoLatinSquare::BRANCH_NONE, "none",
                " ");
  opt.branching(OrthoLatinSquare::BRANCH_CBS_MAX_SD,
                "cbs_max_sd", "maxSD counting base search");
  opt.branching(OrthoLatinSquare::BRANCH_CBS_A_AVG_SD,
                "cbs_a_avg_sd", "aAvgSD counting base search");

  opt.parse(argc,argv);
  Script::run<OrthoLatinSquare,DFS,SizeOptions>(opt);
  return 0;
}

// STATISTICS: example-any

