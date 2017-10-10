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
#include <map>

#include "cbs.hpp"

using namespace Gecode;

using CBSBrancherInt = CBSBrancher<Int::IntView>;

using zVarId = VarId;
using xVarId = VarId;

using IdMap = SharedHashMap<zVarId, std::pair<xVarId, xVarId>>;

std::pair<int, int> vals(int n, int zVal) {
  return std::make_pair(zVal % n + 1, zVal / n + 1 );
};

class OLBrancher : public CBSBrancherInt {
protected:
  const int n;
  ViewArray<Int::IntView> z;
  IdMap zId_to_xs;
public:
  OLBrancher(const Home& home, ViewArray<Int::IntView>& x,
             ViewArray<Int::IntView>& z0,
             IdMap& zId_to_xs0, int n0)
    : CBSBrancherInt(home, x), n(n0), z(z0), zId_to_xs(zId_to_xs0) {}

  OLBrancher(Space& home, bool share, OLBrancher& b)
    : CBSBrancherInt(home, share, b), n(b.n) {
    z.update(home,share,b.z);
    zId_to_xs.update(home,share,b.zId_to_xs);
  }
  Candidate getChoice(Space& home) override {

    std::map<xVarId, std::map<Val, Dens>> densities;
    for_every_log_entry([&](PropId prop_id, SlnCnt slnCnt,
                            VarId var_id, Val val, SlnCnt dens) {

      if (densities.find(var_id) == densities.end())
        densities[var_id] = std::map<Val, Dens>{};

      if (densities[var_id].find(val) == densities[var_id].end()) {
        densities[var_id][val] = dens;
      } else {
        densities[var_id][val] = std::max(densities[var_id][val], dens);
      }
    });

    struct { int pos; Val val; Dens max; } best{-1, 0, 0};

    for (int i=0; i<z.size(); i++) {
      if (z[i].assigned()) continue;
      auto xIds = zId_to_xs[z[i].id()];
      for (Int::ViewValues<Int::IntView> val(z[i]); val(); ++val) {
        auto xVals = vals(n, val.val());
        double dens = 0;
        if (densities.find(xIds.first) != densities.end()) {
          dens = std::max(densities[xIds.first][xVals.first], dens);
        }
        if (densities.find(xIds.second) != densities.end()) {
          dens = std::max(densities[xIds.second][xVals.second], dens);
        }
        if (dens == 0) continue;
        if (dens > best.max)
          best = {i, val.val(), dens};
      }
    }

    assert(best.pos != -1);
    return {(unsigned int)best.pos, best.val};
  }
  ExecStatus commit(Space& home, const Choice& c, unsigned int a) override {
    const auto& pvc = static_cast<const PosValChoice<int>&>(c);
    int pos=pvc.pos().pos, val=pvc.val();
    if (a == 0)
      return me_failed(z[pos].eq(home,val)) ? ES_FAILED : ES_OK;
    else
      return me_failed(z[pos].nq(home,val)) ? ES_FAILED : ES_OK;
  }
  Brancher* copy(Space& home, bool share) override  {
    return new (home) OLBrancher(home,share,*this);
  }
  static void post(Home home, ViewArray<Int::IntView>& x,
                   ViewArray<Int::IntView>& z,
                   IdMap zId_to_xs, int n) {
    (void) new (home) OLBrancher(home, x, z, zId_to_xs, n);
  }
};

void olbrancher(Home home, const IntVarArgs& z, const IntVarArgs& x1,
                const IntVarArgs& x2, int n) {

  ViewArray<Int::IntView> y;
  {
    auto s = x1.size() + x2.size();
    y = ViewArray<Int::IntView>(home, s);

    int i = 0;
    for (int j=0; j<x1.size(); j++) y[i++] = x1[j];
    for (int j=0; j<x2.size(); j++) y[i++] = x2[j];
    assert(i == s);
  }

  ViewArray<Int::IntView> z0(home,z);

  IdMap zId_to_xs;
  {
    zId_to_xs.init();
    auto id = [](const decltype(z[0])& x) { return x.varimp()->id(); };
    for (int i = 0; i < z.size(); i++)
      zId_to_xs[id(z[i])] = std::make_pair(id(x1[i]), id(x2[i]));
  }

  OLBrancher::post(home, y, z0, zId_to_xs, n);
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

