#ifndef __CBS_HPP__
#define __CBS_HPP__

#include <gecode/int.hh>
#include <ext/hash_map>

using namespace Gecode;

class BranchingHeuristic : public CBS {
public:
  struct Candidate {
    int var_id;
    int val;
  };
  // TODO: Quoi faire si multithread? Est-ce que je fais juste clear
  // mon data global....????
//  static virtual void dispose() = 0;
  virtual Candidate get(void) = 0;
};

class MaxSD : public BranchingHeuristic {
public:
  MaxSD() {
    best_candidate.c.var_id = -1;
    best_candidate.density = 0;
  }
  virtual void set(unsigned int var_id, int val, double density) {
    if (density > best_candidate.density) {
      best_candidate.c.var_id = var_id;
      best_candidate.c.val = val;
      best_candidate.density = density;
    }
  }
  virtual Candidate get(void) {
    return best_candidate.c;
  }
private:
  struct {
    Candidate c;
    double density;
  } best_candidate;
};

class aAvgSD : public BranchingHeuristic {
public:
  aAvgSD(ViewArray<Int::IntView>& x) {
    if (densities_sum == NULL) {
      minVal = INT_MAX;
      int maxVal = INT_MIN;
      for (unsigned int i=0; i<x.size(); i++) {
        pos[x[i].id()] = i;
        if (x[i].min() < minVal) minVal = x[i].min();
        if (x[i].max() > maxVal) maxVal = x[i].max();
      }
      assert(minVal != INT_MAX && maxVal != INT_MIN);

      width = maxVal - minVal + 1;
      assert(width > 1);

      int size = x.size() * width;
      densities_sum = heap.alloc<double>(size);
      count = heap.alloc<int>(size);
    }

    best_candidate.c.var_id = -1;
    best_candidate.density_moy = 0;

    for (int i=0; i<pos.size()*width; i++) {
      densities_sum[i] = 0;
      count[i] = 0;
    }
  }
  virtual void set(unsigned int var_id, int val, double density) {
    assert(densities_sum != NULL);
    unsigned int i = pos[var_id] * width + val - minVal;
    densities_sum[i] += density;
    count[i] += 1;

    double density_moy = densities_sum[i] / (double)count[i];
    if (density_moy > best_candidate.density_moy) {
      best_candidate.c.var_id = var_id;
      best_candidate.c.val = val;
      best_candidate.density_moy = density_moy;
    }
  }
  virtual Candidate get(void) {
    return best_candidate.c;
  }
private:
  static int minVal;
  static int width;
  static double *densities_sum;
  static int *count;
  static __gnu_cxx::hash_map<unsigned int, unsigned int> pos;
private:
  struct {
    Candidate c;
    double density_moy;
  } best_candidate;
};

int aAvgSD::minVal = INT_MAX;
int aAvgSD::width = -1;
double *aAvgSD::densities_sum = NULL;
int *aAvgSD::count = NULL;
__gnu_cxx::hash_map<unsigned int, unsigned int> aAvgSD::pos;

class CBSBrancher : public Brancher {
protected:
  ViewArray<Int::IntView> x;
public:
  CBSBrancher(Home home, ViewArray<Int::IntView>& x0)
    : Brancher(home), x(x0) {}
  static void post(Home home, ViewArray<Int::IntView>& x) {
    (void) new (home) CBSBrancher(home,x);
  }
  virtual size_t dispose(Space& home) {
    (void) Brancher::dispose(home);
    return sizeof(*this);
  }
  CBSBrancher(Space& home, bool share, CBSBrancher& b)
    : Brancher(home,share,b) {
    x.update(home,share,b.x);
  }
  virtual Brancher* copy(Space& home, bool share) {
    return new (home) CBSBrancher(home,share,*this);
  }
  // status
  virtual bool status(const Space& home) const {
    Space& h = const_cast<Space&>(home);
    for (Propagators p(h, PropagatorGroup::all); p(); ++p)
      if (p.propagator().cbs(h, NULL))
        return true;
    return false;
  }
  // choice
  virtual Choice* choice(Space& home) {
    aAvgSD heur(x);
    for (Propagators p(home, PropagatorGroup::all); p(); ++p)
      p.propagator().cbs(home, &heur);

    for (int i=0; i<x.size(); i++)
      if (x[i].id() == heur.get().var_id)
        return new PosValChoice<int>(*this,2,i,heur.get().val);

    GECODE_NEVER;
    return NULL;
  }
  virtual Choice* choice(const Space&, Archive& e) {
    int pos, val;
    e >> pos >> val;
    return new PosValChoice<int>(*this,2,pos,val);
  }
  // commit
  virtual ExecStatus commit(Space& home, const Choice& c, unsigned int a) {
    const PosValChoice<int>& pvc = static_cast<const PosValChoice<int>&>(c);
    int pos=pvc.pos().pos, val=pvc.val();
    if (a == 0)
      return me_failed(x[pos].eq(home,val)) ? ES_FAILED : ES_OK;
    else
      return me_failed(x[pos].nq(home,val)) ? ES_FAILED : ES_OK;
  }
  // print
  virtual void print(const Space& home, const Choice& c, unsigned int a,
                     std::ostream& o) const {
    const PosValChoice<int>& pvc = static_cast<const PosValChoice<int>&>(c);
    int pos=pvc.pos().pos, val=pvc.val();
    if (a == 0)
      o << "x[" << pos << "] = " << val;
    else
      o << "x[" << pos << "] != " << val;
  }
};

enum CBSStrategy {
 MAX_SD,
 A_AVG_SD
};

void cbsbranch(Home home, const IntVarArgs& x) {
  if (home.failed()) return;
  ViewArray<Int::IntView> y(home,x);
  CBSBrancher::post(home,y);
}

#endif //__CBS_HPP__
