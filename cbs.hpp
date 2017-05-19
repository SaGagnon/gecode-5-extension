#ifndef __CBS_HPP__
#define __CBS_HPP__

#include <gecode/int.hh>
#include <ext/hash_map>

using namespace Gecode;

typedef __gnu_cxx::hash_map<unsigned int, unsigned int> HashMap;

class SharedHashMap : public SharedHandle {
protected:
  class SharedHashMapObject : public SharedHandle::Object {
  public:
    HashMap hash_map;
  public:
    SharedHashMapObject(void) {std::cout << "1: only one time" << std::endl;}
    SharedHashMapObject(const SharedHashMapObject& shmo)
      : hash_map(shmo.hash_map) {std::cout << "2: never" << std::endl;}
    virtual Object* copy(void) const {
      std::cout << "3: never" << std::endl;
      return new SharedHashMapObject(*this);
    }
    virtual ~SharedHashMapObject(void) {}
  };
public:
  SharedHashMap(void) {}
  SharedHashMap(const SharedHashMap& shm)
    : SharedHandle(shm) {}
  void init(void) {
    assert(object() == NULL);
    std::cout << "4: only one time" << std::endl;
    object(new SharedHashMapObject());
  }
  HashMap* get(void) const {
    return &static_cast<SharedHashMapObject*>(object())->hash_map;
  }
  // some inherited members
  void update(Space& home, bool share, SharedHandle& sh) {
    SharedHandle::update(home,share,sh);
  }
};

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

template<class View>
class aAvgSD : public BranchingHeuristic {
private:
  SharedHashMap positions;
public:
  aAvgSD(Space& home, const ViewArray<View>& x) {
    assert(densities_sum == NULL);
    assert(count == NULL);
    positions.init();

    minVal = INT_MAX;
    int maxVal = INT_MIN;
    for (unsigned int i=0; i<x.size(); i++) {
      (*positions.get())[x[i].id()] = i;
      if (x[i].min() < minVal) minVal = x[i].min();
      if (x[i].max() > maxVal) maxVal = x[i].max();
    }
    assert(minVal != INT_MAX && maxVal != INT_MIN);

    width = maxVal - minVal + 1;
    assert(width > 1);

    int size = x.size() * width;
    densities_sum = home.alloc<double>(size);
    count = home.alloc<int>(size);
    for (int i=0; i<size; i++) {
      densities_sum[i] = 0;
      count[i] = 0;
    }
  }
  aAvgSD(Space& home, const aAvgSD& a)
    : positions(a.positions), minVal(a.minVal), width(a.width) {
    int size = (unsigned int)(*positions.get()).size() * width;
    densities_sum = home.alloc<double>(size);
    count = home.alloc<int>(size);
    memcpy(densities_sum, a.densities_sum, size * sizeof(double));
    memcpy(count, a.count, size * sizeof(int));
  }
  void clear() {
    best_candidate.c.var_id = -1;
    best_candidate.density_moy = 0;
    for (int i=0; i<(*positions.get()).size()*width; i++) {
      densities_sum[i] = 0;
      count[i] = 0;
    }

  }
  virtual void set(unsigned int var_id, int val, double density) {
    assert(densities_sum != NULL);
    unsigned int i = (*positions.get())[var_id] * width + val - minVal;
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
  ~aAvgSD() {
    positions.~SharedHashMap();
  }
private:
  int minVal;
  int width;
  double *densities_sum;
  int *count;
private:
  struct {
    Candidate c;
    double density_moy;
  } best_candidate;
};


template<class View>
class CBSBrancher : public Brancher {
protected:
  ViewArray<View> x;
  aAvgSD<View> heur;
public:
  CBSBrancher(Home home, ViewArray<View>& x0)
    : Brancher(home), x(x0), heur(home,x0) {
    // Because we use heur and must the desctructor of the SharedHashMap
    home.notice(*this,AP_DISPOSE);
  }
  static void post(Home home, ViewArray<View>& x) {
    (void) new (home) CBSBrancher(home,x);
  }
  virtual size_t dispose(Space& home) {
    home.ignore(*this, AP_DISPOSE);
    heur.~aAvgSD();
    (void) Brancher::dispose(home);
    return sizeof(*this);
  }
  CBSBrancher(Space& home, bool share, CBSBrancher& b)
    : Brancher(home,share,b), heur(home, b.heur) {
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
    heur.clear();
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
  CBSBrancher<Int::IntView>::post(home,y);
}

void cbsbranch(Home home, const BoolVarArgs& x) {
  if (home.failed()) return;
  ViewArray<Int::BoolView> y(home,x);
  CBSBrancher<Int::BoolView>::post(home,y);
}

#endif //__CBS_HPP__
