#ifndef __CBS_HPP__
#define __CBS_HPP__

#include <gecode/int.hh>
#include <ext/hash_map>
#include <ext/hash_set>
#include <cstring>
#include <tuple>
#include <functional>
#include <map>
#include <unordered_map>

#include <fstream>
#include <iostream>
#include <set>


#ifdef SQL
#include <sql-interface.hh>
#endif

using namespace Gecode;

using PropId  = unsigned int;
using VarId   = unsigned int;
using Val     = int;
using Dens    = double;
using SlnCnt  = double;

class PropInfo {
public:
  struct Record {
    VarId var_id;
    Val val;
    Dens dens;
  };
private:
  struct Records {
    // During the insertion in *records, we need to keep track of the current
    // position. At the end, curr_pos == records_size
    unsigned int pos{0};

    // Density information for every (var,val) pair
    Record *x{nullptr};

    // Total number of (var,val) pair the propagator shares with the brancher
    // and needs to communicate.
    size_t size{0};
  } records;

  SlnCnt slnCount{0};

  // Sum of all domains of non assigned variable in the propagator. We need this
  // information along with records.size because we can't rely only on
  // records.size to know how many changes occured in the propagator (if we
  // want to trigger recomputation with a changes pourcentage cutoff).
  size_t domAggr{0};

public:
  PropInfo() = default;

  PropInfo(Space& home, size_t domAggr0, size_t domAggrB)
    : domAggr(domAggr0) {
    records.x = home.alloc<Record>(domAggrB);
    records.size = domAggrB;
  }

  PropInfo(Space& home, const PropInfo& o)
    : domAggr(o.domAggr), records(o.records), slnCount(o.slnCount) {
    records.x = home.alloc<Record>(records.size);
    memcpy(records.x, o.records.x, records.pos * sizeof(Record));
  }

  Record& operator[](unsigned int i) {
    assert(records.x != nullptr);
    assert(i < records.size);
    return records.x[i];
  }

  Record operator[](unsigned int i) const {
    return const_cast<PropInfo*>(this)->operator[](i);
  }

  size_t getDomAggr() const { return domAggr; }
  size_t getDomAggrB() const { return records.size; }
  unsigned int getPosRec() const { return records.pos; }

  SlnCnt getSlnCnt() const { return slnCount; }
  void setSlnCnt(SlnCnt cnt) {
    assert(cnt > 1);
    slnCount = cnt;
  }

  void insert_record(const Record&& r) {
    assert(records.x != nullptr);
    assert(r.dens>0 && r.dens<1);
    records.x[records.pos] = r;
    records.pos++;
    assert(records.size >= records.pos);
  }

  void reuse_space(size_t s) {
    // We discard the previous entries by setting the count to 0 (we
    // thus reuse previous allocated memory. The number of records can't
    // grow).
    records.pos = 0;
    records.size = s;
  }
};


typedef __gnu_cxx::hash_map<PropId, PropInfo,
  __gnu_cxx::hash<PropId>, __gnu_cxx::equal_to<PropId>,
  Gecode::space_allocator<PropInfo> >
  LogProp;

class VarIdToPos : public SharedHandle {
protected:
  class VarIdToPosO : public SharedHandle::Object {
  public:
    typedef __gnu_cxx::hash_map<unsigned int, unsigned int> HashMap;
    HashMap hash_map;
  public:
    VarIdToPosO(void) {}
    VarIdToPosO(const VarIdToPosO& o)
      : hash_map(o.hash_map) {}
    virtual Object* copy(void) const {
      return new VarIdToPosO(*this);
    }
    virtual ~VarIdToPosO(void) {}
  };
public:
  VarIdToPos(void) {}
  void init(void) {
    assert(object() == NULL);
    object(new VarIdToPosO());
  }
  bool isIn(unsigned int i) const {
    VarIdToPosO::HashMap *hm = &static_cast<VarIdToPosO*>(object())->hash_map;
    return hm->find(i) != hm->end();
  }
  unsigned int& operator[](unsigned int i) {
    return static_cast<VarIdToPosO*>(object())->hash_map[i];
  }
  unsigned int operator[](unsigned int i) const {
    assert(isIn(i));
    return static_cast<VarIdToPosO*>(object())->hash_map[i];
  }
};

// View array description
struct VADesc {
  template<class View>
  VADesc(const ViewArray<View> x) {
    // The VarIdToPos object is first implicitly constructed with the default
    // constructor and its shared hashmap is not yet allocated. init() thus
    // allocate memory for the hashmap
    positions.init();

    size = x.size();

    // We can assign an index for each variable id
    for (unsigned int i=0; i<x.size(); i++)
      positions[x[i].id()] = i;

    minVal = std::numeric_limits<int>::max();
    int maxVal = std::numeric_limits<int>::min();
    for (unsigned int i=0; i<x.size(); i++) {
      if (x[i].min() < minVal) minVal = x[i].min();
      if (x[i].max() > maxVal) maxVal = x[i].max();
    }
    assert(minVal != INT_MAX && maxVal != INT_MIN);

    width = maxVal - minVal + 1;
    assert(width > 1);
  }
  VADesc(Space& home, bool share, VADesc& xD)
    : size(xD.size), minVal(xD.minVal), width(xD.width) {
    // We tell that we have a subscription to x the new constructed space
    // The hashmap of the VarIdPos object is shared between all spaces. We must
    // tell here that we want to access the same memory that was allocated in
    // the default constructor for the hashmap even if we change space. The
    // only exception is when we use multithreading; the hash map will be copied
    // and shared amongs the spaces in the new thread
    positions.update(home,share,xD.positions);
  }
  // Hash map that assign an index for each variable in x given its id. It is
  // useful if we want to store computations in a continuous memory space
  VarIdToPos positions;
  // Number of variables
  int size;
  // The minimum value in all the variable in x
  int minVal;
  // The width of the union of all variable domains in x
  int width;
};

unsigned int varvalpos(const VADesc& xD, VarId var_id, Val val) {
  return xD.positions[var_id] * xD.width + val - xD.minVal;
}
//
unsigned int varpos(const VADesc& xD, VarId var_id) {
  return xD.positions[var_id];
}

/**
 * \brief Base class for collecting densities from propagators
 *
 * Before beginning to explain what this class does, we must talk about the
 * slndist() method in Gecode::Propagator and the CBS class from which we inherit.
 *
 * Gecode::Propagator is the base class for every constraint in Gecode.
 * The algorithm to compute solution densities for a given (variable,value) pair
 * depends on the constraints in which the variable is involved (its
 * propagators). For this reason, there's a virtual method in Gecode::Propagator
 * that any given concrete propagator can overload to tell how it compute
 * densities for its variables. Here is the signature of the method:
 *
 *   Gecode::Propagator::slndist(Space& home, CBS* densities) const;
 *
 * As the time of writting this comment, the distinct propagator (see
 * gecode/int/distinct.hh) and regular propagator (see
 * gecode/int/extensional.hh) specialise the slndist method to specify its
 * algorithm to compute solutions densities.
 *
 * Of interest to us is the second argument of this method: CBS* densities.
 * The class CBS is an interface (virtual pure class) that contains only one
 * method:
 *
 * virtual void Gecode::CBS::set(unsigned int var_id, int val, double density) = 0;
 *
 * When a propagator overloads its slndist() method, all it knows is he has access
 * to a object CBS with a set() method that enables him to communicate the
 * results of its computation (i.e. densities for each of its
 * (variable,value) pair).
 *
 * When braching, CBSBrancher must call the slndist() method on each active
 * propagator that overloads the method and pass an object that inherits the
 * class CBS to store densities and compute a choice for branching.
 *
 * This class is a specialisation of CBS that we use for caching densities
 * between computation in a Log and reuse them if possible.
 */

class CBSPos : public SharedHandle {
protected:
  class CBSPosO : public SharedHandle::Object {
  public:
    class BoolArray : public SolnDistributionSize {
    private:
      VarId min_id{0};
      VarId max_id{0};
      bool *inBrancher{nullptr};
    private:
      int size() const { return max_id - min_id + 1; }
      bool& get(VarId var_id) { return inBrancher[var_id - min_id]; }
      bool get(VarId var_id) const {
        return const_cast<BoolArray*>(this)->get(var_id);
      }
    public:
//      BoolArray() = default;
      template<class View>
      explicit BoolArray(const ViewArray<View>& x) {
        auto minmax = std::minmax_element(x.begin(), x.end(), [](View a, View b) {
          return a.id() < b.id();
        });
        min_id = minmax.first->id();
        max_id = minmax.second->id();

        inBrancher = heap.alloc<bool>(size());
        for (int i=0; i<size(); i++)
          inBrancher[i] = false;

        for (auto var : x)
          get(var.id()) = true;
      }
      BoolArray(const BoolArray& o)
        : min_id(o.min_id), max_id(o.max_id) {
        inBrancher = heap.alloc<bool>(size());
        for (int i=0; i<size(); i++)
          inBrancher[i] = o.inBrancher[i];
      }
      ~BoolArray() {
        std::cout << "CECI DOIT APPARAITRE" << std::endl;
        heap.free(inBrancher, size());
      }
      virtual bool varInBrancher(VarId var_id) const {
        if (var_id-min_id < 0) return false;
        if (var_id-min_id >= size()) return false;
        return get(var_id);
      }
    };
    BoolArray bool_array;
  public:
//    CBSPosO() = default;
    template<class View>
    explicit CBSPosO(const ViewArray<View>& x) : bool_array(x) {}
    CBSPosO(const CBSPosO& o) : bool_array(o.bool_array) {}
    virtual Object* copy() const { return new CBSPosO(*this); }
    ~CBSPosO() override = default;
  };
public:
  template<class View>
  explicit CBSPos(const ViewArray<View>& x) {
    assert(object() == nullptr);
    object(new CBSPosO(x));
  }
  CBSPosO::BoolArray& getObject() const {
    return static_cast<CBSPosO*>(object())->bool_array;
  }
};




template<class View, typename Derived>
class CBSBrancher : public Brancher, public SolnDistribution {
public:
  // A choice for branching
  struct Candidate {
    unsigned int idx; // Index of the variable in x
    Val val; // Value in the domain of the variable
  };
protected:
  const double recomputation_ratio;
  ViewArray<View> x;
  VADesc xD;
  CBSPos cbs_pos;
//  BranchingHeur<View> heur;
  LogProp logProp;
public:
  CBSBrancher(Home home, ViewArray<View>& x0, double recomputation_ratio0)
    : Brancher(home), recomputation_ratio(recomputation_ratio0),
      x(x0), xD(x), cbs_pos(x0),
      logProp(LogProp::size_type(), LogProp::hasher(),
                LogProp::key_equal(), LogProp::allocator_type(home)) {
    // Because we must call the destructor of aAvgSD
    home.notice(*this,AP_DISPOSE);
  }
  static void post(Home home, ViewArray<View>& x,
                   double recomputation_ratio=1) {
    assert(recomputation_ratio > 0 && recomputation_ratio <= 1);
    (void) new (home) Derived(home,x,recomputation_ratio);
  }
  virtual size_t dispose(Space& home) {
    home.ignore(*this, AP_DISPOSE);
    // ~aAvgSD() calls ~SharedHashMap() which calls ~SharedHashMapObject() to
    // deallocate the hash map when the refcount of SharedHashMapObject is 0
    // TODO: Est-ce que la déallocation fonctionne encore ????
//    heur.~BranchingHeur<View>();
    (void) Brancher::dispose(home);
    return sizeof(*this);
  }
  CBSBrancher(Space& home, bool share, CBSBrancher& b)
    : Brancher(home,share,b), recomputation_ratio(b.recomputation_ratio),
      xD(home,share,b.xD), cbs_pos(b.cbs_pos),
      logProp(b.logProp.begin(), b.logProp.end(),
                 LogProp::size_type(), LogProp::hasher(),
                 LogProp::key_equal(), LogProp::allocator_type(home)) {
    x.update(home,share,b.x);
    for (auto& elem : logProp)
      elem.second = PropInfo(home, b.logProp[elem.first]);
  }
  virtual Brancher* copy(Space& home, bool share) {
//    CBSBrancher *ret = home.alloc<CBSBrancher>(1);
    return new (home) Derived(home,share,*this);
  }
  virtual bool status(const Space& home) const {
    Space& h = const_cast<Space&>(home);
    for (Propagators p(h, PropagatorGroup::all); p(); ++p) {
      // Sum of domains of all variable in propagator
      unsigned int domAggr;
      // Same, but for variables that are also in this brancher.
      unsigned int domAggrB;
      p.propagator().slndistsize(&cbs_pos.getObject(), domAggr, domAggrB);
      // If there's still a propagator that has an unassigned variable that is
      // also in the brancher, we tell our brancher has still work to do.
      if (domAggrB > 0)
        return true;
    }
    return false;
  }

  virtual Candidate getChoice(Space& home) = 0;

  virtual bool compute(VarId var_id) const {
    return true;
  }
  // Method used by all propagators for communicating calculated densities for
  // each of its (variable,value) pair.
  virtual void setMarginalDistribution(PropId prop_id, VarId var_id, Val val,
                                       Dens density) {
    assert(var_id != 0);
    if (!xD.positions.isIn(var_id)) return;
    assert(!x[xD.positions[var_id]].assigned());
    assert(x[xD.positions[var_id]].in(val));
    assert(logProp.find(prop_id) != logProp.end());

    logProp[prop_id].insert_record(PropInfo::Record{var_id, val, density});
  }
  virtual void setSupportSize(PropId prop_id, SlnCnt count) {
    logProp[prop_id].setSlnCnt(count);
  }



  virtual const Choice* choice(Space& home) {
    // Active propagators and the size we need for their log.
    // We considere a propagator as "active" only if
    // - it supports slndist()
    // - it has unassigned variables that are also in the brancher
    struct Psize { size_t domAggr; size_t domAggrB; };
    __gnu_cxx::hash_map<PropId, Psize > activeProps;
    for (Propagators p(home, PropagatorGroup::all); p(); ++p) {
      unsigned int domAggr, domAggrB;
      p.propagator().slndistsize(&cbs_pos.getObject(), domAggr, domAggrB);
      if (domAggrB != 0) {
        activeProps[p.propagator().id()] = {domAggr, domAggrB};
      }
    }

    // We delete log elements corresponding to non active propagators
    {
      std::vector<PropId> propsToDelete;
      for(const auto& kv : logProp)
        if (activeProps.find(kv.first) == activeProps.end())
          propsToDelete.push_back(kv.first);

      for (auto prop_id : propsToDelete)
        logProp.erase(prop_id);
    }

    for (Propagators p(home, PropagatorGroup::all); p(); ++p) {
      auto prop_id = p.propagator().id();
      if (activeProps.find(prop_id) == activeProps.end()) continue;
      auto aProp = &activeProps[prop_id];

      // Propagator already in the log?
      bool in_log = logProp.find(prop_id) != logProp.end();
      bool changed = true;

      if (in_log) {
        auto prop = &logProp[prop_id];
        changed = prop->getDomAggr()*recomputation_ratio > aProp->domAggr;
        if (changed)
          prop->reuse_space(aProp->domAggrB);
      } else {
        // We create a new propagator
        logProp[prop_id] = PropInfo(home, aProp->domAggr, aProp->domAggrB);
      }

      if (!in_log || changed) {
//        p.propagator().slndist(home,&heur,SolnDistribution::MAX_PER_PROP);
        p.propagator().slndist(home,&*this);
      }
    }


    // We find the choice.
    Candidate c = getChoice(home);
    assert(!x[c.idx].assigned());
    assert(x[c.idx].in(c.val));

    return new PosValChoice<int>(*this,2,c.idx,c.val);
  }
  virtual const Choice* choice(const Space&, Archive& e) {
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
public:
  void for_every_log_entry(const std::function<void(PropId,SlnCnt,VarId,Val,Dens)>& f) {
    for (const auto& elem : logProp) {
      auto prop_id = elem.first;
      auto prop = &elem.second;
      for (int i=0; i<prop->getPosRec(); i++) {
        auto r = (*prop)[i];
        f(prop_id, prop->getSlnCnt(), r.var_id, r.val, r.dens);
      }
    }
  }
};

template<class View>
class MAXSD : public CBSBrancher<View,MAXSD<View>> {
  using CBSBrancher<View,MAXSD<View>>::x;
  using CBSBrancher<View,MAXSD<View>>::xD;
  using CBSBrancher<View,MAXSD<View>>::for_every_log_entry;
  typedef typename CBSBrancher<View,MAXSD<View>>::Candidate Candidate;
public:
  MAXSD(const Home& home, ViewArray<View>& x0, double recomputation_ratio0)
    : CBSBrancher<View,MAXSD<View>>(home, x0, recomputation_ratio0) {}

  MAXSD(Space& home, bool share, CBSBrancher<View,MAXSD<View>>& b)
    : CBSBrancher<View,MAXSD<View>>(home, share, b) {}

  Candidate getChoice(Space& home) override {
    PropInfo::Record best{0,0,0};
    for_every_log_entry([&](PropId prop_id, SlnCnt slnCnt,
                            VarId var_id, Val val, SlnCnt dens) {
      unsigned int pos = varpos(xD, var_id);
      if (!x[pos].assigned() && x[pos].in(val))
        if (dens > best.dens || (dens == best.dens && var_id < best.var_id))
          best = {var_id, val, dens};
    });
    assert(best.var_id != 0);
    return {xD.positions[best.var_id],best.val};
  }
};

enum CBSBranchingHeuristic {
  MAX_SD,
  MAX_REL_SD,
  MAX_REL_RATIO,
  A_AVG_SD,
  W_SC_AVG,
  AI
};


template<class View, class T>
void _cbsbranch(Home home, const T& x, CBSBranchingHeuristic s) {
  if (home.failed()) return;
  ViewArray<View> y(home,x);
  switch(s) {
    case MAX_SD:
      MAXSD<View>::post(home,y,1);
      break;
    default:
      GECODE_NEVER;
  }
}

void cbsbranch(Home home, const IntVarArgs& x, CBSBranchingHeuristic s) {
  _cbsbranch<Int::IntView,IntVarArgs>(home,x,s);
}

void cbsbranch(Home home, const BoolVarArgs& x, CBSBranchingHeuristic s) {
  _cbsbranch<Int::BoolView,BoolVarArgs>(home,x,s);
}



#endif //__CBS_HPP__
