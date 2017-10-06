#ifndef __CBS_HPP__
#define __CBS_HPP__

#include <iostream>
#include <ext/hash_map>

#include <gecode/int.hh>

#ifdef SQL
#include <sql-interface.hh>
#endif

using namespace Gecode;

using PropId  = unsigned int;
using VarId   = unsigned int;
using Val     = int;
using Dens    = double;
using SlnCnt  = double;

/**
 * \brief Density records for a propagator
 */
class PropInfo {
public:
  struct Record {
    VarId var_id;
    Val val;
    Dens dens;
  };
private:
  struct Records {
    /// Position for inserting the next record
    unsigned int pos{0};
    /// List of records
    Record *x{nullptr};
    /// Number of records allocated
    size_t size{0};
  } records;

  SlnCnt slnCount{0};
  /// Sum of all domains of non assigned variable
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
  void reuse_mem(size_t s) {
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

/**
 * \brief Maps variables ids to indexes
 */
class VarIdToPos : public SharedHandle {
protected:
  class VarIdToPosO : public SharedHandle::Object {
  public:
    typedef __gnu_cxx::hash_map<unsigned int, unsigned int> HashMap;
    HashMap hash_map;
  public:
    VarIdToPosO() = default;
    VarIdToPosO(const VarIdToPosO& o)
      : hash_map(o.hash_map) {}
    Object* copy() const override {
      return new VarIdToPosO(*this);
    }
    ~VarIdToPosO() override = default;
  };
public:
  VarIdToPos() = default;
  void init() {
    assert(object() == nullptr);
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

class VarInBrancher : public SharedHandle {
protected:
  class VarInBrancherO : public SharedHandle::Object {
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
      bool varInBrancher(VarId var_id) const override {
        if (var_id-min_id < 0) return false;
        if (var_id-min_id >= size()) return false;
        return get(var_id);
      }
    };
    BoolArray bool_array;
  public:
    template<class View>
    explicit VarInBrancherO(const ViewArray<View>& x) : bool_array(x) {}
    VarInBrancherO(const VarInBrancherO& o) : bool_array(o.bool_array) {}
    Object* copy() const override { return new VarInBrancherO(*this); }
    ~VarInBrancherO() override = default;
  };
public:
  template<class View>
  explicit VarInBrancher(const ViewArray<View>& x) {
    assert(object() == nullptr);
    object(new VarInBrancherO(x));
  }
  VarInBrancherO::BoolArray& getObject() const {
    return static_cast<VarInBrancherO*>(object())->bool_array;
  }
};

template<class View>
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
  VarIdToPos varpos;
  VarInBrancher varInBrancher;
  LogProp logProp;
public:
  CBSBrancher(Home home, ViewArray<View>& x0, double recomputation_ratio0=1)
    : Brancher(home), recomputation_ratio(recomputation_ratio0),
      x(x0), varInBrancher(x0),
      logProp(LogProp::size_type(), LogProp::hasher(),
              LogProp::key_equal(), LogProp::allocator_type(home)) {
    assert(recomputation_ratio0 > 0 && recomputation_ratio0 <= 1);
    // Because we must call the destructor of aAvgSD
    home.notice(*this,AP_DISPOSE);
    // The VarIdToPos object is first implicitly constructed with the default
    // constructor and its shared hashmap is not yet allocated. init() thus
    // allocate memory for the hashmap
    varpos.init();

    // We assign an index for each variable id
    for (unsigned int i=0; i<x.size(); i++)
      varpos[x[i].id()] = i;
  }
  size_t dispose(Space& home) override {
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
      varInBrancher(b.varInBrancher),
      logProp(b.logProp.begin(), b.logProp.end(),
                 LogProp::size_type(), LogProp::hasher(),
                 LogProp::key_equal(), LogProp::allocator_type(home)) {
    // We tell that we have a subscription to x the new constructed space
    // The hashmap of the VarIdPos object is shared between all spaces. We must
    // tell here that we want to access the same memory that was allocated in
    // the default constructor for the hashmap even if we change space. The
    // only exception is when we use multithreading; the hash map will be copied
    // and shared amongs the spaces in the new thread
    varpos.update(home,share,b.varpos);

    x.update(home,share,b.x);
    for (auto& elem : logProp)
      elem.second = PropInfo(home, b.logProp[elem.first]);
  }
  bool status(const Space& home) const override {
    auto& h = const_cast<Space&>(home);
    for (Propagators p(h, PropagatorGroup::all); p(); ++p) {
      // Sum of domains of all variable in propagator
      unsigned int domAggr;
      // Same, but for variables that are also in this brancher.
      unsigned int domAggrB;
      p.propagator().slndistsize(&varInBrancher.getObject(), domAggr, domAggrB);
      // If there's still a propagator that has an unassigned variable that is
      // also in the brancher, we tell our brancher has still work to do.
      if (domAggrB > 0)
        return true;
    }
    return false;
  }
  virtual Candidate getChoice(Space& home) = 0;

  bool compute(VarId var_id) const override {
    return true;
  }
  // Method used by all propagators for communicating calculated densities for
  // each of its (variable,value) pair.
  void setMarginalDistribution(PropId prop_id, VarId var_id, Val val,
                                       Dens density) override {
    assert(var_id != 0);
    if (!varpos.isIn(var_id)) return;
    assert(!x[varpos[var_id]].assigned());
    assert(x[varpos[var_id]].in(val));
    assert(logProp.find(prop_id) != logProp.end());

    logProp[prop_id].insert_record(PropInfo::Record{var_id, val, density});
  }
  void setSupportSize(PropId prop_id, SlnCnt count) override {
    logProp[prop_id].setSlnCnt(count);
  }
  const Choice* choice(Space& home) override {
    // Active propagators and the size we need for their log.
    // We considere a propagator as "active" only if
    // - it supports slndist()
    // - it has unassigned variables that are also in the brancher
    struct Psize { size_t domAggr; size_t domAggrB; };
    __gnu_cxx::hash_map<PropId, Psize > activeProps;
    for (Propagators p(home, PropagatorGroup::all); p(); ++p) {
      unsigned int domAggr, domAggrB;
      p.propagator().slndistsize(&varInBrancher.getObject(), domAggr, domAggrB);
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
          prop->reuse_mem(aProp->domAggrB);
      } else {
        // We create a new propagator
        logProp[prop_id] = PropInfo(home, aProp->domAggr, aProp->domAggrB);
      }

      if (!in_log || changed) {
        p.propagator().slndist(home,this);
      }
    }

    // We find the choice.
    Candidate c = getChoice(home);
    assert(!x[c.idx].assigned());
    assert(x[c.idx].in(c.val));

    return new PosValChoice<int>(*this,2,c.idx,c.val);
  }
  const Choice* choice(const Space&, Archive& e) override {
    int pos, val;
    e >> pos >> val;
    return new PosValChoice<int>(*this,2,pos,val);
  }
  ExecStatus commit(Space& home, const Choice& c, unsigned int a) override {
    const auto& pvc = static_cast<const PosValChoice<int>&>(c);
    int pos=pvc.pos().pos, val=pvc.val();
    if (a == 0)
      return me_failed(x[pos].eq(home,val)) ? ES_FAILED : ES_OK;
    else
      return me_failed(x[pos].nq(home,val)) ? ES_FAILED : ES_OK;
  }
  void print(const Space& home, const Choice& c, unsigned int a,
                     std::ostream& o) const override {
    const auto& pvc = static_cast<const PosValChoice<int>&>(c);
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
class MAXSD : public CBSBrancher<View> {
  using CBSBrancher<View>::x;
  using CBSBrancher<View>::varpos;
  using CBSBrancher<View>::for_every_log_entry;
  typedef typename CBSBrancher<View>::Candidate Candidate;
public:
  MAXSD(const Home& home, ViewArray<View>& x0, double recomputation_ratio0)
    : CBSBrancher<View>(home, x0, recomputation_ratio0) {}

  MAXSD(Space& home, bool share, CBSBrancher<View>& b)
    : CBSBrancher<View>(home, share, b) {}

  Candidate getChoice(Space& home) override {
    PropInfo::Record best{0,0,0};
    for_every_log_entry([&](PropId prop_id, SlnCnt slnCnt,
                            VarId var_id, Val val, SlnCnt dens) {
      unsigned int pos = varpos[var_id];
      if (!x[pos].assigned() && x[pos].in(val))
        if (dens > best.dens || (dens == best.dens && var_id < best.var_id))
          best = {var_id, val, dens};
    });
    assert(best.var_id != 0);
    return {varpos[best.var_id],best.val};
  }
  Brancher* copy(Space& home, bool share) override {
    return new (home) MAXSD(home,share,*this);
  }
  static void post(Home home, ViewArray<View>& x,
                   double recomputation_ratio=1) {
    (void) new (home) MAXSD(home,x,recomputation_ratio);
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
      assert(false);
  }
}

void cbsbranch(const Home& home, const IntVarArgs& x, CBSBranchingHeuristic s) {
  _cbsbranch<Int::IntView,IntVarArgs>(home,x,s);
}

void cbsbranch(const Home& home, const BoolVarArgs& x, CBSBranchingHeuristic s) {
  _cbsbranch<Int::BoolView,BoolVarArgs>(home,x,s);
}

#endif //__CBS_HPP__