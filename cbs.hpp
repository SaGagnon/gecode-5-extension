#ifndef __CBS_HPP__
#define __CBS_HPP__

#include <iostream>
#include <ext/hash_map>

#include <gecode/int.hh>
#include <map>

#ifdef SQL
#include <sql-interface.hh>
#endif

#define MD

const int COUSSIN_MD = 0;


// TODO: Ne pas déclarer using namespace dans un .hpp...
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
  size_t domsum{0};

  unsigned int mindom{0};
public:
  PropInfo() = default;
  PropInfo(Space& home, size_t domsum0, size_t domsum_b, unsigned int md)
    : domsum(domsum0), mindom(md) {
    records.x = home.alloc<Record>(domsum_b);
    records.size = domsum_b;
  }
  PropInfo(Space& home, const PropInfo& o)
    : domsum(o.domsum), records(o.records), slnCount(o.slnCount),
      mindom(o.mindom) {
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
  size_t getDomSum() const { return domsum; }
  size_t getDomSumB() const { return records.size; }
  unsigned int getPosRec() const { return records.pos; }

  SlnCnt getSlnCnt() const { return slnCount; }
  void setSlnCnt(SlnCnt cnt) {
    assert(cnt > 1);
    slnCount = cnt;
  }
  void insert_record(const Record&& r) {
    assert(records.x != nullptr);
    assert(r.dens>=0 && r.dens<=1.001);
    records.x[records.pos] = r;
    records.pos++;
    assert(records.size >= records.pos);
  }
  void reuse_mem(size_t domsum0, size_t domsum_b) {
    // We discard the previous entries by setting the count to 0 (we
    // thus reuse previous allocated memory. The number of records can't
    // grow).
    records.pos = 0;
    records.size = domsum_b;
    domsum = domsum0;
  }
  void setMinDom(unsigned int md) {
    mindom = md;
  }
  unsigned int getMinDom() {
    return mindom;
  }
};

using LogProp = __gnu_cxx::hash_map<PropId, PropInfo, __gnu_cxx::hash<PropId>,
  __gnu_cxx::equal_to<PropId>, Gecode::space_allocator<PropInfo>>;

/**
 * \brief Maps variables ids to indexes
 */
template<class Key, class Val>
class SharedHashMap : public SharedHandle {
  using HashMap =  __gnu_cxx::hash_map<Key,Val>;
protected:
  class SharedHashMapO : public SharedHandle::Object {
  public:
    HashMap hash_map;
  public:
    SharedHashMapO() = default;
    SharedHashMapO(const SharedHashMapO& o)
      : hash_map(o.hash_map) {}
    Object* copy() const override {
      return new SharedHashMapO(*this);
    }
    ~SharedHashMapO() override = default;
  };
public:
  SharedHashMap() = default;
  void init() {
    assert(object() == nullptr);
    object(new SharedHashMapO());
  }
  bool isIn(unsigned int i) const {
    auto *hm = &static_cast<SharedHashMapO*>(object())->hash_map;
    return hm->find(i) != hm->end();
  }
//  Val& operator[](unsigned int i) {
//    return static_cast<SharedHashMapO*>(object())->hash_map.at(i);
//  }
  Val operator[](unsigned int i) const {
//    assert(isIn(i));
    return static_cast<SharedHashMapO*>(object())->hash_map.at(i);
  }
  void insert(unsigned int k, unsigned int v) {
    static_cast<SharedHashMapO*>(object())->hash_map.insert(std::make_pair(k,v));
  }
};

class VarInBrancher : public SharedHandle {
protected:
  class VarInBrancherO : public SharedHandle::Object {
  public:
    class BoolArray {
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
      bool inbrancher(VarId var_id) const {
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
  bool inbrancher(VarId var_id) const {
    return getObject().inbrancher(var_id);
  }
};

template<class View>
class CBSBrancher : public Brancher {
public:
  // A choice for branching
  struct Candidate {
    unsigned int idx; // Index of the variable in x
    Val val; // Value in the domain of the variable
  };
protected:
  mutable bool naif_branching;
  const double recomputation_ratio;
  ViewArray<View> x;
  SharedHashMap<VarId, unsigned int> varpos;
  VarInBrancher varInBrancher;
  LogProp logProp;
  unsigned int global_mindom;
public:
  CBSBrancher(Home home, ViewArray<View>& x0, double recomputation_ratio0=1)
    : Brancher(home), naif_branching(false),
      recomputation_ratio(recomputation_ratio0),
      x(x0), varInBrancher(x0),
      logProp(LogProp::size_type(), LogProp::hasher(),
              LogProp::key_equal(), LogProp::allocator_type(home)),
      global_mindom(0) {
    assert(recomputation_ratio0 >= 0 && recomputation_ratio0 <= 1);
    // Because we must call the destructor of aAvgSD
    home.notice(*this,AP_DISPOSE);
    // The VarIdToPos object is first implicitly constructed with the default
    // constructor and its shared hashmap is not yet allocated. init() thus
    // allocate memory for the hashmap
    varpos.init();

    // We assign an index for each variable id
    for (unsigned int i=0; i<x.size(); i++)
      varpos.insert(x[i].id(), i);
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
    : Brancher(home,share,b), naif_branching(b.naif_branching),
      recomputation_ratio(b.recomputation_ratio),
      varInBrancher(b.varInBrancher),
      logProp(b.logProp.begin(), b.logProp.end(),
                 LogProp::size_type(), LogProp::hasher(),
                 LogProp::key_equal(), LogProp::allocator_type(home)),
      global_mindom(b.global_mindom) {
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
      unsigned int domsum;
      // Same, but for variables that are also in this brancher.
      unsigned int domsum_b;
      p.propagator().domainsizesum(
        // (STD::BIND) Passing VarInBrancher as a standart function pointer: we
        // need to bind it with an instance of the object (the "this" argument)
        std::bind(&VarInBrancher::inbrancher, &varInBrancher, std::placeholders::_1),
        domsum, domsum_b);
      // If there's still a propagator that has an unassigned variable that is
      // also in the brancher, we tell our brancher has still work to do.
      if (domsum_b > 0)
        return true;
    }

    naif_branching = true;
    for (int i =0; i<x.size(); i++)
      if (!x[i].assigned())
        return true;

    return false;
  }
  virtual Candidate getChoice(Space& home) = 0;

//  bool compute(VarId var_id) const override {
//    return true;
//  }
//  Type type() const override {
//    return SolnDistrib::ALL;
//  }
  // Method used by all propagators for communicating calculated densities for
  // each of its (variable,value) pair.
  void marginaldistrib(PropId prop_id, VarId var_id, Val val, Dens density) {
    assert(var_id != 0);
    if (!varpos.isIn(var_id)) return;
    assert(!x[varpos[var_id]].assigned());
    assert(x[varpos[var_id]].in(val));
    assert(logProp.find(prop_id) != logProp.end());

    logProp[prop_id].insert_record(PropInfo::Record{var_id, val, density});
  }
//  void supportsize(PropId prop_id, SlnCnt count) override {
//    logProp[prop_id].setSlnCnt(count);
//  }
  const Choice* choice(Space& home) override {

    if (!naif_branching) {
      // Active propagators and the size we need for their log.
      // We considere a propagator as "active" only if
      // - it supports slndist()
      // - it has unassigned variables that are also in the brancher
      struct Psize {
        size_t domsum;
        size_t domsum_b;
        unsigned int mindom;
      };

      global_mindom = INT_MAX;

      __gnu_cxx::hash_map<PropId, Psize> activeProps;
      for (Propagators p(home, PropagatorGroup::all); p(); ++p) {
        unsigned int domsum, domsum_b, mindom;
        p.propagator().domainsizesum(
          // See "STD::BIND" comment
          std::bind(&VarInBrancher::inbrancher, &varInBrancher, std::placeholders::_1),
          domsum, domsum_b);
        if (domsum_b != 0) {
          p.propagator().mindom(
            std::bind(&VarInBrancher::inbrancher, &varInBrancher, std::placeholders::_1),
            mindom
          );
          if (mindom < global_mindom)
            global_mindom = mindom;

          activeProps[p.propagator().id()] = {domsum, domsum_b, mindom};
        }
      }

      // We delete log elements corresponding to non active propagators
      {
        std::vector<PropId> propsToDelete;
        for (const auto& kv : logProp)
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
          changed = prop->getDomSum() * recomputation_ratio > aProp->domsum;
          #ifdef MD
//          if (changed && aProp->mindom > global_mindom + COUSSIN_MD) {
//            static int cnt = 0;
//            std::cout << aProp->mindom << " " << global_mindom << " " << ++cnt << std::endl;
//          }
          changed = changed && aProp->mindom <= global_mindom + COUSSIN_MD;
          #endif
          if (changed) {
            prop->reuse_mem(aProp->domsum, aProp->domsum_b);
            prop->setMinDom(aProp->mindom);
          }
        } else {
          // We create a new propagator
          logProp[prop_id] = PropInfo(home, aProp->domsum, aProp->domsum_b, aProp->mindom);
        }



        using namespace std::placeholders;

        if (!in_log || changed) {
          p.propagator().solndistrib(
            home,
            // See "STD::BIND" comment
            std::bind(&CBSBrancher::marginaldistrib, this, _1, _2, _3, _4)
//          ,  Propagator::SolnDistribCalc::MAX_PER_PROP
          );

//          static int cnt = 0;
//          std::cout << "calc: " << ++cnt << std::endl;
        }

      }

      // We find the choice.
//      static int cnt = 0;
//      std::cout << "getChoice " << ++cnt << std::endl;
      Candidate c = getChoice(home);
    assert(!x[c.idx].assigned());
    assert(x[c.idx].in(c.val));
      return new PosValChoice<int>(*this, 2, c.idx, c.val);
    } else {
      Candidate c;
      int mindom = INT_MAX;
      for (int i=0; i<x.size(); i++) {
        if (x[i].assigned()) continue;
        if (x[i].size() < mindom) {
          mindom = x[i].size();
          int minval = INT_MAX;
          for (Int::ViewValues<View> val(x[i]); val(); ++val) {
            if (val.val() < minval) {
              c = {i, val.val()};
              minval = val.val();
            }
          }
        }
      }
      return new PosValChoice<int>(*this, 2, c.idx, c.val);
    }
  }
  const Choice* choice(const Space&, Archive& e) override {
    int pos, val;
    e >> pos >> val;
    return new PosValChoice<int>(*this,2,pos,val);
  }
  ExecStatus commit(Space& home, const Choice& c, unsigned int a) override {
    const auto& pvc = static_cast<const PosValChoice<int>&>(c);
    int pos=pvc.pos().pos, val=pvc.val();
    if (a == 0) {
//      std::cout << "g " << pos << " " << val << " ";
//      print_stats();
      return me_failed(x[pos].eq(home, val)) ? ES_FAILED : ES_OK;
    }
    else {
//      std::cout << "b " << pos << " " << val << " ";
//      print_stats();
      return me_failed(x[pos].nq(home, val)) ? ES_FAILED : ES_OK;
    }
  }
  void print_stats() {
    if (naif_branching) return;

    std::map<VarId, std::map<Val, std::vector<Dens>>> densities;
    for_every_log_entry([&](PropId, SlnCnt, VarId varId, Val val, Dens dens) {
      densities[varId][val].push_back(dens);
    });

    double sum = 0;
    double nb = 0;

    for (auto var : densities) {
      for (auto val : var.second) {
        auto n = val.second.size();
        if (n != 1) {

          double tot = 0;
          for (auto dens : val.second)
            tot += dens;

          double mean = tot / n;

          double std = 0;
          for (auto dens : val.second)
            std += std::pow(dens - mean, 2);

          std /= n - 1;
          std = std::sqrt(std);

          sum += std;
        }
        nb++;
      }
    }

    std::cout << sum / nb << " " << nb << std::endl;
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
  void for_every_log_entry(const std::function<void(PropId,SlnCnt,VarId,Val,Dens)>& f) const {
    for (const auto& elem : logProp) {
      auto prop_id = elem.first;
      auto prop = &elem.second;
      for (int i=0; i<prop->getPosRec(); i++) {
        auto pos = varpos[(*prop)[i].var_id];
        // With recomputation_ratio, it is possible some records are no more good.
        if (!x[pos].assigned() && x[pos].in((*prop)[i].val))
          f(prop_id, prop->getSlnCnt(), (*prop)[i].var_id, (*prop)[i].val, (*prop)[i].dens);
      }
    }
  }
};


template<class View>
void print_densities_first_node(CBSBrancher<View>* brancher) {
  static bool _print = true;
  if (_print) {
//    brancher->for_every_log_entry([](PropId prop_id, SlnCnt slnCnt,
//                                     VarId var_id, Val val, Dens dens) {
//      std::cout << prop_id << ',' << var_id << ',' << val << ',' << dens
//                << std::endl;
//    });
//    brancher->print_stats();
    _print = false;
  }
}


template<class View>
class maxSD : public CBSBrancher<View> {
  using CBSBrancher<View>::x;
  using CBSBrancher<View>::varpos;
  using CBSBrancher<View>::for_every_log_entry;
  typedef typename CBSBrancher<View>::Candidate Candidate;
public:
  maxSD(const Home& home, ViewArray<View>& x0, double recomputation_ratio0)
    : CBSBrancher<View>(home, x0, recomputation_ratio0) {}

  maxSD(Space& home, bool share, CBSBrancher<View>& b)
    : CBSBrancher<View>(home, share, b) {}

  Candidate getChoice(Space& home) override {
    print_densities_first_node(this);

    PropInfo::Record best{0,0,0};
    for_every_log_entry([&](PropId prop_id, SlnCnt slnCnt,
                            VarId var_id, Val val, Dens dens) {

      unsigned int pos = varpos[var_id];
//      double diff = dens - best.dens;
//      double precision = 0.000001;
//      if (diff > precision || (std::abs(diff) < precision && var_id < best.var_id))
      if (dens > best.dens)
        best = {var_id, val, dens};
    });
    assert(best.var_id != 0);
    return {varpos[best.var_id],best.val};
  }
  Brancher* copy(Space& home, bool share) override {
    return new (home) maxSD(home,share,*this);
  }
  static void post(Home home, ViewArray<View>& x,
                   double recomputation_ratio=1) {
    (void) new (home) maxSD(home,x,recomputation_ratio);
  }
};

template<class View>
class maxSDMD : public CBSBrancher<View> {
  using CBSBrancher<View>::x;
  using CBSBrancher<View>::varpos;
  using CBSBrancher<View>::for_every_log_entry;
  using CBSBrancher<View>::logProp;
  using CBSBrancher<View>::global_mindom;
  typedef typename CBSBrancher<View>::Candidate Candidate;
public:
  maxSDMD(const Home& home, ViewArray<View>& x0, double recomputation_ratio0)
    : CBSBrancher<View>(home, x0, recomputation_ratio0) {}

  maxSDMD(Space& home, bool share, CBSBrancher<View>& b)
    : CBSBrancher<View>(home, share, b) {}

  Candidate getChoice(Space& home) override {
    print_densities_first_node(this);

    PropInfo::Record best{0,0,0};
    for_every_log_entry([&](PropId prop_id, SlnCnt slnCnt,
                            VarId var_id, Val val, Dens dens) {
      if (logProp[prop_id].getMinDom() <= global_mindom + COUSSIN_MD) {
        unsigned int pos = varpos[var_id];
        double diff = dens - best.dens;
        double precision = 0.000001;
        if (diff > precision ||
            (std::abs(diff) < precision && var_id < best.var_id))
          best = {var_id, val, dens};
      }
    });
    assert(best.var_id != 0);
    return {varpos[best.var_id],best.val};
  }
  Brancher* copy(Space& home, bool share) override {
    return new (home) maxSDMD(home,share,*this);
  }
  static void post(Home home, ViewArray<View>& x,
                   double recomputation_ratio=1) {
    (void) new (home) maxSDMD(home,x,recomputation_ratio);
  }
};

template<class View>
class avgSD : public CBSBrancher<View> {
  using CBSBrancher<View>::x;
  using CBSBrancher<View>::varpos;
  using CBSBrancher<View>::for_every_log_entry;
  typedef typename CBSBrancher<View>::Candidate Candidate;
public:
  avgSD(const Home& home, ViewArray<View>& x0, double recomputation_ratio0)
    : CBSBrancher<View>(home, x0, recomputation_ratio0) {}

  avgSD(Space& home, bool share, CBSBrancher<View>& b)
    : CBSBrancher<View>(home, share, b) {}

  Candidate getChoice(Space& home) override {
    print_densities_first_node(this);

    // Aggregation of densities
    using TotDens = Dens;
    using NbProp = unsigned int;
    std::map<VarId, std::map<Val, std::pair<TotDens, NbProp>>> aggr;

    for_every_log_entry([&](PropId prop_id, SlnCnt slnCnt,
                            VarId var_id, Val val, Dens dens) {
      if (aggr[var_id].find(val) == aggr[var_id].end())
        aggr[var_id][val] = std::make_pair(0,0);
      aggr[var_id][val].first += dens;
      aggr[var_id][val].second += 1;
    });

    PropInfo::Record best{0,0,0};
    for (auto varid : aggr) {
      for (auto val : varid.second) {
        double moy = val.second.first / val.second.second;
//        double moy = val.second.first / 2;
        if (moy > best.dens) {
          best = {varid.first, val.first, moy};
        }
      }
    }

    assert(best.var_id != 0);
    return {varpos[best.var_id], best.val};
  }
  Brancher* copy(Space& home, bool share) override {
    return new (home) avgSD(home,share,*this);
  }
  static void post(Home home, ViewArray<View>& x,
                   double recomputation_ratio=1) {
    (void) new (home) avgSD(home,x,recomputation_ratio);
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
void _cbsbranch(Home home, const T& x, CBSBranchingHeuristic s,
                double recomputation_ratio = 1) {
  if (home.failed()) return;
  ViewArray<View> y(home,x);
  switch(s) {
    case MAX_SD:
//      #ifdef MD
//      maxSDMD<View>::post(home,y,recomputation_ratio);
//      #else
      maxSD<View>::post(home,y,recomputation_ratio);
//      #endif

      break;
    case A_AVG_SD:
      avgSD<View>::post(home,y,recomputation_ratio);
      break;
    default:
      assert(false);
  }
}

void cbsbranch(const Home& home, const IntVarArgs& x, CBSBranchingHeuristic s,
               double recomputation_ratio = 1) {
  _cbsbranch<Int::IntView,IntVarArgs>(home,x,s,recomputation_ratio);
}

void cbsbranch(const Home& home, const BoolVarArgs& x, CBSBranchingHeuristic s,
               double recomputation_ratio = 1) {
  _cbsbranch<Int::BoolView,BoolVarArgs>(home,x,s,recomputation_ratio);
}

#endif //__CBS_HPP__