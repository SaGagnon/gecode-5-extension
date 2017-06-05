#ifndef __CBS_HPP__
#define __CBS_HPP__

#include <gecode/int.hh>
#include <ext/hash_map>
#include <ext/hash_set>
#include <cstring>
#include <tuple>

using namespace Gecode;

struct Record { unsigned int var_id; int val; double density; };
typedef __gnu_cxx::hash_map< unsigned int, std::pair<size_t,Record*>,
  __gnu_cxx::hash<unsigned int>, __gnu_cxx::equal_to<unsigned int>,
  Gecode::space_allocator<std::pair<size_t,Record*>> > LogDensity;

typedef __gnu_cxx::hash_map< unsigned int, std::pair<size_t, double>,
  __gnu_cxx::hash<unsigned int>, __gnu_cxx::equal_to<unsigned int>,
  Gecode::space_allocator<std::pair<size_t, double> > > LogProp;


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
  unsigned int& operator[](unsigned int i) {
    return static_cast<VarIdToPosO*>(object())->hash_map[i];
  }
  unsigned int operator[](unsigned int i) const {
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

    minVal = INT_MAX;
    int maxVal = INT_MIN;
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

unsigned int varvalpos(const VADesc& xD, unsigned int var_id, int val) {
  return xD.positions[var_id] * xD.width + val - xD.minVal;
}

unsigned int varpos(const VADesc& xD, unsigned int var_id) {
  return xD.positions[var_id];
}

class AbstractFeature {
public:
  virtual AbstractFeature* copy(Space& home, bool share) = 0;
  virtual double get(const VADesc& xD, unsigned int var_id, int val) = 0;
  virtual void aggregate(const VADesc& xD, unsigned int prop_id, double slnCnt,
                         const Record& r) = 0;
  virtual void clear(void) = 0;
};

template<class T>
class Feature : public AbstractFeature {
protected:
  SharedArray<T> arr;
public:
  Feature() {}
  Feature(int size) {
    arr.init(size);
    for (int i=0; i<size; i++)
      arr[i] = 0;
  }
  Feature(Space& home, bool share, Feature& f) {
    arr.update(home,share,f.arr);
  }
  void clear(void) {
    for (int i=0; i<arr.size(); i++)
      arr[i] = 0;
  }
};

template<class T, class Concrete>
class FeatureVarVal : public Feature<T> {
protected:
  using Feature<T>::arr;
public:
  FeatureVarVal() : Feature<T>() {}
  FeatureVarVal(VADesc& xD)
    : Feature<T>(xD.size * xD.width) {}
  FeatureVarVal(Space& home, bool share, FeatureVarVal& f)
    : Feature<T>(home,share,f) {}
  virtual AbstractFeature* copy(Space& home, bool share) {
    Concrete *ret = home.alloc<Concrete>(1);
    *ret = Concrete(home,share,*this);
    return ret;
  }
  virtual double get(const VADesc& xD, unsigned int var_id, int val) {
    return arr[varvalpos(xD,var_id,val)];
  }
};

// This
#define CLASS_FEATURE(Name,Type) \
  class Name : public FeatureVarVal<Type,Name> { \
    using FeatureVarVal<Type,Name>::arr; \
    using FeatureVarVal<Type,Name>::FeatureVarVal; \
  public:


CLASS_FEATURE(MaxDensity, double)
  virtual void aggregate(const VADesc& xD, unsigned int prop_id, double slnCnt,
                         const Record& r) {
    double *dens = &arr[varvalpos(xD, r.var_id, r.val)];
    if (*dens < r.density)
      *dens = r.density;
  }
};

CLASS_FEATURE(DensitySum, double)
  virtual void aggregate(const VADesc& xD, unsigned int prop_id, double slnCnt,
                         const Record& r) {
    arr[varvalpos(xD, r.var_id, r.val)] += r.density;
  }
};

CLASS_FEATURE(VarPropCount, double)
  virtual void aggregate(const VADesc& xD, unsigned int prop_id, double slnCnt,
                         const Record& r) {
    arr[varvalpos(xD, r.var_id, r.val)] += 1;
  }
};

CLASS_FEATURE(SlnCntProdDens, double)
  virtual void aggregate(const VADesc& xD, unsigned int prop_id, double slnCnt,
                         const Record& r) {
    arr[varvalpos(xD, r.var_id, r.val)] += slnCnt * r.density;
  }
};

CLASS_FEATURE(SlnCntSum, double)
  virtual void aggregate(const VADesc& xD, unsigned int prop_id, double slnCnt,
                         const Record& r) {
    arr[varvalpos(xD, r.var_id, r.val)] += slnCnt;
  }
};

/**
 * \brief Base class for collecting densities from propagators
 *
 * Before beginning to explain what this class does, we must talk about the
 * cbs() method in Gecode::Propagator and the CBS class from which we inherit.
 *
 * Gecode::Propagator is the base class for every constraint in Gecode.
 * The algorithm to compute solution densities for a given (variable,value) pair
 * depends on the constraints in which the variable is involved (its
 * propagators). For this reason, there's a virtual method in Gecode::Propagator
 * that any given concrete propagator can overload to tell how it compute
 * densities for its variables. Here is the signature of the method:
 *
 *   Gecode::Propagator::cbs(Space& home, CBS* densities) const;
 *
 * As the time of writting this comment, the distinct propagator (see
 * gecode/int/distinct.hh) and regular propagator (see
 * gecode/int/extensional.hh) specialise the cbs method to specify its algorithm
 * to compute solutions densities.
 *
 * Of interest to us is the second argument of this method: CBS* densities.
 * The class CBS is an interface (virtual pure class) that contains only one
 * method:
 *
 * virtual void Gecode::CBS::set(unsigned int var_id, int val, double density) = 0;
 *
 * When a propagator overloads its cbs() method, all it knows is he has access
 * to a object CBS with a set() method that enables him to communicate the
 * results of its computation (i.e. densities for each of its
 * (variable,value) pair).
 *
 * When braching, CBSBrancher must call the cbs() method on each active
 * propagator that overloads the method and pass an object that inherits the
 * class CBS to store densities and compute a choice for branching.
 *
 * This class is a specialisation of CBS that we use for caching densities
 * between computation in a Log and reuse them if possible.
 */
template<class View>
class BranchingHeuristic : public CBS {
public:
  // A choice for branching
  struct Candidate {
    int idx; // Index of the variable in x
    int val; // Value in the domain of the variable
  };
private:
  // Propagator that is currently using the set() method.
  int current_prop;
protected:
  AbstractFeature **features;
  unsigned int n_features;

  // Array of variables we are using for branching
  ViewArray<View> x;
  VADesc xD;
  // The log in which the set method is currently inserting
  LogDensity *logDensity;
  LogProp *logProp;
public:
  /**
   * Constructor for posting.
   *
   * This constructor will be called only once at the beginning of the problem
   * when we post the brancher.
   *
   * @param home Space in which we construct this object
   * @param x0 Variables concerning the branching heuristic
   */
  BranchingHeuristic(Space& home, const ViewArray<View>& x0,
                     unsigned int n_features0)
    : x(x0), xD(x0), n_features(n_features0) {
    features = home.alloc<AbstractFeature*>(n_features);
  }
  /**
   * Constructor for cloning spaces
   *
   * @param home New constructed space
   * @param share Can we share internal data with the new clone
   * @param bh Object that is being cloned from the parent space
   */
  BranchingHeuristic(Space& home, bool share, BranchingHeuristic& bh)
    : xD(home,share,bh.xD), n_features(bh.n_features) {
    features = home.alloc<AbstractFeature*>(n_features);
    for (int i=0; i<n_features; i++)
      features[i] = bh.features[i]->copy(home,share);
    x.update(home,share,bh.x);
  }
  // Method for specifying the logDensity the set() method uses
  void set_log_density(LogDensity *_logDensity) {
    assert(_logDensity != NULL);
    logDensity = _logDensity;
  }
  void set_log_sln_cnt(LogProp *_logSlnCnt) {
    assert(_logSlnCnt != NULL);
    logProp = _logSlnCnt;
  }
  // Method for specifying the propagator that is currently using the set
  // method of this class
  void set_current_prop(unsigned int prop_id) {
    current_prop = prop_id;
  }
  // Method used by all propagators for communicating calculated densities for
  // each of its (variable,value) pair.
  virtual void set(unsigned int var_id, int val, double density) {
    assert(current_prop != -1);
    assert(density>0 && density<1);
    Record r; r.var_id=var_id; r.val=val; r.density=density;
    // Number of records for the current propagator
    size_t *nb_record = &(*logDensity)[current_prop].first;
    (*logDensity)[current_prop].second[(*nb_record)++] = r;
  }
  virtual void setSlnCnt(double slnCnt) {
    assert(current_prop != -1);
    (*logProp)[current_prop].second = slnCnt;
  }

  virtual Candidate getChoice(Space& home) {
    // Some memory area for computing features may be shared between all
    // the spaces. We ask the clear thosees areas before computing.
    for (int f=0; f<n_features; f++) {
      features[f]->clear();
    }

    // isVarSub is true if it is inside at least one propagator that supports
    // cbs. If we have no record in the log for this variable, it will be false
    // and we don't want to consider it in our choice recommendation
    __gnu_cxx::hash_set<unsigned int> isVarSub;

    // For every propagator in the log
    for (LogDensity::iterator it=logDensity->begin(); it!=logDensity->end(); ++it) {
      unsigned int prop_id = it->first;
      size_t nb_records = it->second.first;
      // Aggregation for every record concerning the propagator (i.e. every
      // (variable,value) pair and their corresponding density)
      for (unsigned int i=0; i<nb_records; ++i) {
        Record r = it->second.second[i];
        for (int f=0; f<n_features; f++) {
          features[f]->aggregate(xD,prop_id,(*logProp)[prop_id].second,r);
        }
      }
    }

    Candidate c;
    c.idx = -1;
    double best_score = INT_MIN;

    for (int i=0; i<x.size(); i++) {

      bool instrumented = false;
      for (SubscribedPropagators sp(x[i]); sp(); ++sp) {
        if (sp.propagator().cbs(home,NULL)) {
          instrumented = true;
          break;
        }
      }
      if (!instrumented) continue;

      if (x[i].assigned()) continue;
      for (Int::ViewValues<View> val(x[i]); val(); ++val) {
        double score = getScore(xD,x[i].id(),val.val());
        if (score > best_score) {
          c.idx = i;
          c.val = val.val();
          best_score = score;
        }
      }
    }
    return c;
  }

protected:
  virtual double getScore(const VADesc& xD, unsigned int var_id, int val) = 0;
};

template<typename F>
void assign(int i, Space& home, AbstractFeature** af, VADesc& xD) {
  af[i] = home.alloc<F>(1);
  *static_cast<F*>(af[i]) = F(xD);
}

#define USING_BRANCHING_HEUR \
  using BranchingHeuristic<View>::BranchingHeuristic; \
  using BranchingHeuristic<View>::n_features; \
  using BranchingHeuristic<View>::features; \
  using BranchingHeuristic<View>::x; \
  using BranchingHeuristic<View>::xD;



template<class View>
class maxSD : public BranchingHeuristic<View> {
  USING_BRANCHING_HEUR
public:
  maxSD(Space& home, const ViewArray<View>& x)
    : BranchingHeuristic<View>(home,x,1) {
    assign<MaxDensity>(0,home,features,xD);
  }
  virtual double getScore(const VADesc& xD, unsigned int var_id, int val)  {
    return features[0]->get(xD,var_id,val);
  }
};

template<class View>
class aAvgSD : public BranchingHeuristic<View> {
  USING_BRANCHING_HEUR
public:
  aAvgSD(Space& home, const ViewArray<View>& x)
    : BranchingHeuristic<View>(home,x,2) {
    assign<DensitySum>(0,home,features,xD);
    assign<VarPropCount>(1,home,features,xD);
  }
  virtual double getScore(const VADesc& xD, unsigned int var_id, int val)  {
    return features[0]->get(xD,var_id,val) / features[1]->get(xD,var_id,val);
  }
};


template<class View>
class maxRelSD : public BranchingHeuristic<View> {
  USING_BRANCHING_HEUR
public:
  maxRelSD(Space& home, const ViewArray<View>& x)
    : BranchingHeuristic<View>(home,x,1) {
    assign<MaxDensity>(0,home,features,xD);
  }
  virtual double getScore(const VADesc& xD, unsigned int var_id, int val)  {
    return features[0]->get(xD,var_id,val) - 1.0/x[varpos(xD,var_id)].size();
  }
};


template<class View>
class maxRelRatio : public BranchingHeuristic<View> {
  USING_BRANCHING_HEUR
public:
  maxRelRatio(Space& home, const ViewArray<View>& x)
    : BranchingHeuristic<View>(home,x,1) {
    assign<MaxDensity>(0,home,features,xD);
  }
  virtual double getScore(const VADesc& xD, unsigned int var_id, int val)  {
    return features[0]->get(xD,var_id,val) * x[varpos(xD,var_id)].size();
  }
};


template<class View>
class wcSCAvg : public BranchingHeuristic<View> {
  USING_BRANCHING_HEUR
public:
  wcSCAvg(Space& home, const ViewArray<View>& x)
    : BranchingHeuristic<View>(home,x,2) {
    assign<SlnCntProdDens>(0,home,features,xD);
    assign<SlnCntSum>(1,home,features,xD);
  }
  virtual double getScore(const VADesc& xD, unsigned int var_id, int val)  {
    return features[0]->get(xD,var_id,val) / features[1]->get(xD,var_id,val);
  }
};



template<class View, template<class> class BranchingHeur>
class CBSBrancher : public Brancher {
  typedef typename BranchingHeuristic<View>::Candidate Candidate;
protected:
  ViewArray<View> x;
  BranchingHeur<View> heur;
  LogDensity logDensity;
  LogProp logProp;
public:
  CBSBrancher(Home home, ViewArray<View>& x0)
    : Brancher(home), x(x0), heur(home,x0),
      logDensity(LogDensity::size_type(), LogDensity::hasher(),
                 LogDensity::key_equal(), LogDensity::allocator_type(home)),
      logProp(LogProp::size_type(), LogProp::hasher(),
                LogProp::key_equal(), LogProp::allocator_type(home))
  {
    // Because we must call the destructor of aAvgSD
    home.notice(*this,AP_DISPOSE);
  }
  static void post(Home home, ViewArray<View>& x) {
    (void) new (home) CBSBrancher(home,x);
  }
  virtual size_t dispose(Space& home) {
    home.ignore(*this, AP_DISPOSE);
    // ~aAvgSD() calls ~SharedHashMap() which calls ~SharedHashMapObject() to
    // deallocate the hash map when the refcount of SharedHashMapObject is 0
    heur.~BranchingHeur<View>();
    (void) Brancher::dispose(home);
    return sizeof(*this);
  }
  CBSBrancher(Space& home, bool share, CBSBrancher& b)
    : Brancher(home,share,b), heur(home,share,b.heur),
      logDensity(b.logDensity.begin(), b.logDensity.end(),
                 LogDensity::size_type(), LogDensity::hasher(),
                 LogDensity::key_equal(), LogDensity::allocator_type(home)),
      logProp(b.logProp.begin(), b.logProp.end(),
                 LogProp::size_type(), LogProp::hasher(),
                 LogProp::key_equal(), LogProp::allocator_type(home)) {
    x.update(home,share,b.x);

    for (LogDensity::iterator it=logDensity.begin(); it!=logDensity.end(); ++it) {
      unsigned int prop_id = it->first;
      size_t count = b.logDensity[prop_id].first;
      it->second.second = home.alloc<Record>(count);
      memcpy(it->second.second, b.logDensity[prop_id].second, count*sizeof(Record));
    }
  }
  virtual Brancher* copy(Space& home, bool share) {
//    CBSBrancher *ret = home.alloc<CBSBrancher>(1);
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
  virtual const Choice* choice(Space& home) {
    // Active propagators and the size we need for their log.
    // We considere a propagator as "active" only if he supports cbs().
    // TODO: Commentaire
    __gnu_cxx::hash_map<unsigned int, int> activeProps;
    for (Propagators p(home, PropagatorGroup::all); p(); ++p) {
      int log_size = p.propagator().cbs(home,NULL);
      // If the propagator supports cbs and has records
      if (log_size != 0) {
        assert(log_size > 1);
        activeProps[p.propagator().id()] = log_size;
      }
    }

    // We delete log elements corresponding to non active propagators
    {
      std::vector<unsigned int> propsToDelete;

      // Propagators to delete (not active and in log)
      for (LogDensity::iterator it = logDensity.begin(); it != logDensity.end(); ++it)
        if (activeProps.find(it->first) == activeProps.end())
          propsToDelete.push_back(it->first);

      // We delete them from the log
      for (std::vector<unsigned int>::iterator it = propsToDelete.begin();
           it != propsToDelete.end(); ++it) {
        unsigned int prop_id = *it;
        // TODO: Trouver une manière de free comme il le faut (avec vrai grandeur, pas nb element en ce moment)...
//        home.free(log[prop_id].second, log[prop_id].first);
        logDensity.erase(prop_id);
        logProp.erase(prop_id);
      }
    }

    // We specify the log that will be modified when the propagators use the
    // CBS::set()
    heur.set_log_density(&logDensity);
    heur.set_log_sln_cnt(&logProp);
    for (Propagators p(home, PropagatorGroup::all); p(); ++p) {
      if (!p.propagator().cbs(home,NULL)) continue;
      unsigned int prop_id = p.propagator().id();
      // Propagator already in the log?
      bool in_log = logDensity.find(prop_id) != logDensity.end();

      if (in_log) {
        bool changed = logProp[prop_id].first != activeProps[prop_id];
        if (changed) {
          // We discard the previous entries by setting the count to 0 (we
          // thus reuse previous allocated memory. The number of records can't
          // grow).
          logDensity[prop_id].first = 0;
          // TODO: Comment
          logProp[prop_id].first = (size_t)activeProps[prop_id];
        } else {
          // We continue, meaning we will use the previous entries
          continue;
        }
      } else {
        // Otherwise, we need to allocate space for the records of the
        // new propagator
        logDensity[prop_id] =
          std::make_pair(0, home.alloc<Record>(activeProps[prop_id]));
        logProp[prop_id] = std::make_pair(activeProps[prop_id], 0);
      }
      heur.set_current_prop(prop_id);
      p.propagator().cbs(home,&heur);
    }

    // We find the choice.
    Candidate c = heur.getChoice(home);
    static int count=0;
    assert(!x[c.idx].assigned());
    assert(x[c.idx].in(c.val));

    for (Propagators p(home, PropagatorGroup::all); p(); ++p) {
      if (!p.propagator().cbs(home, NULL)) continue;
      unsigned int prop_id = p.propagator().id();
      assert(logProp[prop_id].first == logDensity[prop_id].first);
    }

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
};

enum CBSBranchingHeuristic {
  MAX_SD,
  MAX_REL_SD,
  MAX_REL_RATIO,
  A_AVG_SD,
  W_SC_AVG
};


template<class View, class T>
void _cbsbranch(Home home, const T& x, CBSBranchingHeuristic s) {
  if (home.failed()) return;
  ViewArray<View> y(home,x);

  switch(s) {
    case MAX_SD:
      CBSBrancher<View,maxSD>::post(home,y);
      break;
    case MAX_REL_SD:
      CBSBrancher<View,maxRelSD>::post(home,y);
      break;
    case MAX_REL_RATIO:
      CBSBrancher<View,maxRelRatio>::post(home,y);
      break;
    case A_AVG_SD:
      CBSBrancher<View,aAvgSD>::post(home,y);
      break;
    case W_SC_AVG:
      CBSBrancher<View,wcSCAvg>::post(home,y);
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
