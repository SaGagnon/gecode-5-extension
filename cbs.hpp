#ifndef __CBS_HPP__
#define __CBS_HPP__

#include <gecode/int.hh>
#include <ext/hash_map>
#include <ext/hash_set>
#include <cstring>

using namespace Gecode;

class ChangedPropagators : public LocalHandle {
protected:
  class ChangedPropagatorO : public LocalObject {
  public:
    typedef __gnu_cxx::hash_set<unsigned int, __gnu_cxx::hash<unsigned int>,
      __gnu_cxx::equal_to<unsigned int>,
      Gecode::space_allocator<unsigned int> > Set;
    Set changed;
  public:
    ChangedPropagatorO(Space& home)
      : LocalObject(home),
        changed(Set::size_type(), Set::hasher(), Set::key_equal(),
                Set::allocator_type(home)) {}
    ChangedPropagatorO(Space& home, bool share,
                            ChangedPropagatorO& o)
      : LocalObject(home,share,o),
        changed(o.changed.begin(), o.changed.end(),
                Set::size_type(), Set::hasher(), Set::key_equal(),
                Set::allocator_type(home)) {}
    virtual LocalObject* copy(Space& home, bool share) {
      return new (home) ChangedPropagatorO(home,share,*this);
    }
  };
public:
  ChangedPropagators(void)
    : LocalHandle() {}
  ChangedPropagators(Space& home)
    : LocalHandle(new (home) ChangedPropagatorO(home)) {}
  void insert(unsigned int prop_id) {
    ChangedPropagatorO* o =
      static_cast<ChangedPropagatorO*>(object());
    o->changed.insert(prop_id);
  }
  bool contains(unsigned int prop_id) const {
    ChangedPropagatorO* o =
      static_cast<ChangedPropagatorO*>(object());
    return o->changed.find(prop_id) != o->changed.end();
  }
  void clear() {
    ChangedPropagatorO* o =
      static_cast<ChangedPropagatorO*>(object());
    o->changed.clear();
  }
  ChangedPropagatorO::Set* get(void) {
    ChangedPropagatorO* o =
      static_cast<ChangedPropagatorO*>(object());
    return &o->changed;
  }
};

/**
 * \brief Class for tracking changes in variable domains during propagation
 *
 * ViewUpdateLooker is a propagator whose only goal is to track domain changes
 * in its variables via advisors; it does not do any propagation.
 *
 * Each time a variable is modified, ViewUpdateLooker receives a notification
 * via its advise method. By modifying changedProp accordingly, we can transfer
 * this information to the BranchingHeuristic before making a choice.
 *
 * PC_GEN_NONE means that the propagator is not scheduled for propagation when
 * one of its variable is modified.
 */
template<class View>
class ViewUpdateLooker : public NaryPropagator<View,PC_GEN_NONE> {
  using NaryPropagator<View,PC_GEN_NONE>::x;
protected:
  // An advisor only concern a single variable. For this reason, we must have
  // a council to manage every advisor for each of our variables.
  Council<ViewAdvisor<View>> c;
  // Shared object with BranchingHeuristic to track changes in variable domains.
  ChangedPropagators changedProps;

  // Constructor for posting
  ViewUpdateLooker(Home home, ViewArray<View>& x, const ChangedPropagators& cp)
    : NaryPropagator<View,PC_GEN_NONE>(home,x), c(home), changedProps(cp) {
    for (int i=0; i<x.size(); i++)
      (void) new (home) ViewAdvisor<View>(home, *this, c, x[i]);
  }
  // Constructor for cloning
  ViewUpdateLooker(Space& home, bool share, ViewUpdateLooker<View>& p)
    : NaryPropagator<View,PC_GEN_NONE>(home,share,p) {
    c.update(home,share,p.c);
    changedProps.update(home,share,p.changedProps);
  }
public:
  static ExecStatus post(Home home, ViewArray<View>& x,
                         const ChangedPropagators & cp) {
    (void) new (home) ViewUpdateLooker<View>(home,x,cp);
    return ES_OK;
  }
  virtual ExecStatus advise(Space& home, Advisor &_a, const Delta& d) {
    ViewAdvisor<View>& a(static_cast<ViewAdvisor<View>&>(_a));
    View v(a.view());
    for (SubscribedPropagators sp(v); sp(); ++sp)
      changedProps.insert(sp.propagator().id());

    // If the view is assigned, we don't need its advisor anymore.
    if (a.view().assigned())
      a.dispose(home,c);

    // We can delete the propagator if we have no advisor left
    return c.empty() ? home.ES_SUBSUMED(*this) : ES_FIX;
  }
  virtual ExecStatus propagate(Space&, const ModEventDelta&) {
    // The propagator is not subscribed to any of its variable. Thus, this
    // method won't be called.
    GECODE_NEVER;
    return ES_FIX;
  }
  virtual Actor* copy(Space& home, bool share) {
    return new (home) ViewUpdateLooker<View>(home,share,*this);
  }
  virtual size_t dispose(Space& home) {
    for (Advisors<ViewAdvisor<View>> va(c); va(); ++va)
      va.advisor().dispose(home,c);
    c.dispose(home);
    (void) NaryPropagator<View,PC_GEN_NONE>::dispose(home);
    return sizeof(*this);
  }
};

struct Record { unsigned int var_id; int val; double density; };
typedef __gnu_cxx::hash_map< unsigned int, std::pair<size_t,Record*>,
  __gnu_cxx::hash<unsigned int>, __gnu_cxx::equal_to<unsigned int>,
  Gecode::space_allocator<unsigned int> > Log;

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
  VarIdToPosO::HashMap* get(void) const {
    return &static_cast<VarIdToPosO*>(object())->hash_map;
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
  // Array of variables we are using for branching
  ViewArray<View> x;
  // Hash map that assign an index for each variable in x given its id. It is
  // useful if we want to store computations in a continuous memory space
  VarIdToPos positions;
  // The log in which the set method is currently inserting
  Log *log;
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
  BranchingHeuristic(Space& home, const ViewArray<View>& x0)
    : x(x0) {
    // The VarIdToPos object is first implicitly constructed with the default
    // constructor and its shared hashmap is not yet allocated. init() thus
    // allocate memory for the hashmap
    positions.init();
    // We can assign an index for each variable id
    for (unsigned int i=0; i<x.size(); i++)
      (*positions.get())[x[i].id()] = i;
  }
  /**
   * Constructor for cloning spaces
   *
   * @param home New constructed space
   * @param share Can we share internal data with the new clone
   * @param bh Object that is being cloned from the parent space
   */
  BranchingHeuristic(Space& home, bool share, BranchingHeuristic& bh) {
    // We tell that we have a subscription to x the new constructed space
    x.update(home,share,bh.x);
    // The hashmap of the VarIdPos object is shared between all spaces. We must
    // tell here that we want to access the same memory that was allocated in
    // the default constructor for the hashmap even if we change space. The
    // only exception is when we use multithreading; the hash map will be copied
    // and shared amongs the spaces in the new thread
    positions.update(home,share,bh.positions);
  }
  // Method for specifying the log the set() method uses
  void set_log(Log *_log) {
    assert(_log != NULL);
    log = _log;
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
    Record r; r.var_id=var_id; r.val=val; r.density=density;
    // Number of records for the current propagator
    size_t *nb_record = &(*log)[current_prop].first;
    (*log)[current_prop].second[(*nb_record)++] = r;
  }

  virtual Candidate get(Space& home) = 0;
};

template<class View>
class aAvgSD : public BranchingHeuristic<View> {
  typedef typename BranchingHeuristic<View>::Candidate Candidate;
  using BranchingHeuristic<View>::x;
  using BranchingHeuristic<View>::positions;
  using BranchingHeuristic<View>::log;
private:
  int minVal;
  int width;
  SharedArray<double> densities_sum;
  SharedArray<int> count;
public:
  aAvgSD(Space& home, const ViewArray<View>& x)
    : BranchingHeuristic<View>(home,x) {
    minVal = INT_MAX;
    int maxVal = INT_MIN;
    for (unsigned int i=0; i<x.size(); i++) {
      if (x[i].min() < minVal) minVal = x[i].min();
      if (x[i].max() > maxVal) maxVal = x[i].max();
    }
    assert(minVal != INT_MAX && maxVal != INT_MIN);

    width = maxVal - minVal + 1;
    assert(width > 1);

    int size = x.size() * width;

    densities_sum.init(size);
    count.init(size);
    for (int i=0; i<size; i++) {
      densities_sum[i] = 0;
      count[i] = 0;
    }
  }
  aAvgSD(Space& home, bool share, aAvgSD& a)
    : BranchingHeuristic<View>(home,share,a), minVal(a.minVal), width(a.width) {
    densities_sum.update(home,share,a.densities_sum);
    count.update(home,share,a.count);
  }
  virtual Candidate get(Space& home) {
    // densities_sum and count are shared between all the spaces. This way, we
    // don't need to allocate memory each time we want to calculate a new
    // choice. This also mean that we have to reset everything to 0
    for (int i=0; i<(*positions.get()).size()*width; i++) {
      densities_sum[i] = 0;
      count[i] = 0;
    }

    // isVarSub is true if it is inside at least one propagator that supports
    // cbs. If we have no record in the log for this variable, it will be false
    // and we don't want to consider it in our choice recommendation
    __gnu_cxx::hash_set<unsigned int> isVarSub;

    // For every propagator in the log
    for (Log::iterator it=log->begin(); it!=log->end(); ++it) {
      unsigned int prop_id = it->first;
      size_t nb_records = it->second.first;
      // Aggregation for every record concerning the propagator (i.e. every
      // (variable,value) pair and their corresponding density)
      for (unsigned int i=0; i<nb_records; ++i) {
        Record r = it->second.second[i];
        unsigned int pos = (*positions.get())[r.var_id] * width + r.val - minVal;
        densities_sum[pos] += r.density;
        count[pos] += 1;
      }
    }

    Candidate c;
    c.idx = -1;
    double best_density_moy = 0;

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
        int idx = i*width + val.val() - minVal;
//        if (count[idx]==0) {
//          printf("COUNT=0. prop_id=%i\n",);
//        }
        assert(count[idx] != 0);
        double dens_moy = densities_sum[idx] / (double)count[idx];
        if (dens_moy > best_density_moy) {
          c.idx = i;
          c.val = val.val();
          best_density_moy = dens_moy;
        }
      }
    }
    return c;
  }
};

template<class View>
class CBSBrancher : public Brancher {
  typedef typename BranchingHeuristic<View>::Candidate Candidate;
protected:
  ViewArray<View> x;
  aAvgSD<View> heur;
  ChangedPropagators changedProps;
  Log log;
public:
  CBSBrancher(Home home, ViewArray<View>& x0, const ChangedPropagators& spc)
    : Brancher(home), x(x0), heur(home,x0), changedProps(spc),
      log(Log::size_type(), Log::hasher(), Log::key_equal(),
          Log::allocator_type(home)) {
    // Because we must call the destructor of aAvgSD
    home.notice(*this,AP_DISPOSE);
  }
  static void post(Home home, ViewArray<View>& x,
                   const ChangedPropagators& spc) {
    (void) new (home) CBSBrancher(home,x,spc);
  }
  virtual size_t dispose(Space& home) {
    home.ignore(*this, AP_DISPOSE);
    // ~aAvgSD() calls ~SharedHashMap() which calls ~SharedHashMapObject() to
    // deallocate the hash map when the refcount of SharedHashMapObject is 0
    heur.~aAvgSD();
    (void) Brancher::dispose(home);
    return sizeof(*this);
  }
  CBSBrancher(Space& home, bool share, CBSBrancher& b)
    : Brancher(home,share,b), heur(home, share, b.heur),
      log(b.log.begin(), b.log.end(), Log::size_type(), Log::hasher(),
          Log::key_equal(), Log::allocator_type(home)) {
    x.update(home,share,b.x);
    changedProps.update(home,share,b.changedProps);

    for (Log::iterator it=log.begin(); it!=log.end(); ++it) {
      unsigned int prop_id = it->first;
      size_t count = b.log[prop_id].first;
      it->second.second = home.alloc<Record>(count);
      memcpy(it->second.second, b.log[prop_id].second, count*sizeof(Record));
    }
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
  virtual const Choice* choice(Space& home) {
    // Active propagators and the size we need for their log.
    // TODO: Est-ce que je peux avoir les active props de cette manière?
    // TODO: Je pense que oui. Les propagateurs disabled sont seulement mis
    // TODO: dans la queue idle sans faire de propagations, ce qui ne change
    // TODO: rien pour nous. Il faut cependant savoir si on utilise les
    // TODO: propagateurs idle.
    //
    // We considere a propagator as "active" only if he supports cbs().
    __gnu_cxx::hash_map<unsigned int, size_t> activeProps;
    for (int i=0; i<x.size(); i++) {
      for (SubscribedPropagators sp(x[i]); sp(); ++sp) {
        if (!sp.propagator().cbs(home,NULL)) continue;
        unsigned int prop_id = sp.propagator().id();
        bool contains = activeProps.find(prop_id) != activeProps.end();
        if (contains)
          activeProps[prop_id] += x[i].size();
        else
          activeProps[prop_id] = x[i].size();
      }
    }

    // We delete log elements corresponding to non active propagators
    {
      std::vector<unsigned int> propsToDelete;

      // Propagators to delete (not active and in log)
      for (Log::iterator it = log.begin(); it != log.end(); ++it)
        if (activeProps.find(it->first) == activeProps.end())
          propsToDelete.push_back(it->first);

      // We delete them from the log
      for (std::vector<unsigned int>::iterator it = propsToDelete.begin();
           it != propsToDelete.end(); ++it) {
        unsigned int prop_id = *it;
        // TODO: Trouver une manière de free comme il le faut (avec vrai grandeur, pas nb element en ce moment)...
//        home.free(log[prop_id].second, log[prop_id].first);
        log.erase(prop_id);
      }
    }

    // We specify the log that will be modified when the propagators use the
    // CBS::set()
    heur.set_log(&log);
    for (Propagators p(home, PropagatorGroup::all); p(); ++p) {
      if (!p.propagator().cbs(home,NULL)) continue;
      unsigned int prop_id = p.propagator().id();
      // Activity in propagator since last branching
      bool changed = changedProps.contains(prop_id);
      // Propagator already in the log?
      bool in_log = log.find(prop_id) != log.end();

      if (in_log) {
        if (changed)
          // We discard the previous entries by setting the count to 0 (we
          // thus reuse previous allocated memory. The number of records can't
          // grow).
          log[prop_id].first = 0;
        else
          // We continue, meaning we will use the previous entries
          continue;
      } else {
        // Otherwise, we need to allocate space for the records of the
        // new propagator
        log[prop_id] =
          std::make_pair(0, home.alloc<Record>(activeProps[prop_id]));
      }
      heur.set_current_prop(prop_id);
      p.propagator().cbs(home,&heur);
    }

    // Clear for the next round of propagations (use of advisor in our case)
    changedProps.clear();

    // We find the choice.
    Candidate c = heur.get(home);
    static int count=0;
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
};

enum CBSStrategy {
 MAX_SD,
 A_AVG_SD
};

void cbsbranch(Home home, const IntVarArgs& x) {
  if (home.failed()) return;
  ViewArray<Int::IntView> y(home,x);
  ChangedPropagators spc(home);
  ViewUpdateLooker<Int::IntView>::post(home,y,spc);
  CBSBrancher<Int::IntView>::post(home,y,spc);
}

void cbsbranch(Home home, const BoolVarArgs& x) {
  if (home.failed()) return;
  ViewArray<Int::BoolView> y(home,x);
  ChangedPropagators spc(home);
  ViewUpdateLooker<Int::BoolView>::post(home,y,spc);
  CBSBrancher<Int::BoolView>::post(home,y,spc);
}

#endif //__CBS_HPP__
