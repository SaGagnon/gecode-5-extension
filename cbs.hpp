#ifndef __CBS_HPP__
#define __CBS_HPP__

#include <gecode/int.hh>
#include <ext/hash_map>
#include <ext/hash_set>
#include <cstring>

using namespace Gecode;


class SharedPropChanged : public LocalHandle {
protected:
  class SharedPropChangedObject : public LocalObject {
  public:
    typedef __gnu_cxx::hash_set<unsigned int, __gnu_cxx::hash<unsigned int>,
      __gnu_cxx::equal_to<unsigned int>,
      Gecode::space_allocator<unsigned int> > Set;
    Set changed;
  public:
    SharedPropChangedObject(Space& home)
      : LocalObject(home),
        changed(Set::size_type(), Set::hasher(), Set::key_equal(),
                Set::allocator_type(home)) {
    }
    SharedPropChangedObject(Space& home, bool share,
                            SharedPropChangedObject& o)
      : LocalObject(home,share,o),
        changed(o.changed.begin(), o.changed.end(),
                Set::size_type(), Set::hasher(), Set::key_equal(),
                Set::allocator_type(home)) {}
    virtual LocalObject* copy(Space& home, bool share) {
      return new (home) SharedPropChangedObject(home,share,*this);
    }
    // TODO: Essayer de delete ca. Ca devrait rien faire.
    virtual size_t dispose(Space& home) {
      return sizeof(*this);
    }
  };
public:
  SharedPropChanged(void)
    : LocalHandle() {}
  SharedPropChanged(Space& home)
    : LocalHandle(new (home) SharedPropChangedObject(home)) {}
  SharedPropChanged(const SharedPropChanged& spc)
    : LocalHandle(spc) {}
  void insert(unsigned int prop_id) {
    SharedPropChangedObject* o =
      static_cast<SharedPropChangedObject*>(object());
    o->changed.insert(prop_id);
  }
  bool contains(unsigned int prop_id) const {
    SharedPropChangedObject* o =
      static_cast<SharedPropChangedObject*>(object());
    return o->changed.find(prop_id) != o->changed.end();
  }
  void clear() {
    SharedPropChangedObject* o =
      static_cast<SharedPropChangedObject*>(object());
    o->changed.clear();
  }
  SharedPropChangedObject::Set* get(void) {
    SharedPropChangedObject* o =
      static_cast<SharedPropChangedObject*>(object());
    return &o->changed;
  }
};

template<class View>
class ViewUpdateLooker : public NaryPropagator<View,PC_GEN_NONE> {
protected:
  using NaryPropagator<View,PC_GEN_NONE>::x;

  Council<ViewAdvisor<View>> c;
  SharedPropChanged changedProps;

  /// Constructor for posting
  ViewUpdateLooker(Home home, ViewArray<View>& x, const SharedPropChanged& spc)
    : NaryPropagator<View,PC_GEN_NONE>(home,x), c(home), changedProps(spc) {
    for (int i=0; i<x.size(); i++) {
      (void) new (home) ViewAdvisor<View>(home, *this, c, x[i]);
    }
  }
  /// Constructor for cloning \a p
  ViewUpdateLooker(Space& home, bool share, ViewUpdateLooker<View>& p)
    : NaryPropagator<View,PC_GEN_NONE>(home,share,p) {
    c.update(home,share,p.c);
    changedProps.update(home,share,p.changedProps);
  }
public:
  static ExecStatus post(Home home, ViewArray<View>& x,
                         const SharedPropChanged & spc) {
    (void) new (home) ViewUpdateLooker<View>(home,x,spc);
    return ES_OK;
  }
  virtual ExecStatus advise(Space& home, Advisor &_a, const Delta& d) {
    ViewAdvisor<View>& a(static_cast<ViewAdvisor<View>&>(_a));
    View v(a.view());
//    std::cout << "advise ";
    for (SubscribedPropagators sp(v); sp(); ++sp) {
      changedProps.insert(sp.propagator().id());
    }
  }
  virtual ExecStatus propagate(Space&, const ModEventDelta&) {
    return ES_FIX;
  }
  virtual PropCost cost(const Space&, const ModEventDelta&) const {
    return PropCost::record();
  }
  virtual Actor* copy(Space& home, bool share) {
    return new (home) ViewUpdateLooker<View>(home,share,*this);
  }
  virtual size_t dispose(Space& home) {
    c.dispose(home);
    (void) NaryPropagator<View,PC_GEN_NONE>::dispose(home);
    return sizeof(*this);
  }
};


//******************************************************************************

struct Record { unsigned int var_id; int val; double density; };
typedef __gnu_cxx::hash_map< unsigned int, std::pair<size_t,Record*>,
  __gnu_cxx::hash<unsigned int>, __gnu_cxx::equal_to<unsigned int>,
  Gecode::space_allocator<unsigned int> > Log;

class SharedHashMap : public SharedHandle {
protected:
  class SharedHashMapObject : public SharedHandle::Object {
  public:
    typedef __gnu_cxx::hash_map<unsigned int, unsigned int> HashMap;
    HashMap hash_map;
  public:
    SharedHashMapObject(void) {}
    SharedHashMapObject(const SharedHashMapObject& shmo)
      : hash_map(shmo.hash_map) {}
    virtual Object* copy(void) const {
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
    object(new SharedHashMapObject());
  }
  SharedHashMapObject::HashMap* get(void) const {
    return &static_cast<SharedHashMapObject*>(object())->hash_map;
  }
};

template<typename T>
class SharedTable : public SharedHandle {
protected:
  class SharedTableObject : public SharedHandle::Object {
  public:
    typedef __gnu_cxx::hash_map<unsigned int, unsigned int> HashMap;
    T *table;
  public:
    SharedTableObject(size_t size) {
      table = new T[size];
    }
    SharedTableObject(const SharedTableObject& sto)
      : table(sto.table) {}
    virtual Object* copy(void) const {
      return new SharedTableObject(*this);
    }
    virtual ~SharedTableObject(void) {
      delete[] table;
    }
  };
public:
  SharedTable(void) {}
  SharedTable(const SharedTable& st)
    : SharedHandle(st) {}
  void init(size_t size) {
    assert(object() == NULL);
    object(new SharedTableObject(size));
  }
  T& operator[](int i) {
    return static_cast<SharedTableObject*>(object())->table[i];
  }
};

class BranchingHeuristic : public CBS {
public:
  struct Candidate {
    int idx;
    int val;
  };
  // TODO: Quoi faire si multithread? Est-ce que je fais juste clear
  // mon data global....????
//  static virtual void dispose() = 0;
  virtual Candidate get(void) = 0;
};

//class MaxSD : public BranchingHeuristic {
//public:
//  MaxSD() {
//    best_candidate.c.idx = -1;
//    best_candidate.density = 0;
//  }
//  virtual void set(unsigned int var_id, int val, double density) {
//    if (density > best_candidate.density) {
//      best_candidate.c.idx = var_id;
//      best_candidate.c.val = val;
//      best_candidate.density = density;
//    }
//  }
//  virtual Candidate get(void) {
//    return best_candidate.c;
//  }
//private:
//  struct {
//    Candidate c;
//    double density;
//  } best_candidate;
//};

template<class View>
class aAvgSD : public BranchingHeuristic {
private:
  SharedHashMap positions;
  Log *log;
  int current_prop;
public:
  aAvgSD(Space& home, const ViewArray<View>& x) {
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

    size_t size = (size_t)x.size() * width;

    densities_sum.init(size);
    count.init(size);
//    densities_sum = new double[size];
//    count = new int[size];
    for (int i=0; i<size; i++) {
      densities_sum[i] = 0;
      count[i] = 0;
    }
  }
  aAvgSD(Space& home, const aAvgSD& a)
    : positions(a.positions), minVal(a.minVal), width(a.width),
      densities_sum(a.densities_sum), count(a.count) {
//    int size = (unsigned int)(*positions.get()).size() * width;
//    densities_sum = home.alloc<double>(size);
//    count = home.alloc<int>(size);
//    memcpy(densities_sum, a.densities_sum, size * sizeof(double));
//    memcpy(count, a.count, size * sizeof(int));
  }
  void setup(Log *_log) {
    log = _log;
//    current_prop = -1;
//
//    best_candidate.c.var_id = -1;
//    best_candidate.density_moy = 0;

  }
  void set_current_prop(int _current_prop) {
    current_prop = _current_prop;
  }
  virtual void set(unsigned int var_id, int val, double density) {
    assert(current_prop != -1);
    //TODO changer ça dans une interface.
    Record r; r.var_id=var_id; r.val=val; r.density=density;
    size_t *count = &(*log)[current_prop].first;
    (*log)[current_prop].second[*count] = r;
    // TODO: essayer de mettre ça plus clean ici. (*count)++ ???
    *count = *count + 1;
  }
  virtual Candidate get(void) {
    // clearing
    for (int i=0; i<(*positions.get()).size()*width; i++) {
      densities_sum[i] = 0;
      count[i] = 0;
    }

    // For every propagator in the log
    for (Log::iterator it=log->begin(); it!=log->end(); ++it) {
      // For every record concerning the propagator
      unsigned int prop_id = it->first;
      size_t nb_records = it->second.first;
      for (unsigned int i=0; i<nb_records; ++i) {
        Record r = it->second.second[i];

        unsigned int idx = (*positions.get())[r.var_id];
//        if (idx==552 && r.val==18) {
//          printf("Log: prop=%i, pos=%i, val=%i, dens=%f\n",prop_id,idx,r.val,r.density);
//        }

        unsigned int pos = (*positions.get())[r.var_id] * width + r.val - minVal;
        densities_sum[pos] += r.density;
        count[pos] += 1;
      }
    }

    Candidate c;
    c.idx = -1;
    double best_density_moy = 0;

    // For each variable
    for (int i=0; i<(*positions.get()).size(); i++) {
      // For each value
      for (int j=0; j<width; j++) {
        int idx = i*width + j;
        double dens_moy = densities_sum[idx] / (double)count[idx];
        if (dens_moy > best_density_moy) {
          c.idx = i;
          c.val = j - minVal;
          best_density_moy = dens_moy;
        }
      }
    }


//    printf("dens:%f, ",best_density_moy);
    return c;
  }
private:
  int minVal;
  int width;
  SharedTable<double> densities_sum;
  SharedTable<int> count;
};

//******************************************************************************

template<class View>
class CBSBrancher : public Brancher {
protected:
  ViewArray<View> x;
  aAvgSD<View> heur;
  SharedPropChanged changedProps;
  Log log;
public:
  CBSBrancher(Home home, ViewArray<View>& x0, const SharedPropChanged& spc)
    : Brancher(home), x(x0), heur(home,x0), changedProps(spc),
      log(Log::size_type(), Log::hasher(), Log::key_equal(),
          Log::allocator_type(home)) {
    // Because we must call the destructor of aAvgSD
    home.notice(*this,AP_DISPOSE);
  }
  static void post(Home home, ViewArray<View>& x,
                   const SharedPropChanged& spc) {
    (void) new (home) CBSBrancher(home,x,spc);
  }
  virtual size_t dispose(Space& home) {
    home.ignore(*this, AP_DISPOSE);
    // ~aAvgSD() calls ~SharedHashMap() which calls ~SharedHashMapObject() to
    // deallocate the hash map when the refcount of SharedHashMapObject is 0
    heur.~aAvgSD();
//    log.~Log();
    (void) Brancher::dispose(home);
    return sizeof(*this);
  }
  CBSBrancher(Space& home, bool share, CBSBrancher& b)
    : Brancher(home,share,b), heur(home, b.heur),
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
  virtual Choice* choice(Space& home) {
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

    //--------------------------------------------------------------------------
//    std::vector<unsigned int> __active_props;
//    for (auto x : activeProps)
//      __active_props.push_back(x.first);
//    std::sort(__active_props.begin(), __active_props.end());
//
//    std::vector<unsigned int> __space_active_props;
//      for (Propagators p(home, PropagatorGroup::all); p(); ++p)
//        __space_active_props.push_back(p.propagator().id());
//    std::sort(__space_active_props.begin(), __space_active_props.end());
//    assert(__active_props.size() == __space_active_props.size());

//    for (int i=0; i<__active_props.size(); i++)
//      assert(__active_props[i] == __space_active_props[i]);
    //--------------------------------------------------------------------------
//
//    printf("in log                   :");
//    for (auto x : log)
//      printf("%i, ",x.first);
//    printf("\n");
//
//    printf("active props (subscribed):");
//    for (auto x : activeProps)
//      printf("%i, ",x.first);
//    printf("\n");
//
//    printf("active props (space)     :");
//    for (auto x : __space_active_props)
//      printf("%i, ",x);
//    printf("\n");
//
//    printf("changed props            :");
//    for (auto x : *changedProps.get())
//      printf("%i, ",x);
//    printf("\n");

//    Log backup = log;

    // We delete log elements corresponding to non active propagators
    {
//      __gnu_cxx::hash_set<unsigned int> activeProps;
      std::vector<unsigned int> propsToDelete;

      // Active propagators
//      for (Propagators p(home, PropagatorGroup::all); p(); ++p)
//        activeProps.insert(p.propagator().id());

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
    heur.setup(&log);
    for (Propagators p(home, PropagatorGroup::all); p(); ++p) {
//      if (p.propagator().id()==36) printf("BOUYA\n");
      if (!p.propagator().cbs(home,NULL)) continue;
      unsigned int prop_id = p.propagator().id();
      // Activity in propagator since last branching
      bool changed = changedProps.contains(prop_id);
      // Propagator already in the log?
      bool in_log = log.find(prop_id) != log.end();

//      if (prop_id==36) printf("changed=%i, in_log=%i\n",changed,in_log);

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
    BranchingHeuristic::Candidate c = heur.get();
//    log=backup;


//    printf("idx:%i, val:%i\n",c.idx,c.val);

//    exit(-1);
//    for (int i=0; i<x.size(); i++)
//      if (x[i].id() == c.var_id)
//        return new PosValChoice<int>(*this,2,i,c.val);
    assert(!x[c.idx].assigned());
    assert(x[c.idx].in(c.val));
    return new PosValChoice<int>(*this,2,c.idx,c.val);


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
    if (a == 0) {
//      printf("x[%i]=%i\n\n",pos,val);
      return me_failed(x[pos].eq(home,val)) ? ES_FAILED : ES_OK;
    }
    else {
//      printf("x[%i]!=%i\n\n",pos,val);
      return me_failed(x[pos].nq(home,val)) ? ES_FAILED : ES_OK;
    }
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
  SharedPropChanged spc(home);
  ViewUpdateLooker<Int::IntView>::post(home,y,spc);
  CBSBrancher<Int::IntView>::post(home,y,spc);
}

void cbsbranch(Home home, const BoolVarArgs& x) {
  if (home.failed()) return;
  ViewArray<Int::BoolView> y(home,x);
  SharedPropChanged spc(home);
  ViewUpdateLooker<Int::BoolView>::post(home,y,spc);
  CBSBrancher<Int::BoolView>::post(home,y,spc);
}

#endif //__CBS_HPP__
