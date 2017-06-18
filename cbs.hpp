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

//#define SQL

#ifdef SQL
#include "sql-interface.hh"
#endif

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
//
unsigned int varpos(const VADesc& xD, unsigned int var_id) {
  return xD.positions[var_id];
}

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
class BranchingHeuristic : public SolnDistribution {
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
  BranchingHeuristic(Space& home, const ViewArray<View>& x0)
    : x(x0), xD(x0) {}
  /**
   * Constructor for cloning spaces
   *
   * @param home New constructed space
   * @param share Can we share internal data with the new clone
   * @param bh Object that is being cloned from the parent space
   */
  BranchingHeuristic(Space& home, bool share, BranchingHeuristic& bh)
    : xD(home,share,bh.xD) {
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
  virtual void setMarginalDistribution(unsigned int var_id, int val,
                                       double density) {
    assert(current_prop != -1);
    assert(density>0 && density<1);
    assert(x[xD.positions[var_id]].in(val));
    Record r; r.var_id=var_id; r.val=val; r.density=density;
    // Number of records for the current propagator
    size_t *nb_record = &(*logDensity)[current_prop].first;
    (*logDensity)[current_prop].second[(*nb_record)++] = r;
  }
  virtual void setSupportSize(double count) {
    assert(current_prop != -1);
    (*logProp)[current_prop].second = count;
  }

public:
  void for_every_log_entry(
    std::function<void(unsigned int,double,unsigned int, int, double)> f) {
    // For every propagator in the log
    for (auto prop : *logDensity) {
      unsigned int prop_id = prop.first;
      double slnCnt = (*logProp)[prop_id].second;
      size_t nb_records = prop.second.first;
      Record *records = prop.second.second;
      for (unsigned int i=0; i<nb_records; ++i) {
        Record *r = &records[i];
        f(prop_id, slnCnt, r->var_id, r->val, r->density);
      }
    }
  }

  void for_every_varIdx_val(Space& home,
                            std::function<void(unsigned int, int)> f) {
    for (unsigned int i=0; i<x.size(); i++) {
      bool instrumented = false;
      for (SubscribedPropagators sp(x[i]); sp(); ++sp) {
        if (sp.propagator().cbs(home,NULL)) {
          instrumented = true;
          break;
        }
      }
      if (!instrumented) continue;

      if (x[i].assigned()) continue;
      for (Int::ViewValues<View> val(x[i]); val(); ++val)
        f(x[i].id(),val.val());
    }
  }

  virtual Candidate getChoice(Space& home) = 0;
};

#define USING_BH \
  using BranchingHeuristic<View>::BranchingHeuristic; \
  using BranchingHeuristic<View>::x; \
  using BranchingHeuristic<View>::xD; \
  using BranchingHeuristic<View>::for_every_log_entry; \
  using BranchingHeuristic<View>::for_every_varIdx_val; \
  using BranchingHeuristic<View>::logDensity; \
  using BranchingHeuristic<View>::logProp; \
  typedef typename BranchingHeuristic<View>::Candidate Candidate;


template<class View>
class maxSD : public BranchingHeuristic<View> {
  USING_BH
public:
  virtual Candidate getChoice(Space& home) {
    struct Best {int var_id; int val; double dens; }
      best_candidate{-1,0,0};

    for_every_log_entry([&](unsigned int prop_id, double slnCnt,
                            unsigned int var_id, int val, double dens) {
      if (dens>best_candidate.dens)
          best_candidate = Best{var_id, val, dens};
    });
    return Candidate{xD.positions[best_candidate.var_id],best_candidate.val};
  }
};


template<class View>
class aAvgSD : public BranchingHeuristic<View> {
  USING_BH
protected:
  SharedArray<double> tot_dens;
  SharedArray<int> prop_count;
public:
  aAvgSD(Space& home, const ViewArray<View>& x)
    : BranchingHeuristic<View>(home,x) {
    int size = xD.size * xD.width;
    tot_dens.init(size);
    prop_count.init(size);
    for (int i=0; i<size; i++) {
      tot_dens[i] = 0;
      prop_count[i] = 0;
    }
  }
  aAvgSD(Space& home, bool share, aAvgSD& h)
    : BranchingHeuristic<View>(home,share,h) {
    tot_dens.update(home,share,h.tot_dens);
    prop_count.update(home,share,h.prop_count);
  }

  virtual Candidate getChoice(Space& home) {
    for (int i=0; i<tot_dens.size(); i++) {
      tot_dens[i] = 0;
      prop_count[i] = 0;
    }

    for_every_log_entry([&](unsigned int prop_id, double slnCnt,
                            unsigned int var_id, int val, double dens) {
      unsigned int idx = varvalpos(xD,var_id,val);
      tot_dens[idx] += dens;
      prop_count[idx] += 1;
    });

    struct Best { int var_id; int val; double dens_moy;
    } best_candidate{-1,0,0};

    for_every_varIdx_val(home, [&](unsigned var_id, int val) {
      unsigned int idx = varvalpos(xD,var_id,val);
      double dens_moy = tot_dens[idx] / (double)prop_count[idx];
      if (dens_moy > best_candidate.dens_moy)
        best_candidate = Best{var_id,val,dens_moy};
    });

    assert(best_candidate.var_id != -1);
    return Candidate{xD.positions[best_candidate.var_id],best_candidate.val};
  }
};


const size_t SIZE = 30*30*30;

template<class View>
class ai : public BranchingHeuristic<View> {
  USING_BH
protected:
  static double sum_dens[SIZE];
  static int nb_prop[SIZE];
  static int var_dom_size[SIZE];

  static double sum_slnCnt_x_dens[SIZE];
  static double sum_slnCnt[SIZE];
  static double slnCnt[SIZE];

  static double maxsd[SIZE];
  static double aAvgSD[SIZE];
  static double maxRelSD[SIZE];
  static double maxRelRatio[SIZE];
  static double wSCAvg[SIZE];
  static double wAntiSCAvg[SIZE];

  static double ctrTightness_x_dens[SIZE];
  static double sum_ctrTightness[SIZE];
  static double wTAvg[SIZE];
  static double wAntiTAvg[SIZE];

  static double sum_prop_card_prod_x_density[SIZE];
  static double sum_prop_card_prod[SIZE];
  static double wDAvg[SIZE];
public:
  virtual Candidate getChoice(Space& home) {
    size_t size = (size_t)(xD.size * xD.width);

    /// CLEAR
    for (int i=0; i<SIZE; i++) {
      sum_dens[i] = 0;
      nb_prop[i] = 0;
      var_dom_size[i] = 0;
      sum_slnCnt_x_dens[i] = 0;
      sum_slnCnt[i] = 0;
      slnCnt[i] = 0;
      maxsd[i] = 0;
      aAvgSD[i] = 0;
      maxRelSD[i] = 0;
      maxRelRatio[i] = 0;
      wSCAvg[i] = 0;
      wAntiSCAvg[i] = 0;
      wTAvg[i] = 0;
      wAntiTAvg[i] = 0;

      ctrTightness_x_dens[i] = 0;
      sum_ctrTightness[i] = 0;

      sum_prop_card_prod_x_density[i] = 0;
      sum_prop_card_prod[i] = 0;
      wDAvg[i] = 0;
    }
    /// CLEAR



    std::map<std::pair<unsigned int, unsigned int>, double> var_dens_entropy;

    std::unordered_map<unsigned int, double> prop_card_prod;
    std::unordered_map<unsigned int, double> prop_proj_tightness;

    /**
     * Computation
     */

    // Product of domain of every variables in each proapgators
    for (int i=0; i<x.size(); i++) {
      for (SubscribedPropagators sp(x[i]); sp(); ++sp) {
        if (sp.propagator().cbs(home, NULL)) {
          unsigned int prop_id = sp.propagator().id();
          if (prop_card_prod.find(prop_id) == prop_card_prod.end())
            prop_card_prod[prop_id] = 1;
          prop_card_prod[prop_id] *= x[i].size();
        }
      }
    }

    // Tightness of each propagator
    for (Propagators p(home, PropagatorGroup::all); p(); ++p) {
      if (p.propagator().cbs(home, NULL)) {
        unsigned int prop_id = p.propagator().id();
        double slnCnt = (*logProp)[prop_id].second;
        prop_proj_tightness[prop_id] = slnCnt / prop_card_prod[prop_id];
      }
    }

    for_every_log_entry([&](unsigned int prop_id, double slnCnt,
                            unsigned int var_id, int val, double density) {
      unsigned int idx = varvalpos(xD,var_id,val);
      maxsd[idx] = std::max(maxsd[idx], density);
      sum_dens[idx] += density;
      nb_prop[idx] += 1;
      var_dom_size[idx] = x[varpos(xD,var_id)].size();

      sum_slnCnt_x_dens[idx] += slnCnt * density;
      sum_slnCnt[idx] += slnCnt;

      auto key = std::make_pair(prop_id, var_id);
      if (var_dens_entropy.find(key) == var_dens_entropy.end())
        var_dens_entropy[key] = 0;
      var_dens_entropy[key] -= density*log(density) / log(var_dom_size[idx]);

      ctrTightness_x_dens[idx] = prop_proj_tightness[prop_id] * density;
      sum_ctrTightness[idx] += prop_proj_tightness[prop_id];

      sum_prop_card_prod_x_density[idx] += prop_card_prod[prop_id] * density;
      sum_prop_card_prod[idx] += prop_card_prod[prop_id];
    });


    for_every_log_entry([&](unsigned int prop_id, double slnCnt,
                            unsigned int var_id, int val, double density) {
      unsigned int idx = varvalpos(xD,var_id,val);
      wAntiSCAvg[idx] += (sum_slnCnt[idx] - slnCnt) * density / sum_slnCnt[idx];
      wAntiTAvg[idx] += (sum_ctrTightness[idx] - prop_proj_tightness[prop_id])
                        * density / sum_ctrTightness[idx];
    });

    for_every_varIdx_val(home, [&](unsigned var_id, int val) {
      unsigned int idx = varvalpos(xD,var_id,val);
      aAvgSD[idx] = sum_dens[idx] / (double)nb_prop[idx];
      maxRelSD[idx] = maxsd[idx] - (1.0/(double)var_dom_size[idx]);
      maxRelRatio[idx] = maxsd[idx] / (1.0/(double)var_dom_size[idx]);
      wSCAvg[idx] = sum_slnCnt_x_dens[idx] / sum_slnCnt[idx];
      wTAvg[idx] = ctrTightness_x_dens[idx] / sum_ctrTightness[idx];
      wDAvg[idx] = sum_prop_card_prod_x_density[idx] / sum_prop_card_prod[idx];
    });

    struct Best { int var_id; int val; double score;
    } best_candidate{-1,0,0};


    #ifdef SQL
    // TODO: Mettre un flag qui fait ça ou non ici...
    for (auto prop : (*logDensity)) {
      unsigned int prop_id = prop.first;
      double slnCnt = (*logProp)[prop_id].second;
      size_t nb_records = prop.second.first;
      Record *records = prop.second.second;
      for (unsigned int i=0; i<nb_records; ++i) {
        Record *r = &records[i];
        unsigned int idx = varvalpos(xD,r->var_id,r->val);
        unsigned int var_idx = varpos(xD,r->var_id);

        CBSDB::insert_varval_density_features(
          prop_id, r->var_id, r->val, r->density, slnCnt, sum_slnCnt[idx],
          aAvgSD[idx], var_dom_size[idx],
          var_dens_entropy[std::make_pair(prop_id, r->var_id)], maxRelSD[idx],
          maxRelRatio[idx], wSCAvg[idx], wAntiSCAvg[idx], wTAvg[idx],
          wAntiTAvg[idx], wDAvg[idx]);
      }
    }
    #endif


    for (auto prop : (*logDensity)) {
      unsigned int prop_id = prop.first;
      double slnCnt = (*logProp)[prop_id].second;
      size_t nb_records = prop.second.first;
      Record *records = prop.second.second;
      for (unsigned int i=0; i<nb_records; ++i) {
        Record *r = &records[i];
        unsigned int idx = varvalpos(xD,r->var_id,r->val);
        unsigned int var_idx = varpos(xD,r->var_id);


        // cbs_max_sd, 0-20, ECH10, reglog_full
//        double x = 0;
//        x += 0.222850436077 * maxsd[idx];
//        x += 4.60691374795 * aAvgSD[idx];
//        x += 0.0368418587122 * var_dom_size[idx];
//        x += 4.2498849779 *
//          var_dens_entropy[std::make_pair(prop_id,r->var_id)];
//        x += 2.29865995159 * maxRelSD[idx];
//        x += 1.26150264369 * maxRelRatio[idx];
//        x += 1.76381743389 * wSCAvg[idx];
//        x += 1.26487771943 * wAntiSCAvg[idx];
//        x += 0.0337647421938 * wTAvg[idx];
//        x += -0.877825562551 * wAntiTAvg[idx];
//        x += -0.430734215691 * wDAvg[idx];
//
//        double intercept = -8.65939782;
//        x += intercept;
//
//        double score = 1.0/(1.0 + exp(-x));

        double score = maxsd[idx];

        if (score > best_candidate.score)
          best_candidate = Best{r->var_id, r->val, score};
      }
    }


//    for_every_varIdx_val(home, [&](unsigned var_id, int val) {
//      unsigned int idx = varvalpos(xD,var_id,val);
//      unsigned int var_idx = varpos(xD,var_id);
//
////      0.790347429025
////      dens: 0.00692986720185
////      a_avg_sd: 5.12728647027
////      max_rel_sd: 8.52939316547
////      [-2.72276799]
////      double x = 0;
////      x += 0.189388381864 * maxsd[idx];
////      x += 6.1533440685 * aAvgSD[idx];
////      x += 0.0366048523955 * var_dom_size[idx];
////      x += 4.09477607086 * var_dens_entropy[std::make_pair(prop_id,var_id)];
////      x += 2.29259310496 * maxRelSD[idx];
////      x += 1.26064249745 * maxRelRatio[idx];
////      x += -0.93 * wSCAvg[idx];
////      x += 0.32 * wAntiSCAvg[idx];
//
////      double intercept = -2.72;
////      x += intercept;
//
////      double score = 1.0/(1.0 + exp(-x));
////      printf("%f\n",score);
//
//      double score = maxsd[idx];
//
//      if (score > best_candidate.score)
//        best_candidate = Best{var_id,val,score};
//    });

//    printf("score=%f\n",best_candidate.score);
    assert(best_candidate.var_id != -1);
    return Candidate{xD.positions[best_candidate.var_id],best_candidate.val};
  }

};

template<class View> double ai<View>::sum_dens[SIZE]{};
template<class View> int ai<View>::nb_prop[SIZE]{};
template<class View> int ai<View>::var_dom_size[SIZE]{};

template<class View> double ai<View>::sum_slnCnt_x_dens[SIZE]{};
template<class View> double ai<View>::sum_slnCnt[SIZE]{};
template<class View> double ai<View>::slnCnt[SIZE]{};

template<class View> double ai<View>::maxsd[SIZE]{};
template<class View> double ai<View>::aAvgSD[SIZE]{};
template<class View> double ai<View>::maxRelSD[SIZE]{};
template<class View> double ai<View>::maxRelRatio[SIZE]{};
template<class View> double ai<View>::wSCAvg[SIZE]{};
template<class View> double ai<View>::wAntiSCAvg[SIZE]{};

template<class View> double ai<View>::ctrTightness_x_dens[SIZE]{};
template<class View> double ai<View>::sum_ctrTightness[SIZE]{};
template<class View> double ai<View>::wTAvg[SIZE]{};
template<class View> double ai<View>::wAntiTAvg[SIZE]{};

template<class View> double ai<View>::sum_prop_card_prod_x_density[SIZE]{};
template<class View> double ai<View>::sum_prop_card_prod[SIZE]{};
template<class View> double ai<View>::wDAvg[SIZE]{};


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
                LogProp::key_equal(), LogProp::allocator_type(home)) {
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
  virtual bool status(const Space& home) const {
    Space& h = const_cast<Space&>(home);
    for (Propagators p(h, PropagatorGroup::all); p(); ++p)
      if (p.propagator().cbs(h, NULL))
        return true;
    return false;
  }
  virtual const Choice* choice(Space& home) {
    #ifdef SQL
    CBSDB::new_node();
    #endif

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
      bool changed = true;

      if (in_log) {
        changed = logProp[prop_id].first != activeProps[prop_id];
        if (changed) {
          // We discard the previous entries by setting the count to 0 (we
          // thus reuse previous allocated memory. The number of records can't
          // grow).
          logDensity[prop_id].first = 0;
          // TODO: Comment
          logProp[prop_id].first = (size_t)activeProps[prop_id];
        }
      } else {
        // Otherwise, we need to allocate space for the records of the
        // new propagator
        logDensity[prop_id] =
          std::make_pair(0, home.alloc<Record>(activeProps[prop_id]));
        logProp[prop_id] = std::make_pair(activeProps[prop_id], 0);
      }

      if (!in_log || changed) {
        heur.set_current_prop(prop_id);
        p.propagator().cbs(home,&heur);
      }
      #ifdef SQL
      CBSDB::new_propagator("", prop_id, "", logProp[prop_id].second);
      #endif
    }

//    #ifdef SQL
//    // TODO: Mettre un flag qui fait ça ou non ici...
//    for (auto prop : logDensity) {
//      unsigned int prop_id = prop.first;
//      double slnCnt = logProp[prop_id].second;
//      size_t nb_records = prop.second.first;
//      Record *records = prop.second.second;
//      for (unsigned int i=0; i<nb_records; ++i) {
//        Record *r = &records[i];
//        CBSDB::insert_varval_density(prop_id, r->var_id, r->val, r->density);
//      }
//    }
//    #endif

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
  W_SC_AVG,
  AI
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
      GECODE_NEVER;
//      CBSBrancher<View,maxRelSD>::post(home,y);
      break;
    case MAX_REL_RATIO:
      GECODE_NEVER;
//      CBSBrancher<View,maxRelRatio>::post(home,y);
      break;
    case A_AVG_SD:
      CBSBrancher<View,aAvgSD>::post(home,y);
      break;
    case W_SC_AVG:
      GECODE_NEVER;
//      CBSBrancher<View,wcSCAvg>::post(home,y);
      break;
    case AI:
      CBSBrancher<View,ai>::post(home,y);
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
