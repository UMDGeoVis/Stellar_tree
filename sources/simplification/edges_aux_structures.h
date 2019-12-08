#ifndef EDGES_AUX_STRUCTURES_H
#define EDGES_AUX_STRUCTURES_H

#include <list>

#include "utilities/basic_wrappers.h"

// this structure is used for sorting the priority_queue
// and defines the edge indicies + the edge weight
struct edge_weight
{
    ivect e;
    double weight;

    edge_weight() { weight = 0; }
    edge_weight(ivect &v, double w)
    {
        e = v;
        weight = w;
    }
    edge_weight(const edge_weight &orig)
    {
        e = orig.e;
        weight = orig.weight;
    }

    inline bool operator <(const edge_weight& ew) const { return (weight < ew.weight); }
    inline bool operator ==(const edge_weight& ew) const { return (e == ew.e && weight == ew.weight); }
    inline bool operator !=(const edge_weight& ew) const { return !(*this==ew); }
    inline friend std::ostream& operator<<(std::ostream& out, const edge_weight& p)
    {
        out << "[";
        for(unsigned i=0; i<p.e.size(); i++)
        {
            out << p.e[i] << " ";
        }
        out << "] -> "<<p.weight<<" ";
        return out;
    }

    inline ivect& get_edge() { return e; }
    inline double get_weight() { return weight; }
};

// this structure keeps the history of an edge
struct edge_history
{
   double initial_weight;
   double weight_summation;
   double last_added_weight;

   edge_history() { initial_weight = weight_summation = last_added_weight = 0; }
   edge_history(double iw)
   {
       initial_weight = iw;
       weight_summation = initial_weight;
       last_added_weight = 0;
   }
   inline friend std::ostream& operator<<(std::ostream& out, const edge_history& p)
   {
       out << "(";
       out << p.initial_weight << " ";
       out << p.weight_summation << " ";
       out << p.last_added_weight << ")";
       return out;
   }

//   edge_history(edge_history &orig)
//   {
//       initial_weight = orig.initial_weight;
//       weight_summation = orig.weight_summation;
//       last_added_weight = orig.last_added_weight;
//   }
};

typedef std::queue<ivect> edge_queue;
typedef std::priority_queue<edge_weight> edge_pqueue;
typedef std::map<ivect,edge_history> edge_history_map;
typedef std::list<edge_weight> edge_list;
typedef std::map<int,std::set<std::pair<int,short> > > rejected_edge_map;

#endif // EDGES_AUX_STRUCTURES_H
