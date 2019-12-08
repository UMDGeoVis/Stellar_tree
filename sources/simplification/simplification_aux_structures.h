/*
    This file is part of the Stellar library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)

    This project has been supported by the Italian Ministry of Education and
    Research under the PRIN 2009 program, and by the National Science Foundation
    under grant number IIS-1116747.

    The Stellar library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The Stellar library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Stellar library.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SIMPLIFICATION_AUX_STRUCTURES_H
#define SIMPLIFICATION_AUX_STRUCTURES_H

#include <unordered_map>
#include <bm/bm.h>

#include <stellar_tree/stellar_tree.h>
#include <validator/connectedness_validator_adjs.h>
#include <utilities/string_management.h>
#include <utilities/usage.h>

#include "edges_aux_structures.h"
#include "quadric_aux_structures.h"

class simplifier_parameters
{
public:
    inline void increment_counter(int d) { simplification_counter[d]++;}
    inline void increment_counter(int d, int v) { simplification_counter[d] += v;}
    inline int get_counter(int d) { return  simplification_counter[d]; }
    inline ivect& get_counters() { return simplification_counter; }    

    inline Stellar_Tree& get_tree() { return *tree; }

    inline bool keep_simplifying(Simplicial_Mesh &mesh)
    {
        int sum=0;
        for(unsigned i=0; i<simplification_counter.size(); i++)
            sum += simplification_counter[i];
        return (sum>0 && mesh.count_top_cells_num()>0);
    }

    inline void enable_debug_prints() { debug_print = true; } /// ONLY FOR DEBUG
    inline bool is_enable_debug_prints() { return debug_print; } /// ONLY FOR DEBUG
    inline void disable_debug_prints() { debug_print = false; } /// ONLY FOR DEBUG
    inline void enable_debug_mode() { debug_mode = true; } /// ONLY FOR DEBUG
    inline bool is_enable_debug_mode() { return debug_mode; } /// ONLY FOR DEBUG
    inline void disable_debug_mode() { debug_mode = false; } /// ONLY FOR DEBUG

    inline void add_checked_edge(ivect &e, bool is_contracted)
    {
        pair<ivect,bool> p;
        p.first = e;
        p.second = is_contracted;
        checked_edges.push_back(p);
    }
    inline void write_checked_edges(string path, ivect &v_orig)
    {
        stringstream ss, ssc;
        ss << path << "_check_order.edges";
        ssc << path << "_contracted.edges";

        FILE *output, *outputC;
        output = fopen(ss.str().c_str(),"w");
        outputC = fopen(ssc.str().c_str(),"w");
        fprintf(output,"%d \n",(int)checked_edges.size());
        while(!checked_edges.empty())
        {
            pair<ivect,bool> p = checked_edges.front();
            checked_edges.pop_front();
            /// here the range is from 0 to n-1 (for skeleton blockers)

            ///debug
//            if(((v_orig[p.first[0]-1]-1)==0 && (v_orig[p.first[1]-1]-1)==603) /*||
//                    ((v_orig[p.first[0]-1]-1)==1316 && (v_orig[p.first[1]-1]-1)==5232)*/)
//            {
//                cout<<p.first[0]<<" "<<p.first[1]<<" --> "<<p.second<<endl;
//                cout<<"original ids: "<<v_orig[p.first[0]-1]<<" "<<v_orig[p.first[1]-1]<<endl;
//            }

            fprintf(output,"%d %d \n",v_orig[p.first[0]-1]-1,v_orig[p.first[1]-1]-1);
            if(p.second) /// if the edge has been contracted
                /// conversely here the range is from 1 to n -> to be consistent with the rest of the visualization output
                fprintf(outputC,"%d %d \n",v_orig[p.first[0]-1],v_orig[p.first[1]-1]);
        }
        fclose(output);
        fclose(outputC);
    }

protected:
    ivect simplification_counter;
    list<pair<ivect,bool> > checked_edges;

    Stellar_Tree* tree;
    bool debug_print; /// ONLY FOR DEBUG
    bool debug_mode; /// ONLY FOR DEBUG

    simplifier_parameters() {}
    simplifier_parameters(Stellar_Tree *t, int num_top_types)
    {
        simplification_counter.clear();
        simplification_counter.assign(num_top_types,0);
        tree = t;
        debug_print = false; /// ONLY FOR DEBUG
        debug_mode = false; /// ONLY FOR DEBUG
    }
};

class contraction_parameters : public simplifier_parameters
{
public:
    contraction_parameters(int num_top_types) { simplification_counter.assign(num_top_types,0); } // dummy constructor
    contraction_parameters(Stellar_Tree *t, int num_top_types) : simplifier_parameters(t,num_top_types)
    {
        checked_edges = 0;

        vt_gen = get_et = clean_t_lists = get_vts = get_links = check_link = 0;
        cache_insert = contract_edge =0;
        get_edges = 0;
        check_simplex_in_cob = 0;
//        get_skeleton_simplices = check_skeleton_condition = update_skeleton = gen_skeleton = 0;

        edge_q_max_size = local_vt_max_references = cache_max_size = cache_max_references = 0;

        failed_d_simpl_lc_num = failed_vertices_lc_num = ok_lc = 0;

        failed_sub_simpl_lc_num.clear();
        failed_sub_simpl_lc_num.assign(num_top_types,0);

        failed_top_lc_num = 0;
//        ram_peak = MemoryUsage().getValue_in_MB(false);;

        min_err = INT_MAX; max_err = INT_MIN; avg_err = 0;
        num_edges = 0; min_c = 0; max_c = 0;
        histogram.assign(4,0);
    }
    inline void increment_contracted_edges_counter() { checked_edges++; }
    inline int get_contracted_edges_num() { return checked_edges; }

    inline void print_simplification_counters()
    {
//        cout<<"[STATS] simplified simplices"<<endl;
//        for(unsigned i=0; i<simplification_counter.size(); i++)
//        {
//            if(simplification_counter[i]>0)
//                cout<<"   "<<simplification_counter[i]<<" "<<i+1<<"-simplices simplified"<<endl;
//        }

        if(checked_edges >0)
            cerr<<"[STATS] contracted edges "<<checked_edges<<endl;
    }

    inline void print_simplification_partial_timings()
    {
        double tot = vt_gen+get_et+clean_t_lists+get_vts+cache_insert+get_links+cache_insert+check_link;
        tot += /*get_skeleton_simplices + check_skeleton_condition + update_skeleton + gen_skeleton+*/contract_edge+get_edges+check_simplex_in_cob;

        cerr<<"=== === ****** === ==="<<endl;
        cerr<<"[PARTIAL TIMES]"<<endl;
        cerr<<"       get TARGET edges: "<<get_edges<<" ("<<get_edges*100/tot<<"%)"<<endl;
        cerr<<"       VT generation: "<<vt_gen<<" ("<<vt_gen*100/tot<<"%)"<<endl;
        cerr<<"       get ET: "<<get_et<<" ("<<get_et*100/tot<<"%)"<<endl;
        cerr<<"       get VT: "<<get_vts<<" ("<<get_vts*100/tot<<"%)"<<endl;
        cerr<<"       TOT expand leaf representation: "<<(get_edges+vt_gen+get_et+get_vts)<<" ("<<(get_edges+vt_gen+get_et+get_vts)*100/tot<<"%)"<<endl;
        cerr<<"=== === === ==="<<endl;
        cerr<<"       get links (or target top simplices): "<<get_links<<" ("<<get_links*100/tot<<"%)"<<endl;
        cerr<<"       check simplex in ETop: "<<check_simplex_in_cob<<" ("<<check_simplex_in_cob*100/tot<<"%)"<<endl;
        cerr<<"       check contraction condition: "<<check_link<<" ("<<check_link*100/tot<<"%)"<<endl;
        cerr<<"       TOT check contraction condition: "<<(check_simplex_in_cob+check_link+get_links)<<" ("<<(check_simplex_in_cob+check_link+get_links)*100/tot<<"%)"<<endl;
        cerr<<"=== === === ==="<<endl;
        cerr<<"       TOT contract edges: "<<contract_edge<<" ("<<contract_edge*100/tot<<"%)"<<endl;
        cerr<<"=== === === ==="<<endl;
        cerr<<"       clean top lists: "<<clean_t_lists<<" ("<<clean_t_lists*100/tot<<"%)"<<endl;
        cerr<<"       VT cache insert: "<<cache_insert<<" ("<<cache_insert*100/tot<<"%)"<<endl;
        cerr<<"       TOT finalization leaf step: "<<(clean_t_lists+cache_insert)<<" ("<<(clean_t_lists+cache_insert)*100/tot<<"%)"<<endl;
        cerr<<"=== === === ==="<<endl;
//        cerr<<"       gen skeleton blocker: "<<gen_skeleton<<" ("<<gen_skeleton*100/tot<<"%)"<<endl;
//        cerr<<"       get skeleton simplices: "<<get_skeleton_simplices<<" ("<<get_skeleton_simplices*100/tot<<"%)"<<endl;
//        cerr<<"       check skeleton link condition: "<<check_skeleton_condition<<" ("<<check_skeleton_condition*100/tot<<"%)"<<endl;
//        cerr<<"       update skeleton blocker: "<<update_skeleton<<" ("<<update_skeleton*100/tot<<"%)"<<endl;
//        cerr<<"=== === === ==="<<endl;
        cerr<<"    PARTIAL TOT: "<<tot<<endl;
        cerr<<"[AUX STRUCTURES STATS/STORAGE]"<<endl;
        cerr<<"       MAX EDGE QUEUE SIZE/STORAGE: "<<edge_q_max_size<<" ("<<2*sizeof(int)*edge_q_max_size/(1024.0*1024.0)<<" MBs)"<<endl;
        cerr<<"       MAX LOCAL-VT SIZE/STORAGE: "<<local_vt_max_references<<" ("<<sizeof(int)*local_vt_max_references/(1024.0*1024.0)<<" MBs)"<<endl;
        cerr<<"       MAX CACHE SIZE/STORAGE: max_size -> "<<cache_max_size/*<<" refs -> "<<cache_max_references*/<<" ("<<(sizeof(int)*cache_max_size+cache_max_size*sizeof(int)*local_vt_max_references)/(1024.0*1024.0)<<" MBs)"<<endl;
        /// actually we cannot estimate the "real cache size" thus we pick the biggest local vt encoded and we estimate the cache requirements by multiplying for the max cache size
        cerr<<"=== === ****** === ==="<<endl;
        cerr<<"[SIMPLIFICATION STATISTICS]"<<endl;
        cerr<<"       VERIFIED LINK CONDITIONS: "<<ok_lc<<endl;
        if(this->failed_vertices_lc_num > 0)
            cerr<<"       FAILED LINK CONDITIONS ON VERTICES: "<<failed_vertices_lc_num<<endl;
        long tot_sub_failed_lc = 0;
        for(auto it=failed_sub_simpl_lc_num.begin(); it!=failed_sub_simpl_lc_num.end(); ++it)
            tot_sub_failed_lc += *it;
        if(tot_sub_failed_lc > 0)
        {
            cerr<<"       FAILED LINK CONDITIONS ON SUB-SIMPLICES: ";
            for(auto it=failed_sub_simpl_lc_num.begin(); it!=failed_sub_simpl_lc_num.end(); ++it)
                cerr<<*it<<" ";
            cerr<<endl;
        }
        if(this->failed_top_lc_num>0)
            cerr<<"       FAILED LINK CONDITIONS ON D-SIMPLICES: "<<failed_d_simpl_lc_num<<endl;
        cerr<<"       FAILED CHECKS ON TOP-SIMPLICES CONDITION: "<<failed_top_lc_num<<endl;
        cerr<<"=== === ****** === ==="<<endl;
    }

public: /// FOR DEBUG
    /// timings variables
    double vt_gen, get_et, clean_t_lists, get_vts, get_links, check_link;
    double cache_insert, contract_edge;
    double get_edges;
    double check_simplex_in_cob;
//    double get_skeleton_simplices, check_skeleton_condition, update_skeleton, gen_skeleton;

    /// storage variables
    long edge_q_max_size, local_vt_max_references, cache_max_size, cache_max_references;

    /// simplification stats
    long failed_vertices_lc_num, failed_d_simpl_lc_num, ok_lc;
    vector<long> failed_sub_simpl_lc_num;
    long failed_top_lc_num;

    float ram_peak;

//    vector<Matrix4x4> initialQuadric;
    double min_err, max_err, avg_err;
    int num_edges, min_c, max_c;
    double threshold;
    vector<int> histogram;

protected:
    int checked_edges;
};

class w_contraction_parameters : public contraction_parameters
{
public:
    w_contraction_parameters(Stellar_Tree *t, int num_top_types/*, bool no_homology*/) : contraction_parameters(t,num_top_types)
    {
        swap_contraction_order = false;
        change_state_edges = 0;
//        do_not_preserve_homology = no_homology;
    }

    inline void init_v_degrees_array(int num_v) { v_degrees.assign(num_v,0); }
    inline ivect& get_v_degrees_array() { return v_degrees; }
    inline int get_v_degree(int v) { return v_degrees[v-1]; } /// NOTA!

    void init_history_map(string path, ivect &new_v_pos);
    inline edge_history& get_edge_history(ivect &e) { return e_history[e]; }
    inline edge_history_map::iterator find(ivect &e) { return e_history.find(e); }
    inline edge_history_map::iterator begin() { return e_history.begin(); }
    inline edge_history_map::iterator end() { return e_history.end(); }
    inline void remove_edge(ivect &e) { e_history.erase(e); }
    inline void add_edge(ivect &e, edge_history &eh) { e_history.insert(make_pair(e,eh)); }

    inline void set_contraction_order(bool b) { swap_contraction_order = b; }
    inline bool swap_order() { return swap_contraction_order; }

    inline void increment_change_state_edges() { change_state_edges++; }
    inline void print_change_state_edges_num()
    {
        cerr<<"       EDGES THAT CHANGE STATE DURING SIMPLIFICATION PROCEDURE: "<<change_state_edges<<endl;
        cerr<<"=== === ****** === ==="<<endl;
    }

    inline void init_v_communities_array(int num_v) { v_communities.assign(num_v,list<int>()); }
    inline vector<list<int> >& get_v_communities_array() { return v_communities; }
    inline int get_v_community_size(int v) { return v_communities[v-1].size(); } /// NOTA!
    inline void add_v_to_community(int v, int v_to_add) { v_communities[v-1].push_front(v_to_add); }

    template<class C, class T> void compute_v_community_stats(Mesh<C,T> &mesh);
    template<class C, class T> void write_v_community(string path, ivect &v_orig, Mesh<C, T> &mesh);
    template<class C, class T> void compute_v_degree_stats(Mesh<C, T> &mesh);

    void compute_summation_history_map_stats();

//    inline bool is_do_not_preserve_homology() { return do_not_preserve_homology; }
    inline bool is_weight_below_avg(ivect &e) { return e_history[e].initial_weight < avg_initial_weight; }

protected:
    edge_history_map e_history;
    double min_initial_weight, max_initial_weight, avg_initial_weight;

    ivect v_degrees;
    bool swap_contraction_order;
    int change_state_edges;
    vector<list<int> > v_communities;
//    bool do_not_preserve_homology;
};

template<class C, class T> void w_contraction_parameters::compute_v_community_stats(Mesh<C,T> &mesh)
{
    int min_v_pushed=INT_MAX, max_v_pushed=0;
    double avg_v_pushed=0;
    int min_v_adopted=INT_MAX, max_v_adopted=0;
    double avg_v_adopted=0;
    int min_community=INT_MAX, max_community=0;
    double avg_community=0;
    double num_vertices = 0;

    for(unsigned v=0; v<v_communities.size(); v++)
    {
        if(!mesh.is_vertex_removed(v+1))
        {
            int local_community_summation=0;

            list<int> &community = v_communities[v];
            if(min_v_pushed > community.size())
                min_v_pushed = community.size();
            if(max_v_pushed < community.size())
                max_v_pushed = community.size();
            avg_v_pushed += community.size();
            local_community_summation +=community.size();

            for(list<int>::iterator it=community.begin(); it!=community.end(); ++it)
            {
                int c_size = v_communities[*it-1].size();
                if(min_v_adopted > c_size)
                    min_v_adopted = c_size;
                if(max_v_adopted < c_size)
                    max_v_adopted = c_size;
                avg_v_adopted += c_size;
                local_community_summation += c_size;
            }

            if(min_community > local_community_summation)
                min_community = local_community_summation;
            if(max_community < local_community_summation)
                max_community = local_community_summation;
            avg_community += local_community_summation;

            num_vertices++;
        }
    }

    avg_v_pushed /= num_vertices;
    avg_v_adopted /= num_vertices;
    avg_community /= num_vertices;

    cout<<"[STATS] "<<endl;
    cout<<"   [v_pushed] min: "<<min_v_pushed<<" avg: "<<avg_v_pushed<<" max:"<<max_v_pushed<<endl;
    cout<<"   [v_adopted] min: "<<min_v_adopted<<" avg: "<<avg_v_adopted<<" max:"<<max_v_adopted<<endl;
    cout<<"   [community_size] min: "<<min_community<<" avg: "<<avg_community<<" max:"<<max_community<<endl;
}

template<class C, class T> void w_contraction_parameters::write_v_community(string path, ivect &v_orig, Mesh<C, T> &mesh)
{
    stringstream ss;
    ss << path << ".v_communities";

    ofstream output(ss.str());

    for(unsigned v=0; v<v_communities.size(); v++)
    {
        if(!mesh.is_vertex_removed(v+1))
        {
            list<int> &community = v_communities[v];
            list<int> adopted;

            output << v_orig[v] << endl;
            output << "   ";

            for(list<int>::iterator it=community.begin(); it!=community.end(); ++it)
            {
                output << v_orig[*it-1] << " ";

                list<int> &a_community = v_communities[*it-1];
                for(list<int>::reverse_iterator ait=a_community.rbegin(); ait!=a_community.rend(); ++ait)
                {
                    adopted.push_front(*ait);
                }
            }
            output << endl << "   ";
            for(list<int>::iterator it=adopted.begin(); it!=adopted.end(); ++it)
            {
                output << v_orig[*it-1] << " ";
            }
            output << endl;
        }
    }
}

template<class C, class T> void w_contraction_parameters::compute_v_degree_stats(Mesh<C, T> &mesh)
{
    int min=INT_MAX, max=0;
    double avg=0, num=0;

    for(unsigned v=0; v<v_degrees.size(); v++)
    {
        if(!mesh.is_vertex_removed(v+1))
        {
            if(min > v_degrees[v])
                min = v_degrees[v];
            if(max < v_degrees[v])
                max = v_degrees[v];
            avg += v_degrees[v];

            num++;

//                if(v_degrees[v]==0)
//                    cerr<<v+1<<" degree nullo"<<endl;
        }
    }

    avg /= num;

    cout<<"[STATS] v_degree min: "<<min<<" avg: "<<avg<<" max: "<<max<<endl;
}

///////////////

class simplification_stats
{
public:
    simplification_stats(int num_entries)
    {
        time_simpl = time_up_mesh = time_up_tree = 0;

        num_top_d_simplices.assign(num_entries,0);
        tot_simplification.assign(num_entries,0);
    }

    inline void start_timer() { time.start(); }
    inline void stop_timer() { time.stop(); }

    inline void increment_simplification_timer() { time_simpl += time.get_elapsed_time(); }
    inline void increment_updating_mesh_timer() { time_up_mesh += time.get_elapsed_time(); }
    inline void increment_updating_tree_timer() { time_up_tree += time.get_elapsed_time(); }

    inline void set_num_simplices(int d, int num) { this->num_top_d_simplices[d] = num; }

    inline void increment_simplification_counter(int d, int num) { this->tot_simplification[d] += num; }
    inline void reset_simplification_counter()
    {
        for(unsigned d=0; d<tot_simplification.size(); d++)
            tot_simplification[d] = 0;
    }
    inline int get_simplification_counter(int d) { return this->tot_simplification[d]; }

    inline bool all_top_cells_simplified()
    {
        for(unsigned d=0; d<tot_simplification.size(); d++)
        {
            if(this->tot_simplification[d] != this->num_top_d_simplices[d])
                return false;
        }
        return true;
    }
    inline bool all_top_d_cells_simplified(int d) { return (this->tot_simplification[d]==this->num_top_d_simplices[d]); }

    inline void print_partial_simplification_stats()
    {
        for(unsigned d=0; d<tot_simplification.size(); d++)
        {
            if(tot_simplification.size() > 1)
                cout<<"[STAT] simplified "<<d+1<<"-simplices"<<endl;
            cout<<"[STAT] removed: "<<tot_simplification[d]<<" of "<<num_top_d_simplices[d];
            int perc = (tot_simplification[d]/(double)num_top_d_simplices[d])*100;
            cout<<" ("<<perc<<"%)"<<endl;
        }
    }

    inline void print_simplification_timings()
    {
        double tot = time_simpl+time_up_mesh+time_up_tree;
        cerr<<"[SUMMARY TIMES]"<<endl;
        cerr<<"       simplification: "<<time_simpl<<" ("<<time_simpl*100/tot<<"%)"<<endl;
        cerr<<"       updating mesh: "<<time_up_mesh<<" ("<<time_up_mesh*100/tot<<"%)"<<endl;
        cerr<<"       updating tree: "<<time_up_tree<<" ("<<time_up_tree*100/tot<<"%)"<<endl;
        cerr<<"    TOTAL: "<<tot<<endl;
    }

private:
    Timer time;

    double time_simpl;
    double time_up_tree;
    double time_up_mesh;

    ivect num_top_d_simplices;
    ivect tot_simplification;
};

#endif // SIMPLIFICATION_AUX_STRUCTURES_H
