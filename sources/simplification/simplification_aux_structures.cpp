#include "simplification_aux_structures.h"

void w_contraction_parameters::init_history_map(string path, ivect &new_v_pos)
{
    min_initial_weight=INT_MAX, max_initial_weight=0, avg_initial_weight=0;

    string line;
    int v1, v2;
    ifstream input(path.c_str());
    while( input.good() )
    {
        getline(input,line);

        if(input.eof())
            break;

        vector<string> tokens;
        string_management::tokenize(line,tokens," "); /// check the separator for each file type

        ivect e;
        v1=atoi(tokens[0].c_str());
        v2=atoi(tokens[1].c_str());
        e.push_back(new_v_pos[v1]);
        e.push_back(new_v_pos[v2]);
        sort(e.begin(),e.end());

        double w = atof(tokens[2].c_str());
        edge_history eh(w);
        e_history.insert(make_pair(e,eh));

        if(min_initial_weight > w)
            min_initial_weight = w;
        if(max_initial_weight < w)
            max_initial_weight = w;
        avg_initial_weight += w;
    }

    avg_initial_weight /= (double)e_history.size();

    cout<<"[STATS] initial edge-weights"<<endl;
    cout<<"     min: "<<min_initial_weight<<" avg: "<<avg_initial_weight<<" max: "<<max_initial_weight<<endl;

    double min_h_avg = (avg_initial_weight-min_initial_weight) / 2.0;
    double max_h_avg = (max_initial_weight-avg_initial_weight) / 2.0;

    int num_edges_below_avg = 0;

    ivect e_histogram;
    e_histogram.assign(6,0);
    for(edge_history_map::iterator it=e_history.begin(); it!= e_history.end(); ++it)
    {
        edge_history &eh = it->second;
        if(eh.initial_weight==min_initial_weight)
            e_histogram[0]++;
        else if(min_initial_weight < eh.initial_weight && eh.initial_weight <= (avg_initial_weight-min_h_avg))
            e_histogram[1]++;
        else if((avg_initial_weight-min_h_avg) < eh.initial_weight && eh.initial_weight <= avg_initial_weight)
            e_histogram[2]++;
        else if(avg_initial_weight < eh.initial_weight && eh.initial_weight <= (avg_initial_weight+max_h_avg))
            e_histogram[3]++;
        else if((avg_initial_weight+max_h_avg) < eh.initial_weight && eh.initial_weight < max_initial_weight)
            e_histogram[4]++;
        else if(eh.initial_weight==max_initial_weight)
            e_histogram[5]++;

        if(eh.initial_weight < avg_initial_weight)
            num_edges_below_avg++;
    }
    cout<<"[STATS] edge-weights distribution (tot = "<<e_history.size()<<")"<<endl;
    cout<<"     ";
    for(unsigned i=0; i<e_histogram.size(); i++)
    {
        cout<<e_histogram[i]<<" ("<<(e_histogram[i]/(double)e_history.size())<<"%)";
        if(i==2)
            cout<<endl<<"     ";
        else
            cout<<" - ";
    }
    cout<<endl;
    cout<<"     edges below avg: "<<num_edges_below_avg<<endl;
}

void w_contraction_parameters::compute_summation_history_map_stats()
{
//    double min_initial_weight=INT_MAX, max_initial_weight=0, avg_initial_weight=0;
    double min_sum_weight=INT_MAX, max_sum_weight=0, avg_sum_weight=0;

    for(edge_history_map::iterator it=e_history.begin(); it!= e_history.end(); ++it)
    {
        edge_history &eh = it->second;

//        if(min_initial_weight > eh.initial_weight)
//            min_initial_weight = eh.initial_weight;
//        if(max_initial_weight < eh.initial_weight)
//            max_initial_weight = eh.initial_weight;
//        avg_initial_weight += eh.initial_weight;

        if(min_sum_weight > eh.weight_summation)
            min_sum_weight = eh.weight_summation;
        if(max_sum_weight < eh.weight_summation)
            max_sum_weight = eh.weight_summation;
        avg_sum_weight += eh.weight_summation;
    }

//    avg_initial_weight /= (double)e_history.size();
    avg_sum_weight /= (double)e_history.size();

    cout<<"[STATS] summation edge-weights"<<endl;
    cout<<"        min: "<<min_sum_weight<<" avg: "<<avg_sum_weight<<" max:"<<max_sum_weight<<endl;

    double min_h_avg = (avg_sum_weight-min_sum_weight) / 2.0;
    double max_h_avg = (avg_sum_weight-max_sum_weight) / 2.0;

    ivect e_histogram;
    e_histogram.assign(6,0);
    for(edge_history_map::iterator it=e_history.begin(); it!= e_history.end(); ++it)
    {
        edge_history &eh = it->second;
        if(eh.weight_summation==min_sum_weight)
            e_histogram[0]++;
        else if(min_sum_weight < eh.weight_summation && eh.weight_summation <= (avg_sum_weight-min_h_avg))
            e_histogram[1]++;
        else if((avg_sum_weight-min_h_avg) < eh.weight_summation && eh.weight_summation <= avg_sum_weight)
            e_histogram[2]++;
        else if(avg_sum_weight < eh.weight_summation && eh.weight_summation <= (avg_sum_weight+max_h_avg))
            e_histogram[3]++;
        else if((avg_sum_weight+max_h_avg) < eh.weight_summation && eh.weight_summation < max_sum_weight)
            e_histogram[4]++;
        else if(eh.weight_summation==max_sum_weight)
            e_histogram[5]++;
    }
    cout<<"[STATS] summation edge-weights distribution (tot = "<<e_history.size()<<")"<<endl;
    cout<<"        ";
    for(auto i=e_histogram.begin(); i!=e_histogram.end(); ++i)
        cout<<*i<<" ("<<(*i/(double)e_history.size())<<"%) - ";
    cout<<endl;
}
