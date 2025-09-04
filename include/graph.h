#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <algorithm>
#include <sstream>
#include <string>

class GraphNetwork
{
    private:
        struct Edge
        {
            int node1;
            int node2;
            double weight;
        };
        std::vector<Edge> edges;
        std::vector<int> node_indices;

    public:
        GraphNetwork();
        ~GraphNetwork();
        void AddNode(int index);
        void AddEdge(int node1, int node2, double weight);
        double GetEdgeWeight(int node1, int node2);
        void UpdateEdge(int node1, int node2, double weight);
        void RemoveNode(int index);
        void RemoveEdge(int node1, int node2);
        std::string PrintGraph();
};

#endif
