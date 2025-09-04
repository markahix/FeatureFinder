#include "graph.h"

// #include <iostream>

// Constructor
GraphNetwork::GraphNetwork()
{
    edges = {};
    node_indices = {};
}

// Destructor
GraphNetwork::~GraphNetwork()
{

}

void GraphNetwork::AddNode(int index)
{
    if (std::find(node_indices.begin(),node_indices.end(), index) == node_indices.end())
    {
        node_indices.push_back(index);
    }
}

void GraphNetwork::UpdateEdge(int node1, int node2, double weight)
{
    for (Edge &tmp_edge : edges)
    {
        if ((tmp_edge.node1 == node1) && (tmp_edge.node2 == node2))
        {
            tmp_edge.weight = weight;
        }
    }
}

void GraphNetwork::AddEdge(int node1, int node2, double weight)
{
    AddNode(node1);
    AddNode(node2);
    int idx_a;
    int idx_b;
    if (node1 < node2)
    {
        idx_a = node1;
        idx_b = node2;
    }
    else
    {
        idx_a = node2;
        idx_b = node1;
    }
    Edge tmp_edge = {idx_a, idx_b, weight};
    for (Edge ch_edge : edges)
    {
        if ((ch_edge.node1 == idx_a) && (ch_edge.node2 == idx_b))
        {
            UpdateEdge(idx_a, idx_b, weight);
            return;
        }
    }
    edges.push_back(tmp_edge);
}

void GraphNetwork::RemoveNode(int index)
{
    for (int i=0; i < node_indices.size(); i++)
    {
        if (node_indices[i] == index)
        {
            node_indices.erase(node_indices.begin()+i);
            break;
        }
    }
    std::vector<std::vector<int>> edge_pairs = {};
    for (int i = 0; i < edges.size(); i++)
    {
        if ((edges[i].node1 == index) || (edges[i].node2 == index))
        {
            edge_pairs.push_back({edges[i].node1, edges[i].node2});
        }
    }
    for (std::vector<int> pair : edge_pairs)
    {
        RemoveEdge(pair[0],pair[1]);
    }
}

void GraphNetwork::RemoveEdge(int node1, int node2)
{
    int idx_a;
    int idx_b;
    if (node1 < node2)
    {
        idx_a = node1;
        idx_b = node2;
    }
    else
    {
        idx_a = node2;
        idx_b = node1;
    }
    for (int i = 0; i < edges.size(); i++)
    {
        if ((edges[i].node1 == idx_a) && (edges[i].node2 == idx_b))
        {
            edges.erase(edges.begin() + i);
            break;
        }
    }
}

double GraphNetwork::GetEdgeWeight(int node1, int node2)
{
    int idx_a;
    int idx_b;
    if (node1 < node2)
    {
        idx_a = node1;
        idx_b = node2;
    }
    else
    {
        idx_a = node2;
        idx_b = node1;
    }
    for (int i = 0; i < edges.size(); i++)
    {
        if ((edges[i].node1 == idx_a) && (edges[i].node2 == idx_b))
        {
            return edges[i].weight;
        }
    }
    return -1.0;
}

std::string GraphNetwork::PrintGraph()
{
    std::stringstream buffer;
    buffer.str("");
    buffer << "Graph contains " << node_indices.size() << " nodes and " << edges.size() << " edges." << std::endl << std::endl;
    buffer << "Nodes: " << std::endl;
    for (int i = 0; i < node_indices.size(); i++)
    {
        buffer << node_indices[i] << ", ";
        if (i%10 == 9)
        {
            buffer << std::endl;
        }
    }
    buffer << std::endl << "Edges: " << std::endl;
    for (Edge tmp_edge : edges)
    {
        buffer << tmp_edge.node1 << "," << tmp_edge.node2 << " --> "<< tmp_edge.weight << std::endl;
    }
    buffer << std::endl;
    return buffer.str();
}