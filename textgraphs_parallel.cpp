
#include <algorithm>
#include <limits>
#include <stack>
#include <iostream>
#include <vector>
#include <fstream>
#include <mpi.h>
#include <omp.h>
#include <time.h>
#include <string>
#include <sstream>
#include <queue>
using namespace std;
string file = "Email-Enron.txt";
int total_nodes;
struct Matrix {  // basically it is a sparse matrix
    int row;
    int col;
    int value;
};

class Node {
private:
    int name;

public:
    vector<Node*> edge;

    // Constructor
    Node(int name = -1);
    const vector<Node*>& getEdges() const;
int get_name() const
{
    return this->name;
}

void addEdge(Node* node) 
{
    this->edge.push_back(node);
}
const vector<Node*>& getEdges() const 
{
    return this->edge;
}

    // Destructor
    ~Node(){

    }
};

class Graph 
{
private:
    int total_no_ofNodes;
    vector<Node*> nodes;
    int** adjacency_matrix;
    vector<Matrix> sparse_Matrix_Entries;

public:
    
    Graph(int total_no_ofNodes = 0);

    // member Functions
    vector<Node*>& get_no_Nodes();
    int find_number_of_neighbours(int start);
    const Node* getNode(int start);
int get_Total_noNodes() 
{
    return this->total_no_ofNodes;
}

int read_from_theFile()
 {
    ifstream obj("../" + file);
    if (!obj.is_open()) {
        cerr << "Error opening the file!" << endl;
        return -1;
    }
    string line;
    // Skipping initial lines
    for (int i = 0; i < 3; ++i) {
        getline(obj, line);
    }

    int nodes, edges;
    string messg;
    stringstream ss(line);
    ss >> messg >> messg >> nodes >> messg >> edges;

    // Update the graph
    this->total_no_ofNodes = nodes;

    getline(myFile, line);

    while (getline(myFile, line)) {
        int node1, node2;
        stringstream ss2(line);
        ss2 >> node1 >> node2;

        this->nodes[node1]->addEdge(this->nodes[node2]);
    }

    myFile.close();
    return 0;
}

void print_edges() 
{
    for (const Node* node : nodes) 
    {
        cout << "Node " << node->get_name() << " edges: ";
        for (const Node* edge : node->getEdges()) {
            cout << edge->get_name() << " ";
        }
        cout << endl;
    }
}

int initialize_AdjacencyMatrix() {

    this->adjacency_matrix = new int*[this->total_no_ofNodes];

    for (int i = 0; i < this->total_no_ofNodes; i++) {
        this->adjacency_matrix[i] = new int[this->total_no_ofNodes]{0};
    }

    for (const Node* node : nodes) {
        for (const Node* edge : node->getEdges()) {
            this->adjacency_matrix[node->get_name()][edge->get_name()] = 1;
        }
    }

    return 0;
}

int** get_adjacency_List() {
    return this->adjacency_matrix;
}

void delete_adjacency_List() {
    for (int i = 0; i < total_no_ofNodes; ++i) {
        delete[] this->adjacency_matrix[i];
    }
    delete[] this->adjacency_matrix;
}

vector<Node*>& get_no_Nodes() {
    return nodes;
}

Graph(int total_no_ofNodes) : total_no_ofNodes(total_no_ofNodes) {
    nodes.resize(total_no_ofNodes, NULL);
    for (int i = 0; i < total_no_ofNodes; i++) {
        this->nodes[i] = new Node(i);
    }
}
int find_number_of_neighbours(int start) {
    for (const Node* node : this->get_no_Nodes()) {
        if (node->get_name() == start) {
            return node->getEdges().size();
        }
    }
    return -1;      // Return -1 if the node with the given name is not found

}
    // Destructor
    ~Graph() {
    this->total_no_ofNodes = 0;
    for (auto node : this->nodes) {
        delete node;
    }
    for (int i = 0; i < this->total_no_ofNodes; i++)    // Deallocate memory for AdjacencyMatrix

     {
        delete[] this->adjacency_matrix[i];
    }
    delete[] this->adjacency_matrix;
}
};

Node::Node(int name)
{
   this->name = name;
}
int* findKShortestpath(Graph* graph, int start_node, int end_node, int k) {
    int n = graph->get_Total_noNodes();
    vector<Node*>& nodes = graph->get_no_Nodes(); 

    vector<vector<int>> dis(n + 1, vector<int>(k, numeric_limits<int>::max()));

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    pq.push({0, start_node});
    dis[start_node][0] = 0;

    while (!pq.empty()) {
        int u = pq.top().second;
        int d = pq.top().first;
        pq.pop();
        if (dis[u][k - 1] < d)
            continue;

        const vector<Node*>& edges = nodes[u]->getEdges();

        #pragma omp for
        for (const Node* edge : edges) {
            int v = edge->get_name();
            int cost = 1; // Assuming all edges have a cost of 1

            if (d + cost < dis[v][k - 1]) {
                #pragma omp critical
                {
                    dis[v][k - 1] = d + cost;
                    sort(dis[v].begin(), dis[v].end());
                    pq.push({d + cost, v});
                }
            }
        }
    }

    int* distances = new int[k];
    #pragma omp for
    for (int i = 0; i < k; i++) {
        distances[i] = dis[end_node][i] + 1;
        if (distances[i] <= 0) {
            #pragma omp critical
            {
                distances[i] = 2147483647;
            }
        }
    }
    return distances;
}

const Node* Graph::get_Node(int start) {
    for (const Node* node : this->get_no_Nodes()) {
        if (node->get_name() == start) {
            return node;
        }
    }
    return nullptr;
}

void quick_sort(int* arr, int left, int right) {
    if (left < right) {
        int pivot = arr[(left + right) / 2];
        int i = left, j = right;
        while (i <= j) {
            while (arr[i] < pivot)
                i++;
            while (arr[j] > pivot)
                j--;
            if (i <= j) {
                swap(arr[i], arr[j]);
                i++;
                j--;
            }
        }
        quick_sort(arr, left, j);
        quick_sort(arr, i, right);
    }
}

//main function
int main(int argc, char** argv) {
    ofstream obj("../parallelcodeTimes.txt", ios::app);
    if (!obj.is_open()) {
        cerr << "Error opening times file!" << endl;
        return -1;
    } else {
        // Setting up nodes based on file
        if (file == "Email-EuAll.txt") {
            total_nodes = 265214;
        } else if (file == "Email-Enron.txt") {
            total_nodes = 36692;
        }

        MPI_Init(&argc, &argv);

        int rank, size, k = 3;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        srand(time(0)); 
        int start_node_id = 0;

        srand(time(0)); 
        double startTime, endTime, runtime;

        startTime = MPI_Wtime(); 
        int endNodeId = 10;

        Graph* graph = new Graph(total_nodes);
        graph->readFromFile();
        int find_number_of_neighbours = graph->find_number_of_neighbours(start_node_id);
        int* recv_shortest_distances = new int[k * find_number_of_neighbours];

#pragma omp parallel for
        for (int i = 0; i < k; i++) {
            recv_shortest_distances[i] = 999999999;
        }

        if (size > find_number_of_neighbours) {
            if (rank == 0) {
                cout << "More processes than neighbors.\nYou should have " << find_number_of_neighbours << " processes for ideal computation.\nTerminating." << endl;
            }
            MPI_Finalize();
            return 0;
        }

        const Node* start_node = graph->getNode(start_node_id);
        int* k_shortest_Distance = findKShortestpath(graph, start_node->edge[rank]->get_name(), endNodeId, k);

        MPI_Gather(k_shortest_Distance, k, MPI_INT, recv_shortest_distances, k, MPI_INT, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            quick_sort(recv_shortest_distances, 0, k * size - 1);

            int* shortestK = new int[k];
            cout << "\nthe k shortest paths are from " << start_node_id << " to " << endNodeId << ": ";

#pragma omp parallel for
            for (int i = 0; i < k; i++) {
                shortestK[i] = recv_shortest_distances[i];
                cout << shortestK[i] << " ";
            }
            cout << endl;

            endTime = MPI_Wtime(); // Record the end time
            runtime = endTime - startTime; // Calculate the runtime

            cout << "Runtime on process : " << runtime << " seconds" << endl;

            // Adding to times File
            timeFile << runtime << endl;
        }
    }

    timeFile.close();
    MPI_Finalize();
    return 0;
}
