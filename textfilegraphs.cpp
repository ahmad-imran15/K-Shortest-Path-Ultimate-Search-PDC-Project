
#include <fstream>
#include <string>
#include <sstream>
#include <queue>
#include <chrono>
#include <omp.h>
using namespace std;

int total_nodes;
string FileName = "Email-Enron.txt";
struct Matrix {
    int row;
    int col;
    int value;
};
class Node{
private:
    int name;
public:
    vector<Node*> edge;

    // Constructor
    Node(int name = -1);

    // member functions

  Node(int name): name(name) 
  {

  }
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

    //Destructor
   ~Node() {

   }
};

class Graph{
private:
    int totalNodes;
    vector<Node*> nodes;
    int** adjacency_Matrix;
    vector<Matrix> sparseMatrixEntries;

public:
    // Constructor
    Graph(int totalNodes = 0);

    // member Functions
    const Node * getNode(int start);

    Graph Graph(int total_no_ofNodes) : total_no_ofNodes(total_no_ofNodes) {
    nodes.resize(total_no_ofNodes, NULL);
    for (int i = 0; i < total_no_ofNodes; i++) {
        this->nodes[i] = new Node(i);
    }
}
}

int getTotalNodes(){ 
    return this->totalNodes;
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
    ~Graph(){
    this->totalNodes = 0;
    for(auto node: this->nodes){
        delete node;
    }
    // Deallocate memory for AdjacencyMatrix
    for(int i=0; i<this->totalNodes; i++){
        delete[] this->AdjacencyMatrix[i];
    }
    delete[] this->AdjacencyMatrix;
}
const Node* getNode(int start) {

    for (const Node* node : this->getNodes()) {
        if (node->getName() == start) {
            return node;
        }
    }

    return nullptr;
}

};

int* findKShortestpath(Graph* graph, int startNode, int endNode, int k) {
    int n = graph->get_Total_noNodes();
    vector<Node*>& nodes = graph->getNodes(); 

    vector<vector<int>> dis(n + 1, vector<int>(k, numeric_limits<int>::max()));

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    pq.push({0, startNode});
    dis[startNode][0] = 0;

    while (!pq.empty()) {
        int u = pq.top().second;
        int d = pq.top().first;
        pq.pop();
        if (dis[u][k - 1] < d)
            continue;

        const vector<Node*>& edges = nodes[u]->getEdges();

        #pragma omp for
        for (const Node* edge : edges) {
            int v = edge->getName();
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
        distances[i] = dis[endNode][i] + 1;
        if (distances[i] <= 0) {
            #pragma omp critical
            {
                distances[i] = 2147483647;
            }
        }
    }
    return distances;
}

void quicksort(int* arr, int left, int right) {
    if (left < right) {
        int pivot = arr[(left + right) / 2];
        int i = left, j = right;
        while (i <= j) {
            while (arr[i] < pivot) i++;
            while (arr[j] > pivot) j--;
            if (i <= j) {
                swap(arr[i], arr[j]);
                i++;
                j--;
            }
        }
        quicksort(arr, left, j);
        quicksort(arr, i, right);
    }
}

int main(int argc, char** argv) {
    ofstream obj("../serialcodeTimes.txt", ios::app);  //obj=timefile
    if (!obj.is_open()) {
        cerr << "Error opening time's file!" << endl;
        return -1;
    }
    else{
        auto start = std::chrono::high_resolution_clock::now();  // Start measuring the time
        if(FileName == "Email-EuAll.txt") // Set up nodes based on file
        {
            total_nodes = 265214;
        }
        else if(FileName == "Email-Enron.txt"){
            total_nodes = 36692;
        }

        // Graph Instance
        Graph* graph = new Graph(total_nodes);
        graph->readFromFile();
        int findNumberOfNeighbours = graph->findNumberOfNeighbours(0);
        const Node* startNode = graph->getNode(0);
        int i =0,k=3;

        int * recieveBuffer = new int [k*findNumberOfNeighbours];

        for (Node * edges : startNode->getEdges())
        { 
            int * distances=findKShortest(graph, edges->getName(), 10, k);

            for (int j=0;j<k;i++,j++)
                recieveBuffer[i]=distances[j];
        }
        quicksort(recieveBuffer, 0, k * findNumberOfNeighbours -1 );

        int * shortestK = new int [k];
        cout<<"\n the k shortest paths are from "<<0<<" to "<<10<<": ";

        for(int i = 0; i<k;i++){
            shortestK[i]= recieveBuffer[i];
            cout<<shortestK[i]<<" ";
        }
        cout<<endl;

        // End measuring time 
        auto end = std::chrono::high_resolution_clock::now();

        // Calculate duration of time
        std::chrono::duration<double> duration = end - start;
        std::cout << "Execution time is: " << duration.count() << " seconds" << std::endl;

       
        obj << duration.count() << endl;   // Adding to times File

        // Clean up
        delete graph;

    }

    obj.close();
    return 0;
}