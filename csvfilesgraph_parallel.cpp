#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <mpi.h>
#include <omp.h>
#include <ctime>


using namespace std;


//it is a Class representing a Vertex in the Network
class Vertex {
private:
    string VertexName; 
    // Name of the Vertex
    vector<Vertex*> edges; 
    //it is edges connecting to other Vertexs
    vector<int> weights; 
    // this is weights corresponding to the edges

public:
    Vertex(string name = ""); 
    // Constructor
    ~Vertex(); 
    // Destructor
    string getName() const; 
    // Getter for Vertex name
    void addEdge(Vertex* node, int weight); 
    // Method to add an edge to the Vertex
    const vector<Vertex*>& getConnectedVertexs() const; 
    const vector<int>& getEdgeWeights() const; 
    //this is Getter for edge weights
};

// Class representing the Network of Vertexs
class Network {
private:
    vector<int> indices; 
   
    int totalVertexs; 
    // Total number of Vertexs
    vector<Vertex*> Vertexs; 
    int** AdjMatrix; 
    //this is  Adjacency matrix for the Network

public:
    Network(int totalVertexs = 0); // Constructor
    ~Network(); // Destructor
    int countTotalVertexs(); // Method to count the total Vertexs
    int retrieveNetworkData(); 
    void displayConnections(); 
    int setupAdjacencyMatrix(); 
    // Method to setup the adjacency matrix
    void exhibitAdjacencyMatrix(); 
    // Method to exhibit the adjacency matrix
    const vector<Vertex*>& getVertexs() const; 
    //this is basically Getter for Vertexs
    int getWeightBetweenVertexs(const Vertex* Vertex1, const Vertex* Vertex2); 
    // Method to get weight between two Vertexs
};

//this is teh Constructor for Vertex class
Vertex::Vertex(string name): VertexName(name) {}

//this is the  Destructor for Vertex class
Vertex::~Vertex() {}

// Getter for Vertex name
string Vertex::getName() const {
    return this->VertexName;
}

//here is  method to add an edge to the Vertex
void Vertex::addEdge(Vertex* node, int weight){
    this->edges.push_back(node);
    this->weights.push_back(weight);
}

// Getter for connected Vertexs
const vector<Vertex*>& Vertex::getConnectedVertexs() const {
    return this->edges;
}

// Getter for edge weights
const vector<int>& Vertex::getEdgeWeights() const {
    return this->weights;
}

// thi isc constructor for Network class
Network::Network(int totalVertexs): totalVertexs(totalVertexs) {}

// thi is the destructor for Network class
Network::~Network(){
    //hrere deleto=ing allocated memory for Vertexs
    for(auto Vertex: this->Vertexs){
        delete Vertex;
    }
    // using loop to delete allocated memory for adjacency matrix
    for(int i=0; i<this->totalVertexs; i++){
        delete[] this->AdjMatrix[i];
    }
    delete[] this->AdjMatrix;
    //heer edelleting the memory
}

// method to count total Vertexs
int Network::countTotalVertexs(){ 
    return this->totalVertexs;
}

// this is the getter for Vertexs
const vector<Vertex*>& Network::getVertexs() const {
    return this->Vertexs;
}



// method to get weight between two Vertexs
int Network::getWeightBetweenVertexs(const Vertex* Vertex1, const Vertex* Vertex2) {
    for(size_t i = 0; i < Vertex1->getConnectedVertexs().size(); ++i) {
        if(Vertex1->getConnectedVertexs()[i]->getName() == Vertex2->getName()) {
            return Vertex1->getEdgeWeights()[i];
        }
    }
    return -1;
}

// thi sis teh method to retrieve data from a file and populate the Network
int Network::retrieveNetworkData(){
    // Open data file
    ifstream dataFile("classic-who.csv");
    // Check if file is open or ot
    if (!dataFile.is_open()) {
        cerr << "Error opening file!" << endl;
        return -1;
    }
    string currentLine;
    getline(dataFile, currentLine); // Skipping initial line
    
    //her reading the data line by line
    while(getline(dataFile, currentLine)){
        int weight;
        string Vertex1, Vertex2;
        stringstream ss2(currentLine);
        getline(ss2, Vertex1, ',');
        getline(ss2, Vertex2, ',');
        ss2 >> weight;
        
        //here finding or create Vertexs and add edges
        int index1 = -1;
        for (int i = 0; i < Vertexs.size(); ++i) {
            if (Vertexs[i]->getName() == Vertex1) {
                index1 = i;
                break;
            }
        }
        if (index1 == -1) {
            index1 = Vertexs.size();
            Vertexs.push_back(new Vertex(Vertex1));
        }
        
        int index2 = -1;
        for (int i = 0; i < Vertexs.size(); ++i) {
            if (Vertexs[i]->getName() == Vertex2) {
                index2 = i;
                break;
            }
        }
        if (index2 == -1) {
            index2 = Vertexs.size();
            Vertexs.push_back(new Vertex(Vertex2));
        }
        //here x1 and x3 are indexes 1 and indes2
        
        Vertexs[index1]->addEdge(Vertexs[index2], weight);
        Vertexs[index2]->addEdge(Vertexs[index1], weight);
    }

    dataFile.close();
    return 0;
}

// thi si the method to display connections between Vertexs
void Network::displayConnections() {
    int count = 0;
    for(const Vertex* vertex : Vertexs) { // Corrected variable name to Vertexs
        cout << "<<Vertex "<< ++count <<">>: " << vertex->getName() << " <<connections>>: ";
        for (const Vertex* edge : vertex->getConnectedVertexs()) { // Corrected function name to getConnectedVertexs()
            cout << edge->getName() << ", ";
        }
        cout << endl;
    }
}


// Method to setup the adjacency matrix
int Network::setupAdjacencyMatrix(){
    // Allocate memory for the adjacency matrix
    this->AdjMatrix = new int*[this->totalVertexs];
    for(int i = 0; i < this->totalVertexs; ++i){
        this->AdjMatrix[i] = new int[this->totalVertexs]{0};
    }

    // Populate the adjacency matrix
    for(int i = 0; i < this->totalVertexs; ++i) {
        // Iterate over connected vertices
        for(const Vertex* edge : this->Vertexs[i]->getConnectedVertexs()) {
            int index = -1;
            // Find the index of the connected vertex
            for(int j = 0; j < this->totalVertexs; ++j) {
                if(Vertexs[j]->getName() == edge->getName()) {
                    index = j;
                    break;
                }
            }
            // If the connected vertex is found, set the weight in the adjacency matrix
            if(index != -1) {
                this->AdjMatrix[i][index] = this->getWeightBetweenVertexs(Vertexs[i], Vertexs[index]);
            }
        }
    }

    return 0;
}

// Method to exhibit the adjacency matrix
void Network::exhibitAdjacencyMatrix(){
    cout << "AdjacencyMatrix: " << endl;
    for(int r=0; r<this->totalVertexs; r++){
        for(int c=0; c<this->totalVertexs; c++){
            cout << this->AdjMatrix[r][c] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

// Function to get the index of the Vertex given its name
int VertexIndex(const string& VertexName, const vector<Vertex*>& Vertexs) {
    for (int i = 0; i < Vertexs.size(); ++i) {
        if (Vertexs[i]->getName() == VertexName) {
            return i;
        }
    }
    // return -1 if Vertex name is not found
    return -1;
}

// Function to get the name of the Vertex given its index
string VertexName(int index, const vector<Vertex*>& Vertexs) {
    if (index >= 0 && index < Vertexs.size()) {
        return Vertexs[index]->getName();
    }
    // return an empty string if index is out of range
    return "";
}

// Function to find shortest paths in parallel
int* findShortestPathsParallel(Network* Network, const string& startVertexName, const string& endVertexName, int k) {
    const vector<Vertex*>& Vertexs = Network->getVertexs();
    int startIndex = -1, endIndex = -1;

    // Find the indices of the start and end Vertexs in parallel
    #pragma omp parallel for
    for (int i = 0; i < Vertexs.size(); ++i) {
        if (Vertexs[i]->getName() == startVertexName) {
            startIndex = i;
        }
        if (Vertexs[i]->getName() == endVertexName) {
            endIndex = i;
        }
    }

    // Check if start or end Vertex not found
    if (startIndex == -1 || endIndex == -1) {
        cerr << "Start or end Vertex not found!" << endl;
        return nullptr;
    }

    // Initialize variables
    int n = Vertexs.size();
    vector<vector<int>> distances(n, vector<int>(k, numeric_limits<int>::max()));
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    pq.push(make_pair(0, startIndex));
    distances[startIndex][0] = 0;

    // Dijkstra's algorithm to find k shortest paths
    while (!pq.empty()) {
        int u = pq.top().second;
        int d = pq.top().first;
        pq.pop();

        // If the current distance is greater than the k-th shortest path distance, skip
        if (distances[u][k - 1] < d) {
            continue;
        }

        const vector<Vertex*>& edges = Vertexs[u]->getConnectedVertexs();
        
        // Relaxation step - parallel loop to update distances
        #pragma omp parallel for
        for (size_t i = 0; i < edges.size(); ++i) {
            int v = -1;
            // Find the index of the connected Vertex
            for (size_t j = 0; j < Vertexs.size(); ++j) {
                if (Vertexs[j]->getName() == edges[i]->getName()) {
                    v = j;
                    break;
                }
            }
            // If the connected Vertex is not found, skip
            if (v == -1) {
                continue;
            }
            // Get the weight of the edge between u and v
            int weight = Network->getWeightBetweenVertexs(Vertexs[u], edges[i]);
            
            // Update distances if a shorter path is found
            if (d + weight < distances[v][k - 1]) {
                // Use a critical section to update distances safely
                #pragma omp critical
                {
                    distances[v][k - 1] = d + weight;
                    sort(distances[v].begin(), distances[v].end()); // Sort the distances
                    pq.push(make_pair(d + weight, v)); // Push the new distance and Vertex to priority queue
                }
            }
        }
    }

    // Copy the k shortest paths to the result array
    int* result = new int[k];
    for (int i = 0; i < k; ++i) {
        result[i] = distances[endIndex][i];
    }

    return result;
}

// Main function
int main(int argc, char** argv) {
    // Open times file for writing
    ofstream timeFile("parallelTimes.txt", ios::app);
    // Check if times file is opened successfully
    if (!timeFile.is_open()) {
        cerr << "Error opening times file!" << endl;
        return -1;
    }
    else {
        // Initialize MPI
        MPI_Init(&argc, &argv);

        int rank, size, k = 3;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        srand(time(0));
        double startTime, endTime, runtime;

        startTime = MPI_Wtime();
        
        // Create a Network object and retrieve data
      	Network* network = new Network();
    	network->retrieveNetworkData();

        string startVertexName = "A. H. Millington";
        string endVertexName = "Adric";

        // Find k shortest paths in parallel
        int* shortestPaths = findShortestPathsParallel(network, startVertexName, endVertexName, k); // Corrected variable name to network

        // Synchronize all processes
        MPI_Barrier(MPI_COMM_WORLD);

        // Display results and calculate runtime in rank 0
        if (rank == 0) {
            cout << "K Shortest Paths from " << startVertexName << " to " << endVertexName << ":" << endl;
            for (int i = 0; i < k; ++i) {
                cout << shortestPaths[i] << " ";
            }
            cout << endl;

            endTime = MPI_Wtime();
            runtime = endTime - startTime;

            cout << "Runtime on process : " << runtime << " seconds" << endl;

            // Write runtime to times file
            timeFile << runtime << endl;
        }
    }

    // Close times file and finalize MPI
    timeFile.close();
    MPI_Finalize();
    return 0;
}