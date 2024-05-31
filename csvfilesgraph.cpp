#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <chrono>

using namespace std;

class Vertex{
private:
    string identifier;
    //Here i initialzed a identifier for the vertex
    vector<Vertex*> links;
    vector<int> distances;
    //also the  vector to store distances to neighboring vertices

public:
    //here is the  constructor
    Vertex(string identifier = "");

    //below are the member functions
    string getIdentifier() const;
    void addLink(Vertex* vertex, int distance);
    const vector<Vertex*>& getLinks() const;
    const vector<int>& getDistances() const;
    //above are the function to get the identifier , neighboring vertex of the vertex , and distances to neighboring vertices 

    //this below is teh destructor
    ~Vertex();
};

class Network{
private:
    vector<int> indices;
    int totalVertices;
    // this is the total number of vertices in the network
    vector<Vertex*> vertices;
    //this is the vector to veector to store pointers to vertices
    int** addj_of_matri;
    //this is the adjacency matrix representation of the network

public:
    //this below is the constructor
    Network(int totalVertices = 0);
   ~Network();
    //and these are the member Functions
    int getTotalVertices(); // Function to get the total number of vertices
    int parseFile(); // Function to read data from a file and create the network

    void displayConnections(); 
    //this function is to display connections between vertices
    int initializeaddj_of_matri();
    void printAdjMatrix();
     //abov function is for print the adjacency matrix
    const vector<Vertex*>& getVertices() const;
    int calculateDistance(const Vertex* vertex1, const Vertex* vertex2); 
    //below function to get the identifier of a vertex given its index
    string vertexIdentifier(int index, const vector<Vertex*>& vertices);
    


    //and for every class this is the destructor
  //  ~Network();
};

//this codee is for selecting the file
string name_of_file = "classic-who.csv";
int TOTAL_VERTICES = 377;
//giving the total number of vertices in the network

//this is the constructor
Vertex::Vertex(string identifier): identifier(identifier) {}

//and same this is the destructor
Vertex::~Vertex() {}

//it will return the identifier of the vertex
string Vertex::getIdentifier() const {
    return this->identifier;
}


void Vertex::addLink(Vertex* vertex, int distance){
    this->links.push_back(vertex);
    //this push back will add the neighboring vertex to the vector
    this->distances.push_back(distance);
    //and for thsi push back for distance will add the distance rto the neifghbour
}

const vector<Vertex*>& Vertex::getLinks() const {
    return this->links;
    //it will return the neighboring vertices

}

const vector<int>& Vertex::getDistances() const {
    return this->distances;
    //her eit will return the distances to neighboring vertices
}


Network::Network(int totalVertices): totalVertices(totalVertices) {} //this is the constructor 

Network::~Network(){
    this->totalVertices = 0;
     //it will reset the total number of vertices
    for(auto vertex: this->vertices){
        delete vertex; 
        // it will dlete each vertex
    }
    // Deallocate memory for addj_of_matri
    for(int i=0; i<this->totalVertices; i++){
        delete[] this->addj_of_matri[i]; 
        //this abive delete will delete the each row of the adjacency matrix
    }
    delete[] this->addj_of_matri;
     //similarly this delete will delete the adjacency matrix
}

int Network::getTotalVertices(){ 
    return this->totalVertices;
    //this return will return the total number of vertices
}

const vector<Vertex*>& Network::getVertices() const {
    return this->vertices;
    //ye wala return lkakrde ga  vertices of the network
}
int Network::parseFile(){
    ifstream file("classic-who.csv"); // Open the file
    if (!file.is_open()) {
        cerr << "Error..... Opening the file!" << endl; // Print error message if file cannot be opened
        return -1; // Return -1 to indicate failure
    }
    string currentLine;
    
    // Skipping initial line
    getline(file, currentLine); // Skip the first line of the file
    
    while(getline(file, currentLine)){
        int distance;
        string vertex1, vertex2;

        //here Creating stringstream instance to extract meaningful information effectively
        stringstream dataStream(currentLine);
        getline(dataStream, vertex1, ','); 
        getline(dataStream, vertex2, ','); 
        // Extract the first and the second vertex
        dataStream >> distance; // Extract the distance
        
        // Check if vertex1 already exists, if not, create and store it
        int x1 = -1;
        for (int i = 0; i < vertices.size(); ++i) {
            if (vertices[i]->getIdentifier() == vertex1) {
                x1 = i;
                break;
            }
        }
        if (x1 == -1) {
            x1 = vertices.size();
            vertices.push_back(new Vertex(vertex1)); // Create and store the vertex
        }
        
        // Check if vertex2 already exists, if not, it will create and than store it
        int x2 = -1;
        for (int i = 0; i < vertices.size(); ++i) {
            if (vertices[i]->getIdentifier() == vertex2) {
                x2 = i;
                break;
            }
        }
        if (x2 == -1) {
            x2 = vertices.size();
            vertices.push_back(new Vertex(vertex2)); 
            // here creating and store the vertex
        }
        
        // Add link between vertices
        vertices[x1]->addLink(vertices[x2], distance);
        //here adding a link from vertex1 to vertex2
        
        
        vertices[x2]->addLink(vertices[x1], distance);
         //like vice versa it will add a link from vertex2 to vertex1
    }

    file.close(); // Close the file
    return 0; // Return 0 to indicate success
}


void Network::displayConnections(){
    int count = 0;
    for(const Vertex* vertex : vertices) {
        cout << "<<Vertex "<< ++count <<">>: " << vertex->getIdentifier() << " <<connection>>: ";
        for (const Vertex* neighbor : vertex->getLinks()) {
            cout << neighbor->getIdentifier() << ", "; 
            //here printing the identifiers of neighboring vertices
        }
        cout <<"\n\n";
    }
}


int Network::initializeaddj_of_matri(){
    //here i am  creating the space for teh matrix
    this->addj_of_matri = new int*[this->totalVertices];
    for(int i = 0; i < this->totalVertices; ++i){
        this->addj_of_matri[i] = new int[this->totalVertices]{0}; 
        //also here i am  initializing  the adjacency matrix and giving values as  zeros
    }

    // this below for loop will iterate over vertices to find indices
    for(int i = 0; i < this->totalVertices; ++i) {
        for(const Vertex* neighbor : this->vertices[i]->getLinks()) {
            //and will find the index of the neighbor vertex
            int index = -1;
            for(int j = 0; j < this->totalVertices; ++j) {
                if(vertices[j]->getIdentifier() == neighbor->getIdentifier()) {
                    index = j;
                    break;
                }
            }
            if(index != -1) {
                // it will set adjacency matrix value to distance if link exists
                this->addj_of_matri[i][index] = this->calculateDistance(vertices[i], vertices[index]);
            }
        }
    }

    return 0; 
    
}
//function to print the adjoint of patrix
void Network::printAdjMatrix(){
    cout << "The addjoint of of matrix is: " << endl;
    for(int r=0; r<this->totalVertices; r++){
        for(int c=0; c<this->totalVertices; c++){
            cout << this->addj_of_matri[r][c] << " "; 
            //it will print each element of the adjacency matrix
        }
        cout << endl;
    }
    cout << endl;
}

// Function to get the identifier of the vertex given its index
string Network::vertexIdentifier(int index, const vector<Vertex*>& vertices) {
    if (index >= 0 && index < vertices.size()) {
        return vertices[index]->getIdentifier(); 
        //it will return the identifier of the vertex at the given index
    }
    //here below it will return an empty string if index out of rangee hoya
    return "";
}

//rhis is the function to get the index of the vertex given its identifier
int vertexIndex(const string& vertexIdentifier, const vector<Vertex*>& vertices) {
    for (int i = 0; i < vertices.size(); ++i) {
        if (vertices[i]->getIdentifier() == vertexIdentifier) {
            return i; 
            //it will return the index of the vertex with the given identifier
        }
    }
    //if not found that ye return -1 arde ga
    return -1;
}


int Network::calculateDistance(const Vertex* vertex1, const Vertex* vertex2) {
    //here it will find the index of vertex2 in the links of vertex1
    for(size_t i = 0; i < vertex1->getLinks().size(); ++i) {
        if(vertex1->getLinks()[i]->getIdentifier() == vertex2->getIdentifier()) {
            // and return the distance associated with the link
            return vertex1->getDistances()[i];
        }
    }
    // If link not found, return a default distance or throw an exception
    return -1;
     //this return will  will adjust the value according to your needs
}

int* findKShortestPaths(Network* network, const string& startVertexIdentifier, const string& endVertexIdentifier, int k) {
    
    //it will retrieve the vertices from the network
    const vector<Vertex*>& vertices = network->getVertices();

    //here in this berlwo part of code it wil  find the start and end vertices
    int start = vertexIndex(startVertexIdentifier, vertices);
    int end = vertexIndex(endVertexIdentifier, vertices);

    //here start and aend  verticles are intiliaed as start and end

    if (start == -1 || end == -1) {
        cerr << "Sorry....The Start or end vertex not found!" << endl;
         //it will dispaly the  error message if start or end vertex not found
        return nullptr;
         //and than it will return the  nullptr 
         //that shows ke kuch nhi mila
    }

    int n = vertices.size();

    // yahan ,ene 2d vector intiliaze kiya
    vector<vector<int>> distances(n, vector<int>(k, numeric_limits<int>::max()));

    //setting the  Priority queue to store vertices based on distance
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    pq.push(make_pair(0, start));
    distances[start][0] = 0;

    while (!pq.empty()) {
        int u = pq.top().second;
        int d = pq.top().first;
        pq.pop();

        if (distances[u][k - 1] < d) {
            continue;
        }

        const vector<Vertex*>& neighbors = vertices[u]->getLinks();
        
        for (const Vertex* neighbor : neighbors) {
            int v = vertexIndex(neighbor->getIdentifier(), vertices);
            int weight = network->calculateDistance(vertices[u], neighbor);
              
            if (v == -1) {
                continue; 
                //it will skip if the vertex index is not found
            }
            
            if (d + weight < distances[v][k - 1]) {
                distances[v][k - 1] = d + weight;
                sort(distances[v].begin(), distances[v].end());
                pq.push(make_pair(d + weight, v));
            }
        }
    }

    // it will copy the k shortest distances into an array
    int* result = new int[k];
    for (int i = 0; i < k; ++i) {
        result[i] = distances[end][i];
    }

    return result; 
    // in end will return the array of k shortest distances
}

//Main function
int main() {
    ofstream outputFile("serialTimes.txt", ios::app); 
//her ioepning the output file 

    //checking the if file will not open that it will print the error message if file cannot be opened
    if (!outputFile.is_open()) {
        cerr << "Error.... Opening the times file!" << endl; 
        return -1; 
        // and the an rreturn -1 to indicate failure
    }
    else{
        // Start calculating the time for the code woerking
        auto start = std::chrono::high_resolution_clock::now();

        // iyt will craetye the  network and read ythe all data from file
        Network* network = new Network(TOTAL_VERTICES);
        network->parseFile();

        // Find k shortest paths
        string startVertexIdentifier = "A. H. Millington"; 
        //giving teh  start vertex
        string endVertexIdentifier = "Ace"; 
        //and similary the end vertix end vertex
        int k = 3; 
        // Number of shortest paths to find

        int* shortestPaths = findKShortestPaths(network, startVertexIdentifier, endVertexIdentifier, k);

        // Print the k shortest paths
        cout << "The K Shortest Paths from the " << startVertexIdentifier << " to the " << endVertexIdentifier << "is:" << endl;
        for (int i = 0; i < k; ++i) {
            cout << shortestPaths[i] << " "; 
            // it will print each shortest path
        }
        cout << endl;

        // End measuring time
        auto end = std::chrono::high_resolution_clock::now();

        // Calculate duration
        std::chrono::duration<double> duration = end - start;
        std::cout << "The Execution time is: " << duration.count() << " seconds" << std::endl;

        // Adding to times File
        outputFile << duration.count() << endl; 
        //it will write the execution time to the output file

        delete[] shortestPaths;
        delete network; 
        // Deallocate memory 
    }
    outputFile.close(); 
    // Closing the output file
    return 0;
}