//This program written on 4/10/2015 by Ram Rohan Dantu
//The copyrights to this program belong to the course Advanced Computer Networks CE6390

#include <iostream>
#include <iomanip>
#include <vector>
//#include <algorithm>
using namespace std;


//Functions
void printgraph(vector<double> & graph, size_t siz);
void printList(vector<double> & stNodes, size_t sizStNodes);
void makeStNodesGraph(vector<double> & stNodesGraph, size_t sizStNodes);
void printTree(vector<double> & tree, size_t siz);
int findIndex(vector<double> & arr, int start, int end, double value);
void bubbleSort(vector<double> & edges, size_t siz, int option);
int findEdgeIndex(vector<double> & arr, int start, int end, double edge, double node1, double node2);
void deleteEdge(vector<double> & edges, size_t & siz, int index);
void minSpanningTree(vector<double> & graph, vector<double> & MST, size_t sizG, size_t & sizMST, size_t numEdges);
void dijkstra(vector<double> & graph, vector<double> & SPT, vector<double> & MST, size_t sizG, size_t sizSPT, size_t sizMST, int start, int end);




int main(void){
	
	cout << "\n\t\tProgram to calculate the steiner tree by Ram Rohan Dantu" << endl;
	cout << "\n\t\tWritten for the course CE6390 Spring 2015" << endl;
	cout << "\n\n\t\t\"Freedom is the goal of all nature\" - Vivekananda in San Francisco,1900" << endl;
	
	//Declarations
	
	//Adjacency matrix
	size_t size = 10;
	vector<double> graph(size*size); // adjacency matrix of graph
	
	//Set of steiner nodes
	size_t sizeStNodes = 5;
	vector<double> stNodes(sizeStNodes);

	//Complete graph of steiner nodes
	vector<double> stNodesGraph(sizeStNodes*sizeStNodes);

	//list of edge valus and their corresponding end nodes
	size_t numEdges = 12;
	size_t sizeEdges = numEdges*3+1;
	vector<double> edges(sizeEdges);
	
	//declare a location to store the MST
	size_t sizeMST_T1 = 1;
	vector<double> MST_T1(sizeMST_T1);
	

	//declare a location to store Gs
	size_t sizeSPT = size;
	vector<double> SPT(sizeSPT*sizeSPT);

	//Code
	for(int i=0;i<size;i++){
		for(int j=0;j<size;j++)
			if(j==0){
				graph[size*i+j] = static_cast<double>(i);
			}else{
				if(i==0){
					graph[size*i+j] = static_cast<double>(j);
				}else{
					graph[size*i +j] = 0;
				}
			}
	}
	//printgraph(graph,size);
	
	cout << "\nSTART" << endl;
	cout << "\nPopulating Adjacency matrix of graph with nodes and Edges" << endl;
	//populate graph with the problem to the steiner tree
	//connections of node1
	//graph[10] = 1;
	graph[12] = 10;
	graph[19] = 1;

	//connections of node2
	graph[21] = 10;
	graph[23] = 8;
	graph[26] = 1;

	//connections of node3
	graph[32] = 8;
	graph[35] = 2;
	graph[34] = 9;

	//connections of node4
	graph[43] = 9;
	graph[45] = 2;
	
	//connections of node5
	graph[53] = 2;
	graph[54] = 2;
	graph[56] = 1;
	graph[59] = 1;

	//connections of node6
	graph[62] = 1;
	graph[65] = 1;
	graph[67] = 1;

	//connections of node7
	graph[76] = 1;
	graph[78] = 0.5;

	//connections of node8
	graph[87] = 0.5;
	graph[89] = 0.5;

	//connections of node9
	graph[91] = 1;
	graph[95] = 1;
	graph[98] = 0.5;

	
	//Populate the list of steiner nodes
	stNodes[1] = 1;
	stNodes[2] = 2;
	stNodes[3] = 3;
	stNodes[4] = 4;
	
	cout << "\nAdjacency matrix of graph \'G\'" << endl;
	printgraph(graph,size);
	
	cout << "\nList of Steiner nodes \'S\' in graph" << endl;
	printList(stNodes,sizeStNodes);
	
	cout <<"\nStep1: make a complete graph \'G1\' of steiner nodes \'S\' with equal edge weights" << endl;
	makeStNodesGraph(stNodesGraph,sizeStNodes);
	printgraph(stNodesGraph,sizeStNodes);
	
	
	//add the edges in the graph to the list of edges
	int count =1;
	for(int i=1;i<size;i++){
		for(int j=i;j<size;j++){
			if(graph[size*i+j]>0){
				edges[count] = graph[size*i+j];
				edges[count+1] = i;
				edges[count+2] = j;
				count+=3;
			}	
		}
	}
	
	
	numEdges = 6;
	cout <<"\nStep2: Find the Minimal Spanning Tree \'T1\' of \'G1\'\n"<< endl;
	minSpanningTree(stNodesGraph, MST_T1, sizeStNodes, sizeMST_T1, numEdges);
	printTree(MST_T1,sizeMST_T1);
	
	cout <<"\nStep3: Construct the subgraph \'Gs\' of \'G\' by replacing each edge in \'T1\' by its orresponding shortest path in \'G\' \n"<<endl;
	dijkstra(graph, SPT, MST_T1, size, sizeSPT, sizeMST_T1, 1,1);
	printgraph(SPT,sizeSPT);
	
	cout<<"\nStep4,5: Find the Steiner tree \'Th\' from \'Gs\'\n"<<endl;
	minSpanningTree(SPT, MST_T1, sizeSPT, sizeMST_T1, numEdges);	
	printTree(MST_T1,sizeMST_T1);
		
	cout<<"\n\n-----------SUCCESS---------------"<<endl;

	

	return 0;
}


void printgraph(vector<double> & graph, size_t siz){

	for(int i=0;i<siz;i++){
		for(int j=0;j<siz;j++){
			cout << setw(8) << graph[siz*i + j];
		}
		cout << endl;
	}
	cout << endl;
}


void printList(vector<double> & stNodes, size_t sizStNodes){

	for(int i=0;i<sizStNodes;i++){
		cout << setw(8) << stNodes[i];
	}
	cout << endl;
}

void makeStNodesGraph(vector<double> & stNodesGraph, size_t sizStNodes){

	for(int i=0;i<sizStNodes;i++){
		for(int j=0;j<sizStNodes;j++){
			if(i==0){
				stNodesGraph[sizStNodes*i+j] = j;
			}else{
				if(j==0){
					stNodesGraph[sizStNodes*i+j] = i;
				}else{
					if(i!=j){
						stNodesGraph[sizStNodes*i+j] = sizStNodes-1;
					}
				}
			}
		}
	}
}

int findIndex(vector<double> & arr, int start, int end, double value){
	for(int i=start;i<end;i++){
		//cout << i << endl;
		if(arr[i]==value){
			//cout << "A";
			return i;
		}else{
			if(i==end-1 && arr[i]!=value){
				//cout << "B";
				return -1;
			}
		}
	}
}


int findEdgeIndex(vector<double> & arr, int start, int end, double edge, double node1, double node2){
	for(int i=start;i<end;i++){
		//cout << i << endl;
		if(arr[i]==edge && arr[i+1]==node1 && arr[i+2]==node2){
			//cout << "A";
			return i;
		}else{
			if(i==end-3 && arr[i]!=edge){
				//cout << "B";
				return -1;
			}
		}
	}
}




void minSpanningTree(vector<double> & graph, vector<double> & MST, size_t sizG, size_t & sizMST, size_t numEdges){
	//Prim's algorithm
	
	//declaring temporary graph for storage
	vector<double> tempGraph(sizG*sizG);
	vector<double> closedSet(sizG);
	
	//vector of edges
	size_t sizEdges = numEdges*3+1;
	vector<double> edges(sizEdges);
	
		//resizeMST 
	sizMST = sizG;
	MST.resize(sizMST*sizMST);

	for(int i=0;i<sizG;i++){
		for(int j=0;j<sizG;j++){
			tempGraph[sizG*i +j] = graph[sizG*i +j];
			if(i==0 || j==0){
				MST[sizMST*i+j] = graph[sizG*i +j];
			}else{
				MST[sizMST*i+j] = 0;
			}
		}
	}

	for(int i=1;i<sizG;i++){
		closedSet[i] = i;
	}
	
	
	//add the edges in the graph to the list of edges
	int count =1;
	for(int i=1;i<sizG;i++){
		for(int j=i;j<sizG;j++){
			if(tempGraph[sizG*i+j]>0){
				edges[count] = tempGraph[sizG*i+j];
				edges[count+1] = i;
				edges[count+2] = j;
				count+=3;
			}	
		}
	}


	bubbleSort(edges,sizEdges,0);
	while(edges.size() > 1){
		if(closedSet[edges[2]] != closedSet[edges[3]]){ // find, enter condition if loop is not formed
			MST[sizMST*edges[2] + edges[3]] = tempGraph[sizG*edges[2] + edges[3]]; //add edge to MST
			
			//make the union
			int p = closedSet[edges[2]];
			int q = closedSet[edges[3]];
			for(int i=1;i<sizG;i++){
				if(closedSet[i] == q){
					closedSet[i] = p;
				}
			}
			deleteEdge(edges, sizEdges, 1); //delete 1 edge from sorted list
		}else{ // condition if loop is formed
			deleteEdge(edges, sizEdges, 1); // delete the edge from the list of edges as it forms a loop
		}	
	}
	
}


void printTree(vector<double> & tree, size_t siz){

	for(int i=0;i<siz;i++){
		for(int j=0;j<siz;j++){
			if(i != 0){
				if(j==0){
					cout << "\n" << setw(4) << i ;
				}else{
					if(tree[siz*i+j]>0){
						cout << setw(4) <<"->" << setw(4) << j;
					}
					
				}
			}
		}
	}
}


void bubbleSort(vector<double> & edges, size_t siz, int option){
	
	double temp;
	double tempa;
	double tempb;
	for(int k=1;k<siz-3;k+=3){
		for(int i=1;i<siz-3;i+=3){
			if(edges[i] > edges[i+3] && option == 0){
				temp = edges[i+3];
				tempa = edges[i+4];
				tempb = edges[i+5];
				edges[i+3] = edges[i];
				edges[i+4] = edges[i+1];
				edges[i+5] = edges[i+2];
				edges[i] = temp;
				edges[i+1] = tempa;
				edges[i+2] = tempb;
			}
			if(edges[i] < edges[i+3] && option == 1){
				temp = edges[i+3];
				tempa = edges[i+4];
				tempb = edges[i+5];
				edges[i+3] = edges[i];
				edges[i+4] = edges[i+1];
				edges[i+5] = edges[i+2];
				edges[i] = temp;
				edges[i+1] = tempa;
				edges[i+2] = tempb;
			}
		}
	}	
}

void deleteEdge(vector<double> & edges, size_t & siz, int index){
	
	edges[1] = 0;
	edges[2] = 0;
	edges[3] = 0;
	if(edges.size()>4){
		bubbleSort(edges,siz,1);
		for(int i=0;i<index;i++){
			siz -= 3;
		}
		edges.resize(siz);
		bubbleSort(edges,siz,0);
	}else{
		siz -=3;
		edges.resize(siz);
	}

}


	

//Algorithm to find the shortest path between two nodes
void dijkstra(vector<double> & graph, vector<double> & SPT, vector<double> & MST, size_t sizG, size_t sizSPT, size_t sizMST, int start, int end){
	
	//set holdidng the list of nodes added to the shortest path tree set
	vector<double> sptSET(sizG);

	//vector holding the distances of all the nodes
	vector<double> distances(sizG);

	//vector holding the parents of all the nodes
	vector<double> parents(sizG);

	//vector holding the spanning tree
	vector<double> minSpanTree(sizMST);

	
	//initialize min spanning tree vector
	for(int i=1;i<sizMST;i++){
		for(int j=1;j<sizMST;j++){
			if(MST[sizMST*i+j]>0){
				minSpanTree[j] = i;
			}
		}
	}
	
	
	//Initialize Gs
	for(int i=0;i<sizSPT;i++){
		for(int j=0;j<sizSPT;j++){
			if(i==0){
				SPT[sizSPT*i+j] = j;
			}
			if(j==0){
				SPT[sizSPT*i+j] = i;
			}
		}
	}

	//Iterate through minSpanTree vector to replace corresponding edges with shortest paths in MST
	for(int k=1;k<sizMST;k++){
		if(minSpanTree[k]>0){
			start = minSpanTree[k];
			end = k;
			
			//Initialize the shortest path tree set
			//initialize the distances of all the nodes
			//initialize the distance of the starting node to 0
			//initialize all parents to -1
			for(int i=1;i<sizG;i++){
				sptSET[i] = -1;
				distances[i] = 1000;
				parents[i] = -1;
			}
			distances[start] = 0;
			sptSET[start] = 0;
			int current = start;
			
			//DIJKSTRA's ALGORITHM
			//while loop iterates till end is added to the shortest path tree
			while(sptSET[end] == -1){
				for(int j=1;j<sizG;j++){//update the distances of adjacent nodes
					if(graph[sizG*current+j]>0){//edge is present in graph
						if(sptSET[j] == -1){//vertex is not added to shortest path tree set
							if(distances[current] + graph[sizG*current+j] < distances[j]){
								distances[j] = distances[current] + graph[sizG*current+j];
								parents[j] = current;
							}
						}
					}
				}

				//Add the vertex with min distance to the shortest path tree set
				int min = 1;
				double min_distance = 1001;
				
				for(int j=1;j<sizG;j++){
					if(sptSET[j] == -1 && distances[j] < min_distance){
						min = j;
						min_distance = distances[j];
					}
				}
				if(sptSET[min] == -1){
					current = min;
					sptSET[current] = 0;

				}

			}

						
			//generate Gs
			int s = start;
			int m = end;
			while(m != s){
				SPT[sizSPT*m + parents[m]] = graph[sizG*m + parents[m]];
				SPT[sizSPT*parents[m] + m] = graph[sizG*parents[m] + m];
				m = parents[m];
			}



		}

	}

}

