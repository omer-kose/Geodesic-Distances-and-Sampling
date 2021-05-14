#include <iostream>

#include <GLFW/glfw3.h>
#include <boost/heap/priority_queue.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <stack>
#include <algorithm>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh.h"


#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "imgui.h"

#include <iomanip>


using namespace geometrycentral;
using namespace geometrycentral::surface;

//Mesh and Geometry
std::unique_ptr<SurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;



//Geodesic Properties (Might be stored in a class, or struct as a better and safer design)
//For now making these global to easilty access them from IMGUI
std::vector<std::set<std::pair<size_t, double>>> graph;
std::vector<double> distances;
std::vector<int> parent;
std::list<Vector3> path;
std::vector<double> vertexAVGS;
std::vector<double> vertexMGDS;
int source;
int target;
double dijkstraRuntime;
double pathLength;
int nFPSSamples;
//Tau for elimianting samples in Symmetric Invariant Sampling
double tau;


class VertexComparator
{
public:
	//Comparator, that will form min heap wrt eucledian distance
	bool operator() (const std::pair<size_t, double>& lhs, const std::pair<size_t,
		double>& rhs) const
	{
		return lhs.second > rhs.second;
	}
};


void createGraph(const std::unique_ptr<SurfaceMesh>& mesh, const std::unique_ptr<VertexPositionGeometry>& geometry)
{
	graph.resize(mesh->nVertices());
	//Traverse all the faces
	for (const Vertex& v : mesh->vertices())
	{
		for (const Vertex& n : v.adjacentVertices())
		{
			graph[v.getIndex()].insert({ n.getIndex(), (geometry->inputVertexPositions[v] - geometry->inputVertexPositions[n]).norm() });
		}
	}
}

std::vector<double> dijkstraPQ(int source)
{
	//Could be initalized outside
	std::vector<double> dist(graph.size(), std::numeric_limits<double>::max());
	parent.resize(graph.size(), -1);
	boost::heap::priority_queue<std::pair<size_t, double>, boost::heap::compare<VertexComparator>> pq;
	//FOR BENCHMARK (initializations are skipped for benchmark)
	double t1 = glfwGetTime();


	pq.push({ source, 0.0 });
	dist[source] = 0.0;
	parent[source] = source;

	while (!pq.empty())
	{
		std::pair<size_t, double> vertex = pq.top();
		pq.pop();

		//Get the neighbours
		const std::set<std::pair<size_t, double>>& neighbours = graph[vertex.first];
		//Explore the neighbours
		for (const auto& n : neighbours)
		{
			//Relax condition
			if (dist[vertex.first] + n.second < dist[n.first])
			{
				dist[n.first] = dist[vertex.first] + n.second;
				parent[n.first] = vertex.first;
				pq.push({ n.first, dist[n.first] });
			}
		}
	}

	dijkstraRuntime = glfwGetTime() - t1;

	return dist;
}


std::vector<double> dijkstraFH(int source)
{
	//Could be initalized outside
	std::vector<double> dist(graph.size(), std::numeric_limits<double>::max());
	parent.resize(graph.size(), -1);
	std::vector<boost::heap::fibonacci_heap<std::pair<size_t, double>, boost::heap::compare<VertexComparator>>::handle_type> handles(graph.size());
	boost::heap::fibonacci_heap<std::pair<size_t, double>, boost::heap::compare<VertexComparator>> fh;
	//FOR BENCHMARK (initializations are skipped for benchmark)
	double t1 = glfwGetTime();

	//Fill the FH first
	for (size_t i = 0; i < graph.size(); ++i)
	{
		handles[i] = fh.emplace(i, std::numeric_limits<double>::max());
	}

	fh.decrease(handles[source], {source, 0.0});
	dist[source] = 0.0;
	parent[source] = source;

	while (!fh.empty())
	{
		std::pair<size_t, double> vertex = fh.top();
		fh.pop();

		//Get the neighbours
		const std::set<std::pair<size_t, double>>& neighbours = graph[vertex.first];
		//Explore the neighbours
		for (const auto& n : neighbours)
		{
			//Relax condition
			if (dist[vertex.first] + n.second < dist[n.first])
			{
				dist[n.first] = dist[vertex.first] + n.second;
				fh.decrease(handles[n.first], {n.first, dist[n.first]});
			}
			
		}
	}

	dijkstraRuntime = glfwGetTime() - t1;

	return dist;
}


//This is the better optimized version of Dijkstra using Fibo Heap (Still not fast as PQ)
std::vector<double> dijkstraFHOptimized(int source)
{
	//Could be initalized outside
	std::vector<double> dist(graph.size(), std::numeric_limits<double>::max());
	parent.resize(graph.size(), -1);
	std::vector<boost::heap::fibonacci_heap<std::pair<size_t, double>, boost::heap::compare<VertexComparator>>::handle_type> handles(graph.size());
	boost::heap::fibonacci_heap<std::pair<size_t, double>, boost::heap::compare<VertexComparator>> fh;
	//This is needed to determine whether to use decrease key or just to push to the fh
	std::vector<bool> visited(graph.size(), false);
	//FOR BENCHMARK (initializations are skipped for benchmark)
	double t1 = glfwGetTime();

	
	handles[source] = fh.emplace(source, 0.0);
	dist[source] = 0.0;
	parent[source] = source;
	visited[source] = true;

	while (!fh.empty())
	{
		std::pair<size_t, double> vertex = fh.top();
		fh.pop();

		//Get the neighbours
		const std::set<std::pair<size_t, double>>& neighbours = graph[vertex.first];
		//Explore the neighbours
		for (const auto& n : neighbours)
		{
			//Relax condition
			if (dist[vertex.first] + n.second < dist[n.first])
			{
				if (visited[n.first] == false)
				{
					dist[n.first] = dist[vertex.first] + n.second;
					handles[n.first] = fh.emplace(n.first, dist[n.first]);
					visited[n.first] = true;
				}
				else
				{
					dist[n.first] = dist[vertex.first] + n.second;
					fh.decrease(handles[n.first], { n.first, dist[n.first] });
				}

			}
		}
	}

	dijkstraRuntime = glfwGetTime() - t1;

	return dist;
}


void retrievePath()
{
	pathLength = distances[target];
	std::stack<int> st;
	int u = target;
	while (u != source)
	{
		st.push(u);
		u = parent[u];
	}
	st.push(source);

	while (!st.empty())
	{
		path.push_back({ geometry->inputVertexPositions[st.top()] });
		st.pop();
	}
}


std::vector<double> dijkstraArray(int source)
{
	//Could be initalized outside
	std::vector<double> dist(graph.size(), std::numeric_limits<double>::max());
	parent.resize(graph.size(), -1);
	std::vector<std::pair<size_t, double>> vertices;
	//FOR BENCHMARK (initializations are skipped for benchmark)
	double t1 = glfwGetTime();

	
	vertices.emplace_back(source, 0.0);
	dist[source] = 0.0;
	parent[source] = source;

	while (!vertices.empty())
	{
		//Find min element. This is the point where algorithm loses a lot of time
		double minVal = std::numeric_limits<double>::max();
		size_t minIndex;
		for (size_t i = 0; i < vertices.size(); ++i)
		{
			if (vertices[i].second < minVal)
			{
				minVal = vertices[i].second;
				minIndex = i;
			}
		}
		std::pair<size_t, double> vertex = vertices[minIndex];
		//Swap the last element with minIndex to pop_back 
		std::swap(vertices[minIndex], vertices[vertices.size() - 1]);
		vertices.pop_back();

		//Get the neighbours
		const std::set<std::pair<size_t, double>>& neighbours = graph[vertex.first];
		//Explore the neighbours
		for (const auto& n : neighbours)
		{
			//Relax condition
			if (dist[vertex.first] + n.second < dist[n.first])
			{
				dist[n.first] = dist[vertex.first] + n.second;
				parent[n.first] = vertex.first;
				vertices.emplace_back(n.first, dist[n.first]);
			}
		}
	}

	dijkstraRuntime = glfwGetTime() - t1;

	return dist;
}

void createDistanceMatrixFile()
{
	std::ofstream outfile;
	outfile.open("ManGeoMatrix.txt", std::ofstream::out);
	for (size_t i = 0; i < mesh->nVertices(); ++i)
	{
		const std::vector<double> distances = dijkstraPQ(i);
		//Write the row into the file 
		for (size_t j = 0; j < distances.size(); ++j)
		{
			outfile << std::noshowpoint << distances[j] << " "; 
		}
		
		outfile << "\n";
	}

	outfile.close();
}




double AVG(size_t V)
{
	const std::vector<double>& dists = dijkstraPQ(V);
	//Sum up all the distances
	double sum = 0.0;
	for (double distance : dists)
	{
		sum += distance;
	}
	return sum;
}


std::vector<double> getVertexAVGs()
{
	std::vector<double> avgs(mesh->nVertices());
	for (size_t i = 0; i < avgs.size(); ++i)
	{
		avgs[i] = AVG(i);
	}
	return avgs;
}


double MGD(size_t V, const std::vector<size_t>& coarseSet)
{
	const std::vector<double>& dists = dijkstraPQ(V);
	//From the coarse set take the mininmum dist one. This resembles FPS.
	double min = std::numeric_limits<double>::max();
	for (size_t p : coarseSet)
	{
		if (dists[p] < min)
		{
			min = dists[p];
		}
	}

	return min;
}

std::vector<double> getVertexMGDS(const std::vector<size_t>& coarseSet)
{
	//Coarse Set will be evaluated to 0 but since we are only going to take local extremas they do not matter
	std::vector<double> MGDS(mesh->nVertices());

	for (size_t i = 0; i < mesh->nVertices(); ++i)
	{
		MGDS[i] = MGD(i, coarseSet);
	}

	return MGDS;
}

//The starting vertex, which maximizes the AVG is provided from outside to avoid consecutive
//AVG computations in case of consecutive FPS calls with different amount of nSamples.
//This is for testing purposes. I think FPS also should figure outevaluate the starting point itself
std::vector<size_t> FPS(size_t startingPoint, size_t nSamples)
{
	if (nSamples <= 0 || (startingPoint < 0 || startingPoint > mesh->nVertices()))
	{
		//Could have thrown exception but yeah.
		return std::vector<size_t>();
	}

	std::vector<size_t> samples;
	//Start from the first vertex.
	samples.push_back(startingPoint);
	std::vector<double> distances = dijkstraPQ(startingPoint);
	//Find the next pt by finding the current farthest sample
	size_t farthestPoint = std::distance(distances.begin(), std::max_element(distances.begin(), distances.end()));
	samples.push_back(startingPoint);
	//For every chosen farthest point run a dijkstra and update correspondences.
	for (size_t i = 1; i < nSamples; ++i)
	{
		const std::vector<double> geodesicResult = dijkstraPQ(farthestPoint);
		//Compare this result with distances to update correspondenses.
		for (size_t j = 0; j <= distances.size(); ++j)
		{
			if (geodesicResult[j] < distances[j])
			{
				//If the new result is closer to the j'th point, we need to update its value.
				distances[j] = geodesicResult[j];
			}
		}
		//Correspondenses updated. Now choose the new farthest point sample
		farthestPoint = std::distance(distances.begin(), std::max_element(distances.begin(), distances.end()));
		samples.push_back(farthestPoint);
	}
	
	return samples;

}


std::vector<size_t> coarseSet()
{
	std::vector<size_t> samples;
    //Create a copy of AVG and weight these distances with the 1/3 of Area of the 1-ring since SPS method uses that as AVG
    //I think they apply a 2-Hodge Star operator to retrieve the value of the vertex.
	std::vector<double> AVGS = vertexAVGS;
	for (size_t i = 0; i < mesh->nVertices(); ++i)
	{
		//Find area of 1-ring of vertex 
		double totalArea = 0.0;
		for (const Face& face : mesh->vertex(i).adjacentFaces())
		{
			totalArea += geometry->faceArea(face);
		}

		AVGS[i] *= (totalArea / 3.0);
	}

	//For every vertex check if it is local extrema of AVG
	for (size_t i = 0; i < mesh->nVertices(); ++i)
	{
		double avgVal = AVGS[i];
		double maxOfNeighbours = std::numeric_limits<double>::min();
		double minOfNeighbours = std::numeric_limits<double>::max();
		//Check local neighbourhood (1-ring) of vertex
		for (const Vertex& neighbour : mesh->vertex(i).adjacentVertices())
		{
			double nAvgVal = AVGS[neighbour.getIndex()];
			if (nAvgVal > maxOfNeighbours)
			{
				maxOfNeighbours = nAvgVal;
			}

			if (nAvgVal < minOfNeighbours)
			{
				minOfNeighbours = nAvgVal;
			}
		}

		//Take it if it is a critical point
		if ((avgVal <  minOfNeighbours))
		{
			samples.push_back(i);
		}

	}

	return samples;
}

//Need to comptue AVG's and MGD's beforehand, this avoids duplicative AVG and MGD calls for tests. 
std::vector<size_t> SPS()
{
	std::vector<size_t> samples = coarseSet();
	//Augment S1 to S2 using MGDs
	std::vector<size_t> augmentationSet;
	double globalMax = std::numeric_limits<double>::min();
	//Extract local maximas
	for (size_t i = 0; i < mesh->nVertices(); ++i)
	{
		double mgdVal = vertexMGDS[i];
		double maxOfNeighbours = std::numeric_limits<double>::min();

		for (const Vertex& neighbour : mesh->vertex(i).adjacentVertices())
		{
			double nMgdVal = vertexMGDS[neighbour.getIndex()];
			if (nMgdVal > maxOfNeighbours)
			{
				maxOfNeighbours = nMgdVal;
			}
		}

		//Take it if it is a critical point
		if (mgdVal > maxOfNeighbours)
		{
			augmentationSet.push_back(i);
			if (mgdVal > globalMax)
			{
				globalMax = mgdVal;
			}
		}


	}
	
	//Eliminate some samples to increase the uniformity
	for (size_t sample : augmentationSet)
	{
		if (vertexMGDS[sample] > tau * globalMax)
		{
			samples.push_back(sample); //If passes merge it with the original set
		}
	}



	//polyscope::getSurfaceMesh("mesh")->addVertexScalarQuantity("AVG", AVGS);


	return samples;
}



void callback()
{
	//Graph Construction Runtime
	if (ImGui::Button("Create Graph"))
	{
		graph.clear();
		createGraph(mesh, geometry);
	}

	//Dijkstra Callback
	ImGui::InputInt("Source", &source);
	if (ImGui::Button("DijkstraPQ"))
	{
		distances = dijkstraPQ(source);
		//Show the distance distribution (Not necessary but anyway)
		std::vector<double> scalar(mesh->nVertices());
		for (size_t i = 0; i < mesh->nVertices(); ++i)
		{
			scalar[i] = distances[i];
		}

		polyscope::getSurfaceMesh("mesh")->addVertexScalarQuantity("Distances", scalar);
	}

	if (ImGui::Button("DijkstraFH"))
	{
		distances = dijkstraFH(source);
		//Show the distance distribution (Not necessary but anyway)
		std::vector<double> scalar(mesh->nVertices());
		for (size_t i = 0; i < mesh->nVertices(); ++i)
		{
			scalar[i] = distances[i];
		}

		polyscope::getSurfaceMesh("mesh")->addVertexScalarQuantity("Distances", scalar);
	}

	if (ImGui::Button("DijkstraFHOptimized"))
	{
		distances = dijkstraFHOptimized(source);
		//Show the distance distribution (Not necessary but anyway)
		std::vector<double> scalar(mesh->nVertices());
		for (size_t i = 0; i < mesh->nVertices(); ++i)
		{
			scalar[i] = distances[i];
		}

		polyscope::getSurfaceMesh("mesh")->addVertexScalarQuantity("Distances", scalar);
	}

	if (ImGui::Button("DijkstraArray"))
	{
		distances = dijkstraArray(source);
		//Show the distance distribution (Not necessary but anyway)
		std::vector<double> scalar(mesh->nVertices());
		for (size_t i = 0; i < mesh->nVertices(); ++i)
		{
			scalar[i] = distances[i];
		}

		polyscope::getSurfaceMesh("mesh")->addVertexScalarQuantity("Distances", scalar);
	}

	//Show Dijkstra Runtime 
	ImGui::SameLine(200);
	ImGui::Text("Total Time Taken: %f", dijkstraRuntime);

	//Retrieve Path to Target
	ImGui::InputInt("Target", &target);
	if (ImGui::Button("Retrieve Path"))
	{
		path.clear();
		retrievePath();
		std::vector<std::array<size_t, 2>> ind;
		for (size_t i = 0; i < path.size() - 1; ++i)
		{
			ind.push_back({ i, i + 1 });
		}
		polyscope::getSurfaceMesh("mesh")->addSurfaceGraphQuantity("Path", path, ind);
	}
	//Show PathLength
	ImGui::SameLine(200);
	ImGui::Text("Path Length is: %f", pathLength);
	//Visualize AVGS
	if (ImGui::Button("Retrieve AVGS"))
	{
		vertexAVGS = getVertexAVGs();
		polyscope::getSurfaceMesh("mesh")->addVertexScalarQuantity("AVG", vertexAVGS);
	}
	if (ImGui::Button("Retrieve MGDS"))
	{
		vertexMGDS = getVertexMGDS(coarseSet());
		polyscope::getSurfaceMesh("mesh")->addVertexScalarQuantity("MGD", vertexMGDS);
	}
	//Show FPS results (First need to run AVG function)
	ImGui::InputInt("nFPSSamples", &nFPSSamples);
	if (ImGui::Button("FPS"))
	{
		size_t farthestPoint = std::distance(vertexAVGS.begin(), std::max_element(vertexAVGS.begin(), vertexAVGS.end()));
		const std::vector<size_t> samples = FPS(farthestPoint, nFPSSamples);
		std::vector<Vector3> samplePointPositions(samples.size());
		std::vector<std::array<size_t, 2>> ind;
		for (size_t i = 0; i < samples.size(); ++i)
		{
			samplePointPositions[i] = geometry->inputVertexPositions[samples[i]];
		}
		polyscope::getSurfaceMesh("mesh")->addSurfaceGraphQuantity("FPS Samples", samplePointPositions, ind);
	}
	//Symmetric Invariant Sampling
	ImGui::InputDouble("Tau", &tau);
	if (ImGui::Button("SPS"))
	{
		const std::vector<size_t>& samples = SPS();
		std::vector<Vector3> samplePointPositions(samples.size());
		std::vector<std::array<size_t, 2>> ind;
		for (size_t i = 0; i < samples.size(); ++i)
		{
			samplePointPositions[i] = geometry->inputVertexPositions[samples[i]];
		}
		polyscope::getSurfaceMesh("mesh")->addSurfaceGraphQuantity("SPS Samples", samplePointPositions, ind);

	}


}


int main()
{
	std::tie(mesh, geometry) = readSurfaceMesh("Models/Geodesic/fprintf/man0.off");




	polyscope::init();
	auto* psMesh =
		polyscope::registerSurfaceMesh("mesh",
			geometry->inputVertexPositions,
			mesh->getFaceVertexList());


	polyscope::state::userCallback = callback;


	polyscope::show(); // pass control to the gui until the user exits


	return 0;
}