#include <iostream>

#include <GLFW/glfw3.h>
#include <boost/heap/priority_queue.hpp>
#include <stack>
#include <algorithm>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh.h"


#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "imgui.h"


using namespace geometrycentral;
using namespace geometrycentral::surface;

//Mesh and Geometry
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;



//Geodesic Properties
std::vector<std::set<std::pair<size_t, double>>> graph;
std::vector<double> distances;
std::vector<int> parent;
std::list<Vector3> path;
int source;
int target;
double dijkstraRuntime;
double pathLength;

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


void createGraph(const std::unique_ptr<ManifoldSurfaceMesh>& mesh, const std::unique_ptr<VertexPositionGeometry>& geometry)
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
		const std::vector<double>& avgs = getVertexAVGs();
		polyscope::getSurfaceMesh("mesh")->addVertexScalarQuantity("AVG", avgs);
	}

}



int main()
{
	std::tie(mesh, geometry) = readManifoldSurfaceMesh("Models/Geodesic/centaur.off");


	polyscope::init();
	auto* psMesh =
		polyscope::registerSurfaceMesh("mesh",
			geometry->inputVertexPositions,
			mesh->getFaceVertexList());


	polyscope::state::userCallback = callback;

	

	polyscope::show(); // pass control to the gui until the user exits


	return 0;
}