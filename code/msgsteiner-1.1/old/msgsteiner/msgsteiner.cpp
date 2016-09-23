/*
   Copyright 2009-2010 Alfredo Braunstein and Riccardo Zecchina

   This file is part of MSGSTEINER (Max Sum for generalized steiner problems on graphs).

   MSGSTEINER is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   MSGSTEINER is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with MSGSTEINER; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */





#include <boost/config.hpp>

#include "mes.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/random.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/connected_components.hpp>

#include <fstream>
#include <vector>
#include <utility>
#include <string>
#include <math.h>
#include <iomanip>
#include <boost/limits.hpp>
#include <queue>


using namespace boost;
using namespace std;

bool tree = false;
int decision = 10;
double beta = 0;
double rein = 0;
int maxit = 2000;
int depth = 10;
double tolerance = 1e-10;
double noise;

double mincost = inf;

mt19937 gen;
mt19937 mes_gen;

uniform_real<double> const uni_dist(0,1);
variate_generator<mt19937 &, uniform_real<> > real01(gen, uni_dist);
variate_generator<mt19937 &, uniform_real<> > mes_real01(mes_gen, uni_dist);

struct EdgeProperties {
       EdgeProperties() : ij(depth), ji(depth) {
               randomize(ij, mes_real01);
               randomize(ji, mes_real01);
       }
       Mes ij;
       Mes ji;
       double c;
};


struct VertexProperties  {
       VertexProperties() : G(0), c(0) {}
       VertexProperties(string const & name) : name(name), type(normal), G(0),c(0){}
       string name;
       enum typeEnum {
               root,
               terminal,
               normal
       };
       typeEnum type;
       bool isroot() const { return type == root; }
       double G;
       double c;
};


typedef adjacency_list<vecS, vecS, undirectedS,
                      VertexProperties, EdgeProperties> GraphBase;


typedef graph_traits<GraphBase>::vertex_iterator vertex_iterator;
typedef graph_traits<GraphBase>::out_edge_iterator edge_iterator;
typedef graph_traits<GraphBase>::edge_iterator graph_edge_iterator;
typedef graph_traits<GraphBase>::edge_descriptor Edge;
typedef graph_traits<GraphBase>::vertex_descriptor Vertex;

struct Graph : public GraphBase {
	Vertex rootid;
};


Graph g;


inline Mes & getMes(Edge e) 
{
	return source(e, g) < target(e, g) ? g[e].ji : g[e].ij;
}


inline Mes & getMesInv(Edge e) 
{
	return source(e, g) > target(e, g) ? g[e].ji : g[e].ij;
}


double updateroot(Vertex i)
{
	edge_iterator eit, eend;
	for (tie(eit, eend) = out_edges(i, g); eit != eend; ++eit) {
		Mes & out = getMesInv(*eit);
		out.B = -inf;
		for (int d = 0; d < depth; ++d) {
			out.A[d] = d ? -inf : 0;
			out.E[d] = 0;
		}
		out.D = 0;
	}
	return 0;
}


void chooseroot()
{
	double maxtot = -inf;
	Vertex maxi = g.rootid;

	edge_iterator eit, eend;
	for (tie(eit, eend) = out_edges(g.rootid, g); eit != eend; ++eit) {
		Vertex i = target(*eit, g);
		Mes & in = getMes(*eit);
		if (in.F[1] - g[i].G >  maxtot) {
			maxtot = in.F[1] - g[i].G;
			maxi = i;
		}
		g[i].G = 0;
	}
	assert(maxi != g.rootid);
	
	cerr << "... done. root: " << g[maxi].name << endl;
	clear_vertex(g.rootid, g);
	remove_vertex(g.rootid, g);
	g.rootid = maxi;

	
	graph_edge_iterator geit, geend;
	for (tie(geit, geend) = edges(g); geit != geend; ++geit) {
		g[*geit].ij = Mes(depth);
		g[*geit].ji = Mes(depth);
		randomize(g[*geit].ij, mes_real01);
		randomize(g[*geit].ji, mes_real01);
	}

	g[maxi].type = VertexProperties::root;
}


double update(Vertex i) 
{
	if (i == g.rootid)
		return updateroot(i);

	Proba sumE(depth, 0),  maxcEA(depth, -2*inf), maxcEA2(depth, -2*inf);
	double sumD = 0;
	vector<int> maxcEAidx(depth, -1);
	edge_iterator eit, eend;
	int idx = 0;
	for (tie(eit, eend) = out_edges(i, g); eit != eend; ++eit, ++idx) {
		double const c = g[*eit].c;
		Mes const & in = getMes(*eit);
		Mes const & out = getMesInv(*eit);
		for (int d = 0; d < depth; ++d) {
			sumE[d] += in.E[d];
			double cEA = d ? -c -in.E[d] + in.A[d - 1] : -inf;
			cEA += rein * out.F[d];
			if (maxcEA[d] <= cEA) {
				maxcEA2[d] = maxcEA[d];
				maxcEA[d] = cEA;
				maxcEAidx[d] = idx;
			} else if (maxcEA2[d] < cEA)
				maxcEA2[d] = cEA;
		}
		sumD += in.D;
	}

	static Mes old(depth);
	double eps = 0;
	idx = 0;
	for (tie(eit, eend) = out_edges(i, g); eit != eend; ++eit, ++idx) {
		Mes const & in = getMes(*eit);
		double c = g[*eit].c;
		
		Mes & out = getMesInv(*eit);
		swap(old, out);
		double maxA = -inf;
		double sumDp = sumD - in.D;
		out.B = -g[i].c + sumDp + rein * g[i].G;
		
		for (int d = 0; d < depth; ++d) {
			double sumEp = sumE[d] - in.E[d];
			double maxcEAp = (maxcEAidx[d] == idx) ? maxcEA2[d] : maxcEA[d];
			out.A[d] = sumEp + maxcEAp;
			maxA = max(maxA, out.A[d]);
		}

		out.D = max(out.B, maxA);
		
		double C = -inf;
		for (int d = depth; d--;) {
			out.E[d] = max(out.D, C);
			double sumEp = sumE[d] - in.E[d];
			C = -c + sumEp + rein * old.F[d];
		//	double cEA = d ? -c -in.E[d] + in.A[d - 1] : -inf;
		//	out.F[d] = sumE[d] + cEA + rein * old.F[d];
			out.F[d] = d ? C + in.A[d - 1] : -inf;
		}

		out.reduce();

		eps = max(eps, l8dist(out, old));
	}
	
	g[i].G = -g[i].c + sumD + rein * g[i].G;

	return eps;
}


double iterate()
{
	vector<int> permutation(num_vertices(g));
	for (unsigned j = 0; j < num_vertices(g); ++j) 
		permutation[j] = j;
	
	random_shuffle(permutation.begin(), permutation.end());

	double eps = 0;
	for (unsigned j = 0; j < num_vertices(g); ++j) 
		eps = max(eps, update(permutation[j]));
	
	return eps;
}


pair<bool, Edge> marginal(Vertex i) 
{
	if (g[i].isroot())
		return make_pair(false, Edge());

	double maxtot = -inf;
	Edge maxtotedge;

	edge_iterator eit, eend;
	for (tie(eit, eend) = out_edges(i, g); eit != eend; ++eit) {
		Mes const & out = getMesInv(*eit);
		for (int d = 0; d < depth; ++d) {
			if (out.F[d] >= maxtot) {
				maxtot = out.F[d];
				maxtotedge = *eit;
			}
		}
	}

	if (maxtot > g[i].G)
		return make_pair(true, maxtotedge);

	return make_pair(false, Edge());
}

struct TreeChecker {
	typedef adjacency_list<vecS, vecS, undirectedS> SimpleGraph;
	TreeChecker() {
		SimpleGraph g2(num_vertices(g));
		unconnected = 0;
		dmax = 0;
		nnodes = 1;
		ecost = 0;
		vcost = 0;
		cost = 0;

		vertex_iterator vit, vend;
		for (tie(vit, vend) = vertices(g); vit != vend; ++vit) if (*vit != g.rootid) {
			bool a;
			Edge e;
			tie(a,e) = marginal(*vit);
			if (a) {
				SimpleGraph::vertex_descriptor v1 = target(e, g);
				SimpleGraph::edge_descriptor e2 = add_edge(v1, *vit, g2).first;
				ecost += g[e].c;
			} else
				vcost += g[*vit].c;
		}
		std::vector<int> distances(num_vertices(g2));
		breadth_first_search(g2,
				    g.rootid,
				    visitor(make_bfs_visitor(record_distances(&distances[0], on_tree_edge()))));
		for (SimpleGraph::vertex_descriptor j = 0; j < num_vertices(g2); ++j) {
			if (out_degree(j, g2) > 0 && j != g.rootid) {
				++nnodes;
				if (distances[j] == 0)
					++unconnected;
				else
					dmax = max(dmax, distances[j]);
			}
		}
		cost = vcost + ecost;
	}

	int unconnected ;
	int dmax;
	int nnodes;

	double ecost;
	double vcost;
	double cost;
};


double converge()
{
	int it = 0;
	double err;

	double cost = inf;
	mincost = inf;
	int coincide = 0;
	double allcoincide = 0;

	do {
		rein = beta * it;
		err = iterate();
		double oldcost = cost;
		TreeChecker c;
		cost = c.cost;
		allcoincide *= .9;
		if (! c.unconnected) {
			if (cost < mincost) {
				mincost = cost;
				cerr << it << " " << cost << " " << c.ecost << " " << c.vcost 
					<< " " << (allcoincide/10) << " " <<  c.dmax << endl;
			}
			if (fabs(oldcost - c.cost) < tolerance) {
				allcoincide++;
				coincide++;
			} else
				coincide = 0;
		} else
			coincide = 0;

	} while (coincide < decision && err > tolerance && ++it < maxit);
	
	return mincost;
}


int idx(string const & id)
{
	static map <string, int> idx_map;
	map<string, int>::iterator mit = idx_map.find(id);
	if (mit == idx_map.end()) 
		return idx_map[id] = add_vertex(VertexProperties(id), g);
	return mit->second;
}


void read_graph(istream & file)
{
	string tok, tok2;

	bool rootset = false;
	while (file >> tok) {
		if (tok == "T") {
			file >> tok2;
			int id = idx(tok2);
			g[id].c = 10000;
		} else if (tok == "W") {
			double w;
			file >> tok2 >> w;
			int id = idx(tok2);
			g[id].c = w;
		} else if (tok == "R") {
			file >> tok2;
			int id = idx(tok2);
			g[id].type = VertexProperties::root;
			g.rootid = id;
			rootset = 1;
		} else if (tok == "E") {
			string i, j;
			double w;
			file >> i >> j >> w;
			Edge e = add_edge(vertex(idx(i), g), vertex(idx(j), g), g).first;
			g[e].c = w  + noise * real01();
		}
	}

	cerr << num_edges(g) << " edges, " << num_vertices(g) << " vertices" << endl;
	if (!rootset) {
		cerr << "Choosing root..." << endl;
		g.rootid = idx("_ROOTMSGSTEINER");
		vertex_iterator vit, vend;
		for (tie(vit, vend) = vertices(g); vit != vend; ++vit) if (g.rootid != *vit) {
			Edge e = add_edge(g.rootid, *vit, g).first;
			g[e].c = inf / num_vertices(g) + noise * real01();
		}
		int olddecision = decision;
		decision = 100;
		converge();
		decision = olddecision;
		chooseroot();
	}
}


void output_tree()
{
	vertex_iterator vit, vend;
	for (tie(vit, vend) = vertices(g); vit != vend; ++vit) {
		Vertex i = *vit;
		bool b; Edge e;
		tie(b, e) = marginal(i);
		if (b) {
			Vertex j = target(e, g);
			cout << g[i].name << " " << g[j].name << " " << g[i].c << " " << g[e].c << endl;
		}
	}
}


namespace po = boost::program_options;

po::variables_map parse_command_line(int ac, char ** av)
{
	po::options_description desc("Usage: " + string(av[0]) + " <option> ... \n\twhere <option> is one or more of");
	desc.add_options()
		("help", "produce help message")
		("depth,d", po::value(&depth)->default_value(10), "set maximum depth")
		("maxit,t", po::value(&maxit)->default_value(100), "set maximum number of iterations")
		("tolerance,e", po::value(&tolerance)->default_value(1e-5), "set convergence tolerance")
		("noise,r", po::value<double>(&noise)->default_value(0), "set random factor")
		("tree,o", "outputs final tree to std output")
		("seed,s", po::value<unsigned>(), "sets instance seed")
		("mseed,z", po::value<unsigned>(), "sets messages seed")
		("messages,M", "output messages on convergence")
		("rein,g", po::value<double>(&beta)->default_value(0), "sets reinforcement parameter rein")
		("decision,y", po::value<int>(&decision)->default_value(10), "program converges after this # repeats of the decision variables");

	po::variables_map vm;
	po::store(po::parse_command_line(ac, av, desc), vm);
	po::notify(vm);

	if (vm.count("seed")) {
		unsigned s = vm["seed"].as<unsigned>();
		gen.seed(s);
	}
	if (vm.count("mseed")) {
		unsigned s = vm["mseed"].as<unsigned>();
		mes_gen.seed(s);
	}
	
	if (vm.count("help")) {
		cerr << desc << "\n";
		exit(1);
	}

	if (vm.count("maxit")) {
		if (maxit == 0)
			tree = true;
	}
	return vm;
}


void output_messages()
{
	vertex_iterator vit, vend;
	for (tie(vit, vend) = vertices(g); vit != vend; ++vit) {
		edge_iterator eit, eend;
		tie(eit, eend) = out_edges(*vit, g);

		Vertex v = *vit;
		cout << "----------------------------------" << endl 
			<< "Vertex " << g[*vit].name << endl << "c: " << g[v].c << endl 
			<< "G: " << g[v].G << endl;
		for (; eit != eend; ++eit) {
			Vertex w = target(*eit, g);
			cerr << g[v].name << " -> " <<  g[w].name << endl;
			Mes m = getMes(*eit);
			cerr << m << endl;
		}
	}

}


int main(int ac, char** av)
{
	cout.setf(ios_base::fixed, ios_base::floatfield);
	po::variables_map vm = parse_command_line(ac, av);
	
	read_graph(cin);
	
	converge();

	if (vm.count("messages")) 
		output_messages();

	if (vm.count("tree"))
		output_tree();
		
	return 0;
}


