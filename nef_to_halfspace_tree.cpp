#ifdef NDEBUG
#undef NDEBUG
#endif

#define CGAL_EIGEN3_ENABLED

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_nary_union_3.h>
#include <CGAL/Nef_nary_intersection_3.h>
#include <CGAL/aff_transformation_tags.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Extended_cartesian.h>
#include <CGAL/convex_decomposition_3.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_d.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/graph/copy.hpp>

#include <fstream>
#include <iostream>
#include <cmath>
#include <array>
#include <algorithm>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Extended_cartesian<CGAL::Gmpq> Kernel_b;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel_b> Nef_polyhedron_b;
typedef Kernel::Point_3 Point;
typedef std::array<Point, 4> quad;

typedef Nef_polyhedron::Vertex_const_handle Vertex_const_handle;
typedef Nef_polyhedron::Halfedge_const_handle Halfedge_const_handle;
typedef Nef_polyhedron::Halffacet_const_handle Halffacet_const_handle;
typedef Nef_polyhedron::SHalfedge_const_handle SHalfedge_const_handle;
typedef Nef_polyhedron::SHalfloop_const_handle SHalfloop_const_handle;
typedef Nef_polyhedron::Halffacet_cycle_const_iterator Halffacet_cycle_const_iterator;
typedef Nef_polyhedron::SHalfedge_around_facet_const_circulator SHalfedge_around_facet_const_circulator;
typedef Nef_polyhedron::SFace_const_handle SFace_const_handle;
typedef Nef_polyhedron::Volume_const_iterator Volume_const_iterator;
typedef Nef_polyhedron::Shell_entry_const_iterator Shell_entry_const_iterator;
typedef Nef_polyhedron::SM_const_decorator SM_const_decorator;
typedef Nef_polyhedron::SFace_cycle_const_iterator SFace_cycle_const_iterator;
typedef Nef_polyhedron::SHalfedge_around_sface_const_circulator SHalfedge_around_sface_const_circulator;
typedef Nef_polyhedron::SVertex_const_handle SVertex_const_handle;
typedef Nef_polyhedron::Volume_const_handle Volume_const_handle;

typedef CGAL::Epick_d<CGAL::Dimension_tag<4>> KdKernel;
typedef KdKernel::Point_d Point_d;
typedef CGAL::Search_traits_d<KdKernel> TreeTraits;
typedef CGAL::Kd_tree<TreeTraits> Tree;
typedef CGAL::Fuzzy_sphere<TreeTraits> Fuzzy_sphere;
typedef std::map<Kernel::Plane_3, Kernel::Plane_3, PlaneLess> plane_map_t;

std::string dump_facet(Halffacet_const_handle h) {
	std::ostringstream oss;

	const auto& p = h->plane();
	oss << "F plane=" << p << std::endl;
	h->facet_cycles_begin();

	auto fc = h->facet_cycles_begin();
	auto se = SHalfedge_const_handle(fc);
	CGAL_assertion(se != 0);
	SHalfedge_around_facet_const_circulator hc_start(se);
	SHalfedge_around_facet_const_circulator hc_end(hc_start);
	CGAL_For_all(hc_start, hc_end) {
		CGAL_NEF_TRACEN("   add vertex " << hc_start->source()->center_vertex()->point());
		oss << "  co=" << hc_start->source()->center_vertex()->point() << std::endl;
	}
	
	oss << std::endl;

	return oss.str();
}

#if defined(HALFSPACE_EXTENDED)
#define halfspace_tree_result Nef_polyhedron_b 
#else
#define halfspace_tree_result Nef_polyhedron
#endif

struct PlaneLess {
	bool operator()(const Kernel::Plane_3& lhs, const Kernel::Plane_3& rhs) const {
		auto lhs_a = lhs.a();
		auto lhs_b = lhs.b();
		auto lhs_c = lhs.c();
		auto lhs_d = lhs.d();
		auto rhs_a = rhs.a();
		auto rhs_b = rhs.b();
		auto rhs_c = rhs.c();
		auto rhs_d = rhs.d();
		return std::tie(lhs_a, lhs_b, lhs_c, lhs_d) < std::tie(rhs_a, rhs_b, rhs_c, rhs_d);
	}
};

// Define edge types
enum EdgeType { CONCAVE, CONVEX };

struct VertexProperties {
	Halffacet_const_handle facet;
	size_t original_index;
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
	VertexProperties,
	boost::property<boost::edge_weight_t, EdgeType>> Graph;

typedef boost::filtered_graph<Graph, boost::keep_all, std::function<bool(Graph::vertex_descriptor)>> Filtered;

struct Intersection_visitor {
	typedef void result_type;
	void operator()(const Kernel::Point_3& p) const
	{
		std::cout << p << std::endl;
	}
	void operator()(const Kernel::Segment_3& s) const
	{
		std::cout << s << std::endl;
	}
	void operator()(const Kernel::Triangle_3& s) const
	{
		std::cout << s << std::endl;
	}
};

struct PlaneHash {
	size_t operator()(const CGAL::Plane_3<Kernel>& plane) const
	{
		// @todo why can I only get this to work on double?
		std::hash<double> h;
		std::size_t result = h(CGAL::to_double(plane.a()));
		boost::hash_combine(result, h(CGAL::to_double(plane.b())));
		boost::hash_combine(result, h(CGAL::to_double(plane.c())));
		boost::hash_combine(result, h(CGAL::to_double(plane.d())));
		return result;
	}
};

// Can be used to convert polyhedron from exact to inexact and vice-versa
template <class Polyhedron_input,
	class Polyhedron_output>
	struct Copy_polyhedron_to
	: public CGAL::Modifier_base<typename Polyhedron_output::HalfedgeDS>
{
	// CGAL::Cartesian_converter<typename Polyhedron_input::Traits, typename Polyhedron_output::Traits> converter;

	Copy_polyhedron_to(const Polyhedron_input& in_poly)
		: in_poly(in_poly) {}

	void operator()(typename Polyhedron_output::HalfedgeDS& out_hds)
	{
		typedef typename Polyhedron_output::HalfedgeDS Output_HDS;
		typedef typename Polyhedron_input::HalfedgeDS Input_HDS;

		CGAL::Polyhedron_incremental_builder_3<Output_HDS>
			builder(out_hds);

		typedef typename Polyhedron_input::Vertex_const_iterator
			Vertex_const_iterator;
		typedef typename Polyhedron_input::Facet_const_iterator
			Facet_const_iterator;
		typedef typename
			Polyhedron_input::Halfedge_around_facet_const_circulator
			HFCC;

		builder.begin_surface(in_poly.size_of_vertices(),
			in_poly.size_of_facets(),
			in_poly.size_of_halfedges());

		for (Vertex_const_iterator
			vi = in_poly.vertices_begin(), end =
			in_poly.vertices_end();
			vi != end; ++vi) {

			CGAL_assertion(vi->point().x().degree() == 0);
			CGAL_assertion(vi->point().y().degree() == 0);
			CGAL_assertion(vi->point().z().degree() == 0);

			auto x = CGAL::to_double(vi->point().x()[0]);
			auto y = CGAL::to_double(vi->point().y()[0]);
			auto z = CGAL::to_double(vi->point().z()[0]);
			
			Point xyz(x, y, z);

			builder.add_vertex(xyz);
		}

		typedef CGAL::Inverse_index<Vertex_const_iterator> Index;
		Index index(in_poly.vertices_begin(),
			in_poly.vertices_end());

		for (Facet_const_iterator
			fi = in_poly.facets_begin(), end =
			in_poly.facets_end();
			fi != end; ++fi) {
			HFCC hc = fi->facet_begin();
			HFCC hc_end = hc;
			builder.begin_facet();
			do {

				builder.add_vertex_to_facet(index[hc->vertex()]);
				++hc;
			} while (hc != hc_end);
			builder.end_facet();
		}
		builder.end_surface();
	} // end operator()(..)
private:
	const Polyhedron_input& in_poly;
}; // end Copy_polyhedron_to<>

template <class Poly_B, class Poly_A>
void poly_copy(Poly_B& poly_b, const Poly_A& poly_a)
{
	poly_b.clear();
	Copy_polyhedron_to<Poly_A, Poly_B> modifier(poly_a);
	poly_b.delegate(modifier);
}

int main(int, char**) {

	/*
	Nef_polyhedron a(CGAL::Plane_3<Kernel>(0, 0, 1, 0));
	Nef_polyhedron b(CGAL::Plane_3<Kernel>(0, 1, 0, 0));
	auto result = a + b;
	*/


	Polyhedron cube1, cube2, cube3, cube4, cube5;
	Nef_polyhedron nef_cube1, nef_cube2, nef_cube3, nef_cube4, nef_cube5;

	// Create two cubes as Polyhedra
	createCube(cube1, 1.0);
	createCube(cube2, 1.0);
	createCube(cube3, 0.5);
	createCube(cube4, 1.0);
	createCube(cube5, 0.5);

	scalePolyhedron(cube3, 1, 1, 3);
	scalePolyhedron(cube5, 1, 1, 3);

	translatePolyhedron(cube5, Kernel::Vector_3(0.5, 0.5, 0));


	convertToNefPolyhedron(cube1, nef_cube1);
	convertToNefPolyhedron(cube2, nef_cube2);
	convertToNefPolyhedron(cube3, nef_cube3);
	convertToNefPolyhedron(cube4, nef_cube4);
	convertToNefPolyhedron(cube5, nef_cube5);

	auto double_convex = nef_cube1 - nef_cube3 - nef_cube5;

	{
		Polyhedron result_poly;
		double_convex.convert_to_polyhedron(result_poly);
		std::ofstream fs("double_convex.off");
		fs << result_poly;
	}

	{
		Graph G = build_facet_edge_graph(double_convex);
		auto tree = build_halfspace_tree(G, double_convex);
		auto converted = tree->evaluate();
		std::cout << "same = " << (double_convex == converted);

		Polyhedron result_poly;
		converted.convert_to_polyhedron(result_poly);
		std::ofstream fs("converted_double_convex.off");
		fs << result_poly;
	}

	/*

	translatePolyhedron(cube2, Kernel::Vector_3(0, 0, 1.0));
	translatePolyhedron(cube2, Kernel::Vector_3(0.001, 0.001, 0));
	translatePolyhedron(cube3, Kernel::Vector_3(0, 0, 1.0));
	translatePolyhedron(cube4, Kernel::Vector_3(0, 0, 1.0));

	rotatePolyhedron(cube4, 0.01, Kernel::Vector_3(0, 0, 1.0));

	// Convert Polyhedra to Nef_polyhedron
	convertToNefPolyhedron(cube1, nef_cube1);
	convertToNefPolyhedron(cube2, nef_cube2);
	convertToNefPolyhedron(cube3, nef_cube3);
	convertToNefPolyhedron(cube4, nef_cube4);



	std::list<Nef_polyhedron*> operands = { &nef_cube1, &nef_cube2, &nef_cube4 };

	{
		// standard

		CGAL::Nef_nary_union_3<Nef_polyhedron> builder;
		for (auto& p : operands) {
			builder.add_polyhedron(*p);
		}
		auto result = builder.get_union();
		Polyhedron result_poly;
		result.convert_to_polyhedron(result_poly);
		std::ofstream fs("standard.off");
		fs << result_poly;
		std::cout << "standard " << result_poly.size_of_vertices() << " " << result_poly.size_of_facets() << std::endl;
	}

	{
		// snapped halfspaces

		std::list<std::unique_ptr<halfspace_tree>> trees;
		std::list<Kernel::Plane_3> planes;
		std::vector<Point_d> planes_as_point;
		planes_as_point.reserve(planes.size());

		for (auto& p : operands) {
			Graph G = build_facet_edge_graph(*p);
			trees.emplace_back(std::move(build_halfspace_tree(G, *p)));
			trees.back()->accumulate(planes);
			auto test = trees.back()->evaluate();
			std::cout << "tree result is " << (test == *p ? std::string("equal") : std::string("not equal")) << std::endl;
		}	
		plane_map_t plane_map;

		for (auto& p : planes) {
			Point_d p(CGAL::to_double(p.a()), CGAL::to_double(p.b()), CGAL::to_double(p.c()), CGAL::to_double(p.d()));
			planes_as_point.push_back(p);
		}

		Tree kdtree(planes_as_point.begin(), planes_as_point.end());
		
		auto plit = planes.begin();
		for (size_t i = 0; i < planes.size(); ++i) {
			auto& query = planes_as_point[i];
			Fuzzy_sphere fs(query, 0.01, 0.);
			std::cout << "q " << query << std::endl;

			std::list<Point_d> results;
			kdtree.search(std::back_inserter(results), fs);
			for (auto& r : results) {
				std::cout << " " << r << std::endl;
			}
			auto sum = std::accumulate(++results.begin(), results.end(), results.front(), [](Point_d a, Point_d b) {return Point_d(a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]); });
			int N = results.size();
			results.clear();

			Point_d n(-query[0], -query[1], -query[2], -query[3]);
			Fuzzy_sphere fsn(n, 0.01, 0.);
			kdtree.search(std::back_inserter(results), fsn);
			for (auto& r : results) {
				std::cout << " " << r << std::endl;
			}
			N += results.size();
			auto sum2 = std::accumulate(results.begin(), results.end(), sum, [](Point_d a, Point_d b) {return Point_d(a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3]); });

			auto avg = Kernel::Plane_3(sum[0] / N, sum[1] / N, sum[2] / N, sum[3] / N);

			std::cout << *plit << " -> " << avg << std::endl;
			plane_map.insert({ *plit++, avg });
		}

		decltype(trees) trees_snapped;

		for (auto& x : trees) {
			trees_snapped.push_back(std::move(x->map(plane_map)));
		}

		CGAL::Nef_nary_union_3<Nef_polyhedron> builder;
		for (auto& p : trees_snapped) {
			builder.add_polyhedron(p->evaluate());
		}
		auto result = builder.get_union();
		Polyhedron result_poly;
		result.convert_to_polyhedron(result_poly);
		std::ofstream fs("snapped.off");
		fs << result_poly;
		std::cout << "snapped " << result_poly.size_of_vertices() << " " << result_poly.size_of_facets() << std::endl;
	}
	*/


	/*
	auto result = nef_cube1 + nef_cube2;
	auto subtracted = nef_cube1 - nef_cube3;

	{
		CGAL::Polyhedron_3<Kernel> x;
		subtracted.convert_to_polyhedron(x);
		std::ofstream fs("debugg.off");
		fs << x;
	}

	CGAL::Polyhedron_3<Kernel> presult;
	Graph G = build_facet_edge_graph(subtracted);
	auto tree = build_halfspace_tree(G, subtracted);
	std::list<Point_d> points;
	tree->accumulate(points);

	Tree kdtree(points.begin(), points.end());

	std::map<Point_d, Point_d> to_cluster_average;

	for (auto& query : points) {
		Fuzzy_sphere fs(query, 0.01, 0.);
		std::cout << "q " << query << std::endl;
		
		std::list<Point_d> results;
		kdtree.search(std::back_inserter(results), fs);
		for (auto& r : results) {
			std::cout << " " << r << std::endl;
		}
		auto sum = std::accumulate(++results.begin(), results.end(), results.front(), [](Point_d a, Point_d b) {return Point_d(a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]); });
		results.clear();
		int N = results.size();

		Point_d n(-query[0], -query[1], -query[2], -query[3]);
		Fuzzy_sphere fsn(n, 0.01, 0.);
		kdtree.search(std::back_inserter(results), fsn);
		for (auto& r : results) {
			std::cout << " " << r << std::endl;
		}
		N += results.size();
		auto sum2 = std::accumulate(results.begin(), results.end(), sum, [](Point_d a, Point_d b) {return Point_d(a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3]); });

		
		auto avg = Point_d(sum[0] / N, sum[1] / N, sum[2] / N, sum[3] / N);
		to_cluster_average.insert({ query, avg });
	}

	auto evaluated = tree->evaluate();
	std::cout << (evaluated == subtracted) << std::endl;
	evaluated.convert_to_polyhedron(presult);
	*/


	/*
	CGAL::Polyhedron_3<Kernel_b> presult_extended;
	build_halfspace_tree(nef_cube1)->evaluate().convert_to_polyhedron(presult_extended);

	Copy_polyhedron_to<CGAL::Polyhedron_3<Kernel_b>, CGAL::Polyhedron_3<Kernel>> builder(presult_extended);
	presult.delegate(builder);
	*/

	return 0;
}
