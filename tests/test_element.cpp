#include <gtest/gtest.h>
#include <test_utils.h>
#include <nef_to_halfspace_tree.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

TEST(HalfspaceTreeGeneration, Cube) {
	typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
	typedef CGAL::Point_3<Kernel> Point;
	typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
	typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
	typedef Nef_polyhedron::Vertex_const_iterator Vertex_const_iterator;

	int i = 51;
	Polyhedron P;
	Nef_polyhedron NP;
	{
		std::ostringstream oss;
		oss << std::setw(5) << std::setfill('0') << i << ".off";
		auto fn = oss.str();
		std::ifstream(fn.c_str()) >> P;
		ASSERT_TRUE(P.is_valid());
		ASSERT_TRUE(P.is_closed());
		CGAL::Polygon_mesh_processing::triangulate_faces(P);
		CGAL::Polygon_mesh_processing::remove_degenerate_faces(P);
		ASSERT_FALSE(CGAL::Polygon_mesh_processing::does_self_intersect(P)) << i << " self intersects" << std::endl;
		NP = Nef_polyhedron(P);
	}

	{
		auto G = build_facet_edge_graph(NP);
		auto T = build_halfspace_tree(G, NP);
		auto NP1 = T->evaluate();

		auto make_vertex_point_it = [](Vertex_const_iterator p) {
			return boost::make_transform_iterator(p, [](auto v) { return v.point(); });
		};

		std::set<Point> s1(make_vertex_point_it(NP.vertices_begin()), make_vertex_point_it(NP.vertices_end()));
		std::set<Point> s2(make_vertex_point_it(NP1.vertices_begin()), make_vertex_point_it(NP1.vertices_end()));

		Polyhedron P;
		convert_to_polyhedron(NP1, P);
		std::ofstream ofs("converted.off");
		ofs << P;

		ASSERT_EQ(s1, s2) << "At least the vertices are the same...";
	}
}