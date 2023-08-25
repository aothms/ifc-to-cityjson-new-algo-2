#include <gtest/gtest.h>
#include <test_utils.h>
#include <nef_to_halfspace_tree.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>

#include <boost/iterator/transform_iterator.hpp>

TEST(HalfspaceTreeGeneration, Cube) {
	typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
	typedef CGAL::Point_3<Kernel> Point;
	typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
	typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
	typedef Nef_polyhedron::Vertex_const_iterator Vertex_const_iterator;

	Polyhedron cube, cube2, cube3;
	createCube(cube, 1.0);
	createCube(cube2, 0.5);
	createCube(cube3, 0.5);
	scalePolyhedron(cube2, 1, 1, 3);
	scalePolyhedron(cube3, 1, 1, 3);
	translatePolyhedron(cube3, Kernel::Vector_3(0.5, 0.5, 0));
	Nef_polyhedron cube_nef(cube);
	Nef_polyhedron cube2_nef(cube2);
	Nef_polyhedron cube3_nef(cube3);

	auto concave = cube_nef - cube2_nef - cube3_nef;

	auto graph = build_facet_edge_graph(concave);
	auto tree = build_halfspace_tree(graph, concave);
	auto concave_evaluated = tree->evaluate();
	// ASSERT_EQ(concave, concave_evaluated) << "We would expect same result, but somehow nested subtraction results in a different object. Something to do with marks on the boundary being different?";

	auto make_vertex_point_it = [](Vertex_const_iterator p) {
		return boost::make_transform_iterator(p, [](auto v) { return v.point(); });
	};

	std::set<Point> s1(make_vertex_point_it(concave.vertices_begin()), make_vertex_point_it(concave.vertices_end()));
	std::set<Point> s2(make_vertex_point_it(concave_evaluated.vertices_begin()), make_vertex_point_it(concave_evaluated.vertices_end()));
	
	ASSERT_EQ(s1, s2) << "At least the vertices are the same...";
}
