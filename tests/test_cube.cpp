#include <gtest/gtest.h>
#include <test_utils.h>
#include <nef_to_halfspace_tree.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>

TEST(HalfspaceTreeGeneration, Cube) {
	typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
	typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
	typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;

	Polyhedron cube;
	createCube(cube, 1.0);
	Nef_polyhedron cube_nef(cube);

	auto graph = build_facet_edge_graph(cube_nef);
	auto tree = build_halfspace_tree(graph, cube_nef);
	ASSERT_EQ(cube_nef, tree->evaluate());
}
