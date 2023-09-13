#include <gtest/gtest.h>
#include <test_utils.h>
#include <nef_to_halfspace_tree.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;

TEST(HalfspaceTreeGeneration, Cube) {
	Polyhedron cube, cube2;
	createCube(cube, 1.0);
	createCube(cube2, 1.0);
	translatePolyhedron(cube2, Kernel::Vector_3(0.5, 0.5, 0));
	Nef_polyhedron cube_nef(cube);
	Nef_polyhedron cube2_nef(cube2);

	auto concave = cube_nef - cube2_nef;

	auto graph = build_facet_edge_graph(concave);
	auto tree = build_halfspace_tree<Kernel>(graph, concave);
	ASSERT_EQ(concave, tree->evaluate());
}
