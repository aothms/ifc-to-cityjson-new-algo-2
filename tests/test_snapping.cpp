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

	Polyhedron cube, cube1, cube2;
	createCube(cube, 1.0);
	createCube(cube1, 1.0);
	createCube(cube2, 1.0);

	translatePolyhedron(cube1, Kernel::Vector_3(0.001, 0.001, 0));
	rotatePolyhedron(cube2, 0.01, Kernel::Vector_3(0, 0, 1.0));

	Nef_polyhedron cube_nef(cube);
	Nef_polyhedron cube1_nef(cube1);
	Nef_polyhedron cube2_nef(cube2);

	std::list<Nef_polyhedron*> operands = { &cube_nef, &cube1_nef, &cube2_nef };

	{
		// standard, apparently the result is 40 verts

		CGAL::Nef_nary_union_3<Nef_polyhedron> builder;
		for (auto& p : operands) {
			builder.add_polyhedron(*p);
		}
		auto result = builder.get_union();

		ASSERT_EQ(result.number_of_vertices(), 40);
	}

	{
		std::list<Kernel::Plane_3> planes;
		std::list<std::unique_ptr<halfspace_tree<Kernel>>> trees, snapped;

		for (auto& p : operands) {
			auto G = build_facet_edge_graph(*p);
			trees.emplace_back(std::move(build_halfspace_tree(G, *p)));
			trees.back()->accumulate(planes);
		}

		auto plane_map = snap_halfspaces(planes, 0.5);
		for (auto& x : trees) {
			snapped.push_back(std::move(x->map(plane_map)));
		}

		CGAL::Nef_nary_union_3<Nef_polyhedron> builder;
		for (auto& p : snapped) {
			builder.add_polyhedron(p->evaluate());
		}
		auto result = builder.get_union();

		ASSERT_EQ(result.number_of_vertices(), 8);
	}
}
