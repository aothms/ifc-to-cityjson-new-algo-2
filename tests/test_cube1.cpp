#include <gtest/gtest.h>
#include <test_utils.h>
#include <nef_to_halfspace_tree.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>

#include <CGAL/Gmpq.h>
#include <CGAL/Extended_cartesian.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef CGAL::Extended_cartesian<CGAL::Gmpq> ExKernel;

template <typename Kernel>
class Shell_explorer {
	std::pair<bool, bool> has_plane_;
	CGAL::Plane_3<Kernel> plane_;
public:
	std::set<CGAL::Point_3<Kernel>> points;
	size_t nf = 0;

	typedef typename CGAL::Nef_polyhedron_3<Kernel>::Vertex_const_handle Vertex_const_handle;
	typedef typename CGAL::Nef_polyhedron_3<Kernel>::Halfedge_const_handle Halfedge_const_handle;
	typedef typename CGAL::Nef_polyhedron_3<Kernel>::Halffacet_const_handle Halffacet_const_handle;
	typedef typename CGAL::Nef_polyhedron_3<Kernel>::SHalfedge_const_handle SHalfedge_const_handle;
	typedef typename CGAL::Nef_polyhedron_3<Kernel>::SHalfloop_const_handle SHalfloop_const_handle;
	typedef typename CGAL::Nef_polyhedron_3<Kernel>::SFace_const_handle SFace_const_handle;

	Shell_explorer(const CGAL::Plane_3<Kernel>& plane)
		: plane_(plane) {}

	void visit(Vertex_const_handle h) {
		points.insert(h->point());
	}

	void visit(Halfedge_const_handle) {}
	void visit(SHalfedge_const_handle) {}
	void visit(SHalfloop_const_handle) {}
	void visit(SFace_const_handle) {}

	void visit(Halffacet_const_handle h) {
		++nf;
		if (h->plane() == plane_) {
			has_plane_.first = true;
		} else if (h->twin()->plane() == plane_) {
			has_plane_.second = true;
		}
	}

	const std::pair<bool, bool>& has_plane() const {
		return has_plane_;
	}
};

TEST(HalfspaceTreeGeneration, Cube) {


	CGAL::Plane_3<ExKernel> p2(0, 0, 1, 1);
	CGAL::Nef_polyhedron_3<ExKernel> plane_nef(p2, CGAL::Nef_polyhedron_3<ExKernel>::Boundary::EXCLUDED);
	CGAL::Nef_polyhedron_3<ExKernel> full_nef(CGAL::Nef_polyhedron_3<ExKernel>::COMPLETE);
	auto tv = full_nef - plane_nef;
	std::cout << tv.number_of_volumes() << std::endl;
	for (auto it = tv.volumes_begin(); it != tv.volumes_end(); ++it) {
		std::cout << "m " << it->mark() << std::endl;
		Shell_explorer<ExKernel> exp(p2);
		tv.visit_shell_objects(it->shells_begin(), exp);
		std::cout << exp.has_plane().first << " " << exp.has_plane().second << " " << exp.nf << std::endl;
		for (auto& p : exp.points) {
			std::cout << " " << p << std::endl;
		}
		std::cout << std::endl << std::endl;
	}

	Shell_explorer<ExKernel> exp(p2);
	plane_nef.visit_shell_objects(plane_nef.volumes_begin()->shells_begin(), exp);
	std::cout << exp.has_plane().first << " " << exp.has_plane().second << " " << exp.nf << std::endl;

	/*
	Polyhedron cube, cube2;
	createCube(cube, 1.0);
	createCube(cube2, 1.0);
	translatePolyhedron(cube2, Kernel::Vector_3(0.5, 0.5, 0));
	Nef_polyhedron cube_nef(cube);
	Nef_polyhedron cube2_nef(cube2);

	auto concave = cube_nef - cube2_nef;

	auto graph = build_facet_edge_graph(concave);
	auto tree = build_halfspace_tree<Kernel, ExKernel>(graph, concave);
	ASSERT_EQ(concave, tree->evaluate());
	*/
}
