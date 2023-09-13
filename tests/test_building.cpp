#include <gtest/gtest.h>
#include <test_utils.h>
#include <nef_to_halfspace_tree.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Nef_polyhedron::Vertex_const_handle Vertex_const_handle;
typedef Nef_polyhedron::Halfedge_const_handle Halfedge_const_handle;
typedef Nef_polyhedron::Halffacet_const_handle Halffacet_const_handle;
typedef Nef_polyhedron::SHalfedge_const_handle SHalfedge_const_handle;
typedef Nef_polyhedron::SHalfloop_const_handle SHalfloop_const_handle;
typedef Nef_polyhedron::SFace_const_handle SFace_const_handle;
typedef Nef_polyhedron::Volume_const_iterator Volume_const_iterator;
typedef Nef_polyhedron::Shell_entry_const_iterator Shell_entry_const_iterator;

class Shell_explorer {
public:
	void visit(Vertex_const_handle) {}
	void visit(Halfedge_const_handle) {}
	void visit(Halffacet_const_handle h) {
		std::cout << "facet mark " << h->mark() << " " << h->plane() <<  std::endl;
	}
	void visit(SHalfedge_const_handle) {}
	void visit(SHalfloop_const_handle) {}
	void visit(SFace_const_handle) {}
};


TEST(HalfspaceTreeGeneration, Cube) {
	typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
	typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
	typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;

	std::list<Nef_polyhedron> elements;
	
	// foundation
	// std::vector<int> included = { 54, 97, 13, 126, 1, 117, 12};
	// front wall
	// std::vector<int> included = { 133, 51, 52, 53, 76, 77, 78 };
	// for (auto& i : included) {
	
	// problematic window
	// 30, 31, 32, 33, 89

	// co-planar
	// 47
	
	// for (size_t i = 1; i <= 133; ++i) {

	// std::vector<int> included = { 32, 89 }; // 30, 31, 32, 33,

	// std::vector<int> included = { 47 };

	// for (size_t i = 1; i <= 133; ++i) {
	
	std::vector<int> included = { 32, 89 };
	for (auto& i : included) {
	
		std::ostringstream oss;
		oss << std::setw(5) << std::setfill('0') << i << ".off";
		auto fn = oss.str();
		Polyhedron P;
		std::ifstream(fn.c_str()) >> P;
		ASSERT_TRUE(P.is_valid());
		ASSERT_TRUE(P.is_closed());
		CGAL::Polygon_mesh_processing::triangulate_faces(P);
		CGAL::Polygon_mesh_processing::remove_degenerate_faces(P);
		ASSERT_FALSE(CGAL::Polygon_mesh_processing::does_self_intersect(P)) << i << " self intersects" << std::endl;
		elements.emplace_back(P);
	}

	{
		CGAL::Nef_nary_union_3<Nef_polyhedron> builder;
		for (auto& p : elements) {
			builder.add_polyhedron(p);

			for (auto it = p.volumes_begin(); it != p.volumes_end(); ++it) {
				if (!it->mark()) {
					continue;
				}
				for (auto jt = it->shells_begin(); jt != it->shells_end(); ++jt) {
					Shell_explorer SE;
					p.visit_shell_objects(SFace_const_handle(jt), SE);
				}
			}

		}
		auto result = builder.get_union();

		int i = 0;
		for (auto it = result.volumes_begin(); it != result.volumes_end(); ++it) {
			if (!it->mark()) {
				continue;
			}
			for (auto jt = it->shells_begin(); jt != it->shells_end(); ++jt, ++i) {
				Polyhedron P;
				result.convert_inner_shell_to_polyhedron(jt, P);
				std::string s = "default-shell-" + std::to_string(i) + ".off";
				std::ofstream ofs(s.c_str());
				ofs << P;
			}
		}

		EXPECT_GT(result.number_of_volumes(), 2);
	}

	{
		std::list<Kernel::Plane_3> planes;
		std::list<std::unique_ptr<halfspace_tree<Kernel>>> trees, snapped;

		for (auto& p : elements) {
			std::cout << "== OPERAND ==" << std::endl;

			auto G = build_facet_edge_graph(p);
			std::cout << "Removed " << edge_contract(G) << " facets" << std::endl;

			dump_facets(G);

			trees.emplace_back(std::move(build_halfspace_tree(G, p)));
			trees.back()->accumulate(planes);
		}

		auto plane_map = snap_halfspaces(planes, 0.005);

		for (auto& p0 : plane_map) {
			for (auto& p1 : plane_map) {
				if (p0.second.opposite() == p1.second) {
					std::cout << "YES!" << std::endl;
				}
			}
		}

		for (auto& x : trees) {
			snapped.push_back(std::move(x->map(plane_map)));
		}

		int ii = 0;

		CGAL::Nef_nary_union_3<Nef_polyhedron> builder;

		// Nef_polyhedron result;

		for (auto& t : snapped) {
			// closure necessary?
			auto p = t->evaluate().closure().regularization();

			{
				Polyhedron P;
				p.convert_to_polyhedron(P);

				std::string s = "snapped-operand-" + std::to_string(ii++) + ".off";
				std::ofstream ofs(s.c_str());
				ofs << P;
			}

			for (auto it = p.volumes_begin(); it != p.volumes_end(); ++it) {
				if (!it->mark()) {
					continue;
				}
				for (auto jt = it->shells_begin(); jt != it->shells_end(); ++jt) {
					Shell_explorer SE;
					p.visit_shell_objects(SFace_const_handle(jt), SE);
				}
			}
			
			builder.add_polyhedron(p);
			// result = (result + p).regularization();
		}

		// regularization necessary?
		auto result = builder.get_union().regularization();

		size_t vi = 0;
		Polyhedron P;
		while (convert_to_polyhedron(result, P, vi++)) {
			std::string s = "snapped-shell-" + std::to_string(vi) + ".off";
			std::ofstream ofs(s.c_str());
			ofs << P;
		}

		EXPECT_TRUE(result.is_simple());

		ASSERT_EQ(result.number_of_volumes(), 2);
	}
}
