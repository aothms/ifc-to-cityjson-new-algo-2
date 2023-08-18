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

template <typename LoopType>
void extrude(LoopType bottom, const Kernel::Vector_3& V, Polyhedron& P) {
	std::list<LoopType> face_list = { bottom };

	for (auto current_vertex = bottom.begin(); current_vertex != bottom.end(); ++current_vertex) {

		auto next_vertex = current_vertex;
		++next_vertex;

		if (next_vertex == bottom.end()) {
			next_vertex = bottom.begin();
		}

		LoopType side = { {
			*next_vertex,
			*current_vertex,
			*current_vertex + V,
			*next_vertex + V	  ,
		} };

		face_list.push_back(side);
	}

	auto top = bottom;
	for (auto& v : top) {
		v += V;
	}
	std::reverse(top.begin(), top.end());

	face_list.push_back(top);

	std::vector<Point> unique_points;
	std::vector<std::vector<std::size_t>> facet_vertices;
	std::map<Point, size_t> points;

	for (auto &face : face_list) {
		facet_vertices.emplace_back();
		for (auto &point : face) {
			auto p = points.insert({ point, points.size() });
			if (p.second) {
				unique_points.push_back(point);
			}
			facet_vertices.back().push_back(p.first->second);
		}
	}

	CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(unique_points, facet_vertices, P);
}

void createCube(Polyhedron& P, double d) {
	quad bottom = { {
		Point(-d, -d, -d),
		Point(+d, -d, -d),
		Point(+d, +d, -d),
		Point(-d, +d, -d)
	} };

	Kernel::Vector_3 V(0, 0, d * 2);
	
	extrude(bottom, V, P);
}

void rotatePolyhedron(Polyhedron& P, double angle, Kernel::Vector_3 axis) {
	axis = axis / std::sqrt(CGAL::to_double(axis.squared_length()));

	double c = std::cos(angle);
	double s = std::sin(angle);

	// Compute the 4x4 rotation matrix
	CGAL::Aff_transformation_3<Kernel> rotation_matrix(
		c + (1 - c) * axis.x() * axis.x(),
		(1 - c) * axis.x() * axis.y() - s * axis.z(),
		(1 - c) * axis.x() * axis.z() + s * axis.y(),
		0.0,
		(1 - c) * axis.y() * axis.x() + s * axis.z(),
		c + (1 - c) * axis.y() * axis.y(),
		(1 - c) * axis.y() * axis.z() - s * axis.x(),
		0.0,
		(1 - c) * axis.z() * axis.x() - s * axis.y(),
		(1 - c) * axis.z() * axis.y() + s * axis.x(),
		c + (1 - c) * axis.z() * axis.z(),
		0.0
	);

	for (auto it = P.vertices_begin(); it != P.vertices_end(); ++it) {
		it->point() = rotation_matrix.transform(it->point());
	}
}

template <typename T>
void scalePolyhedron(Polyhedron& P, T sx, T sy, T sz) {
	CGAL::Aff_transformation_3<Kernel> scale(
		sx, 0, 0, 0,
		0, sy, 0, 0,
		0, 0, sz, 0);
	for (auto it = P.vertices_begin(); it != P.vertices_end(); ++it) {
		it->point() = scale.transform(it->point());
	}
}

void translatePolyhedron(Polyhedron& P, Kernel::Vector_3 v) {
	CGAL::Aff_transformation_3<Kernel> rotation_matrix(
		CGAL::TRANSLATION, v
	);

	for (auto it = P.vertices_begin(); it != P.vertices_end(); ++it) {
		it->point() = rotation_matrix.transform(it->point());
	}
}

void convertToNefPolyhedron(Polyhedron& P, Nef_polyhedron& NP) {
	NP = Nef_polyhedron(P);
}

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


template <typename Visitor>
void custom_visit(SFace_const_handle f, Visitor& V)
{
	std::list<std::pair<Halffacet_const_handle, SFace_const_handle>> SFaceCandidates;
	std::list<std::pair<Halffacet_const_handle, Halffacet_const_handle>> FacetCandidates;
	CGAL::Unique_hash_map<SFace_const_handle, bool> DoneSF(false);
	CGAL::Unique_hash_map<Vertex_const_handle, bool> DoneV(false);
	CGAL::Unique_hash_map<SVertex_const_handle, bool> DoneSV(false);
	CGAL::Unique_hash_map<Halffacet_const_handle, bool> DoneF(false);
	SFaceCandidates.push_back({ nullptr, f });  DoneSF[f] = true;
	while (true) {
		if (SFaceCandidates.empty() && FacetCandidates.empty()) break;
		if (!FacetCandidates.empty()) {
			auto pp = *FacetCandidates.begin();
			Halffacet_const_handle f = pp.second;
			FacetCandidates.pop_front();

			std::cout << "facet " << f->plane() << std::endl;

			V.visit_facet(pp.first, f); // report facet
			Halffacet_cycle_const_iterator fc;
			CGAL_forall_facet_cycles_of(fc, f) {
				if (fc.is_shalfedge()) {
					SHalfedge_const_handle e(fc);
					SHalfedge_const_handle she;
					SHalfedge_around_facet_const_circulator ec(e), ee(e);
					CGAL_For_all(ec, ee) {

						// cannot figure out how to compute CCW on sphere points

						// auto p2 = ec->target()->point();
						// auto p1 = ec->source()->point();
						// auto p0 = ec->prev()->source()->point();
						// 
						// auto e1 = p1 - p0;
						// auto e2 = p2 - p1;
						// std::cout << ">> " << e1 << " " << e2 << std::endl;

						she = ec->twin();
						if (DoneSF[she->incident_sface()]) continue;
						SFaceCandidates.push_back({ f, she->incident_sface() });
						DoneSF[she->incident_sface()] = true;
					}
				} else if (fc.is_shalfloop()) {
					SHalfloop_const_handle l(fc);
					SHalfloop_const_handle ll = l->twin();
					if (DoneSF[ll->incident_sface()]) continue;
					SFaceCandidates.push_back({ f, ll->incident_sface() });
					DoneSF[ll->incident_sface()] = true;
				} else CGAL_error_msg("Damn wrong handle.");
			}
		}
		if (!SFaceCandidates.empty()) {
			auto pp = *SFaceCandidates.begin();
			SFace_const_handle sf = pp.second;
			SFaceCandidates.pop_front();
			// V.visit(sf);
			// if (!DoneV[sf->center_vertex()])
			// 	V.visit(sf->center_vertex()); // report vertex
			DoneV[sf->center_vertex()] = true;
			//      SVertex_const_handle sv;
			SM_const_decorator SD(&*sf->center_vertex());
			/*
			CGAL_forall_svertices(sv,SD){
			  if(SD.is_isolated(sv) && !DoneSV[sv])
				V.visit(sv);
			}
			*/
			SFace_cycle_const_iterator fc;
			CGAL_forall_sface_cycles_of(fc, sf) {
				if (fc.is_shalfedge()) {
					SHalfedge_const_handle e(fc);
					SHalfedge_around_sface_const_circulator ec(e), ee(e);
					CGAL_For_all(ec, ee) {
						// V.visit(SHalfedge_const_handle(ec));
						SVertex_const_handle vv = ec->twin()->source();
						if (!SD.is_isolated(vv) && !DoneSV[vv]) {
							// V.visit(vv); // report edge
							DoneSV[vv] = DoneSV[vv->twin()] = true;
						}
						Halffacet_const_handle f = ec->twin()->facet();
						if (DoneF[f]) continue;
						FacetCandidates.push_back({ pp.first, f }); DoneF[f] = true;
					}
				} else if (fc.is_svertex()) {
					SVertex_const_handle v(fc);
					if (DoneSV[v]) continue;
					// V.visit(v); // report edge
					// V.visit(v->twin());
					DoneSV[v] = DoneSV[v->twin()] = true;
					CGAL_assertion(SD.is_isolated(v));
					SFaceCandidates.push_back({ pp.first, v->twin()->incident_sface() });
					DoneSF[v->twin()->incident_sface()] = true;
					// note that v is isolated, thus twin(v) is isolated too
					//          SM_const_decorator SD;
					//          SFace_const_handle fo;
					//          fo = v->twin()->incident_sface();
					/*
					if(SD.is_isolated(v))
					  fo = v->source()->sfaces_begin();
					else
					  fo = v->twin()->incident_sface();
					*/
				} else if (fc.is_shalfloop()) {
					SHalfloop_const_handle l(fc);
					// V.visit(l);
					Halffacet_const_handle f = l->twin()->facet();
					if (DoneF[f]) continue;
					FacetCandidates.push_back({ pp.first, f }); DoneF[f] = true;
				} else CGAL_error_msg("Damn wrong handle.");
			}
		}
	}
}

class Nef_explorer {
public:
	Nef_explorer() {}
	void visit(Vertex_const_handle v) {
		// std::cout << "vertex" << std::endl;
	}
	void visit(Halfedge_const_handle) {
		// std::cout << "half edge" << std::endl;
	}
	void visit(Halffacet_const_handle h) {
		// std::cout << "half facet" << std::endl;
		// const auto& p = h->plane();
		// std::cout << "  " << CGAL::to_double(p.a()) << " " << CGAL::to_double(p.b()) << " " << CGAL::to_double(p.c()) << " " << CGAL::to_double(p.d()) << std::endl;
	}
	void visit(SHalfedge_const_handle) {
		// std::cout << "s half edge" << std::endl;
	}
	void visit(SHalfloop_const_handle) {
		// std::cout << "s half loop" << std::endl;
	}
	void visit(SFace_const_handle) {
		// std::cout << "s face" << std::endl;
	}
};

class capture_shell_vertices {
public:
	std::list<Vertex_const_handle> vertices;
	capture_shell_vertices() {}
	void visit(Vertex_const_handle v) {
		vertices.push_back(v);
	}
	void visit(Halfedge_const_handle) {
	}
	void visit(Halffacet_const_handle h) {
	}
	void visit(SHalfedge_const_handle) {
	}
	void visit(SHalfloop_const_handle) {
	}
	void visit(SFace_const_handle) {
	}
};

class facet_explorer {
public:
	facet_explorer() {}
	void visit_facet(Halffacet_const_handle h1, Halffacet_const_handle h2) {
		{
			const auto& p = h1->plane();
			std::cout << CGAL::to_double(p.a()) << " " << CGAL::to_double(p.b()) << " " << CGAL::to_double(p.c()) << " " << CGAL::to_double(p.d()) << std::endl;
		}
		{
			const auto& p = h2->plane();
			std::cout << " ->  " << CGAL::to_double(p.a()) << " " << CGAL::to_double(p.b()) << " " << CGAL::to_double(p.c()) << " " << CGAL::to_double(p.d()) << std::endl;
		}
	}
};

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

enum halfspace_operation {
	OP_UNION, OP_SUBTRACTION, OP_INTERSECTION
};

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

typedef CGAL::Epick_d<CGAL::Dimension_tag<4>> KdKernel;
typedef KdKernel::Point_d Point_d;
typedef CGAL::Search_traits_d<KdKernel> TreeTraits;
typedef CGAL::Kd_tree<TreeTraits> Tree;
typedef CGAL::Fuzzy_sphere<TreeTraits> Fuzzy_sphere;
typedef std::map<Kernel::Plane_3, Kernel::Plane_3, PlaneLess> plane_map_t;

class halfspace_tree {
public:
	virtual halfspace_tree_result evaluate(int level=0) const = 0;
	virtual void accumulate(std::list<Kernel::Plane_3>&) const = 0;
	virtual std::unique_ptr<halfspace_tree> map(const plane_map_t&) const = 0;
};

class halfspace_tree_nary_branch : public halfspace_tree {
private:
	halfspace_operation operation_;
	std::list<std::unique_ptr<halfspace_tree>> operands_;

public:
	halfspace_tree_nary_branch(halfspace_operation operation, std::list<std::unique_ptr<halfspace_tree>>&& operands)
		: operation_(operation)
		, operands_(std::move(operands))
	{}
	virtual halfspace_tree_result evaluate(int level) const {
		static const char* const ops[] = { "union", "subtraction", "intersection" };
		std::cout << std::string(level * 2, ' ') << ops[operation_] << " (" << std::endl;

		Nef_polyhedron result;

		if (operation_ == OP_SUBTRACTION) {
			if (operands_.size() != 2) {
				throw std::runtime_error("");
			}
			result = operands_.front()->evaluate(level + 1) - operands_.back()->evaluate(level + 1);
		} else if (operation_ == OP_UNION) {
			CGAL::Nef_nary_union_3<halfspace_tree_result> builder;
			for (auto& op : operands_) {
				builder.add_polyhedron(op->evaluate(level + 1));
			}
			result = builder.get_union();
		} else if (operation_ == OP_INTERSECTION) {
			CGAL::Nef_nary_intersection_3<Nef_polyhedron> builder;
			for (auto& op : operands_) {
				builder.add_polyhedron(op->evaluate(level + 1));
			}
			result = builder.get_intersection();

			/*
			halfspace_tree_result r;
			bool first = true;
			for (auto& op : operands_) {
				if (first) {
					r = op->evaluate(level + 1);
					first = false;
					continue;
				}
				r = r * op->evaluate(level + 1);
			}
			return r;
			*/
		}

		std::cout << std::string(level * 2, ' ') << ")" << std::endl;

		return result
	}
	virtual void accumulate(std::list<Kernel::Plane_3>& points) const {
		for (auto& op : operands_) {
			op->accumulate(points);
		}
	}
	virtual std::unique_ptr<halfspace_tree> map(const plane_map_t& m) const {
		decltype(operands_) mapped;
		for (auto& op : operands_) {
			mapped.emplace_back(op->map(m));
		}
		return std::unique_ptr<halfspace_tree>(new halfspace_tree_nary_branch(operation_, std::move(mapped)));
	}
};

class halfspace_tree_plane : public halfspace_tree {
private:
	Kernel::Plane_3 plane_;

public:
	halfspace_tree_plane(const Kernel::Plane_3& plane)
		: plane_(plane)
	{}
	virtual halfspace_tree_result evaluate(int level) const {

#if defined(HALFSPACE_EXTENDED)
		Kernel_b::Plane_3 plane(plane_.a().exact(), plane_.b().exact(), plane_.c().exact(), plane_.d().exact());
		return Nef_polyhedron_b(plane, Nef_polyhedron_b::Boundary::INCLUDED);
#else
		Kernel::Vector_3 ref;

		static int N = 1;

		std::cout << std::string(level * 2, ' ') << "p " << plane_ << std::endl;

		if (CGAL::abs(plane_.orthogonal_direction().dz()) > CGAL::abs(plane_.orthogonal_direction().dx())) {
			ref = Kernel::Vector_3(1, 0, 0);
		} else {
			ref = Kernel::Vector_3(0, 0, 1);
		}
		auto x = CGAL::cross_product(plane_.orthogonal_direction().vector(), ref);
		x /= std::sqrt(CGAL::to_double(x.squared_length()));
		auto y = CGAL::cross_product(plane_.orthogonal_direction().vector(), x);
		y /= std::sqrt(CGAL::to_double(y.squared_length()));

		std::vector<Point> unique_points = {
			plane_.point() + x * 5 - y * 5,
			plane_.point() - x * 5 - y * 5,
			plane_.point() - x * 5 + y * 5,
			plane_.point() + x * 5 + y * 5,
		};

		auto v = plane_.orthogonal_direction().vector();
		v /= std::sqrt(CGAL::to_double(v.squared_length()));
		v *= 1000.;

		/*
		for (auto& p : unique_points) {
			std::cout << "v " << p << std::endl;
		}
		std::cout << "f " << (N + 0) << " " << (N + 1) << " " << (N + 2) << " " << (N + 3) << std::endl;
		
		N += 4;

		std::cout << std::endl << std::endl;
		*/

		Polyhedron P;
#if defined(HALFSPACE_FACET)
		std::vector<std::vector<std::size_t>> facet_vertices = { {0,1,2,3} };
		CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(unique_points, facet_vertices, P);
#else
		extrude(unique_points, v, P);
#endif
		return Nef_polyhedron(P);

#endif
	}
	virtual void accumulate(std::list<Kernel::Plane_3>& points) const {
		points.push_back(plane_);
	}
	virtual std::unique_ptr<halfspace_tree> map(const plane_map_t& m) const {
		auto it = m.find(plane_);
		if (it != m.end()) {
			return std::unique_ptr<halfspace_tree>(new halfspace_tree_plane(it->second));
		} else {
			return std::unique_ptr<halfspace_tree>(new halfspace_tree_plane(plane_));
		}
	}
};

// Define edge types
enum EdgeType { CONCAVE, CONVEX };

struct VertexProperties {
	Halffacet_const_handle facet;
	size_t original_index;
};

std::vector<CGAL::Triangle_3<Kernel>> triangulate_nef_facet(Halffacet_const_handle f) {
	std::vector<Kernel::Point_3> ps;
	SHalfedge_around_facet_const_circulator it(f->facet_cycles_begin());
	SHalfedge_around_facet_const_circulator first(it);
	Vertex_const_handle lastvh;
	CGAL_For_all(it, first) {
		ps.push_back(it->source()->center_vertex()->point());
	}

	typedef std::tuple<int, int, int> indices_triple;
	std::vector<indices_triple> patches;
	CGAL::Polygon_mesh_processing::triangulate_hole_polyline(ps, std::back_inserter(patches));

	std::vector<CGAL::Triangle_3<Kernel>> result;
	result.reserve(patches.size());
	std::transform(patches.begin(), patches.end(), std::back_inserter(result), [&ps](indices_triple& t) {
		return CGAL::Triangle_3<Kernel>(ps[std::get<0>(t)], ps[std::get<1>(t)], ps[std::get<2>(t)]); 
	});

	return result;
}

// Define the graph type with properties for nodes and edges
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


// Custom visitor to filter edges based on their label
template <typename ComponentMap>
class facet_without_reflex_edge_visitor : public boost::default_bfs_visitor {

private:
	EdgeType edgetype_;
	ComponentMap& components_;
	const Graph& graph_;

public:
	facet_without_reflex_edge_visitor(EdgeType edgetype, ComponentMap& component, const Graph& g)
		: edgetype_(edgetype), components_(component), graph_(g) {}

	template <typename Edge, typename Graph>
	void tree_edge(Edge e, const Graph& g) {
		if (boost::get(boost::edge_weight, graph_, e) == edgetype_) {
			auto srcid = boost::source(e, g);
			auto tgtid = boost::target(e, g);

			auto src = &components_[srcid];
			auto tgt = &components_[tgtid];

			if (*src == -1) {
				std::swap(srcid, tgtid);
				std::swap(src, tgt);
			} 
			
			if (*tgt == -1) {
				bool tgt_has_any_reflex_edge = false;

				typename boost::graph_traits<Graph>::out_edge_iterator ei, ei_end;
				for (boost::tie(ei, ei_end) = boost::out_edges(tgtid, g); ei != ei_end; ++ei) {
					if (boost::get(boost::edge_weight, graph_, *ei) != edgetype_) {
						tgt_has_any_reflex_edge = true;
						break;
					}
				}

				if (!tgt_has_any_reflex_edge) {
					// topological check completed, now check geometry, non topologically connected facets should not geometrically intersect
					bool any_intersecting = false;
					auto& plane_tgt = graph_[tgtid].facet->plane();

					std::cout << "plane_tgt " << plane_tgt << std::endl;

					for (size_t i = 0; i < components_.size(); ++i) {
						if (i == tgtid) {
							continue;
						}

						if (components_[i] == *src) {
							// no need to check for intersection when edge exists
							const bool has_edge = boost::edge(i, tgtid, g).second;
							if (!has_edge) {
								auto triangles_i = triangulate_nef_facet(graph_[i].facet);

								std::cout << "i " << i << ": " << std::endl;
								for (auto& t : triangles_i) {
									std::cout << "  " << t << std::endl;
								}
								
								if (std::any_of(triangles_i.begin(), triangles_i.end(), [&plane_tgt](CGAL::Triangle_3<Kernel>& t) {
									auto x = CGAL::intersection(plane_tgt, t);
									if (x) {
										std::cout << "t " << t << " x " <<  std::endl;
										Intersection_visitor v;
										boost::apply_visitor(v)(*x);
									}
									return (bool)x;
								})) {
									any_intersecting = true;
									break;
								}
							}
						}
					}
					if (!any_intersecting) {
						std::cout << "v " << srcid << " -> " << tgtid << std::endl;
						std::cout << "(" << *src << " " << *tgt << ")" << std::endl;

						*tgt = *src;
					}					
				}
			}
		}
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

Graph build_facet_edge_graph(const Nef_polyhedron& poly) {
	Graph G;

	std::map<Halffacet_const_handle, size_t> facet_to_idx;

	for (auto it = poly.vertices_begin(); it != poly.vertices_end(); ++it) {
		for (auto a = it->shalfedges_begin(); a != it->shalfedges_end(); ++a) {
			auto b = a->snext();
			if (a->facet()->incident_volume()->mark() && a->facet()->incident_volume() == b->facet()->incident_volume() && a->facet()->is_valid() && b->facet()->is_valid()) {
				// cross product of facet normals
				auto avec = a->facet()->plane().orthogonal_vector();
				auto bvec = b->facet()->plane().orthogonal_vector();
				auto cross = CGAL::cross_product(avec, bvec);

				size_t aidx, bidx;
				{
					auto p = facet_to_idx.insert({ a->facet(), facet_to_idx.size() });
					aidx = p.first->second;
					if (p.second) {
						auto& v = G[boost::add_vertex(G)];
						v.facet = a->facet();
						v.original_index = aidx;
					}
				}
				{
					auto p = facet_to_idx.insert({ b->facet(), facet_to_idx.size() });
					bidx = p.first->second;
					if (p.second) {
						auto& v = G[boost::add_vertex(G)];
						v.facet = b->facet();
						v.original_index = bidx;
					}
				}

				// half edge direction
				auto ref = a->next()->source()->center_vertex()->point() - a->source()->center_vertex()->point();
				std::cout << ref << " . " << cross << std::endl;

				// check whether half edge direction conforms to facet normal cross
				auto edge_type = CGAL::scalar_product(cross, ref) > 0
					? CONCAVE
					: CONVEX;

				std::cout << aidx << " -> " << bidx << " " << edge_type << std::endl;

				boost::add_edge(aidx, bidx, edge_type, G);
			}
		}
	}

	for (size_t ii = 0; ii < boost::num_vertices(G); ++ii) {
		std::cout << ii << " " << dump_facet(G[ii].facet) << std::endl;
	}

	return G;
}

std::unique_ptr<halfspace_tree> build_halfspace_tree(Graph& G, Nef_polyhedron& poly, bool negate=false) {
	auto edge_trait = negate ? CONCAVE : CONVEX;

	bool all_convex = true;
	boost::graph_traits<Graph>::edge_iterator ei, ei_end;
	for (boost::tie(ei, ei_end) = boost::edges(G); ei != ei_end; ++ei) {
		if (boost::get(boost::edge_weight, G, *ei) != edge_trait) {
			// all_convex = false;
			break;
		}
	}

	if (!all_convex) {
		std::unique_ptr<halfspace_tree> tree;
		std::list<std::unique_ptr<halfspace_tree>> root_expression;

		CGAL::convex_decomposition_3(poly);
		// the first volume is the outer volume, which is
		// ignored in the decomposition
		Volume_const_iterator ci = ++poly.volumes_begin();
		int NN = 0;

		for (; ci != poly.volumes_end(); ++ci, ++NN) {
			std::list<std::unique_ptr<halfspace_tree>> sub_expression;

			if (ci->mark()) {
				// @todo couldn't get it to work with the multiple volumes of a complex decomposition
				// directly, so for now we need to isolate the individual volumes.
				Polyhedron P;
				poly.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);
				Nef_polyhedron Pnef(P);
				auto Pgraph = build_facet_edge_graph(Pnef);
				for (size_t ii = 0; ii < boost::num_vertices(Pgraph); ++ii) {
					sub_expression.emplace_back(new halfspace_tree_plane(Pgraph[ii].facet->plane()));
				}
			}

			root_expression.emplace_back(new halfspace_tree_nary_branch(OP_INTERSECTION, std::move(sub_expression)));
		}

		tree.reset(new halfspace_tree_nary_branch(OP_UNION, std::move(root_expression)));
		return tree;
	}

	boost::write_graphviz(std::cout, G);

	std::unique_ptr<halfspace_tree> tree_0;
	std::list<std::unique_ptr<halfspace_tree>> root_expression_0;

	// First isolate into completely loose components
	std::vector<int> components_0(boost::num_vertices(G), -1);
	boost::connected_components(G, &components_0[0]);

	for (size_t i = 0; i <= *std::max_element(components_0.begin(), components_0.end()); ++i) {
		std::unique_ptr<halfspace_tree> tree;
		std::list<std::unique_ptr<halfspace_tree>> root_expression;

		auto included_in_vertex_subset_0 = [&components_0, &i](Graph::vertex_descriptor vd) {
			return components_0[vd] == i;
		};
		Filtered sub_graph_filtered_0(G, boost::keep_all{}, included_in_vertex_subset_0);
		Graph sub_graph_0;
		boost::copy_graph(sub_graph_filtered_0, sub_graph_0);

		std::vector<int> components(boost::num_vertices(sub_graph_0), -1);
		int largest_component_idx = -1;

		facet_without_reflex_edge_visitor<decltype(components)> visitor(edge_trait, components, sub_graph_0);
		int num_components = 0;
		for (size_t i = 0; i < boost::num_vertices(sub_graph_0); ++i) {
			if (components[i] == -1) {
				components[i] = num_components++;
				boost::breadth_first_search(sub_graph_0, boost::vertex(i, sub_graph_0), boost::visitor(visitor));
			}
		}

		std::ostream_iterator<size_t> output(std::cout, " ");
		std::cout << "components: " << std::endl;
		std::map<size_t, size_t> comp;
		size_t iii = 0;
		for (auto it = components.begin(); it != components.end(); ++it, ++iii) {
			comp[sub_graph_0[iii].original_index] = *it;
		}
		iii = comp.rbegin()->first;
		for (size_t iiii = 0; iiii <= iii; ++iiii) {
			std::cout << std::setw(2) << iiii << " ";
		}
		std::cout << std::endl;
		for (size_t iiii = 0; iiii <= iii; ++iiii) {
			auto it = comp.find(iiii);
			if (it == comp.end()) {
				std::cout << "__ ";
			} else {
				std::cout << std::setw(2) << it->second << " ";
			}
		}
		std::cout << std::endl;

		std::vector<size_t> component_size(num_components, 0);

		for (auto& idx : components) {
			component_size[idx] ++;
		}

		largest_component_idx = std::distance(component_size.begin(), std::max_element(component_size.begin(), component_size.end()));

		size_t vidx = 0;
		// Convex decomposition is not always optimal, sometimes contains coplanar facets
		std::unordered_set<Kernel::Plane_3, PlaneHash> plane_set;
		for (auto it = components.begin(); it != components.end(); ++it, ++vidx) {
			auto fct = sub_graph_0[vidx].facet;
			if (negate) {
				fct = fct->twin();
			}
			if (*it == largest_component_idx && plane_set.find(fct->plane()) == plane_set.end()) {
				root_expression.emplace_back(new halfspace_tree_plane(fct->plane()));
				plane_set.insert(fct->plane());
			}
		}

		tree.reset(new halfspace_tree_nary_branch(OP_INTERSECTION, std::move(root_expression)));

		auto included_in_vertex_subset = [&components, &largest_component_idx](Graph::vertex_descriptor vd) {
			return components[vd] != largest_component_idx;
		};
		Filtered sub_graph_filtered(sub_graph_0, boost::keep_all{}, included_in_vertex_subset);
		Graph sub_graph;
		// @todo can we go without this copy?
		boost::copy_graph(sub_graph_filtered, sub_graph);

		std::cout << "nv " << boost::num_vertices(sub_graph) << std::endl;

		// @nb counting vertices on filtered_graph returns the original amount
		if (boost::num_vertices(sub_graph)) {
			auto remainder = build_halfspace_tree(sub_graph, poly, !negate);

			std::list<std::unique_ptr<halfspace_tree>> sub_expression;
			sub_expression.emplace_back(std::move(tree));
			sub_expression.emplace_back(std::move(remainder));

			tree.reset(new halfspace_tree_nary_branch(OP_SUBTRACTION, std::move(sub_expression)));
		}

		root_expression_0.emplace_back(std::move(tree));
	}

	if (root_expression_0.size() == 1) {
		return std::move(root_expression_0.front());
	}
	tree_0.reset(new halfspace_tree_nary_branch(OP_UNION, std::move(root_expression_0)));
	return std::move(tree_0);
}

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


/*




int main() {
	// Create a graph
	Graph G;

	// Add nodes and edges
	boost::add_edge(0, 1, FRIENDSHIP, G);
	boost::add_edge(1, 2, FRIENDSHIP, G);
	boost::add_edge(2, 3, FRIENDSHIP, G);
	boost::add_edge(3, 4, FOLLOW, G);
	boost::add_edge(4, 5, FOLLOW, G);

	// Create a property map to extract components
	std::vector<int> component(boost::num_vertices(G), -1);

	// Extract connected components using only FRIENDSHIP edges
	int num_components = 0;
	for (size_t i = 0; i < boost::num_vertices(G); ++i) {
		if (component[i] == -1) {
			component[i] = num_components;
			friendship_edge_filter<std::vector<int>> vis(component, G);
			boost::breadth_first_search(G, boost::vertex(i, G), boost::visitor(vis));
			++num_components;
		}
	}

	// Print the connected components using FRIENDSHIP edges
	for (int i = 0; i < num_components; ++i) {
		std::cout << "Connected Component " << i << " (FRIENDSHIP): ";
		for (size_t j = 0; j < boost::num_vertices(G); ++j) {
			if (component[j] == i) {
				std::cout << j << " ";
			}
		}
		std::cout << std::endl;
	}

	return 0;
}
*/