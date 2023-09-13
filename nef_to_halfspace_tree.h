#ifndef NEF_TO_HALFSPACE_TREE_H
#define NEF_TO_HALFSPACE_TREE_H


class halfspace_tree {
public:
	enum halfspace_operation {
		OP_UNION, OP_SUBTRACTION, OP_INTERSECTION
	};

	virtual halfspace_tree_result evaluate(int level = 0) const = 0;
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

		return result;
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


std::unique_ptr<halfspace_tree> build_halfspace_tree(Graph& G, Nef_polyhedron& poly, bool negate = false) {
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

		convex_subcomponent_visitor<decltype(components)> visitor(edge_trait, components, sub_graph_0);
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

#endif