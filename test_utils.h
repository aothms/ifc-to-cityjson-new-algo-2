#ifndef TEST_UTILS_H
#define TEST_UTILS_H

template <typename LoopType, typename Kernel>
void extrude(LoopType bottom, const CGAL::Vector_3<Kernel>& V, CGAL::Polyhedron_3<Kernel>& P) {
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

template <typename Kernel>
void createCube(CGAL::Polyhedron_3<Kernel>& P, double d) {
	quad bottom = { {
		Point(-d, -d, -d),
		Point(+d, -d, -d),
		Point(+d, +d, -d),
		Point(-d, +d, -d)
	} };

	Kernel::Vector_3 V(0, 0, d * 2);

	extrude(bottom, V, P);
}

#endif