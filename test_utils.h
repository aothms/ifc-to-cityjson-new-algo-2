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


#endif