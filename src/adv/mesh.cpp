#include <adv/mesh.hpp>

namespace phys{

	double mesh<triangle>::maxNorm(const std::vector<double>& v)
	{
		double result = 0.;
		for(size_t i = 0; i < 3; ++i)
			if(squaredNorm(item[i]()) > 0.01)
				result = maths::max(result, abs(v[item[i].get_index()])); 
		if(next)
			result = maths::max(result, ((mesh<triangle>*)next)->maxNorm(v));
		return result;
	}

	void mesh<triangle>::refineBoundary(size_t levels)
	{
		mesh<triangle>* runner = this;
		for(size_t i = 1; i < levels; ++i)
			for(int j = 0; j < maths::power(2,int(i)); ++j)
			{
				//find mid_point between first two verticies
				point2 vertex0 = (*runner)()[0]();
				point2 vertex1 = (*runner)()[1]();
				point2 midpoint = (vertex0 + vertex1)/2.;

				//find the angle subtended at the origin of coordinate system
				double angle1 = acos(vertex1[0]);
				double delta_angle = acos(sqrt(squaredNorm(midpoint)));
				if(j >= maths::power(2,int(i-1))){
					angle1 = -angle1;
					delta_angle = -delta_angle;
				}
				//angle increment on the circular boundary				
				double mid_angle = angle1 + delta_angle;

				//compute new vertex position based on the midpoint angle.
				point2 new_point(cos(mid_angle), sin(mid_angle));
				node<point2> new_vertex(new_point);
				//create new triangle element using the new vetex and append to the mesh.
				triangle t1(runner->item(0), new_vertex, runner->item(1));
				append(t1);

				//repeat for second new triangle
				mid_angle = angle1 - delta_angle;
				new_vertex = node<point2>(point2(cos(mid_angle), sin(mid_angle)));
				triangle t2(runner->item(1), new_vertex, runner->item(2));
				append(t2);
				//move runner to next triangle element.
				runner = (mesh<triangle>*)runner->next;
			}
	}

	void mesh<triangle>::refine_neighbour(	node<point2>& nI,
											node<point2>& nJ,
											node<point2>& nIJ )
	{
		int ni = nI < item; //is vertex (nI) in element (item)
		int nj = nJ < item;

		if (ni && nj)
		{
			--ni;
			--nj;
			int nk = 0;
			while ((nk == ni) || (nk == nj))
				++nk;
			triangle t1(nI, nIJ, item(nk));
			triangle t2(nJ, nIJ, item(nk));
			insert_next_item(t2);
			insert_next_item(t1);
			drop_first_item();
		}
		else
		{
			if (next)
				((mesh<triangle>*)next)->refine_neighbour(nI, nJ, nIJ);
			else
			{
				node<point2> new_node((1. / sqrt(squaredNorm(nIJ()))) * nIJ());
				triangle t1(nI, nIJ, new_node);
				triangle t2(nJ, nIJ, t1(2));
				insert_next_item(t2);
				insert_next_item(t1);
			}
		}
	}

	void mesh<triangle>::refine(const std::vector<double>& v, double thresh)
	{
		for (size_t i = 0; i < 3; ++i)
			for (size_t j = i + 1; j < 3; ++j)//j = 2; j > i; --j
				if ((item[i].get_index() >= 0) &&
					(item[j].get_index() >= 0) &&
					(abs(v[item[i].get_index()]
					- v[item[j].get_index()])  > thresh))
				{
					node<point2> itemij = (item[i]() + item[j]()) / 2.;
					int k = 0;
					while ((k == i) || (k == j))
						++k;
					triangle t1(item(i), itemij, item(k));
					triangle t2(item(j), itemij, item(k));
					if (next)
						((mesh<triangle>*)next)->refine_neighbour(item(i), item(j), t1(1));
					insert_next_item(t2);
					insert_next_item(t1);
					drop_first_item();
					refine(v, thresh);
					return;
				}
		if (next)
			((mesh<triangle>*)next)->refine(v, thresh);
	}


}