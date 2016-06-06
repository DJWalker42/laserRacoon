#include <iomanip>
#include <cmath>

#include "Visualise.h"
#include "FormatOutput.h"

enum {
	BLUE, GREEN, RED, CYAN, PURPLE, YELLOW, VIOLET, LIME,
	PINK, GREY, WHITE, BLACK
};

const static cv::Scalar colour[] = { 
	cv::Scalar(255, 0, 0),		//blue
	cv::Scalar(0, 255, 0),		//green
	cv::Scalar(0, 0, 255),		//red
	cv::Scalar(255, 255, 0),	//cyan
	cv::Scalar(255, 0, 255),	//purple
	cv::Scalar(0, 255, 255),	//yellow										
	cv::Scalar(255, 127, 127),	//violet
	cv::Scalar(127, 255, 127),	//lime
	cv::Scalar(127, 127, 255),	//pink
	cv::Scalar(127, 127, 127),	//grey
	cv::Scalar(255, 255, 255),	//white
	cv::Scalar(0, 0, 0)			//black			
};

const static size_t colour_size = sizeof(colour) / sizeof(colour[0]);

namespace phys{
	namespace visual{

		// Helper functions
		//----------------------------------------------------------------------------
		namespace{

			/* Axis origin computed relative to the entire plot window */
			cv::Point compute_origin(	const cv::Rect& plot,
										const axis_props& x_axis,
										const axis_props& y_axis)
			{
				//X-axis
				double originX = x_axis.scale * x_axis.origin;
				if (x_axis.min * x_axis.max < 0.0) // max +ve, min -ve
				{				
					originX += plot.x + x_axis.scale * fabs(x_axis.min);
				}
				else if (fabs(x_axis.max) > fabs(x_axis.min)) //both +ve
				{
					originX = plot.x;
				}
				else //both -ve
				{
					originX = plot.x + plot.width;
				}

				//Y-axis
				double originY = y_axis.scale * y_axis.origin;
				if (y_axis.min * y_axis.max < 0.0)
				{
					originY = (plot.y + y_axis.scale * y_axis.max) - originY;

				}
				else if (fabs(y_axis.max) > fabs(y_axis.min)) //both +ve
				{
					originY = plot.y + plot.height;
				}
				else //both -ve.
				{
					originY = plot.y;
				}
				
				return cv::Point(int(originX), int(originY));
			}

			/*	This version of compute_data interpolates (linear) the data to sit on pixel values in the background matrix
			Use when your function is singular valued i.e. there exists a unique y(x) at all x.*/
			coords compute_data_1(const stdVec_d& x_vals,
				const stdVec_d& y_vals,
				const cv::Point& origin, //pixel location of axis origin
				const std::pair<double, double>& data_origin, //values of the data at the origin
				double m_scale_x,
				double m_scale_y)
			{
				using namespace phys::maths;
				stdVec_i x_axis, y_axis;
				stdVec_d x_interp, y_interp;

				double origin_val_x = data_origin.first;
				double origin_val_y = data_origin.second;

				int x_new, x_old = int(round(m_scale_x*(x_vals[0] - origin_val_x)));
				x_interp.push_back(double(x_old / m_scale_x + origin_val_x));
				x_axis.push_back(x_old + origin.x);
				for (size_t i = 1; i < x_vals.size(); ++i)
				{
					x_new = int(round(m_scale_x * (x_vals[i] - origin_val_x)));
					if (x_new != x_old){
						x_interp.push_back(double(x_new) / m_scale_x + origin_val_x);
						x_axis.push_back(x_new + origin.x);
					}
					x_old = x_new;
				}

				phys::interp::data plot_data(x_vals, y_vals);
				bool sort = true;

				phys::interp::Linear linear_interp(plot_data, sort);

				for (size_t i = 0; i < x_interp.size(); ++i)
					y_interp.push_back(linear_interp.interpolate(x_interp[i]));

				for (size_t i = 0; i < x_interp.size(); ++i){
					int y_temp = origin.y - int(round(m_scale_y*(y_interp[i] - origin_val_y)));
					y_axis.push_back(y_temp);
				}

				assert(x_axis.size() == y_axis.size());

				coords retval;

				for (size_t i = 0; i < x_axis.size(); ++i){
					retval.push_back(std::pair<int, int>(x_axis[i], y_axis[i]));
				}
				return retval;
			}

			/*	This version of compute_data rounds the floating point data into ints to plot on the background matrix.
			Depending on the resolution of the background and your data this may lead to small artifical discontinuities in the plot.
			Use when the function you are ploting is non-sigular valued i.e. there exists multiple y(x) for a given x.
			For example, the phase space plot of an oscillator, or orbital motion, say. */
			coords compute_data_2(const stdVec_d& x_vals,
				const stdVec_d& y_vals,
				const cv::Point& origin, //pixel location of axis origin
				const std::pair<double, double>& data_origin, //data values at origin
				double m_scale_x,
				double m_scale_y)
			{
				using namespace phys::maths;
				stdVec_i x_axis, y_axis;

				double origin_val_x = data_origin.first;
				double origin_val_y = data_origin.second;

				int x_new, x_old = int(round(m_scale_x * (x_vals[0] - origin_val_x)));
				int y_new, y_old = int(round(m_scale_y * (y_vals[0] - origin_val_y)));
				x_axis.push_back(x_old + origin.x);
				y_axis.push_back(origin.y - y_old);
				
				for (size_t i = 1; i < x_vals.size(); ++i)
				{
					x_new = int(round(m_scale_x * (x_vals[i] - origin_val_x)));
					y_new = int(round(m_scale_y * (y_vals[i] - origin_val_y)));
					if (x_new != x_old){
						x_axis.push_back(x_new + origin.x);
						y_axis.push_back(origin.y - y_new);
					}
					else if (y_new != y_old)
					{
						x_axis.push_back(x_new + origin.x);
						y_axis.push_back(origin.y - y_new);
					}
					x_old = x_new;
					y_old = y_new;
				}

				assert(x_axis.size() == y_axis.size());

				coords retval;

				for (size_t i = 0; i < x_axis.size(); ++i){
					retval.push_back(std::pair<int, int>(x_axis[i], y_axis[i]));
				}

				return retval;
			}
		}

		//Viewer functions	
		//--------------------------------------------------------------------------------------------------------
		// ODE storage argument - makes a choice then calls usual plot function
		void Viewer::plot(const phys::storage::ODE_Storage& data,
			graph_choice choice,
			size_t dim1,
			size_t dim2)
		{
			if (data.get_independent().size() < 2){
				throw std::runtime_error("Not enough data to make a plot - check your code\n");
			}

			stdVec_d x_var;
			stdVec_d y_var;
			m_g_choice = choice;
			switch (choice)
			{
			case DEPEN://0
				m_x_axis.title = data.get_x_name();
				m_y_axis.title = data.get_y_name();
				x_var = data.get_independent();
				y_var = data.get_dependent(dim1);
				break;
			case DERIV://1
				m_x_axis.title = data.get_x_name();
				m_y_axis.title = data.get_dy_name();
				x_var = data.get_independent();
				y_var = data.get_first_deriv(dim1);
				break;
			case XY://2
				m_x_axis.title = data.get_y_name() + "_x";
				m_y_axis.title = data.get_y_name() + "_y";
				x_var = data.get_dependent(dim1);
				y_var = data.get_dependent(dim2);
				break;
			case PHASE://3
				m_x_axis.title = data.get_y_name();
				m_y_axis.title = data.get_dy_name();
				x_var = data.get_dependent(dim1);
				y_var = data.get_first_deriv(dim1);
				break;
			}

			plot(x_var, y_var);
		}

		//General Storage<double> argument -- grabs data and calls relevant plot function
		void Viewer::plot(const phys::storage::Storage<double>& data, bool splitAxes)
		{
			m_x_axis.title = data.get_x_name();

			if (data.get_multi().empty())
			{
				m_y_axis.title = data.get_y_name();
				plot(data.get_independent(), data.get_dependent());
			}
			else if (splitAxes)
			{
				m_y_axis.title = data.get_name(0);
				m_y_axis2.title = data.get_name(1);
				plot_split(data.get_independent(), data.get_multi());
			}
			else
			{
				m_key_titles = data.get_key_names();
				plot(data.get_independent(), data.get_multi());
			}
		}


		void Viewer::plot_split(const stdVec_d& x, const std::vector<stdVec_d>& y)
		{
			assert(y.size() == 2);

			plot(x, y[0]);

			processAxis(y[1], Y2);

			do{
				//m_clr_counter++;
				++m_clr_counter %= colour_size;
				m_dt_colour = colour[BLUE + m_clr_counter];
			} while (m_dt_colour == m_bg_colour);

			cv::Point temp_origin = m_axis_origin;
			temp_origin.x += m_plot_area.width;

			double val_space_y = (m_y_axis2.max - m_y_axis2.min) / 10.0;
			draw_labels(Y2, m_val_origin_y, val_space_y, m_y_axis2.scale, temp_origin);
			drawYaxis(temp_origin.x, false); //false == draw axis title on the right

			coords data;

			std::pair<double, double> data_origin(m_val_origin_x, m_val_origin_y);

			if (m_g_choice < 2)
				data = compute_data_1(x, y[1], m_axis_origin, data_origin, m_x_axis.scale, m_y_axis2.scale);
			else
				data = compute_data_2(x, y[1], m_axis_origin, data_origin, m_x_axis.scale, m_y_axis2.scale);

			if (m_using_range_x || m_using_range_y) chk_data(data);

			if (m_animate)
			{
				draw_key(1, m_dt_colour);
				plot_ani_impl(data);
				m_BG = cv::Mat::zeros(m_rows, m_cols, CV_8UC3); //clear display matrix
			}
			else
			{
				draw_key(1, m_dt_colour);
				plot_impl(data);
				m_BG = cv::Mat::zeros(m_rows, m_cols, CV_8UC3);//clear display matrix
			}

		}

		/* Typical plot function - vector of doubles x vs. vector of doubles y(x) */
		void Viewer::plot(const stdVec_d& x, const stdVec_d& y)
		{
			assert(x.size() == y.size());

			m_is_drawn = true;
		
			processAxis(x, X);
			processAxis(y, Y);
			
			m_axis_origin = compute_origin(m_plot_area, m_x_axis, m_y_axis);
			
			double val_space_x = (m_x_axis.max - m_x_axis.min) / 10.0;
			double val_space_y = (m_y_axis.max - m_y_axis.min) / 10.0;

			draw_labels(X, m_val_origin_x, val_space_x, m_x_axis.scale, m_axis_origin);
			draw_labels(Y, m_val_origin_y, val_space_y, m_y_axis.scale, m_axis_origin);

			draw_axes(m_axis_origin);

			coords data;

			std::pair<double, double> data_origin(m_val_origin_x, m_val_origin_y);

			if (m_g_choice < 2)
				data = compute_data_1(x, y, m_axis_origin, data_origin, m_x_axis.scale, m_y_axis.scale);
			else
				data = compute_data_2(x, y, m_axis_origin, data_origin, m_x_axis.scale, m_y_axis.scale);

			if (m_using_range_x || m_using_range_y) chk_data(data);

			cv::namedWindow(m_plot_name, CV_WINDOW_AUTOSIZE);

			if (m_animate)
			{
				draw_key(0, m_dt_colour);
				plot_ani_impl(data);
				//m_BG = cv::Mat::zeros(m_rows, m_cols, CV_8UC3); //clears the display matrix
			}
			else
			{
				draw_key(0, m_dt_colour);
				plot_impl(data);
				//m_BG = cv::Mat::zeros(m_rows, m_cols, CV_8UC3); //clears the display matrix
			}
		}


		void Viewer::plot(const stdVec_d& x, const std::vector<stdVec_d>& y)
		{
			for (size_t i = 0; i < y.size(); ++i)
				assert(x.size() == y[i].size());

			m_is_drawn = true;

			processAxis(x, X);
			processAxis(y);

			m_axis_origin = compute_origin(m_plot_area, m_x_axis, m_y_axis);
			draw_axes(m_axis_origin);

			double val_space_x = (m_x_axis.max - m_x_axis.min) / 10.0;
			double val_space_y = (m_y_axis.max - m_y_axis.min) / 10.0;

			draw_labels(X, m_val_origin_x, val_space_x, m_x_axis.scale, m_axis_origin);
			draw_labels(Y, m_val_origin_y, val_space_y, m_y_axis.scale, m_axis_origin);

			cv::namedWindow(m_plot_name, CV_WINDOW_AUTOSIZE);

			std::pair<double, double> data_origin(m_val_origin_x, m_val_origin_y);

			for (size_t i = 0; i < y.size(); ++i)
			{
				do{
					m_clr_counter = m_clr_counter % colour_size;
					m_dt_colour = colour[BLUE + m_clr_counter++];
				} while (m_dt_colour == m_bg_colour);

				coords data;

				if (m_g_choice < 2)
					data = compute_data_1(x, y[i], m_axis_origin, data_origin, m_x_axis.scale, m_y_axis.scale);
				else
					data = compute_data_2(x, y[i], m_axis_origin, data_origin, m_x_axis.scale, m_y_axis.scale);


				if (m_using_range_x || m_using_range_y) chk_data(data);

				if (m_animate){
					draw_key(i, m_dt_colour);
					plot_ani_impl(data);
				}
				else{
					draw_key(i, m_dt_colour);
					plot_impl(data);
				}
			}

			return;
		}

		void Viewer::draw_key(size_t which, const cv::Scalar& colour)
		{
			if (m_key_titles.empty()) return; //Do nothing.
			int offset = m_rows / 25;
			//compute location of where to put key.
			cv::Point key_loc = cv::Point(m_cols - m_cols / 10, m_rows / 10 + offset*which);
			//draw it
			cv::putText(m_BG, m_key_titles[which], key_loc, CV_FONT_HERSHEY_PLAIN, m_font_scale*1.2, colour);
		}


		void Viewer::plot_impl(const coords& data)
		{
			for (ccoords_it it = data.begin(); it != data.end(); ++it)
				m_pShape->draw(this, cv::Point(it->first, it->second));

			if (m_lines){
				for (ccoords_it it = data.begin(); it != data.end() - 1; ++it)
				{
					if (it->first > 0 && (it + 1)->first > 0)
						cv::line(m_BG, cv::Point(it->first, it->second),
						cv::Point((it + 1)->first, (it + 1)->second), m_dt_colour);
				}
			}
			cv::imshow(m_plot_name, m_BG);
			cv::waitKey(m_pause);
			return;
		}

		void Viewer::plot_ani_impl(const coords& data, int delay)
		{
			cv::putText(m_BG, "Press Esc to continue...", cv::Point(10, 10), CV_FONT_HERSHEY_PLAIN, m_font_scale*0.7, m_dt_colour);
			cv::Mat fixed = m_BG.clone();
			cv::namedWindow(m_plot_name, CV_WINDOW_KEEPRATIO);
			char key = 1;
			do{
				ccoords_it it = data.begin();
				while (it != data.end() && key != 27) //Esc 
				{
					cv::imshow(m_plot_name, m_BG);
					m_pShape->draw(this, cv::Point(it->first, it->second));
					++it;
					key = cv::waitKey(delay);
				}
				m_BG = fixed.clone();
			} while (key != 27);

			//redraw data in case mulitple data being shown.
			for (size_t i = 0; i < data.size(); ++i)
				m_pShape->draw(this, cv::Point(data[i].first, data[i].second));

			return;
		}

		//user to compute the gradient and y intercept
		void Viewer::draw_line(double m, double c)
		{
			using namespace phys::maths;
			cv::Point temp, start(-1, -1), end(-1, -1);

			//compute intercept point
			int intercept = int(round(m_y_axis.scale * c));
			temp = cv::Point(m_axis_origin.x, m_axis_origin.y - intercept);

			/***** We need to maintain dx and dy parity throughout ******/

			// +ve half of the X axis first - dx will be zero if axis is -ve only.			
			int dx = m_plot_area.x + m_plot_area.width - m_axis_origin.x; //dx here will be >= 0
			int dy = int(round(m_y_axis.scale * dx * m / m_x_axis.scale));

			if (dx) //i.e. we have a +ve x axis.
			{
				if (m > 0) //+ve gradient -> +ve dy
				{
					if (temp.y - dy < m_plot_area.y)
					{
						dy = temp.y - m_plot_area.y; //dy +ve
						dx = int(round(dy*m_x_axis.scale / m_y_axis.scale / m)); //dx +ve
					}
				}
				else // -ve gradient -> -ve dy
				{
					if (temp.y - dy > m_plot_area.y + m_plot_area.height)
					{
						dy = -(m_plot_area.y + m_plot_area.height - temp.y); //maintain dy parity
						dx = int(round(dy*m_x_axis.scale / m_y_axis.scale / m)); // dx +ve
					}
				}

				end = temp + cv::Point(dx, -dy);
			}

			//-ve half of the X-axis - dx will be zero if +ve axis only
			dx = m_plot_area.x - m_axis_origin.x; // dx here will be <= 0
			dy = int(round(m_y_axis.scale * dx * m / m_x_axis.scale));

			if (dx) //i.e. we have a -ve x axis.
			{
				if (m > 0) //+ve gradient -> -ve dy
				{
					if (temp.y - dy > m_plot_area.y + m_plot_area.height)
					{
						dy = -(m_plot_area.y + m_plot_area.height - temp.y); //maintain dy parity
						dx = int(round(dy*m_x_axis.scale / m_y_axis.scale / m)); //dx -ve
					}
				}
				else // -ve gradient -> +ve dy
				{
					if (temp.y - dy < m_plot_area.y)
					{
						dy = temp.y - m_plot_area.y; // dy +ve
						dx = int(round(dy*m_x_axis.scale / m_y_axis.scale / m)); //dx -ve
					}
				}

				start = temp + cv::Point(dx, -dy);
			}

			//if either start or end not set then assign temp (if set coordinates are _ALWAYS_ positive)
			if (start.x < 0) start = temp;
			if (end.x   < 0) end = temp;

			//set the colour for the line; 				
			do{
				//ensure we wrap back to the start of the list
				++m_clr_counter %= colour_size; //change the colour of the line from data already plotted.
				m_dt_colour = colour[BLUE + m_clr_counter];
			} while (m_dt_colour == m_bg_colour); //avoid plotting data same colour as background

			//draw the line
			cv::line(m_BG, start, end, m_dt_colour);

			cv::imshow(m_plot_name, m_BG);
			cv::waitKey();
		}

		void Viewer::drawXaxis(int xo)
		{
			//compute axis start and end points ...
			cv::Point xstart(m_plot_area.x, xo);
			cv::Point xend(m_plot_area.x + m_plot_area.width, xo);
			//... and draw it
			cv::line(m_BG, xstart, xend, m_ax_colour, 2);
			//compute axis title position ...
			cv::Point x_pos(int(m_cols / 2), m_plot_area.y + m_plot_area.height + int(m_rows*m_offset_y / 150));
			//... and draw it
			cv::putText(m_BG, m_x_axis.title, x_pos, CV_FONT_HERSHEY_PLAIN, 1.2*m_font_scale, m_dt_colour);
		}

		void Viewer::drawYaxis(int yo, bool onLeft)
		{
			cv::Point ystart(yo, m_plot_area.y + m_plot_area.height);
			cv::Point yend(yo, m_plot_area.y);
			cv::line(m_BG, ystart, yend, m_ax_colour, 2);

			//choice based on position of axis
			int title_x = onLeft ? int(m_cols * m_offset_x / 100 / 2) : m_plot_area.x + m_plot_area.width + int(m_cols * m_offset_x / 100 / 2);

			cv::Point y_pos(title_x, int(m_rows*m_offset_y / 100 / 2));
			cv::putText(m_BG, onLeft ? m_y_axis.title : m_y_axis2.title, y_pos, CV_FONT_HERSHEY_PLAIN, 1.2*m_font_scale, m_dt_colour);
		}

		void Viewer::draw_axes(const cv::Point& origin)
		{
			drawXaxis(origin.y);
			drawYaxis(origin.x);
		}

		void Viewer::draw_labels(	axis AXIS,
									double val_o,
									double val_s,
									double scale,
									const cv::Point& axis_o)
		{
			const cv::Scalar label_clr = m_dt_colour;
			const int tick_space = int(val_s * scale);
			double disp_val = val_o;
			std::ostringstream convert;
			cv::Point tsize;
			int axis_end;
			int tick_loc;
			cv::Point tloc;
			cv::Point label_loc;

			const int X_offset = int(m_rows * 5 / 100);
			int Y_offset = int(m_cols*m_offset_x / 100 / 2);
			if (AXIS == Y2)
				Y_offset = m_cols - int(m_cols * m_offset_x / 100 / 2);

			
			//Note that we start drawing the ticks at the graph origin hence the need to have seperate while loops
			//to cope with potentialy both postive and negative halves of the axes. 
			//- is there a better way of doing this as we have a significant amount of repeated code? - there may not be.
			//There are subtle differences between the X and Y axis tick locations hence the need to switch
			switch (AXIS)
			{
			case X:
				// first do +ve half of axis; it will place origin even if no +ve values in data.
				tsize = cv::Point(0, 10);
				axis_end = m_plot_area.x + m_plot_area.width;
				tick_loc = axis_o.x;
				do{
					tloc = cv::Point(tick_loc, axis_o.y);
					label_loc = cv::Point(tick_loc - 10, axis_o.y + X_offset);
					cv::line(m_BG, tloc, tloc + tsize, m_ax_colour);
					convert << format(disp_val, 3, 1.e-2, 1.e4);
					cv::putText(m_BG, convert.str(), label_loc, CV_FONT_HERSHEY_PLAIN, .8*m_font_scale, label_clr);
					convert.str(""); //clear stringstream
					tick_loc += tick_space;
					disp_val += val_s;
				} while (tick_loc <= axis_end);

				//second do -ve half of axis if required
				axis_end = m_plot_area.x;
				tick_loc = axis_o.x - tick_space; //already drawn tick and label at origin
				disp_val = val_o - val_s;

				while (tick_loc >= axis_end){
					tloc = cv::Point(tick_loc, axis_o.y);
					label_loc = cv::Point(tick_loc - 10, axis_o.y + X_offset);
					cv::line(m_BG, tloc, tloc + tsize, m_ax_colour);
					convert << format(disp_val, 3, 1.e-2, 1.e4);
					cv::putText(m_BG, convert.str(), label_loc, CV_FONT_HERSHEY_PLAIN, .8*m_font_scale, label_clr);
					convert.str("");
					tick_loc -= tick_space;
					disp_val -= val_s;
				}
				break;
			case Y: //fall through by design
			case Y2:
				//first +ve half of axis; it will place origin even if no +ve values in data.
				tsize = (AXIS == Y) ? cv::Point(-10, 0) : cv::Point(10, 0);
				axis_end = m_plot_area.y;
				tick_loc = axis_o.y;
				do{
					tloc = cv::Point(axis_o.x, tick_loc);
					label_loc = cv::Point(Y_offset, tick_loc);
					cv::line(m_BG, tloc, tloc + tsize, m_ax_colour);
					convert << format(disp_val, 3, 1.e-2, 1.e4);
					cv::putText(m_BG, convert.str(), label_loc, CV_FONT_HERSHEY_PLAIN, .8*m_font_scale, label_clr);
					convert.str("");
					tick_loc -= tick_space;
					disp_val += val_s;
				} while (tick_loc >= axis_end);

				//second -ve half of axis if required
				axis_end = m_plot_area.y + m_plot_area.height;
				tick_loc = axis_o.y + tick_space;
				disp_val = val_o - val_s;

				while (tick_loc <= axis_end)
				{
					tloc = cv::Point(axis_o.x, tick_loc);
					label_loc = cv::Point(Y_offset, tick_loc);
					cv::line(m_BG, tloc, tloc + tsize, m_ax_colour);
					convert << format(disp_val, 3, 1.e-2, 1.e4);
					cv::putText(m_BG, convert.str(), label_loc, CV_FONT_HERSHEY_PLAIN, .8*m_font_scale, label_clr);
					convert.str("");
					tick_loc += tick_space;
					disp_val -= val_s;
				}
				break;
			}
		}

		void Viewer::add_data(const phys::storage::ODE_Storage& data,
			size_t dim1,
			size_t dim2)
		{
			switch (m_g_choice)
			{
			case DEPEN:
				add_data(data.get_independent(), data.get_dependent(dim1), true);
				break;
			case DERIV:
				add_data(data.get_independent(), data.get_first_deriv(dim1), true);
				break;
			case XY:
				add_data(data.get_dependent(dim1), data.get_dependent(dim2), false);
				break;
			case PHASE:
				add_data(data.get_dependent(dim1), data.get_first_deriv(dim1), false);
				break;
			}
		}


		void Viewer::add_data(const phys::storage::Storage<double>& data,
			bool is_data_singular)
		{
			add_data(data.get_independent(), data.get_dependent(), is_data_singular);
		}

		void Viewer::add_data(const stdVec_d& x,
			const std::vector<stdVec_d>& y,
			bool is_data_singular)
		{
			for (size_t i = 0; i < y.size(); ++i){
				draw_key(i, colour[BLUE + m_clr_counter]);
				add_data(x, y[i], is_data_singular);
			}
		}


		void Viewer::add_data(const stdVec_d& x,
			const stdVec_d& y,
			bool is_data_singular)
		{
			if (!m_is_drawn)
				throw(std::runtime_error("No plot drawn on which to add data"));

			do{
				m_clr_counter = m_clr_counter % colour_size;
				m_dt_colour = colour[BLUE + m_clr_counter++];
			} while (m_dt_colour == m_bg_colour);

			coords data;
			std::pair<double, double> data_origin(m_val_origin_x, m_val_origin_y);
			if (is_data_singular)
				data = compute_data_1(	x,
										y,
										m_axis_origin,
										data_origin,
										m_x_axis.scale,
										m_y_axis.scale);
			else
				data = compute_data_2(	x,
										y,
										m_axis_origin,
										data_origin,
										m_x_axis.scale,
										m_y_axis.scale);
			chk_data(data);
			if (m_animate)
				plot_ani_impl(data);
			else
				plot_impl(data);
		}

		void Viewer::chk_data(coords& data)
		{
			size_t count = 0;

			//lft_bound always assigned in loop.
			coords_it lft_bound, rht_bound = data.begin();

			//find coordinates outside of the plot area and set thier x values to some arbritrary negative int.
			//store iterators relating to were x values are within the plot bounds to erase those values if needed
			for (coords_it it = data.begin(); it != data.end(); ++it)
			{
				//check x is within bounds
				if (it->first >= m_plot_area.x && it->first <= m_plot_area.x + m_plot_area.width){
					++count;
					if (count == 1) lft_bound = it; //first value within x-bounds
					//if x in bounds check to see if y out of bounds.
					if (it->second < m_plot_area.y || it->second > m_plot_area.y + m_plot_area.height)
					{
						it->first = -10;
					}//else y in bounds; do nothing
				} else { //x out of bounds.
					it->first = -10;
					if (count > 0){
						rht_bound = it; //one past last value within x-bounds
						count = 0;
					}
				}
			}

			//erase data that is out of x bounds for x vs. y(x) type graphs only.
			if (m_g_choice < 2 && rht_bound != data.begin() && lft_bound != data.begin()){
				//do right bound values first to avoid iterator reallocation then erase left bound values
				data.erase(rht_bound, data.end());
				data.erase(data.begin(), lft_bound - 1);
			}//else do nothing -- data that is out-of-bounds is still "plotted" but won't be seen.
		}

		void Viewer::set_shape(shape_choice sc, size_t sz)
		{
			//delete the pointee of the m_pShape pointer first.
			delete m_pShape;
			//reassign m_pShape pointer to the m_pShape wanted
			switch (sc){
			case CIRCLE:
				m_pShape = new Circle(sz, m_dt_colour);
				break;
			case SQUARE:
				m_pShape = new Square(sz, m_dt_colour);
				break;
			case TRIANGLE:
				m_pShape = new Triangle(sz, m_dt_colour);
				break;
			case CROSS:
				m_pShape = new Cross(sz, m_dt_colour);
				break;
			case ARROW:
				m_pShape = new Arrow(sz, sz / 4, 0.0, m_dt_colour);
				break;
			default:
				std::cerr << "Shape not recognised. Using a circle.\n"; 
				m_pShape = new Circle(sz, m_dt_colour);
				break;
			}
		}

		//----------------------------------------------------------------------------------------
		//	Shape functions
		//----------------------------------------------------------------------------------------

		void Circle::draw(Viewer* v, const cv::Point& loc)
		{
			cv::circle(v->m_BG, loc, m_hull_x, v->m_dt_colour);
		}

		void Square::draw(Viewer* v, const cv::Point& loc)
		{
			cv::Rect r_loc(loc.x - m_hull_x / 2, loc.y - m_hull_x / 2, m_hull_x, m_hull_x);
			cv::rectangle(v->m_BG, r_loc, v->m_dt_colour);
		}

		void Triangle::draw(Viewer* v, const cv::Point& loc)
		{
			//compute the three vertices from central location
			cv::Point top(loc.x, loc.y - m_hull_y / 2);
			cv::Point lft(loc.x - m_hull_x / 2, loc.y + m_hull_y / 2);
			cv::Point rht(loc.x + m_hull_x / 2, loc.y + m_hull_y / 2);

			//draw the m_lines top-left, top-right, and left-right.
			cv::line(v->m_BG, top, lft, v->m_dt_colour);
			cv::line(v->m_BG, top, rht, v->m_dt_colour);
			cv::line(v->m_BG, lft, rht, v->m_dt_colour);
		}

		void Cross::draw(Viewer* v, const cv::Point& loc)
		{
			//four points (corners of the hull), two m_lines
			cv::Point top_lft(loc.x - m_hull_x / 2, loc.y - m_hull_y / 2);
			cv::Point top_rht(loc.x + m_hull_x / 2, loc.y - m_hull_y / 2);
			cv::Point btm_lft(loc.x - m_hull_x / 2, loc.y + m_hull_y / 2);
			cv::Point btm_rht(loc.x + m_hull_x / 2, loc.y + m_hull_y / 2);

			//draw m_lines between opposite corners
			cv::line(v->m_BG, top_lft, btm_rht, v->m_dt_colour);
			cv::line(v->m_BG, top_rht, btm_lft, v->m_dt_colour);
		}

		void Arrow::draw(Viewer* v, const cv::Point& loc)
		{
			// define local pi
			const double PI = 3.141592653;

			//compute start and end point of the principal line
			//m_hull_x == arrow length
			double dx = m_hull_x * cos(m_angle*PI / 180) / 2.0;
			double dy = m_hull_x * sin(m_angle*PI / 180) / 2.0;
			cv::Point start(loc.x - int(dx), loc.y + int(dy));
			cv::Point end(loc.x + int(dx), loc.y - int(dy));

			//Draw the principal line
			cv::line(v->m_BG, start, end, v->m_dt_colour);

			//Compute the coordinates of the first line segment
			//m_hull_y == line segment size
			start.x = end.x - int(m_hull_y * cos(m_angle*PI / 180 + PI / 4));
			start.y = end.y - int(m_hull_y * sin(m_angle*PI / 180 + PI / 4));

			//draw line segment
			cv::line(v->m_BG, start, end, v->m_dt_colour);

			//Compute the coordinates of the second line segment
			start.x = end.x - int(m_hull_y * cos(m_angle*PI / 180 - PI / 4));
			start.y = end.y - int(m_hull_y * sin(m_angle*PI / 180 - PI / 4));

			//draw line segment
			cv::line(v->m_BG, start, end, v->m_dt_colour);
		}

	}
}
