#ifndef VISUALISE_HPP
#define VISUALISE_HPP

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include "Interpolation.h"
#include "Storage.h"
#include "Maths.h"

#include <exception> 
#include <memory>

/*
 Common OpenCV types and explanations:

 cv::Mat		--	General array class; can be used to display images in .png, .jpg, .bmp, and so on.
 Each element of the matrix is a pixel and they have dimensions of m_rows, columns,
 and channels (if the matrix is a leaf in a book then the number of channels defines
 the number of leaves). To represent a colour we need three channels; red, green,
 and blue (other colour schemes do exist). A greyscale image only has one channel.

 cv::Scalar	--	Class to define colour (amongst other things). The constructor can take three integer
 parameters between the values of 0 and 255 inclusively, representing the three colour
 channels blue, green, and red respectively (yes reverse order to RGB). Note that
 the max value is not arbitrary but is equal to 2^8 - 1, i.e. an 8-bit unsigned
 integer. Hence, for a greyscale image (one channel) there are actually 256 shades of
 grey -- someone please use that as a title for a book!

 cv::Point	--	Class defining an (x,y) coordinate of a pixel in a cv::Mat. It has members x and y
 that can be accessed directly.

 cv::Rect	--	Class defining a rectangle in terms of pixel coordinates. Can be constructed using
 the top-left (x,y) and width and height parameters, or can be specified by the top-
 left and bottom-right cv::Points. Members can be access directly. Note this is NOT
 a drawing function; it is typically used to specify a region of interest (ROI) in
 an image. To draw a rectangle you have to use the cv::rectangle function.
 */

namespace phys {

using uint = unsigned int;

/* typedefs for shorthand */
typedef std::vector<std::pair<int, int>> coords;
typedef coords::const_iterator ccoords_it;
typedef coords::iterator coords_it;

/*	Enum list for choice of data plot m_shape		*/
enum shape_choice {
	CIRCLE, SQUARE, TRIANGLE, CROSS, ARROW
};

/*	Forward declaration for the Viewer class		*/
class Viewer;

/*	Virtual base class Shape	*/
class Shape {
protected:
	/* Constructor */
	Shape(uint h_x, uint h_y, cv::Scalar clr) :
			m_hull_x(h_x), m_hull_y(h_y), m_shape_colour(clr) {
	}
public:
	virtual ~Shape() {
	}
public:
	virtual void draw(Viewer *v, const cv::Point &loc) = 0;
	uint hull_area() const {
		return m_hull_x * m_hull_y;
	}

protected:
	uint m_hull_x;				//!< (convex) hull x size
	uint m_hull_y;				//!< (convex) hull y size
	cv::Scalar m_shape_colour;//!< m_shape colour e.g. cv::Scalar(255,0,0) is blue.
};

/*	Derived classes of Shape	*/
class Circle: public Shape {
public:
	Circle(uint radius, cv::Scalar clr) :
			Shape(radius, radius, clr) {
	}
	void draw(Viewer *v, const cv::Point &loc);
};

class Square: public Shape {
public:
	Square(uint side, cv::Scalar clr) :
			Shape(side, side, clr) {
	}
	void draw(Viewer *v, const cv::Point &loc);
};

class Triangle: public Shape {
public:
	Triangle(uint height, cv::Scalar clr) :
			Shape(height, height, clr) {
	}
	void draw(Viewer *v, const cv::Point &loc);
};

class Cross: public Shape {
public:
	Cross(uint length, cv::Scalar clr) :
			Shape(length, length, clr) {
	}
	void draw(Viewer *v, const cv::Point &loc);
};

class Arrow: public Shape {
public:
	Arrow(uint length, uint head, double ang, cv::Scalar clr) :
			Shape(length, head, clr), m_angle(ang) {
		m_hull_y = (m_hull_y == 0) ? 2 : m_hull_y;
	}
	void draw(Viewer *v, const cv::Point &loc);
private:
	double m_angle;
};

//structure to help clean up code.
struct axis_props {
	double scale;
	double min;
	double max;
	double origin;
	std::string title;

	axis_props(double scl = 0., double min = 0., double max = 0., double ori =
			0., const std::string &name = "axis") :
			scale(scl), min(min), max(max), origin(ori), title(name) {
	}
};

/*
 *  The Viewer class is a bit of a beast that (probably) indicates it should be redesigned.
 *  Things like the graph choice being coded as an enumeration type might be better off as
 *  a "Graph" class with a hierarchy and polymorphic behaviour; an "Axis" class that is
 *  responsible for its own scale, range, tick marks, labels and so on.
 *
 *  It should only be used to give the user a cursory view of the data and not be used
 *  as a precise plotting tool; there are many professional software suites that
 *  are more than adequate at plotting configurable graphs.
 */

/* Viewer class */
class Viewer {
	/* Allows shape sub-classes to access private Viewer members */
	friend class Circle;
	friend class Square;
	friend class Triangle;
	friend class Cross;
	friend class Arrow;
	/* c++11 using typedef for convenience */
	using UPtrShape = std::unique_ptr<Shape>;

private:
	int m_rows;								//!< height of the Viewer window
	int m_cols;								//!< width of the Viewer window
	double m_offset_x;	//!< percentage of the Viewer window width as border
	double m_offset_y;	//!< percentage of the Viewer window height as border
	axis_props m_x_axis;					//!< x axis properties
	axis_props m_y_axis;					//!< y axis properties
	axis_props m_y_axis2;					//!< second y axis properties
	std::pair<double, double> m_x_range;	//!< optional range in the x axis
	std::pair<double, double> m_y_range;	//!< optional range in the y_axis
	bool m_using_range_x;				//!< using a range on the x axis or not
	bool m_using_range_y;				//!< using a range on the y axis or not
	cv::Rect m_plot_area;//!< plot area (i.e. Viewer window minus the borders).
	cv::Point m_axis_origin;				//!< pixel location of the origin
	cv::Scalar m_bg_colour;					//!< background colour
	cv::Scalar m_ax_colour;					//!< axis colour
	cv::Scalar m_dt_colour;					//!< data colour
	uint m_clr_counter; //!< counter to change data colour if multiple data being plotted
	double m_val_origin_x;//!< optional origin x value (in data units not pixels)
	double m_val_origin_y;//!< optional origin y value (in data units not pixels)
	std::vector<std::string> m_key_titles;//!< key titles if plotting multiple data
	std::string m_plot_name;				//!< Viewer window name
	int m_g_choice;	//!< choice between x vs. y(x), x vs. dy/dx, y1(x) vs. y2(x), y(x) vs. dy/dx
	int m_pause;//!< optional delay makes plot pause; default zero: Viewer hangs until key pressed.
	bool m_lines;						//!< do you want lines with your points?
	bool m_pulses;					//!< pulse lines from x-axis to data point
	bool m_animate;			//!< switches between a static and animated plot.
	bool m_is_drawn;//!< check that a plot has been plotted before trying to add data.
	double m_font_scale;			//!< overall font scale for labels drawn.
	cv::Mat m_BG;//!< opencv matrix - this is the array to which we draw everything.
	UPtrShape m_pShape;	//!< m_shape pointer - choice between circle, square, triangle and cross.
	enum axis {
		X, Y, Y2
	};					//!< enum to select between different axes
public:
	Viewer() :
			m_rows(450), m_cols(800), m_offset_x(15.0), m_offset_y(15.0), m_x_axis(), m_y_axis(), m_y_axis2(), m_x_range(
					0.0, 0.0), m_y_range(0.0, 0.0), m_using_range_x(false), m_using_range_y(
					false), m_plot_area(), m_axis_origin(
					cv::Point(m_rows / 2, m_cols / 2)), m_bg_colour(
					cv::Scalar(0, 0, 0)), m_ax_colour(cv::Scalar(255, 0, 0)), m_dt_colour(
					cv::Scalar(0, 255, 0)), m_clr_counter(1), //default to green
			m_val_origin_x(0.0), m_val_origin_y(0.0), m_key_titles(), m_plot_name(
					"Plot"), m_g_choice(0), m_pause(0), m_lines(false), m_pulses(
					false), m_animate(false), m_is_drawn(false), m_font_scale(
					1.0), m_BG(cv::Mat(m_rows, m_cols, CV_8UC3, m_bg_colour)), m_pShape(
					UPtrShape { new Circle { 1, m_dt_colour } }) {
		cv::Point tl(int(m_cols * m_offset_x / 100),
				int(m_rows * m_offset_y / 100));
		cv::Point br(m_cols - int(m_cols * m_offset_x / 100),
				m_rows - int(m_rows * m_offset_y / 100));
		m_plot_area = cv::Rect(tl, br);
		//m_font_scale = m_rows*m_cols / 1000 / 1000;
	}

	Viewer(const int width, const int height, const cv::Scalar &bg_col =
			cv::Scalar(0, 0, 0),
			const cv::Scalar &ax_col = cv::Scalar(255, 0, 0),
			const cv::Scalar &dt_col = cv::Scalar(0, 255, 0),
			const std::string title = "Plot", const std::string &x_name = "x",
			const std::string &y_name = "y", bool lines_wanted = false,
			bool pulses_wanted = false, bool ani = false,
			const std::vector<std::string> key_names =
					std::vector<std::string>()) :
			m_rows(height), m_cols(width), m_offset_x(15.0), m_offset_y(15.0), m_x_range(
					0.0, 0.0), m_y_range(0.0, 0.0), m_using_range_x(false), m_using_range_y(
					false), m_plot_area(), m_axis_origin(
					cv::Point(m_rows / 2, m_cols / 2)), m_bg_colour(bg_col), m_ax_colour(
					ax_col), m_dt_colour(dt_col), m_clr_counter(1), //default to green
			m_val_origin_x(0.0), m_val_origin_y(0.0), m_key_titles(key_names), m_plot_name(
					title), m_g_choice(0), m_pause(0), m_lines(lines_wanted), m_pulses(
					pulses_wanted), m_animate(ani), m_is_drawn(false), m_font_scale(
					m_rows * m_cols / 1000. / 1000.), m_BG(
					cv::Mat(m_rows, m_cols, CV_8UC3, m_bg_colour)), m_pShape(
					UPtrShape { new Circle { 1, m_dt_colour } }) {
		cv::Point tl(int(m_cols * m_offset_x / 100.),
				int(m_rows * m_offset_y / 100.));
		cv::Point br(m_cols - int(m_cols * m_offset_x / 100.),
				m_rows - int(m_rows * m_offset_y / 100.));
		m_plot_area = cv::Rect(tl, br);
		m_x_axis.title = x_name;
		m_y_axis.title = y_name;
	}

	~Viewer() {
		//delete m_pShape;
	}

	enum graph_choice {
		DEPEN, DERIV, XY, PHASE
	};

	void plot(const phys::ODEStorage &data,
			graph_choice choice = DEPEN, uint dim1 = 1, uint dim2 = 2);

	void plot(const phys::Storage<double> &data,
			bool splitAxes = false);

	void plot(const stdVec_d &x, const stdVec_d &y);

	void plot(const stdVec_d &x, const std::vector<stdVec_d> &y);

	void plot_split(const stdVec_d &x, const std::vector<stdVec_d> &y);

	void add_data(const phys::ODEStorage &data, uint dim1 = 1,
			uint dim2 = 2);

	void add_data(const phys::Storage<double> &data,
			bool is_data_singular = true);

	void add_data(const stdVec_d &x, const stdVec_d &y, bool is_data_singular =
			true);

	void add_data(const stdVec_d &x, const std::vector<stdVec_d> &y,
			bool is_data_singular = true);

	/* Specify a path and filname including .extension (e.g. .png, .jpg, .bmp) */
	void save(const std::string &filename) const {
		cv::imwrite(filename, m_BG);
	}

	/*	Draw a straight line on the Viewer, specified by the y-intercept and the gradient */
	void draw_line(double gradient, double intercept);

	/* Clear the current Viewer including the axes and labels - these are redrawn when new data graphed */
	void clear() {
		set_bg_col(m_bg_colour);
	}

	/*	Warning: calls matrix constructor and assigns */
	void set_resolution(int num_cols, int num_rows, const cv::Scalar &bg_clr =
			cv::Scalar()) {
		m_rows = num_rows;
		m_cols = num_cols;
		m_font_scale = m_rows * m_cols / 1000 / 1000;
		if (bg_clr == cv::Scalar())
			m_BG = cv::Mat(m_rows, m_cols, CV_8UC3, m_bg_colour);
		else
			set_bg_col(bg_clr);
	}
	/*	Warning: calls matrix constructor and assigns
	 This will clear the entire image matrix and set the colour to that passed */
	void set_bg_col(const cv::Scalar &colour) {
		m_bg_colour = colour;
		m_BG = cv::Mat(m_rows, m_cols, CV_8UC3, m_bg_colour);
	}

	void set_shape(shape_choice shp_ch, uint shp_sz);
	void set_ax_col(const cv::Scalar &colour) {
		m_ax_colour = colour;
	}
	void set_dt_col(const cv::Scalar &colour) {
		m_dt_colour = colour;
	}

	void set_x_name(const std::string &name) {
		m_x_axis.title = name;
	}
	void set_y_name(const std::string &name) {
		m_y_axis.title = name;
	}
	void set_y2_name(const std::string &name) {
		m_y_axis2.title = name;
	}
	void set_key_name(const std::vector<std::string> &names) {
		m_key_titles = names;
	}
	void set_plot_name(const std::string &name) {
		m_plot_name = name;
	}

	void set_x_origin(double value) {
		m_val_origin_x = value;
	}
	void set_y_origin(double value) {
		m_val_origin_y = value;
	}

	void set_pause(int new_pause) {
		m_pause = new_pause;
	}

	void set_x_range(double x_min, double x_max) {
		assert(x_max > x_min);
		m_x_range = std::pair<double, double>(x_min, x_max);
		m_using_range_x = true;
	}
	void set_y_range(double y_min, double y_max) {
		assert(y_max > y_min);
		m_y_range = std::pair<double, double>(y_min, y_max);
		m_using_range_y = true;
	}
	void withLines() {
		if (m_animate)
			std::cout
					<< "Warning\nLines cannot be drawn when using animation.\n";
		m_lines = true;
	}
	void withPulses() {
		if (m_animate)
			std::cout
					<< "Warning\nPulses cannot be drawn when using animation.\n";
		m_pulses = true;
	}
	void withPoints() {
		m_lines = false;
	}
	void animation() {
		if (m_lines)
			std::cout
					<< "Warning\nLines cannot be drawn when using animation.\n";
		if (m_pulses)
			std::cout
					<< "Warning\nPulses cannot be drawn when using animation.\n";
		m_animate = true;
	}
	void noAnimation() {
		m_animate = false;
	}
private:

	//sets the axis based on the data values
	template<typename T>
	void processAxis(const std::vector<T> &values, axis which);

	//finds the largest range in the multiple dependent data and sets the yaxis accordingly
	template<typename T>
	void processAxis(const std::vector<std::vector<T>> &multipleValues);

	void draw_axes(const cv::Point &origin);

	void drawXaxis(int x_origin);
	void drawYaxis(int y_origin, bool onLeft = true);

	void draw_labels(axis which, double value_origin, double value_spacing,
			double scale, const cv::Point &m_axis_origin);

	void plot_impl(const coords &data);
	void plot_ani_impl(const coords &data, int delay = 10);

	void draw_key(uint which, const cv::Scalar &colour);

	void chk_data(coords &data);
};

template<typename T>
void Viewer::processAxis(const std::vector<T> &values, axis which) {
	T minVal, maxVal;
	double minScale, maxScale;

	switch (which) {
	case X:

		if (!m_using_range_x) {
			phys::minMax(values, minVal, maxVal);
		} else {
			minVal = T(m_x_range.first);
			maxVal = T(m_x_range.second);
			m_val_origin_x = minVal;
		}

		//If min value is less than the origin then scale using min value. If min value greater than the origin
		//and not using ranges scale using the origin else scale with min value.
		minScale = (minVal < m_val_origin_x) ? minVal : m_val_origin_x;
		//If max value is greater than the origin then scale using max value. If max value less than the origin
		//and not using ranges scale using the origin else scale with max value.
		maxScale = (maxVal > m_val_origin_x) ? maxVal : m_val_origin_x;

		if (minScale == maxScale) {
			throw std::runtime_error(
					"Problem with data: have you passed in a vector of constant values?\n");
		}

		m_x_axis.scale = m_plot_area.width / (maxScale - minScale);
		m_x_axis.min = minVal;
		m_x_axis.max = maxVal;
		m_x_axis.origin = m_val_origin_x;

		break;
	case Y:

		if (!m_using_range_y) {
			phys::minMax(values, minVal, maxVal);
		} else {
			minVal = T(m_y_range.first);
			maxVal = T(m_y_range.second);
			m_val_origin_y = minVal;
		}

		//If min value is less than the origin then scale using min value. If min value greater than the origin
		//and not using ranges scale using the origin else scale with min value.
		minScale = (minVal < m_val_origin_y) ? minVal : m_val_origin_y;
		//If max value is greater than the origin then scale using max value. If max value less than the origin
		//and not using ranges scale using the origin else scale with max value.
		maxScale = (maxVal > m_val_origin_y) ? maxVal : m_val_origin_y;

		if (minScale == maxScale) {
			throw std::runtime_error(
					"Problem with data: have you passed in a vector of constant values?\n");
		}

		m_y_axis.scale = m_plot_area.height / (maxScale - minScale);
		m_y_axis.min = minVal;
		m_y_axis.max = maxVal;
		m_y_axis.origin = m_val_origin_y;

		break;
	case Y2:

		phys::minMax(values, minVal, maxVal);

		//If min value is less than the origin then scale using min value.
		minScale = (minVal < m_val_origin_y) ? minVal : m_val_origin_y;
		//If max value is greater than the origin then scale using max value.
		maxScale = (maxVal > m_val_origin_y) ? maxVal : m_val_origin_y;

		if (minScale == maxScale) {
			throw std::runtime_error(
					"Problem with data: have you passed in a vector of constant values?\n");
		}

		m_y_axis2.scale = m_plot_area.height / (maxScale - minScale);
		m_y_axis2.min = minVal;
		m_y_axis2.max = maxVal;
		m_y_axis2.origin = m_val_origin_y;

		break;
	}
}

//Function to find the overall min and max values in a set of multiple dependent data.
template<typename T>
void find_min_max(const std::vector<std::vector<T>> &values, T &minVal,
		T &maxVal) {
	std::vector<T> min(values.size()), max(values.size());

	//For each y find the max and min value and store to a vector.
	for (uint i = 0; i < values.size(); ++i)
		phys::minMax(values[i], min[i], max[i]);

	T temp;

	//find the smallest of the minimums, and the largest of the maximums
	// here we keep sign convention, e.g. -2 < -1.
	phys::minMax(min, minVal, temp);
	phys::minMax(max, temp, maxVal);
}

template<typename T>
void Viewer::processAxis(const std::vector<std::vector<T>> &values) {
	T minVal, maxVal;
	double minScale, maxScale;

	if (!m_using_range_y) {
		find_min_max(values, minVal, maxVal);
	} else {
		minVal = T(m_y_range.first);
		maxVal = T(m_y_range.second);
		m_val_origin_y = minVal;
	}

	minScale = (minVal < m_val_origin_y) ? minVal : m_val_origin_y;
	maxScale = (maxVal > m_val_origin_y) ? maxVal : m_val_origin_y;

	m_y_axis.scale = m_plot_area.height / (maxScale - minScale);
	m_y_axis.min = minVal;
	m_y_axis.max = maxVal;
	m_y_axis.origin = m_val_origin_y;
}

} //namespace

#endif //header guard
