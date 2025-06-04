#include <windows.h>
#include <commdlg.h>
#include <iostream>
#include <stack>
#include <vector>
#include <list>
#include <fstream>
#include <cmath>

using std::cout, std::endl, std::stack, std::vector, std::swap, std::max, std::min, std::list, std::ofstream, std::ios;

// Command IDs
enum MenuID {
    ID_COLOR_RED = 1,
    ID_COLOR_GREEN,
    ID_COLOR_BLUE,
    ID_CURSOR_ARROW,
    ID_CURSOR_CROSS,
    ID_CLEAR_SCREEN,
    ID_SAVE_IMAGE,
    ID_LOAD_IMAGE,
    ID_DDA_LINE,
    ID_MIDPOINT_LINE,
    ID_PARAMETRIC_LINE,
    ID_DIRECT_CIRCLE,
    ID_POLAR_CIRCLE,
    ID_ITER_POLAR_CIRCLE,
    ID_MIDPOINT_CIRCLE,
    ID_MOD_MIDPOINT_CIRCLE,
    ID_FILL_CIRCLE_LINES,
    ID_FILL_CIRCLE_CIRCLES,
    ID_FILL_SQUARE_HERMIT,
    ID_FILL_RECTANGLE_BEZIER,
    ID_FILL_CONVEX,
    ID_FILL_NONCONVEX,
    ID_FLOOD_RECURSIVE,
    ID_FLOOD_NONRECURSIVE,
    ID_CARDINAL_SPLINE,
    ID_ELLIPSE_DIRECT,
    ID_ELLIPSE_POLAR,
    ID_ELLIPSE_MIDPOINT,
    ID_CLIP_RECT_POINT,
    ID_CLIP_RECT_LINE,
    ID_CLIP_RECT_POLYGON,
    ID_CLIP_SQUARE_POINT,
    ID_CLIP_SQUARE_LINE
};

HBITMAP loaded_bitmap {nullptr};
COLORREF current_color {RGB(255, 0, 0)};
COLORREF bg_color {RGB(255, 255, 255)};
POINT points[200];
int point_count {0};
int current_algorithm {ID_DDA_LINE};
HCURSOR current_cursor {LoadCursor(nullptr, IDC_ARROW)};

typedef struct {
    int xleft, xright;
} EdgeTable[80000];

struct Node {
    double x, minv;
    int ymax;
    explicit Node(const double _x = 0, const int _ymax = 0, const double _minv = 0) : x{_x}, minv{_minv}, ymax{_ymax} {}

    bool operator<(const Node& other) const {
        return x < other.x;
    }
};

// Utility to draw line using DDA
void draw_dda_line(const HDC hdc, POINT p1, POINT p2, const COLORREF color) {
    const int dx = p2.x - p1.x;
    const int dy = p2.y - p1.y;

    SetPixel(hdc, p1.x, p1.y, color);

    if (abs(dx) >= abs(dy)) {
        const double m = static_cast<double>(dy) / dx;
        if (p1.x > p2.x) {
            swap(p1.x, p2.x);
            swap(p1.y, p2.y);
        }
        double y = p1.y;
        for (int x = p1.x; x <= p2.x; ++x) {
            y += m;
            SetPixel(hdc, x, round(y), color);
        }
    } else {
        const double mi {static_cast<double>(dx) / dy};
        if (p1.y > p2.y) {
            swap(p1.x, p2.x);
            swap(p1.y, p2.y);
        }
        double x = p1.x;
        for (int y = p1.y; y <= p2.y; ++y) {
            x += mi;
            SetPixel(hdc, round(x), y, color);
        }
    }
}

// Square clipping window boundaries
struct square_clip_window {
    int left, right, top, bottom;
    int size;
};

struct rect_clip_window {
    int left, right, top, bottom;
    int width, height;
};

square_clip_window square_window = {0, 0, 0, 0, 100}; // Default 100x100 square
rect_clip_window rectangle_window = {0, 0, 0, 0, 150, 100}; // Default 150x100 rectangle

// Union for outcodes (from your algorithm)
union Outcode {
    struct {
        unsigned left: 1;
        unsigned right: 1;
        unsigned top: 1;
        unsigned bottom: 1;
    };
    unsigned all: 4;
};

// Initialize clipping window in the center of the screen
void initialize_square_clip_window(const HWND hwnd) {
    RECT rect;
    GetClientRect(hwnd, &rect);
    const int width = rect.right - rect.left;
    const int height = rect.bottom - rect.top;

    square_window.left = (width - square_window.size) / 2;
    square_window.top = (height - square_window.size) / 2;
    square_window.right = square_window.left + square_window.size;
    square_window.bottom = square_window.top + square_window.size;
}

void initialize_rect_clip_window(const HWND hwnd) {
    RECT rect;
    GetClientRect(hwnd, &rect);
    const int screenWidth = rect.right - rect.left;
    const int screenHeight = rect.bottom - rect.top;

    rectangle_window.left = (screenWidth - rectangle_window.width) / 2;
    rectangle_window.top = (screenHeight - rectangle_window.height) / 2;
    rectangle_window.right = rectangle_window.left + rectangle_window.width;
    rectangle_window.bottom = rectangle_window.top + rectangle_window.height;
}

// Draw the clipping window
void draw_square_clip_window(const HDC hdc) {
    const HPEN pen {CreatePen(PS_SOLID, 2, RGB(255, 0, 0))}; // Red border
    const auto old_pen {static_cast<HPEN>(SelectObject(hdc, pen))};

    Rectangle(hdc, square_window.left, square_window.top, square_window.right + 1, square_window.bottom + 1);

    SelectObject(hdc, old_pen);
    DeleteObject(pen);
}

void draw_rect_clip_window(const HDC hdc) {
    const HPEN pen {CreatePen(PS_SOLID, 2, RGB(0, 0, 255))}; // Blue border to distinguish from square
    const auto old_pen {static_cast<HPEN>(SelectObject(hdc, pen))};

    Rectangle(hdc, rectangle_window.left, rectangle_window.top,
              rectangle_window.right + 1, rectangle_window.bottom + 1);

    SelectObject(hdc, old_pen);
    DeleteObject(pen);
}

// Generate outcode for a point (from your algorithm)
Outcode get_outcode(const POINT p, const int xl, const int xr, const int yb, const int yt) {
    Outcode result = {0};

    if (p.x < xl) {
        result.left = 1;
    }

    if (p.x > xr) {
        result.right = 1;
    }
    if (p.y < yb) {
        result.bottom = 1;
    }
    if (p.y > yt) {
        result.top = 1;
    }

    return result;
}

// Vertical intersection (from your algorithm)
POINT vertical_intersect(const POINT p1, const POINT p2, const int x_edge) {
    POINT result;
    result.x = x_edge;
    result.y = p1.y + (x_edge - p1.x) * static_cast<double>(p2.y - p1.y) / (p2.x - p1.x);
    return result;
}

// Horizontal intersection (from your algorithm)
POINT horizontal_intersect(const POINT p1, const POINT p2, const int y_edge) {
    POINT result;
    result.y = y_edge;
    result.x = p1.x + (y_edge - p1.y) * static_cast<double>(p2.x - p1.x) / (p2.y - p1.y);
    return result;
}

// Point clipping using square window
bool clip_point_square(const HDC hdc, const POINT p, const COLORREF color) {
    // Check if point is inside the clipping window
    if (p.x >= square_window.left && p.x <= square_window.right &&
        p.y >= square_window.top && p.y <= square_window.bottom) {
        // Point is inside, draw it
        SetPixel(hdc, p.x, p.y, color);
        return true;
    }
    return false; // Point is clipped (outside window)
}

bool clip_point_rect(const HDC hdc, const POINT p, const COLORREF color) {
    // Check if point is inside the rectangle clipping window
    if (p.x >= rectangle_window.left && p.x <= rectangle_window.right &&
        p.y >= rectangle_window.top && p.y <= rectangle_window.bottom) {
        // Point is inside, draw it
        SetPixel(hdc, p.x, p.y, color);
        return true;
        }
    return false; // Point is clipped (outside window)
}

// Line clipping using Cohen-Sutherland algorithm with square window (adapted from your code)
void clip_line_square(const HDC hdc, POINT p1, POINT p2, const COLORREF color) {
    Outcode out1 {get_outcode(p1, square_window.left, square_window.right, square_window.top, square_window.bottom)};
    Outcode out2 {get_outcode(p2, square_window.left, square_window.right, square_window.top, square_window.bottom)
};

    while (true) {
        // Both points inside window
        if (out1.all == 0 && out2.all == 0) {
            draw_dda_line(hdc, p1, p2, color);
            return;
        }

        // Both points outside window on same side
        if ((out1.all & out2.all) != 0) {
            return; // Line is completely outside
        }

        // Line crosses window boundary
        Outcode outcode{};
        POINT* point;

        if (out1.all != 0) {
            outcode = out1;
            point = &p1;
        } else {
            outcode = out2;
            point = &p2;
        }

        // Find intersection point
        if (outcode.left) {
            *point = vertical_intersect(p1, p2, square_window.left);
        } else if (outcode.right) {
            *point = vertical_intersect(p1, p2, square_window.right);
        } else if (outcode.bottom) {
            *point = horizontal_intersect(p1, p2, square_window.top);
        } else if (outcode.top) {
            *point = horizontal_intersect(p1, p2, square_window.bottom);
        }

        // Update outcode for the modified point
        if (point == &p1) {
            out1 = get_outcode(p1, square_window.left, square_window.right, square_window.top, square_window.bottom);
        } else {
            out2 = get_outcode(p2, square_window.left, square_window.right, square_window.top, square_window.bottom);
        }
    }
}

void clip_line_rect(const HDC hdc, POINT p1, POINT p2, const COLORREF color) {
    Outcode out1 {get_outcode(p1, rectangle_window.left, rectangle_window.right,
                              rectangle_window.top, rectangle_window.bottom)};
    Outcode out2 {get_outcode(p2, rectangle_window.left, rectangle_window.right,
                              rectangle_window.top, rectangle_window.bottom)};

    while (true) {
        // Both points inside window
        if (out1.all == 0 && out2.all == 0) {
            draw_dda_line(hdc, p1, p2, color);
            return;
        }

        // Both points outside window on same side
        if ((out1.all & out2.all) != 0) {
            return; // Line is completely outside
        }

        // Line crosses window boundary
        Outcode outcode{};
        POINT* point;

        if (out1.all != 0) {
            outcode = out1;
            point = &p1;
        } else {
            outcode = out2;
            point = &p2;
        }

        // Find intersection point
        if (outcode.left) {
            *point = vertical_intersect(p1, p2, rectangle_window.left);
        } else if (outcode.right) {
            *point = vertical_intersect(p1, p2, rectangle_window.right);
        } else if (outcode.bottom) {
            *point = horizontal_intersect(p1, p2, rectangle_window.top);
        } else if (outcode.top) {
            *point = horizontal_intersect(p1, p2, rectangle_window.bottom);
        }

        // Update outcode for the modified point
        if (point == &p1) {
            out1 = get_outcode(p1, rectangle_window.left, rectangle_window.right,
                              rectangle_window.top, rectangle_window.bottom);
        } else {
            out2 = get_outcode(p2, rectangle_window.left, rectangle_window.right,
                              rectangle_window.top, rectangle_window.bottom);
        }
    }
}

// Midpoint Line Algorithm
void draw_midpoint_line(const HDC hdc, const POINT p1, const POINT p2, const COLORREF color) {
    int dx = p2.x - p1.x;
    int dy = p2.y - p1.y;
    int x = p1.x, y = p1.y;

    const int xInc = (dx > 0) ? 1 : -1;
    const int yInc = (dy > 0) ? 1 : -1;
    dx = abs(dx);
    dy = abs(dy);

    if (dx >= dy) { // Slope <= 1
        int d = 2 * dy - dx;
        const int d1 = 2 * dy;
        const int d2 = 2 * (dy - dx);

        SetPixel(hdc, x, y, color);

        while (x != p2.x) {
            if (d < 0) {
                d += d1;
            } else {
                d += d2;
                y += yInc;
            }
            x += xInc;
            SetPixel(hdc, x, y, color);
        }
    } else { // Slope > 1
        int d = 2 * dx - dy;
        const int d1 = 2 * dx;
        const int d2 = 2 * (dx - dy);

        SetPixel(hdc, x, y, color);
        while (y != p2.y) {
            if (d < 0) {
                d += d1;
            } else {
                d += d2;
                x += xInc;
            }
            y += yInc;
            SetPixel(hdc, x, y, color);
        }
    }
}

// Parametric Line Algorithm
void draw_parametric_line(const HDC hdc, const POINT p1, const POINT p2, const COLORREF color) {
    for (float t {0}; t <= 1; t += 0.001) {
        const int x = p1.x + t * (p2.x - p1.x);
        const int y = p1.y + t * (p2.y - p1.y);
        SetPixel(hdc, x, y, color);
    }
}

void draw_8points(const HDC hdc, const POINT center, const int x, const int y, const COLORREF color) {
    SetPixel(hdc, center.x + x, center.y + y, color);
    SetPixel(hdc, center.x - x, center.y + y, color);
    SetPixel(hdc, center.x + x, center.y - y, color);
    SetPixel(hdc, center.x - x, center.y - y, color);
    SetPixel(hdc, center.x + y, center.y + x, color);
    SetPixel(hdc, center.x - y, center.y + x, color);
    SetPixel(hdc, center.x + y, center.y - x, color);
    SetPixel(hdc, center.x - y, center.y - x, color);
}

// Direct Circle Algorithm (using Draw8Points)
void draw_direct_circle(const HDC hdc, const POINT center, const int radius, const COLORREF color) {
    int x {0}, y {radius};
    const int radius_sq {radius * radius};

    while (x <= y) {
        draw_8points(hdc, center, x, y, color);
        ++x;
        y = round(sqrt(radius_sq - x * x));
    }
}

// Polar Circle Algorithm (using Draw8Points)
void draw_polar_circle(const HDC hdc, const POINT center, const int radius, const COLORREF color) {
    const double delta_theta = 1.0 / radius;
    for (double theta {0}; theta <= M_PI/4; theta += delta_theta) {
        const int x = round(radius * cos(theta));
        const int y = round(radius * sin(theta));
        draw_8points(hdc, center, x, y, color);
    }
}

// Iterative Polar Circle Algorithm (using Draw8Points)
void draw_iterative_polar_circle(const HDC hdc, const POINT center, const int radius, const COLORREF color) {
    double x = radius, y {0};
    const double delta_theta = 1.0 / radius;
    const double cd_theta = cos(delta_theta), sd_theta = sin(delta_theta);

    while (x >= y) {
        draw_8points(hdc, center, round(x), round(y), color);

        // Rotate point
        const double x_new = x * cd_theta - y * sd_theta;
        y = x * sd_theta + y * cd_theta;
        x = x_new;
    }
}

// Midpoint Circle Algorithm (using Draw8Points)
void draw_midpoint_circle(const HDC hdc, const POINT center, const int radius, const COLORREF color) {
    int x {0}, y {radius};
    int d {1 - radius};

    while (x <= y) {
        draw_8points(hdc, center, x, y, color);

        if (d < 0) {
            d += 2 * x + 3;
        } else {
            d += 2 * (x - y) + 5;
            --y;
        }
        ++x;
    }
}

// Modified Midpoint Circle Algorithm (using Draw8Points)
void draw_modified_midpoint_circle(const HDC hdc, const POINT center, const int radius, const COLORREF color) {
    int x {0}, y {radius};
    int d {1 - radius};
    int d1 {3};
    int d2 {-2 * radius + 5};

    while (x <= y) {
        draw_8points(hdc, center, x, y, color);
        if (d < 0) {
            d += d1;
            d1 += 2;
            d2 += 2;
        } else {
            d += d2;
            d1 += 2;
            d2 += 4;
            --y;
        }
        ++x;
    }
}

void draw_direct_ellipse(const HDC hdc, const POINT center, const int rx, const int ry, const COLORREF color) {
    if (rx <= 0 || ry <= 0) {
        return;
    }

    const double rx2 {static_cast<double>(rx) * rx};
    const double ry2 {static_cast<double>(ry) * ry};

    for (int xi {0}; xi <= rx; ++xi) {
        const auto xd {static_cast<double>(xi)};
        const double yd {ry * sqrt(1.0 - (xd * xd) / rx2)};
        const int yi {static_cast<int>(round(yd))};

        SetPixel(hdc, center.x + xi, center.y + yi, color);
        SetPixel(hdc, center.x - xi, center.y + yi, color);
        SetPixel(hdc, center.x + xi, center.y - yi, color);
        SetPixel(hdc, center.x - xi, center.y - yi, color);
    }

    for (int yi {0}; yi <= ry; ++yi) {
        const auto yd {static_cast<double>(yi)};
        const double xd {rx * sqrt(1.0 - (yd * yd) / ry2)};
        const int xi {static_cast<int>(round(xd))};

        SetPixel(hdc, center.x + xi, center.y + yi, color);
        SetPixel(hdc, center.x - xi, center.y + yi, color);
        SetPixel(hdc, center.x + xi, center.y - yi, color);
        SetPixel(hdc, center.x - xi, center.y - yi, color);
    }
}

// Polar Ellipse Algorithm
void draw_polar_ellipse(const HDC hdc, const POINT center, const int rx, const int ry, const COLORREF color) {
    if (rx <= 0 || ry <= 0) {
        return;
    }

    const double delta_theta {1.0 / max(rx, ry)};
    for (double theta {0}; theta <= M_PI/2; theta += delta_theta) {
        const int x = round(rx * cos(theta));
        const int y = round(ry * sin(theta));
        SetPixel(hdc, center.x + x, center.y + y, color);
        SetPixel(hdc, center.x - x, center.y + y, color);
        SetPixel(hdc, center.x + x, center.y - y, color);
        SetPixel(hdc, center.x - x, center.y - y, color);
    }
}

// Midpoint Ellipse Algorithm
void draw_midpoint_ellipse(const HDC hdc, const POINT center, const int rx, const int ry, const COLORREF color) {
    if (rx <= 0 || ry <= 0) {
        return;
    }

    const double rx2 = rx * rx;
    const double ry2 = ry * ry;
    double x {0}, y = ry;
    double p = ry2 - rx2 * ry + 0.25 * rx2;

    // Region 1
    while (ry2 * x < rx2 * y) {
        SetPixel(hdc, center.x + x, center.y + y, color);
        SetPixel(hdc, center.x - x, center.y + y, color);
        SetPixel(hdc, center.x + x, center.y - y, color);
        SetPixel(hdc, center.x - x, center.y - y, color);

        if (p < 0) {
            p += 2 * ry2 * x + ry2;
        } else {
            --y;
            p += 2 * ry2 * x - 2 * rx2 * y + ry2;
        }
        ++x;
    }

    // Region 2
    p = ry2 * (x + 0.5) * (x + 0.5) + rx2 * (y - 1) * (y - 1) - rx2 * ry2;
    while (y >= 0) {
        SetPixel(hdc, center.x + x, center.y + y, color);
        SetPixel(hdc, center.x - x, center.y + y, color);
        SetPixel(hdc, center.x + x, center.y - y, color);
        SetPixel(hdc, center.x - x, center.y - y, color);

        if (p > 0) {
            p += -2 * rx2 * y + rx2;
        } else {
            ++x;
            p += 2 * ry2 * x - 2 * rx2 * y + rx2;
        }
        --y;
    }
}

void flood_fill_recursive(const HDC hdc, const int x, const int y, const COLORREF target_color, const COLORREF fill_color) {
    // Get current pixel color

    // Base case: if current color is not target color or already filled
    if (const COLORREF current_color {GetPixel(hdc, x, y)}; current_color != target_color || current_color == fill_color) {
        return;
    }

    // Fill current pixel
    SetPixel(hdc, x, y, fill_color);

    // Recursively fill in 4 directions (4-connected)
    flood_fill_recursive(hdc, x + 1, y, target_color, fill_color);
    flood_fill_recursive(hdc, x - 1, y, target_color, fill_color);
    flood_fill_recursive(hdc, x, y + 1, target_color, fill_color);
    flood_fill_recursive(hdc, x, y - 1, target_color, fill_color);
}

void flood_fill_iterative(const HDC hdc, const int start_x, const int start_y, const COLORREF target_color, const COLORREF fill_color) {
    // Use stack to simulate recursion
    stack<POINT> pixel_stack;

    // Push starting point
    const POINT start_point {start_x, start_y};
    pixel_stack.push(start_point);

    while (!pixel_stack.empty()) {
        const auto [x, y] {pixel_stack.top()};
        pixel_stack.pop();

        if (const COLORREF current_color = GetPixel(hdc, x, y); current_color != target_color || current_color == fill_color) {
            continue;
        }

        SetPixel(hdc, x, y, fill_color);

        // Add neighboring pixels to stack (4-connected)
        POINT neighbors[4] = {
            {x + 1, y},     // Right
            {x - 1, y},     // Left
            {x, y + 1},     // Down
            {x, y - 1}      // Up
        };

        for (auto neighbor : neighbors) {
            pixel_stack.push(neighbor);
        }
    }
}

void draw_hermite_vertical_curve_points(vector<int>& pixelsY, const int y0, const int y1, const int t0y, const int t1y) {
    const double a0 = y0;
    const double a1 = t0y;
    const double a2 = -3 * y0 + 3 * y1 - 2 * t0y - t1y;
    const double a3 = 2 * y0 - 2 * y1 + t0y + t1y;
    pixelsY.clear();

    for (double t {0.0}; t <= 1.0; t += 0.001) {
        const double y = a0 + a1 * t + a2 * t * t + a3 * t * t * t;
        pixelsY.push_back(static_cast<int>(y));
    }
}

void fill_square_with_vertical_hermite(const HDC hdc, const POINT p1, const POINT p2, const COLORREF color) {
    const int width = abs(p2.x - p1.x);
    const int height = abs(p2.y - p1.y);
    const int side = min(width, height);
    const int dir_x = (p2.x > p1.x) ? 1 : -1;
    const int dir_y = (p2.y > p1.y) ? 1 : -1;
    const int x0 = p1.x;
    const int y0 = p1.y;
    const int x1 = x0 + dir_x * side;
    const int y1 = y0 + dir_y * side;
    const int t0y = -side / 4;
    const int t1y = -side / 4;
    const int start_x = min(x0, x1);
    const int end_x = max(x0, x1);
    const int start_y = min(y0, y1);
    const int end_y = max(y0, y1);
    constexpr int total_steps {1000};
    vector<vector<int>> all_y(end_x - start_x + 1);
    for (int xi {0}; xi <= end_x - start_x; ++xi) {
        vector<int>& pixels_y {all_y[xi]};
        draw_hermite_vertical_curve_points(pixels_y, end_y, start_y, t0y, t1y);
    }
    for (int step {0}; step < total_steps; ++step) {
        for (int xi {0}; xi <= end_x - start_x; ++xi) {
            const int x {start_x + xi};
            if (const vector<int>& pixels_y {all_y[xi]}; step < pixels_y.size()) {
                SetPixel(hdc, x, pixels_y[step], color);
            }
        }
    }
}

void draw_bezier_horizontal_curve(const HDC hdc, const int x, const POINT p0, const POINT p1, const POINT p2, const POINT p3) {
    MoveToEx(hdc, x, p0.y, nullptr);
    for (double t {0.0}; t <= 1.0; t += 0.01) {
        const double y {pow(1 - t, 3) * p0.y +
                   3 * pow(1 - t, 2) * t * p1.y +
                   3 * (1 - t) * t * t * p2.y +
                   t * t * t * p3.y};

        LineTo(hdc, x, static_cast<int>(y));
    }
}

void fill_rectangle_with_horizontal_bezier(const HDC hdc, const POINT p1, const POINT p2, const COLORREF color) {
    const int width = abs(p2.x - p1.x);
    const int height = abs(p2.y - p1.y);

    const int dir_x = (p2.x > p1.x) ? 1 : -1;
    const int dir_y = (p2.y > p1.y) ? 1 : -1;

    const int x_start = p1.x;
    const int x_end = p1.x + dir_x * width;

    POINT top    = { 0, p1.y };
    POINT ctrl1  = { 0, p1.y + dir_y * height / 3 };
    POINT ctrl2  = { 0, p1.y + dir_y * 2 * height / 3 };
    POINT bottom = { 0, p1.y + dir_y * height };

    const HPEN pen {CreatePen(PS_SOLID, 1, color)};
    const auto old_pen {static_cast<HPEN>(SelectObject(hdc, pen))};

    for (int x = min(x_start, x_end); x <= max(x_start, x_end); ++x) {
        top.x = ctrl1.x = ctrl2.x = bottom.x = x;
        draw_bezier_horizontal_curve(hdc, x, top, ctrl1, ctrl2, bottom);
    }

    SelectObject(hdc, old_pen);
    DeleteObject(pen);
}

void init(EdgeTable tbl) {
    for (int i {0}; i < 80000; ++i) {
        tbl[i].xleft = 10000;
        tbl[i].xright = -10000;
    }
}

void edge2table(POINT &p1, POINT &p2, EdgeTable tbl){
    if (p1.y == p2.y) {
        return;
    }
    if (p1.y > p2.y) {
        swap(p1, p2);
    }
    int y = p1.y;
    double x = p1.x;
    const double d {static_cast<double>(p2.x - p1.x) / (p2.y - p1.y)};

    while (y < p2.y) {
        if (x < tbl[y].xleft) {
            tbl[y].xleft = static_cast<int>(ceil(x));
        }
        if (x > tbl[y].xright) {
            tbl[y].xright = static_cast<int>(floor(x));
        }

        ++y;
        x += d;
    }
}

void polygon2table(POINT p[], const int n, EdgeTable tbl) {
    POINT point = p[n - 1];
    for (int i {0}; i < n; ++i) {
        POINT next_point {p[i]};
        edge2table(point, next_point, tbl);

        point = p[i];
    }
}

void table2screen(const HDC hdc, EdgeTable tbl, const COLORREF color) {
    for (int y {0}; y < 800; ++y) {
        if (tbl[y].xleft < tbl[y].xright) {
            draw_dda_line(hdc, {tbl[y].xleft, y}, {tbl[y].xright, y}, color);
        }
    }
}

void convex_fill(const HDC hdc, POINT p[], const int n, const COLORREF color) {
    EdgeTable tbl;
    init(tbl);
    polygon2table(p, n, tbl);
    table2screen(hdc, tbl, color);
}

void edge2table_general(POINT &p1, POINT &p2, list<Node> tbl[]) {
    if (p1.y == p2.y) {
        return;
    }
    if (p1.y > p2.y) {
        swap(p1, p2);
    }
    tbl[p1.y].emplace_back(p1.x, p2.y, static_cast<double>(p2.x - p1.x) / (p2.y - p1.y));
}

void polygon2table_general(POINT p[], const int n, list<Node> tbl[]) {
    POINT point {p[n - 1]};
    for (int i {0}; i < n; ++i) {
        POINT next_point {p[i]};
        edge2table_general(point, next_point, tbl);
        point = p[i];
    }
}

void table2screen_general(const HDC hdc, list<Node> tbl[], const COLORREF color) {
    int y {0};
    while (y < 80000 && tbl[y].empty()) {
        ++y;
    }
    list active {tbl[y]};
    while (!active.empty()) {
        active.sort();
        for (auto it {active.begin()}; it != active.end(); ++it) {
            const auto it2 = it++;
            const int x_left = static_cast<int>(ceil(it->x));
            const int x_right = static_cast<int>(floor(it2->x));
            draw_dda_line(hdc, {x_left, y}, {x_right, y}, color);
        }
        ++y;
        for (auto it {active.begin()}; it != active.end();) {
            if (it->ymax == y) {
                it = active.erase(it);
            } else {
                it->x += it->minv;
                ++it;
            }
        }
        for (auto &node : tbl[y]) {
            active.push_back(node);
        }
    }
}

void nonconvex_fill(const HDC hdc, POINT p[], const int n, const COLORREF color) {
    list<Node> tbl[80000];
    polygon2table_general(p, n, tbl);
    table2screen_general(hdc, tbl, color);
}

POINT hermite_interpolate(const POINT p0, const POINT m0, const POINT p1, const POINT m1, const float t) {
    const float t2 {t * t};
    const float t3 {t2 * t};

    const float h1 {2 * t3 - 3 * t2 + 1};
    const float h2 {t3 - 2 * t2 + t};
    const float h3 {-2 * t3 + 3 * t2};
    const float h4 {t3 - t2};

    POINT result;
    result.x = h1 * p0.x + h2 * m0.x + h3 * p1.x + h4 * m1.x;
    result.y = h1 * p0.y + h2 * m0.y + h3 * p1.y + h4 * m1.y;
    return result;
}

void draw_cardinal_spline(const HDC hdc, const vector<POINT>& points, const float c, const COLORREF color) {
    if (points.size() < 4) {
        return;
    }

    vector<POINT> q;
    q.push_back({0, 0});

    for (size_t i = 1; i < points.size() - 1; ++i) {
        const int new_x = (points[i + 1].x - points[i - 1].x) * (c / 2.0f);
        const int new_y = (points[i + 1].x - points[i - 1].x) * (c / 2.0f);
        q.push_back({new_x, new_y});
    }

    for (size_t i {1}; i < points.size() - 2; ++i) {
        for (float t {0}; t < 1.0f; t += 0.0001f) {
            auto [x, y] {hermite_interpolate(points[i], q[i], points[i + 1], q[i + 1], t)};

            SetPixel(hdc, static_cast<int>(x), static_cast<int>(y), color);
        }
    }
}

// Save screen content using BitBlt and file dialog
void save_screen_to_file(HWND hwnd) {
    HDC hdc_screen {GetDC(hwnd)};
    HDC hdc_mem = CreateCompatibleDC(hdc_screen);
    RECT rc;
    GetClientRect(hwnd, &rc);

    HBITMAP h_bitmap = CreateCompatibleBitmap(hdc_screen, rc.right, rc.bottom);
    SelectObject(hdc_mem, h_bitmap);
    BitBlt(hdc_mem, 0, 0, rc.right, rc.bottom, hdc_screen, 0, 0, SRCCOPY);

    OPENFILENAME ofn = { sizeof(ofn) };
    char file_name[MAX_PATH] {"untitled.bmp"};  // Default name with extension

    ofn.lpstrFilter = "Bitmap Files (*.bmp)\0*.bmp\0All Files (*.*)\0*.*\0";
    ofn.lpstrFile = file_name;
    ofn.nMaxFile = MAX_PATH;
    ofn.lpstrDefExt = "bmp";  // Default extension if user doesn't specify
    ofn.lpstrTitle = "Save As Bitmap";
    ofn.Flags = OFN_OVERWRITEPROMPT | OFN_PATHMUSTEXIST;

    if (GetSaveFileName(&ofn)) {
        // Ensure the extension is .bmp
        if (char* ext {strrchr(file_name, '.')}; ext == nullptr || _stricmp(ext, ".bmp") != 0) {
            if (strlen(file_name) + 4 < MAX_PATH) {  // Check buffer space
                strcat(file_name, ".bmp");
            } else {
                MessageBox(hwnd, "Filename too long", "Error", MB_OK | MB_ICONERROR);
                return;
            }
        }

        BITMAPFILEHEADER bmf_header = {};
        BITMAPINFOHEADER bi = {};
        bi.biSize = sizeof(BITMAPINFOHEADER);
        bi.biWidth = rc.right;
        bi.biHeight = -rc.bottom;  // Negative for top-down DIB
        bi.biPlanes = 1;
        bi.biBitCount = 24;
        bi.biCompression = BI_RGB;

        DWORD dw_bmp_size = ((rc.right * 3 + 3) & ~3) * abs(rc.bottom);  // 24-bit aligned row size
        auto lp_bits {new BYTE[dw_bmp_size]};

        GetDIBits(hdc_mem, h_bitmap, 0, rc.bottom, lp_bits, reinterpret_cast<BITMAPINFO *>(&bi), DIB_RGB_COLORS);

        if (ofstream file(file_name, ios::binary); file.is_open()) {
            bmf_header.bfType = 0x4D42;  // 'BM'
            bmf_header.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
            bmf_header.bfSize = bmf_header.bfOffBits + dw_bmp_size;

            file.write(reinterpret_cast<char *>(&bmf_header), sizeof(bmf_header));
            file.write(reinterpret_cast<char *>(&bi), sizeof(bi));
            file.write(reinterpret_cast<char *>(lp_bits), dw_bmp_size);
            file.close();
        } else {
            MessageBox(hwnd, "Failed to create file", "Error", MB_OK | MB_ICONERROR);
        }

        delete[] lp_bits;
    }

    DeleteObject(h_bitmap);
    DeleteDC(hdc_mem);
    ReleaseDC(hwnd, hdc_screen);
}

void load_image_from_file(const HWND hwnd) {
    OPENFILENAME ofn = { sizeof(ofn) };
    char file_name[MAX_PATH] {""};

    ofn.hwndOwner = hwnd;
    ofn.lpstrFile = file_name;
    ofn.nMaxFile = MAX_PATH;
    ofn.lpstrFilter = "Bitmap Files (*.bmp)\0*.bmp\0All Files (*.*)\0*.*\0";
    ofn.lpstrTitle = "Select Bitmap File to Open";
    ofn.Flags = OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST | OFN_HIDEREADONLY;

    if (GetOpenFileName(&ofn)) {
        // Try to open the file to verify it exists
        const HANDLE h_file = CreateFile(
            file_name,
            GENERIC_READ,
            FILE_SHARE_READ,
            nullptr,
            OPEN_EXISTING,
            FILE_ATTRIBUTE_NORMAL,
            nullptr
        );

        if (h_file == INVALID_HANDLE_VALUE) {
            const DWORD err {GetLastError()};
            char msg[256];
            sprintf(msg, "Cannot open file (Error %d)", err);
            MessageBox(hwnd, msg, "Error", MB_OK | MB_ICONERROR);
            return;
        }
        CloseHandle(h_file);

        // Load the bitmap
        const auto h_bitmap {static_cast<HBITMAP>(LoadImage(
            nullptr,
            file_name,
            IMAGE_BITMAP,
            0,
            0,
            LR_LOADFROMFILE | LR_CREATEDIBSECTION | LR_DEFAULTSIZE
        ))};

        if (!h_bitmap) {
            const DWORD err {GetLastError()};
            char msg[256];
            sprintf(msg, "Failed to load image (Error %d)", err);
            MessageBox(hwnd, msg, "Error", MB_OK | MB_ICONERROR);
            return;
        }

        // Clean up previous bitmap if exists
        if (loaded_bitmap) {
            DeleteObject(loaded_bitmap);
        }

        loaded_bitmap = h_bitmap;
        InvalidateRect(hwnd, nullptr, TRUE);
    }
}

void add_menus(const HWND hwnd) {
    const HMENU h_main {CreateMenu()};
    HMENU h_color {CreateMenu()}, h_cursor {CreateMenu()}, h_lines {CreateMenu()}, h_circles {CreateMenu()};
    HMENU h_fills = CreateMenu(), h_clipping = CreateMenu(), h_ellipse = CreateMenu();

    // Color menu
    AppendMenu(h_color, MF_STRING, ID_COLOR_RED, "Red");
    AppendMenu(h_color, MF_STRING, ID_COLOR_GREEN, "Green");
    AppendMenu(h_color, MF_STRING, ID_COLOR_BLUE, "Blue");

    // Cursor menu
    AppendMenu(h_cursor, MF_STRING, ID_CURSOR_ARROW, "Arrow");
    AppendMenu(h_cursor, MF_STRING, ID_CURSOR_CROSS, "Cross");

    // Lines menu
    AppendMenu(h_lines, MF_STRING, ID_DDA_LINE, "DDA Line");
    AppendMenu(h_lines, MF_STRING, ID_MIDPOINT_LINE, "Midpoint Line");
    AppendMenu(h_lines, MF_STRING, ID_PARAMETRIC_LINE, "Parametric Line");

    // Other menus remain the same...
    AppendMenu(h_circles, MF_STRING, ID_DIRECT_CIRCLE, "Direct Circle");
    AppendMenu(h_circles, MF_STRING, ID_POLAR_CIRCLE, "Polar Circle");
    AppendMenu(h_circles, MF_STRING, ID_ITER_POLAR_CIRCLE, "Iterative Polar");
    AppendMenu(h_circles, MF_STRING, ID_MIDPOINT_CIRCLE, "Midpoint Circle");
    AppendMenu(h_circles, MF_STRING, ID_MOD_MIDPOINT_CIRCLE, "Modified Midpoint");

    AppendMenu(h_fills, MF_STRING, ID_FILL_CIRCLE_LINES, "Fill Circle w/ Lines");
    AppendMenu(h_fills, MF_STRING, ID_FILL_CIRCLE_CIRCLES, "Fill Circle w/ Circles");
    AppendMenu(h_fills, MF_STRING, ID_FILL_SQUARE_HERMIT, "Fill Square Hermite");
    AppendMenu(h_fills, MF_STRING, ID_FILL_RECTANGLE_BEZIER, "Fill Rect Bezier");
    AppendMenu(h_fills, MF_STRING, ID_FILL_CONVEX, "Fill Convex");
    AppendMenu(h_fills, MF_STRING, ID_FILL_NONCONVEX, "Fill Non-Convex");
    AppendMenu(h_fills, MF_STRING, ID_FLOOD_RECURSIVE, "Flood Fill Recursive");
    AppendMenu(h_fills, MF_STRING, ID_FLOOD_NONRECURSIVE, "Flood Fill Non-Recursive");
    AppendMenu(h_fills, MF_STRING, ID_CARDINAL_SPLINE, "Cardinal Spline");

    AppendMenu(h_ellipse, MF_STRING, ID_ELLIPSE_DIRECT, "Direct Ellipse");
    AppendMenu(h_ellipse, MF_STRING, ID_ELLIPSE_POLAR, "Polar Ellipse");
    AppendMenu(h_ellipse, MF_STRING, ID_ELLIPSE_MIDPOINT, "Midpoint Ellipse");

    AppendMenu(h_clipping, MF_STRING, ID_CLIP_RECT_POINT, "Clip Point (Rect)");
    AppendMenu(h_clipping, MF_STRING, ID_CLIP_RECT_LINE, "Clip Line (Rect)");
    AppendMenu(h_clipping, MF_STRING, ID_CLIP_RECT_POLYGON, "Clip Polygon (Rect)");
    AppendMenu(h_clipping, MF_STRING, ID_CLIP_SQUARE_POINT, "Clip Point (Square)");
    AppendMenu(h_clipping, MF_STRING, ID_CLIP_SQUARE_LINE, "Clip Line (Square)");

    // Main menu organization
    AppendMenu(h_main, MF_POPUP, reinterpret_cast<UINT_PTR>(h_color), "Color");
    AppendMenu(h_main, MF_POPUP, reinterpret_cast<UINT_PTR>(h_cursor), "Cursor");
    AppendMenu(h_main, MF_POPUP, reinterpret_cast<UINT_PTR>(h_lines), "Line");
    AppendMenu(h_main, MF_POPUP, reinterpret_cast<UINT_PTR>(h_circles), "Circle");
    AppendMenu(h_main, MF_POPUP, reinterpret_cast<UINT_PTR>(h_fills), "Fills");
    AppendMenu(h_main, MF_POPUP, reinterpret_cast<UINT_PTR>(h_ellipse), "Ellipse");
    AppendMenu(h_main, MF_POPUP, reinterpret_cast<UINT_PTR>(h_clipping), "Clipping");

    AppendMenu(h_main, MF_STRING, ID_CLEAR_SCREEN, "Clear");
    AppendMenu(h_main, MF_STRING, ID_SAVE_IMAGE, "Save Screen");
    AppendMenu(h_main, MF_STRING, ID_LOAD_IMAGE, "Load Image");

    SetMenu(hwnd, h_main);
}

LRESULT CALLBACK WndProc(const HWND hwnd, const UINT msg, const WPARAM w_param, const LPARAM l_param) {
    switch (msg) {
    case WM_COMMAND:
        switch (LOWORD(w_param)) {
        case ID_COLOR_RED:
                current_color = RGB(255, 0, 0);
                break;
        case ID_COLOR_GREEN:
                current_color = RGB(0, 255, 0);
                break;
        case ID_COLOR_BLUE:
                current_color = RGB(0, 0, 255);
                break;
        case ID_CURSOR_ARROW:
            current_cursor = LoadCursor(nullptr, IDC_ARROW);
            SetCursor(current_cursor);
            break;

        case ID_CURSOR_CROSS:
            current_cursor = LoadCursor(nullptr, IDC_CROSS);
            SetCursor(current_cursor);
            break;
        case ID_DDA_LINE:
        case ID_MIDPOINT_LINE:
        case ID_PARAMETRIC_LINE:
        case ID_DIRECT_CIRCLE:
        case ID_POLAR_CIRCLE:
        case ID_ITER_POLAR_CIRCLE:
        case ID_MIDPOINT_CIRCLE:
        case ID_MOD_MIDPOINT_CIRCLE:
        case ID_FILL_CIRCLE_LINES:
        case ID_FILL_CIRCLE_CIRCLES:
        case ID_FILL_SQUARE_HERMIT:
        case ID_FILL_RECTANGLE_BEZIER:
        case ID_FILL_CONVEX:
        case ID_FILL_NONCONVEX:
        case ID_FLOOD_RECURSIVE:
        case ID_FLOOD_NONRECURSIVE:
        case ID_CARDINAL_SPLINE:
        case ID_ELLIPSE_DIRECT:
        case ID_ELLIPSE_POLAR:
        case ID_ELLIPSE_MIDPOINT:
        case ID_CLIP_RECT_POINT:
        case ID_CLIP_RECT_LINE:
        case ID_CLIP_RECT_POLYGON:
        case ID_CLIP_SQUARE_POINT:
        case ID_CLIP_SQUARE_LINE:
                point_count = 0;
                current_algorithm = LOWORD(w_param);
                if (current_algorithm == ID_CLIP_SQUARE_POINT || current_algorithm == ID_CLIP_SQUARE_LINE) {
                    InvalidateRect(hwnd, nullptr, TRUE);
                    initialize_square_clip_window(hwnd);
                } else if (current_algorithm == ID_CLIP_RECT_POINT || current_algorithm == ID_CLIP_RECT_LINE) {
                    InvalidateRect(hwnd, nullptr, TRUE);
                    initialize_rect_clip_window(hwnd);
                }
                break;

        case ID_CLEAR_SCREEN:
            point_count = 0;
            if (loaded_bitmap) {
                DeleteObject(loaded_bitmap);
                loaded_bitmap = nullptr;
            }
            InvalidateRect(hwnd, nullptr, TRUE);
            break;

        case ID_SAVE_IMAGE:
            point_count = 0;
            save_screen_to_file(hwnd);
            break;

        case ID_LOAD_IMAGE:
            point_count = 0;
            load_image_from_file(hwnd);
            break;
        default:
                break;
        }
        break;

    case WM_SETCURSOR:
        if (LOWORD(l_param) == HTCLIENT) {
            SetCursor(current_cursor);
            return TRUE;
        }
        break;

    case WM_LBUTTONDOWN: {
        const POINT pt { LOWORD(l_param), HIWORD(l_param) };
        points[point_count] = pt;
        ++point_count;

        // For line algorithms, draw when we have 2 points
        if (point_count == 2 &&
            (current_algorithm == ID_DDA_LINE ||
             current_algorithm == ID_MIDPOINT_LINE ||
             current_algorithm == ID_PARAMETRIC_LINE)) {

            const HDC hdc {GetDC(hwnd)};
            switch (current_algorithm) {
                case ID_DDA_LINE:
                    draw_dda_line(hdc, points[0], points[1], current_color);
                    break;
                case ID_MIDPOINT_LINE:
                    draw_midpoint_line(hdc, points[0], points[1], current_color);
                    break;
                case ID_PARAMETRIC_LINE:
                    draw_parametric_line(hdc, points[0], points[1], current_color);
                    break;
                default:
                    break;
            }
            ReleaseDC(hwnd, hdc);
            point_count = 0; // Reset points for next shape
        } else if (point_count == 2 &&
            (current_algorithm == ID_DIRECT_CIRCLE ||
            current_algorithm == ID_POLAR_CIRCLE ||
            current_algorithm == ID_ITER_POLAR_CIRCLE ||
            current_algorithm == ID_MIDPOINT_CIRCLE ||
            current_algorithm == ID_MOD_MIDPOINT_CIRCLE)) {
                const HDC hdc {GetDC(hwnd)};
                const int radius {static_cast<int>(sqrt(
                    pow(points[1].x - points[0].x, 2) +
                    pow(points[1].y - points[0].y, 2)
                ))};

                switch (current_algorithm) {
                    case ID_DIRECT_CIRCLE:
                        draw_direct_circle(hdc, points[0], radius, current_color);
                        break;
                    case ID_POLAR_CIRCLE:
                        draw_polar_circle(hdc, points[0], radius, current_color);
                        break;
                    case ID_ITER_POLAR_CIRCLE:
                        draw_iterative_polar_circle(hdc, points[0], radius, current_color);
                        break;
                    case ID_MIDPOINT_CIRCLE:
                        draw_midpoint_circle(hdc, points[0], radius, current_color);
                        break;
                    case ID_MOD_MIDPOINT_CIRCLE:
                        draw_modified_midpoint_circle(hdc, points[0], radius, current_color);
                        break;
                    default:
                        break;
                }
                ReleaseDC(hwnd, hdc);
                point_count = 0;
        } else if (point_count == 3 &&
            (current_algorithm == ID_ELLIPSE_DIRECT ||
             current_algorithm == ID_ELLIPSE_POLAR ||
             current_algorithm == ID_ELLIPSE_MIDPOINT))
        {
            const HDC hdc {GetDC(hwnd)};

            // Calculate radii
            const int rx = abs(points[1].x - points[0].x);  // Distance from center to point A (x-radius)
            const int ry = abs(points[2].y - points[0].y);     // Distance from center to point B (y-radius)

            switch (current_algorithm) {
                case ID_ELLIPSE_DIRECT:
                    draw_direct_ellipse(hdc, points[0], rx, ry, current_color);
                    break;
                case ID_ELLIPSE_POLAR:
                    draw_polar_ellipse(hdc, points[0], rx, ry, current_color);
                    break;
                case ID_ELLIPSE_MIDPOINT:
                    draw_midpoint_ellipse(hdc, points[0], rx, ry, current_color);
                    break;
                default:
                    break;
            }
            ReleaseDC(hwnd, hdc);
            point_count = 0;
        } else if (point_count == 1 && (current_algorithm == ID_FLOOD_RECURSIVE || current_algorithm == ID_FLOOD_NONRECURSIVE)) {
            const HDC hdc {GetDC(hwnd)};
            switch (current_algorithm) {
                case ID_FLOOD_RECURSIVE:
                    flood_fill_recursive(hdc, points[0].x, points[0].y, bg_color, current_color);
                    break;
                case ID_FLOOD_NONRECURSIVE:
                    flood_fill_iterative(hdc, points[0].x, points[0].y, bg_color, current_color);
                    break;
                default:
                    break;
            }
        } else if (point_count == 2 && current_algorithm == ID_FILL_SQUARE_HERMIT) {
            const HDC hdc {GetDC(hwnd)};
            fill_square_with_vertical_hermite(hdc, points[0], points[1], current_color);
            ReleaseDC(hwnd, hdc);
            point_count = 0;
        } else if (point_count == 2 && current_algorithm == ID_FILL_RECTANGLE_BEZIER) {
            const HDC hdc {GetDC(hwnd)};
            fill_rectangle_with_horizontal_bezier(hdc, points[0], points[1], current_color);
            ReleaseDC(hwnd, hdc);
            point_count = 0;
        } else if (point_count == 1 && current_algorithm == ID_CLIP_SQUARE_POINT) {
            const HDC hdc {GetDC(hwnd)};
            if (clip_point_square(hdc, points[0], current_color)) {
                cout << "Point (" << points[0].x << ", " << points[0].y << ") is inside the clipping window" << endl;
            } else {
                cout << "Point (" << points[0].x << ", " << points[0].y << ") is clipped (outside window)" << endl;
            }
            ReleaseDC(hwnd, hdc);
            point_count = 0;
        }
        else if (point_count == 2 && current_algorithm == ID_CLIP_SQUARE_LINE) {
            const HDC hdc {GetDC(hwnd)};
            clip_line_square(hdc, points[0], points[1], current_color);
            cout << "Line from (" << points[0].x << ", " << points[0].y << ") to ("
                 << points[1].x << ", " << points[1].y << ") processed for clipping" << endl;
            ReleaseDC(hwnd, hdc);
            point_count = 0;
        } else if (point_count == 1 && current_algorithm == ID_CLIP_RECT_POINT) {
            const HDC hdc {GetDC(hwnd)};
            if (clip_point_rect(hdc, points[0], current_color)) {
                cout << "Point (" << points[0].x << ", " << points[0].y << ") is inside the rectangle clipping window" << endl;
            } else {
                cout << "Point (" << points[0].x << ", " << points[0].y << ") is clipped (outside rectangle window)" << endl;
            }
            ReleaseDC(hwnd, hdc);
            point_count = 0;
        }
        else if (point_count == 2 && current_algorithm == ID_CLIP_RECT_LINE) {
            const HDC hdc {GetDC(hwnd)};
            clip_line_rect(hdc, points[0], points[1], current_color);
            cout << "Line from (" << points[0].x << ", " << points[0].y << ") to ("
                 << points[1].x << ", " << points[1].y << ") processed for rectangle clipping" << endl;
            ReleaseDC(hwnd, hdc);
            point_count = 0;
        }
        //TODO: Add the other algorithms.

        break;
    }
    case WM_RBUTTONDOWN: {
        if (point_count > 2 && current_algorithm == ID_FILL_CONVEX) {
            const HDC hdc {GetDC(hwnd)};
            convex_fill(hdc, points, point_count, current_color);
            ReleaseDC(hwnd, hdc);
            point_count = 0; // Reset points after filling
        }
        else if (point_count > 2 && current_algorithm == ID_FILL_NONCONVEX) {
            const HDC hdc {GetDC(hwnd)};
            nonconvex_fill(hdc, points, point_count, current_color);
            ReleaseDC(hwnd, hdc);
            point_count = 0; // Reset points after filling
        }
        else if (point_count > 0 && current_algorithm == ID_CARDINAL_SPLINE) {
            const HDC hdc {GetDC(hwnd)};
            draw_cardinal_spline(hdc, vector(points, points + point_count), 0.5f, current_color);
            ReleaseDC(hwnd, hdc);
            point_count = 0; // Reset points after drawing
        }
        break;
    }

    case WM_PAINT: {
        PAINTSTRUCT ps;
        const HDC hdc {BeginPaint(hwnd, &ps)};

        if (loaded_bitmap) {
            const HDC hdc_mem = CreateCompatibleDC(hdc);
            const auto old_bitmap = static_cast<HBITMAP>(SelectObject(hdc_mem, loaded_bitmap));
            BITMAP bmp;
            GetObject(loaded_bitmap, sizeof(BITMAP), &bmp);
            BitBlt(hdc, 0, 0, bmp.bmWidth, bmp.bmHeight, hdc_mem, 0, 0, SRCCOPY);
            SelectObject(hdc_mem, old_bitmap);
            DeleteDC(hdc_mem);
        } else {
            RECT rect;
            GetClientRect(hwnd, &rect);
            const auto white_brush = static_cast<HBRUSH>(GetStockObject(WHITE_BRUSH));
            FillRect(hdc, &rect, white_brush);
        }

        if (current_algorithm == ID_CLIP_SQUARE_POINT || current_algorithm == ID_CLIP_SQUARE_LINE) {
            draw_square_clip_window(hdc);
        } else if (current_algorithm == ID_CLIP_RECT_POINT || current_algorithm == ID_CLIP_RECT_LINE) {
            draw_rect_clip_window(hdc);
        }

        EndPaint(hwnd, &ps);
        break;
    }

    case WM_DESTROY:
        if (loaded_bitmap) DeleteObject(loaded_bitmap);
        PostQuitMessage(0);
        break;
    default:
        return DefWindowProc(hwnd, msg, w_param, l_param);
    }
    return 0;
}

int WINAPI WinMain(const HINSTANCE hInstance, HINSTANCE, LPSTR, const int nShowCmd) {
    AllocConsole(); freopen("CONOUT$", "w", stdout);
    cout << "Program started.\n";

    WNDCLASS wc = {};
    wc.lpfnWndProc = WndProc;
    wc.hInstance = hInstance;
    wc.lpszClassName = "GraphicsApp";
    wc.hCursor = current_cursor;
    wc.hbrBackground = reinterpret_cast<HBRUSH>((COLOR_WINDOW + 1));

    RegisterClass(&wc);

    const HWND hwnd {CreateWindow("GraphicsApp", "Graphics Application", WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT, 1000, 700, nullptr, nullptr, hInstance, nullptr)};
    add_menus(hwnd);
    ShowWindow(hwnd, nShowCmd);
    UpdateWindow(hwnd);

    MSG msg;
    while (GetMessage(&msg, nullptr, 0, 0)) {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    return 0;
}
